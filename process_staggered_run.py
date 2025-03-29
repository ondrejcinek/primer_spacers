#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script assesses primers with staggered primers (i.e. inserts in primers, heterogeneity spacers).
Removes the spacers and returns a file with reads, plus some statistics.
"""

#---IMPORTS-----------------------------------------------------------------------------
import os
import re
import argparse
import json
import gzip
import copy
import random
import time
from functools import partial
import subprocess
from Bio.Seq import Seq
from Bio import SeqIO

random.seed(1000) # to ensure a reproducible choice of reads wherever a random selection is needed.



#------------------------- FUNCTIONS --------------------------------------------------------
def decode_file_name(filename_without_path):
    """returns library name, direction of sequencing and whether it has the gz end."""

    reaction_name = None
    direction = None
    gz = False

    if re.search(r"fastq.gz$", filename_without_path, re.I):
        gz = True  # ends with gz, is compressed

    result = re.search(r"([^\/\\]*)_S\d+_L\d+_R([12])_\d+\.fastq", filename_without_path, re.I)
    if result:
        name = result.group(1)
        direction = result.group(2)
        reaction_name = name
        # reaction_name = reaction_name.replace("_", "")
        reaction_name = reaction_name.replace("-", "")
        reaction_name = reaction_name.replace(".", "")
        reaction_name = reaction_name.replace(":", "")

    return reaction_name, direction, gz


#------------------------ CLASSES ----------------------------------------------------------
class SequencingTrackTable:
    """
    Keeps track of the base frequencies at a given position of the sequencing.
    """

    def __init__(self, length):
        self.track_length = length
        self.tracks = {
            "A": [0 for x in range(length)],
            "C": [0 for x in range(length)],
            "G": [0 for x in range(length)],
            "T": [0 for x in range(length)],
            "N": [0 for x in range(length)]
        }
        self.coverage = [0 for x in range(self.track_length)]
        return

    def add_seq(self, sequence):
        """
        Adds a sequence into the track table - the table with frequencies of the bases at given positions of the amplicons.
        """
        seq_id = "NA"
        if isinstance(sequence, Seq):
            seq_id = sequence.id
            sequence = str(sequence)
            #print(sequence)
            #print(len(sequence))

        if len(sequence) > self.track_length:
            raise ValueError("The length of the added sequence %s is larger than the length of the track table: %d > %d" % (seq_id, len(sequence), self.track_length))

        # now it is a text
        allowed_bases = self.tracks.keys()
        for i, base in enumerate(sequence):
            if base in allowed_bases:
                self.tracks[base][i] += 1
                self.coverage[i] += 1

        return

    def __len__(self):
        return self.track_length


    def __add__(self, other):
        if len(self) != len(other):
            raise ValueError("It has no sense to add SequencingTrackTables of unequal lengths")
        result = SequencingTrackTable(len(self))
        for i in range(len(result)):
            for base in ("A", "C", "G", "T", "N"):
                result.tracks[base] = self.tracks[base] + other.tracks[base]
        return result

    def get_proportions(self):
        proportions = copy.deepcopy(self.tracks)
        for i in range(self.track_length):
            for base in ("A", "C", "G", "T", "N"):
                proportions[i][base] = proportions[i][base]/self.coverage[i]
        return proportions

    def get_coverage(self):
        return(copy(self.coverage))

    def tabify(self, name, direction):
        result = []
        for position in range(self.track_length):
            for base in ("A", "C", "G", "T", "N"):
                line = "%s\t%s\t%d\t%s\t%d\t%d" % (name, direction, position, base, self.tracks[base][position], self.coverage[position])
                result.append(line)
        return result


class StaggeredDataset:
    """
    Analysed library set of staggered read files.
    The whole process is launched upon initialization
    """
    def __init__(self, source_dir, target_dir, track_tables_fn = None, primer_description_file = "primers_16SV34_staggered.txt",
                 write_to_files_of_combinations = False, trim_only_spacers = False,  trim_whole_primer = True,
                 trim_last_bases_r1 = 0, trim_last_bases_r2 = 0,
                 merge_pairs_by_usearch = False, usearch = "H:/mikro/programs/nextgen/automating_16S/usearch11.0.667_win32.exe",
                 read_length=300,
                 fastq_maxmergelen = 490, fastq_minmergelen = 350,
                 fastq_maxdiffs = 4, max_reads_per_sample = 0,
                 do_not_overwrite = 0):
        """
        Loads the dictionary with staggered reads and prepares for analysis
        :param source_dir: source directory with the staggered read fastq files
        :param target_dir: directory for the final fastq files
        :param track_tables_fn: output file for storing tables with proportions of signal in libraries
        :param primer_description_file: set of primers.
            One primer on each line. No header.
            name, sequence, f/r, primer mix name, length of the spacer
                columns separated by whitespace
            By default it is this:
            Bact341_0_for5 CCTACGGGAGGCAGCAG f v34 0
            Bact341_2_for5 gaCCTACGGGAGGCAGCAG f v34 2
            Bact341_3_for5 tagCCTACGGGAGGCAGCAG f v34 3
            Bact341_6_for5 agcaattCCTACGGGAGGCAGCAG f v34 7
            Bact806_0_rev7 GGACTACHVGGGTWTCTAAT r v34 0
            Bact806_2_rev7 caGGACTACHVGGGTWTCTAAT r v34 2
            Bact806_3_rev7 tctGGACTACHVGGGTWTCTAAT r v34 3
            Bact806_4_rev7 atctGGACTACHVGGGTWTCTAAT r v34 4
            B806_341_0_for5	GGACTACHVGGGTWTCTAAT f v43 0
            B806_341_2_for5	caGGACTACHVGGGTWTCTAAT f v43 2
            B806_341_3_for5	atctGGACTACHVGGGTWTCTAAT f v43 4
            B806_341_6_for5	tctactGGACTACHVGGGTWTCTAAT f v43 6
            B806_341_0_rev7	CCTACGGGAGGCAGCAG r v43 0
            B806_341_2_rev7	gaCCTACGGGAGGCAGCAG r v43 2
            B806_341_3_rev7	tagCCTACGGGAGGCAGCAG r v43 3
            B806_341_6_rev7	agcaattCCTACGGGAGGCAGCAG r v43 7
        :param write_to_files_of_combinations: if true, makes many small files for specific mix-primer combinations
        :param trim_only_spacers: what to trim - only the heterogeneity spacers, whereas while primers are retained?
        :param trim_whole_primer: What to trim - trim the whole primers?
        :param trim_last_bases_r1: forward read - how many last bases should we trim from R1?
        :param trim_last_bases_r2: reverse read - how many last bases should we trim from R2?

        The init returns nothing. It launches every method needed and produces the result.

        """

        #------------------- Get parameters ---------------------------
        self.source_dir = source_dir
        self.target_dir = target_dir
        self.track_tables_fn = track_tables_fn
        self.primer_description_file = primer_description_file
        self.write_to_files_of_combinations = write_to_files_of_combinations

        self.trim_only_spacers = trim_only_spacers
        self.trim_whole_primer = trim_whole_primer
        self.trim_last_bases_r1 = trim_last_bases_r1
        self.trim_last_bases_r2 = trim_last_bases_r2
        if trim_only_spacers and trim_whole_primer:
            raise ValueError("This combination of parameters is not allowed: trim_only_spacers = True,  trim_whole_primer = True")

        self.merge_pairs_by_usearch = merge_pairs_by_usearch
        if self.merge_pairs_by_usearch == "0":
            self.merge_pairs_by_usearch = False
        self.usearch = usearch

        self.read_length = read_length
        self.fastq_minmergelen = fastq_minmergelen
        self.fastq_maxmergelen = fastq_maxmergelen
        self.fastq_maxdiffs = fastq_maxdiffs
        self.max_reads_per_sample = max_reads_per_sample
        self.do_not_overwrite = do_not_overwrite

        #------------------ parameter check, initialization and preparatory steps ------

        # Make directories
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)
        self.dir_01_r1_r2_trimmed_fastq_by_library = os.path.join(self.target_dir, "01_r1_r2_trimmed_fastq_by_library")
        if not os.path.exists(self.dir_01_r1_r2_trimmed_fastq_by_library):
            os.makedirs(self.dir_01_r1_r2_trimmed_fastq_by_library)

        self.dir_01_r1_r2_fastq_by_combination = os.path.join(self.target_dir, "01_r1_r2_trimmed_fastq_by_library", "by_spacer_combination")
        if write_to_files_of_combinations and not os.path.exists(self.dir_01_r1_r2_fastq_by_combination):
            os.makedirs(self.dir_01_r1_r2_fastq_by_combination)

        # This is only if we merge the two reads by USEARCH
        self.dir_02_merged_by_library = os.path.join(self.target_dir, "02_merged_by_library")
        if self.merge_pairs_by_usearch:
            if not os.path.exists(usearch) or not os.path.isfile(usearch):
                raise ValueError("The USEARCH program has not been found. Please install it from https://www.drive5.com/usearch/download.html \n"
                                 "and supply the path as an usearch argument to this program.\n"
                                 "The free 32-bit version is sufficient.")
            if self.merge_pairs_by_usearch and not os.path.exists(self.dir_02_merged_by_library):
                os.makedirs(os.path.join(self.target_dir, self.dir_02_merged_by_library))
            # a file for usearch protocols
            self.usearch_protocol_file_fn = os.path.join(self.dir_02_merged_by_library, "_usearch_protocol.txt")

        # --------- data on libraries ------------
        self.last_started_library = None
        self.last_processed_library = None
        self.libraries = dict()  # structure for the libraries in the set
        self.track_tables_r1 = {}
        self.track_tables_r2 = {}

        # ------------ Primers -------------------
        if self.primer_description_file == "primers_16SV34_staggered.txt" and (not os.path.exists(self.primer_description_file) or not os.path.isfile(self.primer_description_file)):
            self.write_primers_description_txt(primer_set_name="primer_set_16SV34")
        if self.primer_description_file == "primers_blasto_staggered.txt" and (not os.path.exists(self.primer_description_file) or not os.path.isfile(self.primer_description_file)):
            self.write_primers_description_txt(primer_set_name="blastocystis_subtyping")
        self.primers = []
        # Here we store the signature sequences to tell the "mix and spacer" combinations apart. Is populated by the load_primer_set method.
        self.mix_offset_signature_sequences_f = {}
        self.mix_offset_signature_sequences_f = {}
        # .... {"CCTACGG": "v34_0", "GACCTAC": "v34_2", "TAGCCTA": "v34_3", "AGCAATT": "v34_6",
        #       "GGACTAC": "v43_0", "CAGGACT": "v43_2", "ATCTGGA": "v43_3", "TCTACTG": "v43_6"}

        #This is the overall repertoire what we have. Also populated by the load_primer_set method
        self.mix_offset_combinations = []
        self.overall_count_of_reads_per_mix_offset_combination = {}
        #{"v34_0xv34_0": 0, "v34_0xv34_2": 0, "v34_0xv34_3": 0, "v34_0xv34_4": 0, ......
        # ...... "v43_6xv43_0": 0, "v43_6xv43_2": 0, "v43_6xv43_3": 0, "v43_6xv43_6": 0, "other": 0}
        # the "other" will encompass all that could not be identified.
        self.load_primer_set()  # get the primers

        #------------------------------------------------------------
        # Actual processing of the data
        # ------------------------------------------------------------
        self.usearch_output = ""

        self.find_libraries_and_their_files() #find libraries among the files from the sequencer
        self.process_staggered_primers() #process the staggered primers (including the merging)
        self.save_track_tables()
        return

    def write_primers_description_txt(self, primer_set_name):
        """ This function writes a table with primers, if not found.

        The file describing primers is a header-less text file with columns separated by a whitespace.
        There are 6 columns:
        1) primer name
        2) primer sequence
        3) f = forward, r = reverse
        4) reaction (amplicon) name - for instance, the V3-V4 region is amplified in two reactions,
                    one in the direction of V3 -> V4, the other in the direction of V4 -> V3
        5) count of inserted bases (the length of the spacer)

        """
        primer_sets = dict(
        primer_set_16SV34 = """Bact341_0_for5 CCTACGGGAGGCAGCAG f v34 0
Bact341_2_for5 gaCCTACGGGAGGCAGCAG f v34 2
Bact341_3_for5 tagCCTACGGGAGGCAGCAG f v34 3
Bact341_6_for5 agcaattCCTACGGGAGGCAGCAG f v34 7
Bact806_0_rev7 GGACTACHVGGGTWTCTAAT r v34 0
Bact806_2_rev7 caGGACTACHVGGGTWTCTAAT r v34 2
Bact806_3_rev7 tctGGACTACHVGGGTWTCTAAT r v34 3
Bact806_4_rev7 atctGGACTACHVGGGTWTCTAAT r v34 4
B806_341_0_for5 GGACTACHVGGGTWTCTAAT f v43 0
B806_341_2_for5 caGGACTACHVGGGTWTCTAAT f v43 2
B806_341_3_for5 atctGGACTACHVGGGTWTCTAAT f v43 4
B806_341_6_for5 tctactGGACTACHVGGGTWTCTAAT f v43 6
B806_341_0_rev7 CCTACGGGAGGCAGCAG r v43 0
B806_341_2_rev7 gaCCTACGGGAGGCAGCAG r v43 2
B806_341_3_rev7 tagCCTACGGGAGGCAGCAG r v43 3
B806_341_6_rev7 agcaattCCTACGGGAGGCAGCAG r v43 7
""",
        blastocystis_subtyping = """BT_F0_for5 GGAGGTAGTGACAATAAATC f reaction1_direct 0
BT_F1_for5 aGGAGGTAGTGACAATAAATC f reaction1_direct 1
BT_F2_for5 caGGAGGTAGTGACAATAAATC f reaction1_direct 2
BT_F3_for5 actGGAGGTAGTGACAATAAATC f reaction1_direct 3
BT_R0_rev7 TGCTTTCGCACTTGTTCATC r reaction1_direct 0
BT_R1_rev7 aTGCTTTCGCACTTGTTCATC r reaction1_direct 1
BT_R2_rev7 caTGCTTTCGCACTTGTTCATC r reaction1_direct 2
BT_R3_rev7 actTGCTTTCGCACTTGTTCATC r reaction1_direct 3
BT_R0_for5 TGCTTTCGCACTTGTTCATC f reaction2_reversed 0
BT_R1_for5 aTGCTTTCGCACTTGTTCATC f reaction2_reversed 1
BT_R2_for5 caTGCTTTCGCACTTGTTCATC f reaction2_reversed 2
BT_R3_for5 actTGCTTTCGCACTTGTTCATC f reaction2_reversed 3
BT_F0_rev7 GGAGGTAGTGACAATAAATC r reaction2_reversed 0
BT_F1_rev7 aGGAGGTAGTGACAATAAATC r reaction2_reversed 1
BT_F2_rev7 caGGAGGTAGTGACAATAAATC r reaction2_reversed 2
BT_F3_rev7 actGGAGGTAGTGACAATAAATC r reaction2_reversed 3
"""
        )

        with open(self.primer_description_file, mode = "w") as outhandle:
            print(primer_sets[primer_set_name], file=outhandle)

        return


    def load_primer_set(self):
        """Loads the primer info from a file.
        the primers are one each line,
        name, sequence, f/r, primer mix name, length of the spacer
        columns separated by whitespace

        This is an example of the input primer file for two primer mixes: v34 and v43.
        Bact341_0_for5 CCTACGGGAGGCAGCAG f v34 0
        Bact341_2_for5 gaCCTACGGGAGGCAGCAG f v34 2
        .....
        B806_341_0_rev7	CCTACGGGAGGCAGCAG r v43 0
        B806_341_2_rev7	gaCCTACGGGAGGCAGCAG r v43 2

        The method parses the file, decides on the signature sequences
        This is how the primers look like in the self.primers list of dictionaries
        [{"name": "Bact341_0_for5", "seq_orig": "CCTACGGGAGGCAGCAG", "fr": "f", "mix": "v34", "offset": 0, "sequence": "CCTACGGGAGGCAGCAG", "max_inf_length": 17},
         {"name": "Bact341_2_for5", "seq_orig": "gaCCTACGGGAGGCAGCAG", "fr": "f", "mix": "v34", "offset": 2, "sequence": "GACCTACGGGAGGCAGCAG", "max_inf_length": 19},
         ....
         {"name": "B806_341_6_rev7", "seq_orig": "agcaattCCTACGGGAGGCAGCAG", "fr": "r", "mix": "v43", "offset": 6, "sequence": "AGCAATTCCTACGGGAGGCAGCAG", "max_inf_length": 24}]
        """


        with open(self.primer_description_file) as inhandle:
            lines = inhandle.readlines()
        lines = [line.rstrip("\n") for line in lines]
        lines = [line for line in lines if line != ""]
        self.primers = [dict(zip(["name","seq_orig","fr","mix","offset"], x.split())) for x in lines]
        for primer in self.primers: # no underscores in the mix names!
            primer["mix"] = re.sub("_", "", primer["mix"])
        print(self.primers)

        for primer in self.primers:
            primer["offset"] = int(primer["offset"])
            primer["sequence"] = primer["seq_orig"].upper()
            res1 = re.search(r"[^ACGT]", primer["sequence"], re.I) #find leftmost degenerate base
            if res1:
                primer["max_inf_length"] = res1.span(0)[0] #the maximum informative length is up to that base.
            else:
                primer["max_inf_length"] = len(primer["sequence"]) #otherwise the maximum informative length is the whole primer
        del(primer) #deletes the useless persisting variable.

        #what is the longest common nondegenrate informative part iof primers. In other words, minimum of max length...
        self.max_inf_f = min([x["max_inf_length"] for x in self.primers if x["fr"] == "f"])
        self.max_inf_r = min([x["max_inf_length"] for x in self.primers if x["fr"] == "r"])

        #maximum spacer length
        self.max_spacer_length_f = max([x["offset"] for x in self.primers if x["fr"] == "f"])
        self.max_spacer_length_r = max([x["offset"] for x in self.primers if x["fr"] == "r"])

        #get the sequence of spacers + adjacent primers. The "mix" here means either v34 or v43; mix_offset is v34_0 or v43_7 etc.
        self.mix_offset_signature_sequences_f = {x["sequence"][:self.max_inf_f]:"%s_%d"%(x["mix"], x["offset"]) for x in self.primers if x["fr"] == "f"}
        self.mix_offset_signature_sequences_r = {x["sequence"][:self.max_inf_r]:"%s_%d"%(x["mix"], x["offset"]) for x in self.primers if x["fr"] == "r"}

        self.mix_offset_primer_length_f = {"%s_%d"%(x["mix"], x["offset"]):len(x["sequence"]) for x in self.primers if x["fr"] == "f"}
        self.mix_offset_primer_length_r = {"%s_%d"%(x["mix"], x["offset"]):len(x["sequence"]) for x in self.primers if x["fr"] == "r"}

        if len(set(self.mix_offset_signature_sequences_f)) < len(self.mix_offset_signature_sequences_f):
            raise ValueError("The forward start of mix-spacer combinations cannot be told apart. Either some have the same sequence, or they are too short before the first redundancy.")
        if len(set(self.mix_offset_signature_sequences_r)) < len(self.mix_offset_signature_sequences_r):
            raise ValueError("The reverse start of mix-spacer combinations cannot be told apart. Either some have the same sequence, or they are too short before the first redundancy.")

        # Combination to be reported
        self.mix_offset_combinations = []

        for primer_f in self.primers:
            if primer_f["fr"] != "f":
                continue
            #print(primer_f)
            for primer_r in self.primers:
                if primer_r["fr"] != "r" or primer_f["mix"] != primer_r["mix"]:
                    continue
                #print("...",primer_r)

                #==============================================================================================================================
                mix_offset_combination = "%s_%dx%s_%d" % (primer_f["mix"], primer_f["offset"], primer_r["mix"], primer_r["offset"])
                self.mix_offset_combinations.append(mix_offset_combination)
                self.track_tables_r1[mix_offset_combination] = SequencingTrackTable(305) #to accommodate for the longest possible read.
                self.track_tables_r2[mix_offset_combination] = SequencingTrackTable(305)
                # ==============================================================================================================================
        self.track_tables_r1["other"] = SequencingTrackTable(305)
        self.track_tables_r2["other"] = SequencingTrackTable(305)

        self.overall_count_of_reads_per_mix_offset_combination = dict.fromkeys(self.mix_offset_combinations, 0)
        self.overall_count_of_reads_per_mix_offset_combination["other"] = 0

        #print(json.dumps(self.mix_offset_signature_sequences_f))
        #print(json.dumps(self.overall_count_of_reads_per_mix_offset_combination))
        #print(json.dumps(self.primers))

        return

    def save_track_tables(self):
        if self.track_tables_fn is None or self.track_tables_fn == "":
            return

        tt_dirs = os.path.split(self.track_tables_fn)[0]
        if not os.path.exists(tt_dirs):
            os.makedirs(tt_dirs)

        with open(self.track_tables_fn, mode="w") as fh:
            print("\t".join(["mix_offset_combination", "direction", "pos", "base", "count", "coverage"]), file=fh)
            for mix_offset_combination, table in self.track_tables_r1.items():
                for line in table.tabify(mix_offset_combination, "R1"):
                    print(line, file=fh)
            for mix_offset_combination, table in self.track_tables_r2.items():
                for line in table.tabify(mix_offset_combination, "R2"):
                    print(line, file=fh)

        return

    def find_libraries_and_their_files(self):
        """
        Get the repertoire of sequenced reactions from the sequencer file
        """

        # ---------- Prepare for the processing - make a new directory for split files
        os.chdir(self.source_dir)

        # ----------- Find the files in the hierarchical structure
        for (dirpath, dirnames, filenames) in os.walk(self.source_dir, topdown=True):
            fastq_filenames = [x for x in filenames if re.search(r"(\.fastq|\.fq)(?:\.gz)?$", x, re.I)]
            # print(fastq_filenames)
            for filename in fastq_filenames:
                res0 = re.search(r"^([^\/\\]*)_(S\d+)_L001_R([12])_001\.fastq(?:\.gz)?", filename, re.I)
                if not res0:
                    raise ValueError("Wrong name format of the fastq file: " + filename)
                libname = res0.group(1)
                r1or2 = res0.group(3)              
                
                if libname not in self.libraries: #create an empty structure for read counts per mix_offset_combinations
                    self.libraries[libname] = dict(reads_by_mix_offset_combination = dict.fromkeys(self.mix_offset_combinations, 0))
                    self.libraries[libname]["reads_by_mix_offset_combination"]["other"] = 0
                    self.libraries[libname]["spacer_combination_files"] = dict()
                    self.libraries[libname]["fn_merged_by_library"] = os.path.join(self.dir_02_merged_by_library, "%s_merged.fastq" % libname)
                    
                if filename.endswith('.gz'):
                    out_filename = filename[:-3]
                else:
                    out_filename = filename

                self.libraries[libname]["fn_reads" + r1or2] = dict(
                    fn_source = os.path.join(self.source_dir, dirpath, filename),
                    fn_target = os.path.join(self.dir_01_r1_r2_trimmed_fastq_by_library, out_filename),
                    r1or2 = r1or2,
                    libname = libname,
                    s_number = res0.group(2)
                )


        for libname in self.libraries.keys():
            if "fn_reads1" not in self.libraries[libname] or "fn_reads2" not in self.libraries[libname]:
                raise ValueError("The library does not have both R1 and R2 fastq files " + libname)

        #print(self.libraries)
        return


    def process_staggered_primers(self):
        """
        Process the staggered run.
        Produces the results file.
        :return: nothing
        """

        start_time = time.time()
        previous_library_time = start_time
        libno = 0
        for libname in self.libraries:
            self.last_started_library = libname
            now_time = time.time()
            libno += 1
            print(150*"-")
            print("|   THIS IS LIBRARY %10s, %d of %d (%.1f%%) in this run" % (libname, libno, len(self.libraries), 100*libno/len(self.libraries)))
            print("|   %.1f minutes since the last library and %.1f min since the processing started" % ((now_time-previous_library_time)/60, (now_time-start_time)/60))
            print(150 * "-")
            previous_library_time = now_time

            if libname == "Undetermined":
                print("Skipping the Undetermined reads file")
                previous_library_time = now_time
                continue

            if self.do_not_overwrite:
                if (os.path.exists(self.libraries[libname]["fn_reads1"]["fn_target"]) and os.path.isfile(self.libraries[libname]["fn_reads1"]["fn_target"])
                    and os.path.exists(self.libraries[libname]["fn_reads2"]["fn_target"]) and os.path.isfile(self.libraries[libname]["fn_reads2"]["fn_target"])
                    and (not self.merge_pairs_by_usearch or (os.path.exists(self.libraries[libname]["fn_merged_by_library"]) and os.path.isfile(self.libraries[libname]["fn_merged_by_library"])))
                ):
                    self.last_processed_library = libname
                    print("All files have been found in place.")
                    print("Parameter do_not_overwrite is True. Skipping.")
                    continue

            self.libraries[libname]["report_file"] = os.path.join(self.dir_01_r1_r2_trimmed_fastq_by_library, libname + "_report.txt")
            with open (self.libraries[libname]["report_file"], mode="w") as report_handle:
                print("LIBRARY %s" % libname, file = report_handle)


            # ------ Prepare output filenames ------
            if self.write_to_files_of_combinations:
                # first open combination-specific output files:
                for mix_offset_combination in self.mix_offset_combinations:
                    fname1 = os.path.join(self.dir_01_r1_r2_fastq_by_combination, libname + "_" + mix_offset_combination + "_R1.fastq")
                    fname2 = os.path.join(self.dir_01_r1_r2_fastq_by_combination, libname + "_" + mix_offset_combination + "_R2.fastq")
                    self.libraries[libname]["spacer_combination_files"][mix_offset_combination] = {
                        "r1_fname": fname1, "r1_fh": open(fname1, mode="w"),
                        "r2_fname": fname2, "r2_fh": open(fname2, mode="w")
                    }

                # then also files for the unidentified spacers
                self.libraries[libname]["spacer_combination_files"]["others"] = {
                    "r1_fname": os.path.join(self.dir_01_r1_r2_fastq_by_combination, libname + "_others_R1.fastq"),
                    "r2_fname": os.path.join(self.dir_01_r1_r2_fastq_by_combination, libname + "_others_R2.fastq")
                }
                self.libraries[libname]["spacer_combination_files"]["others"]["r1_fh"] = open(self.libraries[libname]["spacer_combination_files"]["others"]["r1_fname"], mode="w")
                self.libraries[libname]["spacer_combination_files"]["others"]["r2_fh"] = open(self.libraries[libname]["spacer_combination_files"]["others"]["r2_fname"], mode="w")

            #--------------------------------------------------------------------------------
            # Open the source files and parse reads
            #--------------------------------------------------------------------------------

            #partial function used for opening the gzip file. if it is not gzip we use ordinary open
            _open_reads_f = partial(gzip.open, mode="rt") if self.libraries[libname]["fn_reads1"]["fn_source"][-3:] == ".gz" else open
            _open_reads_r = partial(gzip.open, mode="rt") if self.libraries[libname]["fn_reads2"]["fn_source"][-3:] == ".gz" else open
            with  _open_reads_f(self.libraries[libname]["fn_reads1"]["fn_source"]) as fh_reads1,\
                  _open_reads_r(self.libraries[libname]["fn_reads2"]["fn_source"]) as fh_reads2:

                # get the source reads: all, or just a subset
                source_reads1 = []
                source_reads2 = []

                #-------------- ITERATE THROUGH READS ------------------
                total_read_count = 0
                parser_r1 = SeqIO.parse(fh_reads1, "fastq")
                parser_r2 = SeqIO.parse(fh_reads2, "fastq")
                while True:
                    try:
                        read1 = next(parser_r1) #.pop(0)
                    except StopIteration as e:
                        read1 = None
                    try:
                        read2 = next(parser_r2) #.pop(0)
                    except StopIteration as e:
                        read2 = None
                    #if any failed to read from the parser
                    if read1 is None or read1 is None:
                        if read1 is not None or read2 is not None: #verify that this happened at both R1 and R2 ends
                            raise ValueError("Uneven ends of fastq files with R1 and R2 reads in library " + libname)
                        else: #end of iteration in both R1 and R2 read files.
                            break

                    if read1.id != read2.id:
                        raise ValueError("R1 and R2 fastq read names do not correspond in library %s, read number %d" % (libname, total_read_count))

                    # reads OK, append into the list
                    total_read_count += 1
                    source_reads1.append(read1)
                    source_reads2.append(read2)
                # end of reading the source files. Now all are in source_reads1 and source_reads2

            # --------------------------------------------------------------------------------
            # Reduce the reads on input to the required count
            # --------------------------------------------------------------------------------
            reduced_read_count = None
            if self.max_reads_per_sample is not None and total_read_count > self.max_reads_per_sample:
                indices = sorted(random.sample(range(total_read_count), k=self.max_reads_per_sample))
                source_reads1 = [source_reads1[i] for i in indices]
                source_reads2 = [source_reads2[i] for i in indices]
                with open(self.libraries[libname]["report_file"], mode="a") as report_handle:
                    report_line = "The total count of %d reads has been reduced to the requested %d." %(total_read_count, self.max_reads_per_sample)
                    print(report_line, file=report_handle)
                    print(report_line)
                reduced_read_count = self.max_reads_per_sample
            else:
                with open(self.libraries[libname]["report_file"], mode="a") as report_handle:
                    report_line = "All %d reads will be passed downstream" % total_read_count
                    print(report_line, file=report_handle)
                    print(report_line)
                reduced_read_count = total_read_count

            # --------------------------------------------------------------------------------
            # Sort into categories and trim of spacers and primers
            # --------------------------------------------------------------------------------
            collection_of_trimmed_reads1 = []
            collection_of_trimmed_reads2 = []
            for i in range(len(source_reads1)):
                read1 = source_reads1[i]
                read2 = source_reads2[i]

                #what mix and length of offset do we have here?
                f_signature_sequence = str(read1.seq)[:self.max_inf_f]
                r_signature_sequence = str(read2.seq)[:self.max_inf_r]
                f_mix_and_offset = self.mix_offset_signature_sequences_f.get(f_signature_sequence, "other")
                r_mix_and_offset = self.mix_offset_signature_sequences_r.get(r_signature_sequence, "other")

                #Identify the combination of left and right
                identified_mix_offset_combination = f_mix_and_offset + "x" + r_mix_and_offset
                if identified_mix_offset_combination not in self.mix_offset_combinations:
                    identified_mix_offset_combination = "other"
                else:
                    mix_f, f_len_spacer = f_mix_and_offset.split("_")
                    mix_r, r_len_spacer = r_mix_and_offset.split("_")

                self.libraries[libname]["reads_by_mix_offset_combination"][identified_mix_offset_combination] += 1
                self.overall_count_of_reads_per_mix_offset_combination[identified_mix_offset_combination] += 1
                self.track_tables_r1[identified_mix_offset_combination].add_seq(read1)
                self.track_tables_r2[identified_mix_offset_combination].add_seq(read2)

                if identified_mix_offset_combination in self.mix_offset_combinations and mix_f == mix_r: #i.e. if it is not "other" and will be sorted into fastq
                    trim_len_f = 0
                    trim_len_r = 0
                    if self.trim_whole_primer:
                        trim_len_f=self.mix_offset_primer_length_f[f_mix_and_offset]
                        trim_len_r=self.mix_offset_primer_length_r[r_mix_and_offset]
                    elif self.trim_only_spacers:
                        trim_len_f = int(f_len_spacer)
                        trim_len_r = int(r_len_spacer)

                    trimmed_read1 = read1[trim_len_f:]
                    trimmed_read2 = read2[trim_len_r:]

                    if self.trim_last_bases_r1 > 0:
                        trimmed_read1 = trimmed_read1[:-self.trim_last_bases_r1]
                    if self.trim_last_bases_r2 > 0:
                        trimmed_read2 = trimmed_read2[:-self.trim_last_bases_r2]

                    # Writing the reads
                    # For mix v34, the reads will be writted as they are. For mix v43 the reads will be switched
                    # in order to keep the whole resulting fastq file consistent: one fastq starting in v3, the other starting in v4 reverse
                    if (mix_f == mix_r and (mix_f == "v43" or mix_f == "reaction2reversed")):
                        #swap them
                        (trimmed_read2, trimmed_read1) = (trimmed_read1, trimmed_read2)
                    elif (mix_f == mix_r and (mix_f == "v34" or mix_f == "reaction1direct")):
                        pass
                    else:
                        raise ValueError("Unexpected combination of mixes")

                    # Store the trimmed reads to write them later
                    collection_of_trimmed_reads1.append(trimmed_read1)
                    collection_of_trimmed_reads2.append(trimmed_read2)

                    if self.write_to_files_of_combinations:
                        #here it writes into the spacer-combination-specific fastq files.
                        SeqIO.write(
                            trimmed_read1, #trim the spacer
                            self.libraries[libname]["spacer_combination_files"][identified_mix_offset_combination]["r1_fh"],
                            "fastq"
                        )
                        SeqIO.write(
                            trimmed_read2,#trim the spacer
                            self.libraries[libname]["spacer_combination_files"][identified_mix_offset_combination]["r2_fh"],
                            "fastq"
                        )
                else: #not identified, falls into the "other" category
                    if self.write_to_files_of_combinations:
                        # no trimming, primers not identified (not both...)
                        SeqIO.write(read1, self.libraries[libname]["spacer_combination_files"]["others"]["r1_fh"], "fastq")
                        SeqIO.write(read2, self.libraries[libname]["spacer_combination_files"]["others"]["r2_fh"], "fastq")
                # end of the block "read by read"

            # --------------------------------------------------------------------------------
            # Write the trimmed reads
            # --------------------------------------------------------------------------------
            with  open(self.libraries[libname]["fn_reads1"]["fn_target"], mode="w") as fh_trimmed1, \
                  open(self.libraries[libname]["fn_reads2"]["fn_target"], mode="w") as fh_trimmed2:
                    # Write the accumulated reads to the output files
                    SeqIO.write(collection_of_trimmed_reads1, fh_trimmed1, "fastq")
                    SeqIO.write(collection_of_trimmed_reads2, fh_trimmed2, "fastq")


            # --------------------------------------------------------------------------------
            #output of summary statistics:
            # --------------------------------------------------------------------------------
            with open(self.libraries[libname]["report_file"], mode="a") as report_handle:
                report_line = "Absolute read counts - sorted:"  + json.dumps(self.libraries[libname]["reads_by_mix_offset_combination"])
                print(report_line, file=report_handle)
                print(report_line)

                print("Relative read counts", file=report_handle)
                print("Relative read counts")
                for a_key, a_value in self.libraries[libname]["reads_by_mix_offset_combination"].items():
                    percentage = "NA"
                    if reduced_read_count > 0:
                        percentage = "%5.2f%%" % (100 * a_value / reduced_read_count)
                    report_line = "   category %15s: abs = %6d,   %s" % (
                        a_key,
                        a_value,
                        percentage
                    )
                    print(report_line, file=report_handle)
                    print(report_line)

            #close all spacer-combination specific files
            if self.write_to_files_of_combinations:
                for one_combination in self.mix_offset_combinations:
                    self.libraries[libname]["spacer_combination_files"][one_combination]["r1_fh"].close()
                    self.libraries[libname]["spacer_combination_files"][one_combination]["r2_fh"].close()
                    del(self.libraries[libname]["spacer_combination_files"][one_combination]["r1_fh"])
                    del(self.libraries[libname]["spacer_combination_files"][one_combination]["r2_fh"])
                    if os.path.getsize(self.libraries[libname]["spacer_combination_files"][one_combination]["r1_fname"]) == 0:
                        os.remove(self.libraries[libname]["spacer_combination_files"][one_combination]["r1_fname"])
                    if os.path.getsize(self.libraries[libname]["spacer_combination_files"][one_combination]["r2_fname"]) == 0:
                        os.remove(self.libraries[libname]["spacer_combination_files"][one_combination]["r2_fname"])

                self.libraries[libname]["spacer_combination_files"]["others"]["r1_fh"].close()
                self.libraries[libname]["spacer_combination_files"]["others"]["r2_fh"].close()
                del(self.libraries[libname]["spacer_combination_files"]["others"]["r1_fh"])
                del(self.libraries[libname]["spacer_combination_files"]["others"]["r2_fh"])
                if os.path.getsize(self.libraries[libname]["spacer_combination_files"]["others"]["r1_fname"]) == 0:
                    os.remove(self.libraries[libname]["spacer_combination_files"]["others"]["r1_fname"])
                if os.path.getsize(self.libraries[libname]["spacer_combination_files"]["others"]["r2_fname"]) == 0:
                    os.remove(self.libraries[libname]["spacer_combination_files"]["others"]["r2_fname"])

            # --------------------------------------------------------------------------------
            # Merge the left and right output by using USEARCH
            # --------------------------------------------------------------------------------
            if self.merge_pairs_by_usearch:
                command_text = (
                        self.usearch +
                        ' -fastq_mergepairs %s -reverse %s -minhsp 10 -fastq_minovlen 10 -fastq_minmergelen %d  -fastq_maxmergelen %d -fastq_maxdiffs %d -fastqout %s' %
                        (self.libraries[libname]["fn_reads1"]["fn_target"], self.libraries[libname]["fn_reads2"]["fn_target"],
                         self.fastq_minmergelen, self.fastq_maxmergelen, self.fastq_maxdiffs,
                         self.libraries[libname]["fn_merged_by_library"]))
                # print(command_text)
                # RUN!
                this_usearch_command_output = str(subprocess.run(command_text, shell=True, capture_output=True, text=True).stderr)
                self.usearch_output += "\n" + this_usearch_command_output
                print(this_usearch_command_output)
                with open(self.usearch_protocol_file_fn, mode="a") as outhandle:
                    outhandle.write(this_usearch_command_output)
                with open(self.libraries[libname]["report_file"], mode="a") as report_handle:
                    report_handle.write(this_usearch_command_output)

            #-------------------------------------------------------------------------------
            # This is done after each library, so that in case of problems the process can be resumed
            #-------------------------------------------------------------------------------

            # printing + saving provisional results of counting reads - after the last library this will become final...
            #with open(os.path.join(self.dir_01_r1_r2_trimmed_fastq_by_library, "by_files_result.json"), mode="w") as outhandle:
            #    json.dump(self.libraries, outhandle)
            with open(os.path.join(self.dir_01_r1_r2_trimmed_fastq_by_library, "summarized_result.json"), mode="w") as outhandle:
                json.dump(self.overall_count_of_reads_per_mix_offset_combination, outhandle)
            with open(os.path.join(self.dir_01_r1_r2_trimmed_fastq_by_library, "summarized_result.txt"), mode="w") as outhandle:
                    print("combination\tcount", file=outhandle)
                    for the_combination in list(self.overall_count_of_reads_per_mix_offset_combination.keys()) + ["other"]:
                        print(the_combination, self.overall_count_of_reads_per_mix_offset_combination[the_combination], sep="\t", file=outhandle)
            with open(os.path.join(self.dir_01_r1_r2_trimmed_fastq_by_library, "result_of_identifying_staggered_primers.txt"), mode="w") as outhandle:
                print("library\tcombination\tcount", file=outhandle)
                for libname in self.libraries.keys():
                    for the_combination in self.mix_offset_combinations + ["other"]:
                        print(libname, the_combination, self.libraries[libname]["reads_by_mix_offset_combination"][the_combination] , sep="\t", file=outhandle)

            self.last_processed_library = libname
            # end of library-by-library cycle
        return


def process_args_and_run():

    parser = argparse.ArgumentParser(
        description="""
    Processes a sequencing run with staggered primers into a trimmed read set.
    For help type: python assess_staggered_primers.py -h
    """
    )

    parser.add_argument(
        '-source_dir',
        type=str,
        help="The parent directory containing the fastq / fastq.gz files from the sequencing run (required)",
        required=True)

    parser.add_argument(
        '-target_dir',
        type=str,
        help="The target directory for the files of output reads (required)",
        required=True)

    parser.add_argument(
        '-track_tables_fn',
        type=str,
        help="A name for herein generated tables with representation of bases in the sequencing files. Not needed unless you test a new primer sets. (not required)",
        default=None)

    parser.add_argument(
        '-primer_description_file',
        type=str,
        help="A primer description file formatted as specified in the documentation (not required; deafult = primers_16SV34_staggered.txt, for 16S rDNA).",
        default = "primers_16SV34_staggered.txt")

    parser.add_argument(
        '-usearch',
        type=str,
        help="Full path to the USEARCH program. Will be used for merging the paired reads. If not installed, get one at https://www.drive5.com/usearch/download.html - the free 32-bit version is sufficient",
        default=None)

    parser.add_argument(
        '-merge_pairs_by_usearch',
        help="1 / 0 whether we will merge the left and right reads by usearch... default to 0",
        type=int,
        default="0")

    parser.add_argument('-trim_last_bases_r1',
                        type= int,
                        default=0,
                        help = "forward read - how many bases should be trimmed from the end? Some MiSeq runs have the last base of very low quality.  default = 0")
    parser.add_argument('-trim_last_bases_r2',
                        type=int,
                        default=0,
                        help="reverse read - how many bases should be trimmed from the end? Some MiSeq runs have the last base of very low quality.  default = 0")

    parser.add_argument('--trim_only_spacers', action='store_true',
                        help = "What to trim - only the heterogeneity spacers, whereas primers are retained? Default = False")
    parser.add_argument('--trim_whole_primer', action='store_true',
                        help = "What to trim - trim the whole primers? default = True")
    parser.add_argument('--do_not_trim_whole_primer', dest='trim_whole_primer', action='store_false',
                        help = "What to trim - leave the primers on the amplicon? Overrides the above default option.")
    parser.set_defaults(trim_whole_primer=True)

    parser.add_argument('--write_to_files_of_combinations',
                        action='store_true',
                        help="Should fastq files be generated also for all single combination of primers? default = False")

    parser.add_argument('--do_not_overwrite', action='store_true', help="Continue where (if) the previous process stopped? default = False")

    parser.add_argument(
        '-fastq_minmergelen',
        help="minimum length of the merged sequence. Default 350 (because of V3-4 of the 16S rDNA; a bit too short for Blastocystis subtyping)",
        type=int,
        default="350")

    parser.add_argument(
        '-read_length',
        help="The read length in this sequencing. Default 300",
        type=int,
        default="300")

    parser.add_argument(
        '-fastq_maxmergelen',
        help="maximum length of the merged sequence. Default 500 (works well both for the V3-4 of the 16S rDNA, and Blastocystis subtyping)",
        type=int,
        default="500")

    parser.add_argument(
        '-fastq_maxdiffs',
        help="maximum differences between the merged pair memberss . Default 5 (increase for longer overlaps)",
        type=int,
        default="4")

    parser.add_argument(
        '-max_reads_per_sample',
        help="Maximum reads retained per sample. If not given, or 0, all reads will be written to output. Otherwise the given number will be randomly chosen.",
        type=int,
        default="0")


    args = parser.parse_args()

    # now do the job
    st_ds = StaggeredDataset(
        source_dir =                args.source_dir,
        target_dir =                args.target_dir,
        track_tables_fn =           args.track_tables_fn,
        primer_description_file =   args.primer_description_file,
        trim_only_spacers =         args.trim_only_spacers,
        trim_whole_primer =         args.trim_whole_primer,
        trim_last_bases_r1=         args.trim_last_bases_r1,
        trim_last_bases_r2=         args.trim_last_bases_r2,
        write_to_files_of_combinations = args.write_to_files_of_combinations,
        usearch = args.usearch,
        merge_pairs_by_usearch = args.merge_pairs_by_usearch,
        fastq_minmergelen = args.fastq_minmergelen,
        fastq_maxmergelen = args.fastq_maxmergelen,
        fastq_maxdiffs = args.fastq_maxdiffs,
        max_reads_per_sample = args.max_reads_per_sample,
        do_not_overwrite = args.do_not_overwrite,
        read_length = args.read_length
    )

    return


if __name__ == '__main__':
    process_args_and_run()

    # DEBUGGING runs:
    # just note that the script can be also downloaded as a module, and run like this.
    # st_ds = StaggeredDataset(
    #     source_dir="H:/mikro/ngs_16s/41_ibs_po_fmt/run_data_from_sequencer",
    #     target_dir="H:/mikro/ngs_16s/41_ibs_po_fmt/sorted_staggering",
    #     track_tables_fn="H:/mikro/ngs_16s/41_ibs_po_fmt/track_tables.txt",
    #     trim_last_bases_r1=1,
    #     trim_last_bases_r2=1,
    #     fastq_maxdiffs=5,
    #     merge_pairs_by_usearch=1,
    #     max_reads_per_sample = 100000,
    #     do_not_overwrite = True
    #
    # )
    # st_ds = StaggeredDataset(
    #     source_dir=r"H:\mikro\ngs_virome\blasto_202408\run_data_from_sequencer",
    #     target_dir=r"H:\mikro\ngs_virome\blasto_202408\sorted_staggering",
    #     track_tables_fn=r"H:\mikro\ngs_virome\blasto_202408\track_tables.txt",
    #     primer_description_file="primers_blasto_staggered.txt",
    #     trim_last_bases_r1=0,
    #     trim_last_bases_r2=0,
    #     fastq_maxdiffs=5,
    #     merge_pairs_by_usearch=1,
    #     max_reads_per_sample = 10000,
    #     do_not_overwrite = True
    # )
    # st_ds = StaggeredDataset(
    #         source_dir="H:/mikro/ngs_16s/42_lowca/run_data_from_sequencer",
    #         target_dir="H:/mikro/ngs_16s/42_lowca/sorted_staggering",
    #         track_tables_fn="H:/mikro/ngs_16s/42_lowca/track_tables.txt",
    #         trim_last_bases_r1=0,
    #         trim_last_bases_r2=0,
    #         fastq_maxdiffs=5,
    #         merge_pairs_by_usearch=1,
    #         max_reads_per_sample = 500000,
    #         do_not_overwrite = True
    #     )

