#!/usr/env python3

'''
This script should run with the following command line input:
python3 TYTEN_Tracker.py filename_data.txt insertion_data.txt file_list.txt ./bam_dir/ ./fastq_dir/ velvet 21 500 15 ./fasta_dir/ txid559292 saccharomyces_cerevisiae.gff S288C.fsa lax Blast
'''

from TYTEN_functions import *
import sys
import argparse

#improvements to make
#add argument parser
'''
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--runtype", choices = ["Dual", "Blast", "Align"], default = "Blast", help = "chooses the method by which to ID TEs")
parser.add_argument("-f", "--filetype", choices = ["BED", "LIST"], default = "BED", help = "Denotes whether a BED file or LIST file will be used as the initial input" )
parser.add_argument("inputfile", nargs = "*", help = "the name/path to the single input bed file or the two input list files")
'''
#add bedfile start compatibility
#add local blast compatibility 

file_barcode_linkers = sys.argv[1] 
candidate_insertion_events = sys.argv[2] 
file_full_paths = sys.argv[3] 
home_of_BAM_files = sys.argv[4]
home_of_FASTQ_files = sys.argv[5] 
chosen_assembler = sys.argv[6] 
kmer_size = sys.argv[7]  
pairing_distance = sys.argv[8] 
coverage_depth = sys.argv[9] 
home_of_FASTA_files = sys.argv[10] 
yeast_blast_id = sys.argv[11] 
yeast_reference_gff_file = sys.argv[12] 
yeast_reference_genome = sys.argv[13]  
assignment_stringency = sys.argv[14]  
runtype = sys.argv[15] 


if runtype == "Dual":
    create_file_list(home_of_BAM_files)
    build_ty_elements(yeast_reference_gff_file, yeast_reference_genome)
    parse_ty_insertions(file_barcode_linkers, candidate_insertion_events)
    pair_ty_insertions("linked_data.txt", "file_list.txt")
    collect_sequences("paired_data.txt", home_of_BAM_files)
    assemble_insertions(home_of_FASTQ_files, chosen_assembler, kmer_size, pairing_distance, coverage_depth)
    blast_insertions(home_of_FASTA_files, yeast_blast_id, yeast_reference_gff_file, assignment_stringency)
    print("multi-aligning tools coming soon", file = sys.stderr)
    multi_align_insertions(home_of_FASTA_files, "ty_element_scaffolds.fa")
    validate_insertions("TYTEN_calls.txt", "TYTEN_assigns.txt")
    
if runtype == "Blast":
    create_file_list(home_of_BAM_files)
    parse_ty_insertions(file_barcode_linkers, candidate_insertion_events)
    pair_ty_insertions("linked_data.txt", "file_list.txt")
    collect_sequences("paired_data.txt", home_of_BAM_files)
    assemble_insertions(home_of_FASTQ_files, chosen_assembler, kmer_size, pairing_distance, coverage_depth)
    blast_insertions(home_of_FASTA_files, yeast_blast_id, yeast_reference_gff_file,assignment_stringency)
    
if runtype == "Align":
    create_file_list(home_of_BAM_files)
    build_ty_elements(yeast_reference_gff_file, yeast_reference_genome)
    parse_ty_insertions(file_barcode_linkers, candidate_insertion_events)
    pair_ty_insertions("linked_data.txt", "file_list.txt")
    ccollect_sequences("paired_data.txt", home_of_BAM_files)
    assemble_insertions(home_of_FASTQ_files, chosen_assembler, kmer_size, pairing_distance, coverage_depth)
    multi_align_insertions(home_of_FASTA_files, "ty_element_scaffolds.fa")
    
