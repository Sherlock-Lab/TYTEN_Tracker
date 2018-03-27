#!/usr/env python3

'''
'''

import sys
import os
import subprocess
import re
import glob
import time
import gzip
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from subroutines import *

def parse_ty_insertions(ty_event, file_names):

    '''
    This function takes in a list of SNPs/INDEL events and list of barcodes that link SNPs/INDELs to 
    BAM files of orgin in order to isolate insertion events. The function first pairs each unique barcode
    to a give file of origin name. Then the file of origin is assigned to all SNPs and INDELs. Finally, 
    the function determines which entries are insertion events and only prints these to the output file
    generated in the current working directory.
    '''

    BAR_IDX = 1
    F1_CHR_IDX = 9 
    F1_COOR_IDX = 10
    HAND_IDX = 11
    FIR_COOR = 0
    SEC_COOR = 1
    F2_CHR_IDX = 0
    F2_COOR_IDX = 1
    ORG_IDX = 9
    COOR2_IDX = 2
    BAR_DIDX = 0
    HAND_DIDX = 1
    DIV = 2
    SINGLE = 1

    coor1_2_data = {}
    coor2_2_filename = {}
    barcode_2_filename = {}

    #read in candidate site data
    with open(ty_event, "r") as ty_inserts:
        for line in ty_inserts:
            meta_data = line.rstrip("\n").split("\t")
        
            #save identify unique to each sequence yeast strain
            barcode = meta_data[BAR_IDX]
            chromosome = meta_data[F1_CHR_IDX]
            #coordinates may contain a single value or two values separated by '^' and are their for improper to use outright
            improper_coor = meta_data[F1_COOR_IDX]
            #in case of insertions, this value denotes wether the current sample represents the left of the righthand breakpoint location
            handedness = meta_data[HAND_IDX]
            proper_coor = improper_coor.split("^")
        
        
            if len(proper_coor) == DIV:
                #if the coordinate is a breakpont range save both coordinates
                start_coor = proper_coor[FIR_COOR]
                end_coor = proper_coor[SEC_COOR]
            elif len(proper_coor) == SINGLE:
                #if the coordinate is a SNP populate the second coordinate as None
                start_coor = proper_coor[FIR_COOR]
                end_coor = None   
            else:
                #if the coordinate has an irregular length alert the user
                print("coordinate error 1", file = sys.stderr)
        
            #populate a dictionary with a unique location key and other useful info as values 
            key_tuple = (chromosome, start_coor, end_coor)
            value_tuple = (barcode, handedness)
            coor1_2_data[key_tuple] = value_tuple

    #read in sequencing file of origin data
    with open(file_names, "r") as origin:
        for line in origin:
            meta_data = line.rstrip("\n").split("\t")
        
            chromosome = meta_data[F2_CHR_IDX]
            origin_name = meta_data[ORG_IDX]
            #coordinates may contain a single value or two values separated by '^' or '..' and are their for improper to use outright
            improper_coor = meta_data[F2_COOR_IDX]
            #find coordinates with multiple entries, discard value separators, and split up the coordinate values into a list
            coor_edits = re.sub(r"\^", r"\t", improper_coor)
            coor_edits = re.sub(r"\.\.", r"\t", coor_edits)
            proper_coor = coor_edits.split("\t")
        
            if len(proper_coor) == DIV:
                #if the coordinate is a breakpont range save both coordinates
                start_coor = proper_coor[FIR_COOR]
                end_coor = proper_coor[SEC_COOR]
            elif len(proper_coor) == SINGLE:
                #if the coordinate is a SNP populate the second coordinate as None
                start_coor = proper_coor[FIR_COOR]
                end_coor = None   
            else:
                #if the coordinate has an irregular length alert the user
                print("coordinate error 2", file = sys.stderr)
        
            #populate a dictionary with a unique location key and other useful info as values
            key_tuple = (chromosome, start_coor, end_coor)
            coor2_2_filename[key_tuple] = origin_name

    #should any unique location ID from either dictionary above match, populate a dictionary where the barcode is the key and the filename is the value         
    for key, value in coor1_2_data.items():
        if key in coor2_2_filename and value[BAR_DIDX] not in barcode_2_filename:
            barcode_2_filename[value[BAR_DIDX]] = coor2_2_filename[key]

    with open(linked_data.txt, "w") as output_file:
        for key, value in coor1_2_data.items():
            #consider only insertions and not SNPs from this point forward
            if key[COOR2_IDX] is not None:
                #if an insertion's barcode has been linked to a file name, create an output including that filename
                if value[BAR_DIDX] in barcode_2_filename:
                    output_file.write("\t".join(key) + "\t" + value[HAND_DIDX] + "\t" + barcode_2_filename[value[BAR_DIDX]] + "\n")
                #if an insertion's barcode has not been linked to a file name, create an output and note the non-linkage
                else:
                    output_file.write("\t".join(key) + "\t" + value[HAND_DIDX] + "\t" + "No filematch found" + "\n")
   
            
def pair_ty_insertions(input_data, file_names):

    '''
    This function takes in BED-like, candidate insertion list and a list of BAM/compressed BAM
    file names in order to pair each candidate insertion to a BAM file of origin. The function
    compares filenames in the BED file to filenames in the list. If the names match, the full path
    is appended to the BED-like list. Then, since the BED-like list stores each breakpoint created 
    by an insertion event as a different entry, the function compares adjacent entrys and combines
    the left and right breakpoints of insetions. Then the function only outputs candidate insertions
    that have been successful match to a file path and paired to a file generated in the current 
    working directory.
    '''

    input_data = sys.argv[1]
    file_names = sys.argv[2]

    STEP_UP = 1
    LENGTH = 6
    HAND_IDX = 3
    CHR_IDX = 0
    FILE_IDX = 4
    SCOOR_IDX = 1
    ECOOR_IDX = 2
    MIN_DISTANCE = 0
    MAX_DISTANCE = 50
    NAME_IDX = 5
    SINGLE = 1

    line_num = 0
    line_use = 0

    line_2_data = {}
    origin_files = []

    #read in candidate insertion list
    with open(input_data, "r") as insertions:
        for line in insertions:
            #populate a dictionary where line number is the key and insertion data is the value
            data = line.rstrip("\n").split("\t")
            line_2_data[line_num] = data
            line_num = line_num + STEP_UP

    #read in all paths and files associated with sequence data
    with open(file_names, "r") as origin:
        for line in origin:
            data = line.rstrip("\n")
            origin_files.append(data)

    for entry in line_2_data:
        unedited_line = line_2_data[entry]
        unedited_file = unedited_line[FILE_IDX]
    
        #take each filename associated with an insertion candidate and modify it to be a bam or compressed bam filename
        if unedited_file != "No filematch found":
            bam_edit = re.sub(r" \(Variants, KNOWN\)", r".bam", unedited_file)
            gz_edit = re.sub(r" \(Variants, KNOWN\)", r".bam.gz", unedited_file)
        
            #if the bam or compressed bam filename exist as a real file, save it as the new filename for the insertion
            if bam_edit in origin_files:
                line_2_data[entry].append(bam_edit)
            elif gz_edit in origin_files:
                line_2_data[entry].append(gz_edit)
                
    with open("paired_data.txt", "w") as output_file:
        while line_use < len(line_2_data) - STEP_UP:
            #compare each insertion candidate to the candidate appear on the line directly after it within the file
            current_line = line_2_data[line_use]
            comparison_line = line_2_data[(line_use + STEP_UP)]
    
            #should the two lines represent left and right breakpoints of the same insertion, merge them for output
            if len(current_line) == LENGTH and len(comparison_line) == LENGTH:
                if current_line[HAND_IDX] == "Left breakpoint" and comparison_line[HAND_IDX] == "Right breakpoint":
            
                    if current_line[CHR_IDX] == comparison_line[CHR_IDX] and current_line[FILE_IDX] == comparison_line[FILE_IDX]:
                        if MIN_DISTANCE <= (int(comparison_line[SCOOR_IDX]) - int(current_line[ECOOR_IDX])) <= MAX_DISTANCE: 
                    
                            output_file.write(current_line[CHR_IDX] + "\t" + current_line[SCOOR_IDX] + "\t" + comparison_line[ECOOR_IDX] + "\t" + current_line[NAME_IDX] + "\n")  
    
            line_use = line_use + SINGLE


def create_file_list(home_folder):

    '''
    This function takes an input directory containing BAM and compressed BAM files and a output 
    file to pipe results to. The Function searches through the provided directroy for all BAM/
    compressed BAM files and appends each of their names and paths to the output file generated
    in the current working directory. 
    '''
    
    with open("file_list.txt", "w") as output_file:
        #create a list of all bam and compressed bam files
        for root, dirnames, filenames in os.walk(home_folder):
            for filename in filenames:
                if filename.endswith(('.bam', '.gz')):
                    output_file.write(filename)


def collect_sequences(input_file, bam_path):

    '''
    This function completes several tasks, taking a BED-like input list and a home directory containg BAM files.
    The program finds all BAM files within the home directory with names matching feature entries in the provided 
    BED-like file. The BAM files are then indexed and the names of all Bam file reads withing BED file provided 
    coordinates are printed to a temporary list. This list is then used to grab any reads in the BAM whose name
    matches any of those in the list; this proccess will grab reads whose names were originally found in the 
    designated region as well as all their paired-read partners. These reads are stored as temporary SAM files
    which is then restructured into two, paired-end FASTQ data files.
    '''

    file_2_scaffold = {}

    FILE_IDX = 3
    CHR_IDX = 0
    START_IDX = 1
    END_IDX = 2
    SIZE_UP = 530
    SIZE_EX = 30
    RCHR_IDX = 0
    RLB_IDX = 1
    RUB_IDX = 2
    INS_NUM = 0
    INS_END = 2
    TRU_INS = 1
    R1_IDX = 0
    R2_IDX = 1
    INCREASE_STEP = 1
    SEQ_IDX = 0
    QUAL_IDX = 1
    COMPRESS = -3
    DECOMP = "decompressed"
    EMPTY = 0
    DIV = 2
    SINGLE = 1
    STARTER = "@"

    current_output = 1

    #read in candidate insertion site data as a list
    with open(input_file, "r") as regional_data:
        regional_data_list = regional_data.readlines()
        for line in regional_data_list:
            line_data = line.rstrip("\n").split("\t")
            region_file = line_data[FILE_IDX]
            #populate a dictionary with each unique filename as a key with an empty list value
            if region_file not in file_2_scaffold.keys():
                file_2_scaffold[region_file] = []

        for line in regional_data_list:
            #set each item in a candidate insertion site list as a variable
            line_data = line.rstrip("\n").split("\t")
            region_chr = line_data[CHR_IDX]
            region_start = line_data[START_IDX]
            region_end = line_data[END_IDX]
            region_file = line_data[FILE_IDX]
            #create a tuple that holds coordinates for a 500 basepair region up/downstream of the insertion site
            start_tuple = (region_chr, (int(region_start) - SIZE_UP), (int(region_start) + SIZE_EX))
            end_tuple = (region_chr, (int(region_end) - SIZE_EX), (int(region_end) + SIZE_UP))
            #add these tuples to the value-lists of the dictionary keyed by the filename
            file_2_scaffold[region_file].append(start_tuple)
            file_2_scaffold[region_file].append(end_tuple)


    for file_id in file_2_scaffold.keys():
        full_bam_path = None
        gunzip_state = None
        #counts the number of insertion from each file
        insertion_number = INCREASE_STEP
        #collect the full path to each file whose name appears as a key in the dictionary above
        for root, dirnames, filenames in os.walk(bam_path):
            if full_bam_path is None:
            
                for filename in filenames:
                    if filename == file_id:
                        full_bam_path = os.path.abspath(os.path.join(root, filename))
                        if full_bam_path.endswith(".gz"):
                            #if the file is compressed, uncompress it and remove '.gz' from its path name
                            subprocess.call(["gunzip", full_bam_path])
                            full_bam_path = full_bam_path[:COMPRESS]
                            #note that the file was originally compressed 
                            gunzip_state = DECOMP
                        break
                    
            else:
                break
            
        #create an index.bai file for fast access to the alignment .bam file
        subprocess.call(["samtools", "index", full_bam_path])

        for region_tuple in file_2_scaffold[file_id]:
            #for each candidate insertion, grab names of reads aligned both upstream and downstream of the insertion site
            read_names = grab_useful_readnames(full_bam_path, region_tuple)
        
            with open("temp_name_file.txt", "w") as name_file:
                name_file.write("\n".join(read_names) + "\n")
        
            #Using the names of reads near the insertion, collect these reads and their pairs from the bam file 
            subprocess.call(["java", "-jar", "./picard/build/libs/picard.jar", "FilterSamReads", ("I=" + full_bam_path), "O=temp_file.sam", "READ_LIST_FILE=temp_name_file.txt", "FILTER=includeReadList"])
            #remove the no longer needed list of names
            os.remove("temp_name_file.txt")
        
            #grab and order the read sequence/quality for read pair
            r1_r2_reads = create_fastq_data("temp_file.sam")
            r1_data = r1_r2_reads[R1_IDX]
            r2_data = r1_r2_reads[R2_IDX]
        
            #count pairs of upstream/downstream regions of insertion sites
            INS_NUM = INS_NUM + INCREASE_STEP
        
            #create unique fastq file name for each insertion site
            r1_output_name = "Sample" + str(current_output) + "_Insertion" + str(insertion_number) + "_R1_" + str(region_tuple[RCHR_IDX]) + "_" + str(region_tuple[RLB_IDX]) + "_" + str(region_tuple[RUB_IDX]) + ".fastq"
            r2_output_name = "Sample" + str(current_output) + "_Insertion" + str(insertion_number) + "_R2_" + str(region_tuple[RCHR_IDX]) + "_" + str(region_tuple[RLB_IDX]) + "_" + str(region_tuple[RUB_IDX]) + ".fastq"
        
            #pairs upstream and downstream files so each can be counted as one insertion for naming purposes
            if INS_NUM == INS_END:
                INS_NUM = EMPTY
                insertion_number = insertion_number + INCREASE_STEP
        
            #remove no longer needed sam file   
            os.remove("temp_file.sam")
        
            #write fastq file with correct format: read name, sequence, +, quality
            with open(r1_output_name, "w") as new_file1:
                for item in r1_data:
                    #since the list contains the read name as an item followed by a sequence/quality tuple treat even indexes as read names
                    if (int(r1_data.index(item)) % DIV) == EMPTY:
                        new_file1.write(STARTER + str(item) + "\n")
                    #treat odd indexes as sequence/quality tuple
                    elif (int(r1_data.index(item)) % DIV) == SINGLE:
                        new_file1.write(str(item[SEQ_IDX]) + "\n" + "+\n" + str(item[QUAL_IDX]) + "\n")
                    
            #write fastq file with correct format: read name, sequence, +, quality
            with open(r2_output_name, "w") as new_file2:
                for item in r2_data:
                    #since the list contains the read name as an item followed by a sequence/quality tuple treat even indexes as read names
                    if (int(r2_data.index(item)) % DIV) == EMPTY:
                        new_file2.write(STARTER + str(item) + "\n")
                    #treat odd indexes as sequence/quality tuple
                    elif (int(r2_data.index(item)) % DIV) == SINGLE:
                        new_file2.write(str(item[SEQ_IDX]) + "\n" + "+\n" + str(item[QUAL_IDX]) + "\n")          
    
        current_output = current_output + INCREASE_STEP
        #remove unneeded index file for bam no longer in use.
        os.remove((full_bam_path + ".bai"))
        #if a not was made that the current bam file was originally compressed, recompress it.
        if gunzip_state == DECOMP:
            subprocess.call(["gzip", full_bam_path])
   
            
def assemble_insertions(home_folder, assembler, kmer_len, read_length, covergae):

    '''
    This function takes a given directroy containing paired-end fastq files along with a 
    assembler name and its arguements in odrder to create de novo assebled genome FASTA files.  
    '''

    fastq_list = []

    #search recursively through the home directory for all fastq files
    for root, dirnames, filenames in os.walk(home_folder):
        for filename in filenames:
            if filename.endswith(".fastq"):
                #store the name of all fastq files found
                fastq_list.append(filename)
            
    for name in fastq_list:
        #For each R1 file contain the first read of a pair, find its associated R2 file contain the second read of the pair
        if "_R1_" in name:
            pair_one = name
            to_edit = name
            pair_two = re.sub(r"_R1_", r"_R2_", to_edit)
            #create a name for the combine read files to be stored under
            paired_reads = re.sub(r"_R1_", r"_paired_", to_edit)
            edited_pair = re.sub(r"\.", r"_", paired_reads)
        
            #the function currently only supports de novo assembly using the velvet program package
            #later releases will support other assembler packages
            if assembler == "velvet":
                subprocess.call(["perl", "../velvet/shuffleSequences.pl", ("../fastq_data/" + pair_one), ("../fastq_data/" + pair_two), paired_reads])
                output_dir = ("./" + edited_pair + "/")
                subprocess.call(["../velvet/velveth", output_dir, kmer_len, "-fastq", "-shortPaired", paired_reads])
                subprocess.call(["../velvet/velvetg", output_dir , "-ins_length", read_length, "-exp_cov", coverage])
                os.remove(paired_reads)
            
            elif assembler == "xxx":
                pass
            elif assembler == "xxx":
                pass
            elif assembler == "xxx":
                pass

def blast_insertions(scaffold_dir, yeast_strain, genome_feature_file, level):
    
    '''
    This function takes an input directory containg assembled 'synthetic genomes' and assess  whether
    each contig per FASTA file blasts to a location of a provided reference genome ID. These hits are 
    filtered to only include the top hit(s) for each contig and determine if the hits correspond to a 
    provided set of genome features. The list of features associated with hit for the entire FASTA file 
    is then narrowed down to the most hit feature and its ID is assigned to the output file. The output 
    file is then further condenced by comparing linked FASTA file names and only retaining features that
    were the chosen call for both files. This process is optomized to specifically ID transposons in the 
    Saccharoymces cerevisiae reference genome, though could easily be repurposed using different constants.
    '''

    #constants for file.gff parsing
    SEQ_END = 7
    HED_IDX = 0
    FET_IDX = 2
    TY_IDX = 8
    STRS_IDX = 7
    CHR_IDX = 0
    COORS_IDX = 3
    COORE_IDX = 4

    gff_data_list = []

    #constants for file retrieval/blast
    PRESENT = 2
    NAME_TRIM = -11
    HIT_IDX = 0
    NAME_IDX = 0
    NODE = 0
    WAIT = 1
    MAX = 10
    NEWER = 0
    MATCH = 0
    PLURAL = 2
    INCREASE_STEP = 1 

    #constants for transposon hit retrieval
    REF_CHR = 0
    STY = 1
    ETY = 2
    BLEED = 100
    PID = 3
    NPID = 1

    #constants for insertion end-pairing
    EMPTY = 0
    INS_IDX = 0
    HEADER = 'Insertion'
    TRUN_NAME = 20
    BTY = 1
    BSC = 2
    BSI = 3
    BTY_EXP = 0
    BSC_EXP = 1
    BSI_EXP = 2
    AVERAGE = 2
    DATUM = 1
    PAIRED = 4
    DAT_N = 0
    DAT_T = 1
    DAT_S = 2
    DAT_I = 3

    name_2_data = {}
    
    #use the provided level to assign IDs at Family, Sub-family, or copy of origin level
    if level == "lax"
        STRE_IDX = 10
    elif level == "mod"
        STRE_IDX = 12
    
    with open(genome_feature_file, "r") as feature_file:
        for line in feature_file:
            #stop iterating over the gff file if all the lines containing genome features have already been viewed
            if line[:SEQ_END] == "##FASTA":
                break
            else:
                #replace data dividing whitespace with tabs
                line = re.sub("(\s)+", "\t", line)
                #ignore header lines that begin with '#'
                if line[HED_IDX] != "#":
                    genome_feature = line.rstrip('\n').split("\t")
                    feature_id = genome_feature[FET_IDX]
                    element_id = genome_feature[TY_IDX][STRS_IDX:STRE_IDX]
                
                    #grab data for a feature if it is a retrotransposon
                    if feature_id == "LTR_retrotransposon":
                        chromosome = genome_feature[CHR_IDX]
                        start_coor = genome_feature[COORS_IDX]
                        end_coor = genome_feature[COORE_IDX]
                        #append all transposon data as tuples to a list for indexing later
                        gff_data_list.append((chromosome, start_coor, end_coor, element_id))

    #open file to store all inerstion data
    with open("unpaired_hits.txt", "w") as final_data:
        final_data.write("Insertion\tTy Element Type\tTotal Affirming Hits\tAverage Sequence Similiarity\n")

        #locate all synthetic scaffold fasta files
        for root, dirnames, filenames in os.walk(scaffold_dir):
            for filename in filenames:
                if filename.endswith(".fa"):
            
                    #define constants for file parsing
                    contigs = os.path.join(root, filename)
                    output_name = contigs[PRESENT:NAME_TRIM]
                    print(filename, file = sys.stderr)
                    print(output_name, file = sys.stderr)
                    blast_query = list(SeqIO.parse(contigs, "fasta"))
                    single_file_list = []
                    call_types_list = []
                    hsp_list = []
                    current_contig = 0         
                
                    while current_contig < len(blast_query):
                        STEP = NODE
                        percentage = None
                        #Name to save blast_run under TO BE REMOVED IN FINAL CODE
                        out_file = os.path.join(scaffold_dir, "blast_" + str(current_contig) + "_" + output_name + ".xml")
                    
                        #If the fasta file has not been used to create an file.XML, do so now
                        if not os.path.exists(out_file):
                            blast_results_handle = NCBIWWW.qblast("blastn", "nr", blast_query[current_contig].format("fasta"), entrez_query = yeast_strain)
                            #pause for a constant number of seconds between each blast submission
                            time.sleep(WAIT)
                            #Write blast output to an file.XML
                            with open(out_file, "w") as blast_data_file:
                                blast_data_file.write(blast_results_handle.read())
                            
                        #open the previously written file.XML        
                        with open(out_file, "r")  as blast_data_file:      
                            blast_record = NCBIXML.read(blast_data_file)
                    
                            while STEP < MAX:
                                #collect data on HSP if not HSP has been evaluated or if the current HSP is as good as the previous
                                if (percentage is None and STEP == NEWER and len(blast_record.descriptions) > EMPTY) or (total_f >= total and percentage == percentage_f):
                                    name = blast_record.descriptions[STEP].title
                                    #isolate the name of the location
                                    hit_locale = re.findall("chromosome\ [XVI][XVI]?[XVI]?", name)
                                
                                    if len(hit_locale) >= MATCH + INCREASE_STEP:
                                        hit_string = hit_locale[MATCH]
                                        #name the HSP location in the same format as Transposon list data
                                        hit_chromosome = re.sub(r"chromosome\ ", r"chr", hit_string)
                                        start = blast_record.alignments[NODE].hsps[STEP].sbjct_start
                                        stop = blast_record.alignments[NODE].hsps[STEP].sbjct_end
                                        match = blast_record.alignments[NODE].hsps[STEP].identities
                                        total = blast_record.alignments[NODE].hsps[STEP].align_length
                                        percentage = (int(match)/int(total))
                                    
                                        if len(blast_record.alignments[NODE].hsps) > STEP + INCREASE_STEP:
                                            #collect minimal data on the next-in-line HSP for comparison in next while loop iteration
                                            match_f = blast_record.alignments[NODE].hsps[(STEP + INCREASE_STEP)].identities
                                            total_f = blast_record.alignments[NODE].hsps[(STEP + INCREASE_STEP)].align_length
                                            percentage_f = (int(match_f)/int(total_f))
                                    
                                        else:
                                            match_f = NODE
                                            total_f = NODE
                                            percentage = NODE
                                        
                                        #fill list with useful info on selected HSPs
                                        hsp_list.append((hit_chromosome, start, stop, percentage))
                                   
                                    STEP = STEP + INCREASE_STEP
                               
                                #if the current HSP does not meet or exceed quality threshold, exit the while loop    
                                else:
                                    STEP = MAX
                    
                        current_contig = current_contig + INCREASE_STEP
                        
                    if level == "lax" or level == "mod":
                        #compare transposons to all hits kept from blast                               
                        for item in gff_data_list:
                            for hit in hsp_list:
                                if item[REF_CHR] == hit[REF_CHR]:
                        
                                    #if the hit and the transposon intersect save the hit and associated data            
                                    if int(item[STY]) <= (int(hit[STY]) + BLEED) <= int(item[ETY]) or int(item[STY]) <= (int(hit[STY]) - BLEED) <= int(item[ETY]):
                                        single_file_list.append((item[PID], hit[PID]))
                                        call_types_list.append(item[PID])
                            
                                    #if the hit and the transposon intersect save the hit and associated data                   
                                    elif int(item[STY]) <= (int(hit[ETY]) + BLEED) <= int(item[ETY]) or int(item[STY]) <= (int(hit[ETY]) - BLEED) <= int(item[ETY]):
                                        single_file_list.append((item[PID], hit[PID]))                                
                                        call_types_list.append(item[PID])
                                                                                   
                    elif level == "strict":
                        #compare transposons to all hits kept from blast                               
                        for item in gff_data_list:
                            for hit in hsp_list:
                                if item[REF_CHR] == hit[REF_CHR]:
                        
                                    #if the hit and the transposon intersect save the hit and associated data            
                                    if int(item[STY]) <= (int(hit[STY]) + BLEED) <= int(item[ETY]) or int(item[STY]) <= (int(hit[STY]) - BLEED) <= int(item[ETY]):
                                        single_file_list.append((item, hit[PID]))
                                        call_types_list.append(item)
                            
                                    #if the hit and the transposon intersect save the hit and associated data                   
                                    elif int(item[STY]) <= (int(hit[ETY]) + BLEED) <= int(item[ETY]) or int(item[STY]) <= (int(hit[ETY]) - BLEED) <= int(item[ETY]):
                                        single_file_list.append((item, hit[PID]))                                
                                        call_types_list.append(item)                                 
                    
                    #define constants for dominant hit determination    
                    index_of_call = []
                    final_calls = []
                    best_name = None
                
                    if len(call_types_list) != EMPTY:    
                        for item in call_types_list:
                            #determine the number of times each transposon type was defined by a hit and which indexes in the list correspond to them             
                            current_call = list_duplicates_of(call_types_list,item)
                            call_name = item
                            #if a dominat hit is not defined or if the dominat hit is not already the same as the current hit
                            if call_name != best_name:
                    
                                #if the current hit is better than the dominat hit, replace it as the dominant hit
                                if len(current_call) > len(index_of_call):
                                    index_of_call = current_call
                                    best_name = call_name
                        
                                #if the current hit is equal to the dominant hit, do more in-depth comparison
                                elif len(current_call) == len(index_of_call):
                                    challenging_calls = []
                                    encumbent_calls = []
                            
                                    #determine the average sequence similarity between the current transposon and all hits that call it
                                    for call_idx in current_call:
                                        challenging_calls.append(single_file_list[call_idx][NPID])
                                    challenge_q = (sum(challenging_calls) / len(challenging_calls))
                            
                                    #determine the average sequence similarity between the current transposon and all hits that call it
                                    for call_idx in index_of_call:
                                        encumbent_calls.append(single_file_list[call_idx][NPID])
                                    encumbent_q = (sum(encumbent_calls) / len(encumbent_calls))
                            
                                    #if the current hit's transposon call has better average sequence similarity than the dominant hit, replace it
                                    if challenge_q > encumbent_q:
                                        index_of_call = current_call
                                        best_name = call_name
                        
                                #if the dominant hit is better than the current hit, it remains the dominant hit
                                elif len(current_call) < len(index_of_call):
                                    pass
                
                        #Once the best dominant hit has been determined, average its sequence similarity to the transposon and it total # of appearances in the list                                
                        for call_idx in index_of_call:                  
                            final_calls.append(single_file_list[call_idx][NPID])
                        final_q = str(sum(final_calls) / len(final_calls))
                        final_t = str(len(index_of_call) / len(call_types_list))
                        
                        if level == "lax" or level == "mod":
                            #print the dominant hit for each file to the output
                            final_data.write((output_name + "\t" + best_name + "\t" + final_t + "\t" + final_q + "\n"))
                            
                        elif level == "strict":                   
                            best_name = "_".join(best_name)
                            #print the dominant hit for each file to the output
                            final_data.write((output_name + "\t" + best_name + "\t" + final_t + "\t" + final_q + "\n"))                            
                            
    #open insertion data for pairing ends of the insertions                                               
    with open("unpaired_hits.txt", "r") as data_to_pair:
        for line in data_to_pair:
            #replace data dividing whitespace with tabs
            edited_line = re.sub("(\s)+", "\t", line)
            line_data = edited_line.rstrip("\n").split("\t")
        
            #ignore header lines
            if line_data[INS_IDX] != HEADER:
                starting_key = line_data[INS_IDX][:TRUN_NAME]
            
                #If this is first entry for a given insertion add it to the dictionary
                if starting_key not in name_2_data: 
                    name_2_data[starting_key] = (line_data[BTY], line_data[BSC], line_data[BSI])
                
                else:
                    #if the data on the insertion don't agree on the transposon's identity, delete the entry
                    if name_2_data[starting_key][BTY_EXP] != line_data[BTY]:
                        del name_2_data[starting_key]
                
                    #if the data agree on the transposons Identity, average their data together    
                    elif name_2_data[starting_key][BTY_EXP] == line_data[BTY]:
                        scaffold_s1 = float(name_2_data[starting_key][BSC_EXP])
                        scaffold_s2 = float(line_data[BSC])
                        similarity_s1 = float(name_2_data[starting_key][BSI_EXP])
                        similarity_s2 = float(line_data[BSI])
                        finished_scaffold = str((scaffold_s1 + scaffold_s2)/AVERAGE)
                        finished_similarity = str((similarity_s1 + similarity_s2)/AVERAGE)
                        del name_2_data[starting_key]
                    
                        #trim the final output names based on a character overhang of 2
                        if starting_key.endswith("_p"):
                            finished_key = re.sub(r"_p", r"", starting_key)
                            name_2_data[starting_key] = (finished_key, line_data[BTY], finished_scaffold, finished_similarity)
                    
                        #trim the final output names based on a character overhang of 1    
                        elif starting_key.endswith("_"):
                            finished_key = re.sub(r"_", r"", starting_key)
                            name_2_data[starting_key] = (finished_key, line_data[BTY], finished_scaffold, finished_similarity)
                    
                        #trim the final output names based on a character overhang of 0    
                        else:
                            name_2_data[starting_key] = (starting_key, line_data[BTY], finished_scaffold, finished_similarity)

    #Create the final output file of the pipeline                        
    with open("TYTEN_calls.txt", "w") as final_ty:
        #Add a header to the output file
        final_ty.write("Insertion\tTy Element Type\tTotal Affirming Hits\tAverage Sequence Similiarity\n")
        for value in name_2_data.items():
            #if the dict entry was correctly paired and therefor the value a 4 item tuple instead of 3 print it to the output file
            if len(value[DATUM]) == PAIRED:
                final_ty.write(value[DATUM][DAT_N] + "\t" + value[DATUM][DAT_T] + "\t" + value[DATUM][DAT_S] + "\t" + value[DATUM][DAT_I] + "\n")
                 
 
def build_ty_elements(genome_feature_file, reference_genome):

    '''
    This function takes an input FASTA and GFF file and uses feature coordinates pulled from the
    GFF file to grab squences for those regions. The sequences are then stored in an output FASTA
    file where each sequence grabbed is represented as a different contig. This is currently optomized
    for grabbing retrotransposon features from a genome containing low number of retrotransposon 
    copies but the function can be easily repurposed by changing constants.
    '''

    FIRST = 0
    FOOTER = 7
    HED_IDX = 0 
    FET_IDX = 2
    TY_IDX = 8
    STRS_IDX = 7
    STRE_IDX = 12
    CHR_IDX = 0
    COORS_IDX = 3
    COORE_IDX = 4
    TE = "LTR_retrotransposon"
    
    i = 1
    
    unique_id_list = []

    #grab all the S288C refrence chromosome fasta files
    for file in sorted(list(glob.glob("./chr[0-9][0-9].fsa"))):
        with open(file, 'r') as fd:
            contents = fd.read().split("\n")
            
        #replace each file name and fasta header with roman numeral chromosome name
        contents[0] = ">chr" + int_to_roman(i)
        i = i + 1
        with open(file, 'w') as fd:
            fd.write("\n".join(contents))

    #grab all the S288C refrence features from .gff file
    with open(genome_feature_file, "r") as feature_file, open("temp_ty_file.bed", "w") as temp_bed:
        for line in feature_file:
            #stop iterating over the gff file if all the lines containing genome features have already been viewed
            if line[:FOOTER] == "##FASTA":
                break
            else:
                #replace data dividing whitespace with tabs
                line = re.sub("(\s)+", "\t", line)
                #ignore header lines that begin with '#'
                if line[HED_IDX] != "#":
                    genome_feature = line.rstrip("\n").split("\t")
                    feature_id = genome_feature[FET_IDX]
                    element_id = genome_feature[TY_IDX][STRS_IDX:STRE_IDX]
                
                    #grab data for a feature if it is a unique type of retrotransposon
                    if feature_id == TE:
                        if element_id not in unique_id_list:
                            unique_id_list.append(element_id)
                            chromosome = genome_feature[CHR_IDX]
                            start_coor = genome_feature[COORS_IDX]
                            end_coor = genome_feature[COORE_IDX]
                            #output the data in the form of a bed file
                            bed_entry = "\t".join([str(chromosome), str(start_coor), str(end_coor), str(element_id)])
                            temp_bed.write(bed_entry + "\n")

    #grab sequences from the reference genome that correspond to the provided bed file items
    synthetic_genome = open("ty_element_scaffolds.fa", "w")
    subprocess.call(["bedtools", "getfasta", "-name", "-fi", reference_genome, "-bed", "temp_ty_file.bed"], stdout = synthetic_genome)
    #remove intermediate files 
    os.remove("temp_ty_file.bed")
    os.remove("S288C.fsa.fai")              
        
                     
      