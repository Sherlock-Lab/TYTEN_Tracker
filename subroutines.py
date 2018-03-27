#!/usr/env python3

import sys
import subprocess
import os


def grab_useful_readnames(input_bam, candidate_site):
    
    '''
    This function is used to extract the names of all reads that map to a given set of chromosomal coordinates. 
    The function takes the path to a bam file (one which has a indexing .bam.bai file of the same name in the 
    same directory) as well as a tuple containing "(chromosome identity, lower coordinate bound, upper coordinate
    bound) and outputs a list of all read names of interest.
    '''
    
    all_read_names = []
    
    CHR_IDX = 0
    FIR_IDX = 1
    SEC_IDX = 2
    NAM_IDX = 0
    STA_IDX = 6
    READ_P = "="
    
    chromosome = candidate_site[CHR_IDX]
    start_site = candidate_site[FIR_IDX]
    end_site = candidate_site[SEC_IDX]
    
    read_candidates = open("temp_output.sam", "w")
    subprocess.call(["samtools", "view", input_bam, (chromosome + ":" + str(start_site) + "-" + str(end_site))], stdout = read_candidates)
    
    
    with open("temp_output.sam", "r") as candidate_reads:
        for line in candidate_reads:
            read_data = line.rstrip("\n").split("\t")
            read_name = read_data[NAM_IDX]
            read_state = read_data[STA_IDX]
            
            #makes sure every read name is grabbed only once (avoids both reads if paired end data)
            if read_name not in all_read_names:
                if read_state == READ_P:       
                    all_read_names.append(read_name)
                
    #os.remove("temp_output.sam")
    return all_read_names
    


def create_fastq_data(sam_file):

    '''
    This function is used to generate pair-end fastq data sets from a sam file it is passed. The
    function evalutes each read grabbed for errors in the sequence/quality fields since both will 
    appear in the resulting fastq files. 'read error type one' signals the first read in a pair contains
    an error and 'read error type two' means the second read in a pair contains an error. The function 
    also requires that reads be properly paired. 'pair error type one' means the that a read with more
    than one mate was included in the sam file, while 'pair error type two' means the code failed to 
    correctly identify paired reads. finally, the function also make sure the reads in each file appear 
    in the same order. If any reads in list one do not directly align to the order of their pair in list 
    two, a 'file match error' will occur.
    '''
    #list are populated so that reads can be indexed in paired orders in both files regardless of their order in the sam file
    fastq_list1 = ([None] * 25000)
    fastq_list2 = ([None] * 25000)
    
    HED_IDX = 0
    QNAME_IDX = 0
    SEQ_IDX = 9
    QUAL_IDX = 10
    INCREASE_STEP = 1
    WALK = 2
       
    current_step = 0
    current_hit = 0
    current_item = 0
    total_reads = 0
    
    
    with open(sam_file) as isolated_reads:
        for line in isolated_reads:
            if line[HED_IDX] != "@":
                read_data = line.rstrip("\n").split("\t")
                read_name = str(read_data[QNAME_IDX])
                read_seq = str(read_data[SEQ_IDX])
                read_qual = str(read_data[QUAL_IDX])
        
            
                if read_name not in fastq_list1 and read_name not in fastq_list2:        
                    if len(read_seq) == len(read_qual):
                        #if this is the first occurence of a read name and it is error free, replace two 'None's in list one with data
                        fastq_list1[current_hit] = read_name
                        fastq_list1[(current_hit + INCREASE_STEP)] = (read_seq, read_qual)
                        current_hit = current_hit + WALK
                    else:
                        #triggered if there is an error between the sequence and quality fields
                        print("read error type 1", file = sys.stderr)
                        
                elif read_name in fastq_list1 and read_name not in fastq_list2:
                    if len(read_seq) == len(read_qual):
                        #if this is the second occurence of a read name and it is error free, replace two 'None's in list two with data
                        mate_position = fastq_list1.index(read_name)
                        fastq_list2[mate_position] = read_name
                        fastq_list2[(mate_position + INCREASE_STEP)] = (read_seq, read_qual)
                    else:
                        #triggered if there is an error between the sequence and quality fields
                        print("read error type 2", file = sys.stderr)
                
                elif read_name in fastq_list1 and read_name in fastq_list2:
                    #triggers on the third occurence of a read name
                    print("pair error type 1", file = sys.stderr)
        
                elif read_name not in fastq_list1 and read_name in fastq_list2:
                    #should never trigger but completes all possible cases
                    print("pair error type 2", file = sys.stderr)      
        
                current_step = current_step + INCREASE_STEP
    
    while current_item < len(fastq_list1):
        #check that all read pairs have the same indexes in list one and two
        if fastq_list1[current_item] == fastq_list2[current_item]:
            if fastq_list1[current_item] is not None and fastq_list2[current_item] is not None:
                total_reads = total_reads + WALK
        else: 
            print("file match error", file = sys.stderr)
        
        current_item = current_item + WALK
    #remove all excess None entries from finished read pair lists    
    fastq_output1 = fastq_list1[:total_reads]
    fastq_output2 = fastq_list2[:total_reads]    
    
    #return both lists in a tuple
    return (fastq_output1, fastq_output2)


def list_duplicates_of(seq,item):

    '''
    This function was developed by and credited to: 
    PaulMcG via StackOverflow: https://stackoverflow.com/questions/5419204/index-of-duplicates-items-in-a-python-list
    The function takes a list of values and for each unique value returns all the indexes at which it
    exists. For example for the value 'x' in the following list: ['x', 'y', 'x', 'z', 'x'], the function
    would return [0, 1, 3].
    '''
    start_at = -1
    locs = []
    while True:
        try:
            loc = seq.index(item,start_at+1)
        except ValueError:
            break
        else:
            locs.append(loc)
            start_at = loc
    return locs


def int_to_roman(input):

   '''
   Convert an integer to Roman numerals.
   '''
     
   ints = (1000, 900,  500, 400, 100,  90, 50,  40, 10,  9,   5,  4,   1)
   nums = ('M',  'CM', 'D', 'CD','C', 'XC','L','XL','X','IX','V','IV','I')
   result = ""
   for i in range(len(ints)):
      count = int(input / ints[i])
      result += nums[i] * count
      input -= ints[i] * count
   return result

