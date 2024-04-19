#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import csv
import re

from Bio.Seq import Seq
from Bio import SeqIO

def main():
    import argparse

    parser = argparse.ArgumentParser(description="To split the AMR output file (AMRFinder or CARD) contigs detected as plasmids by Platon")
    parser.add_argument("-i", "--input_platon",
                        help="tsv file from Platon", required=True)
    parser.add_argument("-a", "--amr_file",
                        help="TSV file of AMR prediction to split by plamid and chromosome (from AMRFinderPlus or CARD", required=True)
    parser.add_argument("-p", "--output_plasmids",
                        help="output AMR TSV file containing the AMR genes of sequences detected as  \
                        plasmids by Platon", required=True)
    parser.add_argument("-c", "--output_chromosome",
                        help="output AMR TSV file containing the AMR genes of sequences predicted as chromosome", required=True)
    parser.add_argument("-t", "--amr_tool",
                        help="tool used for the prediction: amrfinder or card", \
                        required=True)
 
    args = parser.parse_args()

    plasmids_contigs_id_list = parse_platon_tsv(args.input_platon)
    split_amr_file_by_plasmid(plasmids_contigs_id_list, args.amr_file, args.amr_tool, args.output_plasmids, args.output_chromosome)


def parse_platon_tsv(platon_tsv): 
    plasmid_contigs_lst = []

    tsv_handle = open(platon_tsv, 'r')
    # skip first line
    next(tsv_handle)

    for line in tsv_handle:
        lst_line = line.split('\t')
        plasmid_contigs_lst.append(lst_line[0])

    tsv_handle.close()
    print(plasmid_contigs_lst)
    return(plasmid_contigs_lst)

def split_amr_file_by_plasmid(plas_cont_lst, amr_tsv_in, amr_tool, amr_tsv_out_plas, amr_tsv_out_chr):

    tsv_amr_handle = open(amr_tsv_in, 'r')

    p_csvfile = open(amr_tsv_out_plas, 'w', newline='')
  #  p_samplewriter = csv.writer(p_csvfile, delimiter='\t')

    c_csvfile = open(amr_tsv_out_chr, 'w', newline='')
   # c_samplewriter = csv.writer(p_csvfile, delimiter='\t')

    if(amr_tool == "amrfinder"):
        first_line = tsv_amr_handle.readline()
        # put the header line in both chromosome and plasmid
        p_csvfile.write(first_line)
        c_csvfile.write(first_line)
        for line in tsv_amr_handle:
            lst_line = line.split('\t')
            # check if the contig in is the plasmid contigs list from platon
            if(lst_line[1] in plas_cont_lst):
                print("In the list")
                p_csvfile.write(line)
            else:
                # if not, write in the chromosome output file
                c_csvfile.write(line)

    else:
        # for card output tsv
        if(amr_tool == "card"):
            # here the header looks like this
            # ORF_ID	Contig	Start	Stop	Orientation	Cut_Off	Pass_Bitscore	Best_Hit_Bitscore	Best_Hit_ARO	Best_Identities	ARO	Model_type	SNPs_in_Best_Hit_ARO	Other_SNPs	Drug Class	Resistance Mechanism	AMR Gene Family	Predicted_DNA	Predicted_Protein	CARD_Protein_Sequence	Percentage Length of Reference Sequence	ID	Model_ID	Nudged	Note
            # NODE_1_length_367850_cov_21.288374_25 # 21080 # 22246 # -1 # ID=1_25;partial=00;start_type=ATG;rbs_motif=GGA/GAG/AGG;rbs_spacer=5-10bp;gc_cont=0.429	NODE_1_length_367850_cov_21.288374_25	21080	22246	-	Strict	700	787.719	ugd	98.97	3003577	protein homolog model	n/a	n/a	peptide antibiotic	antibiotic target alteration	pmr phosphoethanolamine transferase	ATGAA...
            # problem is that the ORF numb is stick to the contigID: NODE_1_length_367850_cov_21.288374_25 =>> need to get only NODE_1_length_367850_cov_21.288374
            # remove the last _XX
            first_line = tsv_amr_handle.readline()
            # put the header line in both chromosome and plasmid
            p_csvfile.write(first_line)
            c_csvfile.write(first_line)

            for line in tsv_amr_handle:
                lst_line = line.split('\t')
                orf_num = lst_line[1]
                # remove suffix _XX that gives ORF num at the end of contig_id
                contig_id = re.sub("_[0-9]*$", "", orf_num.rstrip(" ")) 

                if(contig_id in plas_cont_lst):
                    p_csvfile.write(line)
                else:
                    c_csvfile.write(line)
    
        else:
            print("Unkown tool, please state either amrfinder or card")
    
    tsv_amr_handle.close()
    p_csvfile.close()
    c_csvfile.close()


if __name__ == '__main__':
    main()
