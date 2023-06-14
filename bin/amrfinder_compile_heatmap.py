#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import csv

def main():
    import argparse

    parser = argparse.ArgumentParser(description="To compile AMRFinderPlus \
        results for several samples, as table with 0 and 1 for absence/presence")
    parser.add_argument("-i", "--input_dir",
                        help="input directory containing all amrfinder \
                        output files", required=True)
    parser.add_argument("-o", "--output_file",
                        help="output tabular file containing amr genes\
                        for all samples", required=True)
    parser.add_argument("-s", "--files_suffix",
                        help="Suffix of amrfinder ouput files to remove to get \
                        sample_ids", required=False)

    args = parser.parse_args()

    # remove suffix in files names
    suffix="_AMRfinder.txt"
    if(args.files_suffix):
        suffix=args.files_suffix

    dic_genes, lst_samples, dic_info = prep_gene_dic(args.input_dir, suffix)

    write_output_tsv(dic_genes, lst_samples, dic_info, args.output_file)


def prep_gene_dic(in_dir, suf):

    gene_dic = {}
    sample_list = []
    gene_info = {}

    with os.scandir(in_dir) as entries:

        for entry in entries:
        #    print(entry.name)
            if entry.is_file() and entry.name.endswith('.txt'):
                tsv_handle = open(entry, 'r')

                sample=entry.name.rstrip(suf)
                sample_list.append(sample)
                next(tsv_handle)

                for line in tsv_handle:
                    lst_line = line.split('\t')
                    gene_symbol = lst_line[5]
                    # retrieving info about gene, type of AMR
                    gene_info_lst = lst_line[6:12]

                    # if the gene is already in the dic:
                    # add the sample with value 1 for presence
                    if gene_symbol in gene_dic.keys():
                        gene_dic[gene_symbol].append(sample)
                    else:
                        # new gene, add it to the dic first
                        gene_dic[gene_symbol]=[]
                        gene_dic[gene_symbol].append(sample)
                        # adding gene info to specific dic
                        gene_info[gene_symbol]=gene_info_lst

                tsv_handle.close()

    return(gene_dic,sample_list,gene_info)


def write_output_tsv(gene_dic, sample_lst, gene_info_dic, tsv_file):

    with open(tsv_file, 'w', newline='') as csvfile:
        samplewriter = csv.writer(csvfile, delimiter='\t')
        # sort the list of sample_ids to always keep the same order while writing values
        sorted_sample_lst=sorted(sample_lst)
        line1 = ['Gene symbol', 'Seq name', 'Scope', 'Element type', 'Element subtype', 'Class', 'Subclass']

        for sample_id in sorted_sample_lst:
            line1.append(sample_id)
        # write the tsv file header
        samplewriter.writerow(line1)

        for gene, samples in gene_dic.items():
            line = [gene,gene_info_dic[gene][0],gene_info_dic[gene][1], \
                    gene_info_dic[gene][2],gene_info_dic[gene][3], \
                    gene_info_dic[gene][4],gene_info_dic[gene][5]]

            for s_id in sorted_sample_lst:
                pres=0
                if(s_id in samples):
                    pres=1
                line.append(str(pres))

            samplewriter.writerow(line)


if __name__ == '__main__':
    main()
