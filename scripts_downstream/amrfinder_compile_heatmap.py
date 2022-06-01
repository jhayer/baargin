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

    args = parser.parse_args()

    dic_genes, lst_samples = prep_gene_dic(args.input_dir)

    write_output_tsv(dic_genes, lst_samples, args.output_file)


def prep_gene_dic(in_dir):

    gene_dic = {}
    sample_list = []

    with os.scandir(in_dir) as entries:
        for entry in entries:
            print(entry.name)
            if entry.is_file() and entry.name.endswith('.txt'):
                tsv_handle = open(entry, 'r')
                # remove suffix in files names
                # TODO put as options
                suf="_L002_AMRfinder.txt"
                sample=entry.name.rstrip(suf)
                sample_list.append(sample)
                next(tsv_handle)

                for line in tsv_handle:
                    lst_line = line.split('\t')
                    gene_symbol = lst_line[5]
                    # if the gene is already in the dic:
                    # add the sample with value 1 for presence
                    if gene_symbol in gene_dic.keys():
                        gene_dic[gene_symbol].append(sample)
                    else:
                        # new gene, add it to the dic first
                        gene_dic[gene_symbol]=[]
                        gene_dic[gene_symbol].append(sample)

                tsv_handle.close()

    return(gene_dic,sample_list)


def write_output_tsv(gene_dic, sample_lst, tsv_file):

    save_handle = open(tsv_file, 'w')
    # sort the list of sample_ids to always keep the same order while writing values
    sorted_sample_lst=sorted(sample_lst)
    line ="genes"
    for sample_id in sorted_sample_lst:
        line=line+"\t"+sample_id

    save_handle.write(line + os.linesep)

    for gene, samples in gene_dic.items():
        line=gene

        for s_id in sorted_sample_lst:
            pres=0
            if(s_id in samples):
                pres=1
            line=line+"\t"+str(pres)

        save_handle.write(line + os.linesep)

    save_handle.close()


if __name__ == '__main__':
    main()
