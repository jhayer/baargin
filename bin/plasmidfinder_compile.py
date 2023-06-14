#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import csv

def main():
    import argparse

    parser = argparse.ArgumentParser(description="To compile PlasmidFinder \
        results for several samples, as table with 0 and 1 for absence/presence")
    parser.add_argument("-i", "--input_dir",
                        help="input directory containing all PlasmidFinder \
                        output files", required=True)
    parser.add_argument("-o", "--output_file",
                        help="output tabular file containing detected plasmids\
                        for all samples", required=True)
    parser.add_argument("-s", "--files_suffix",
                        help="Suffix of PF ouput files to remove to get \
                        sample_ids", required=False)

    args = parser.parse_args()

    # remove suffix in files names
    suffix="_plasmidfinder_raw_results.tsv"
    if(args.files_suffix):
        suffix=args.files_suffix

    dic_plasmids, lst_samples = prep_plasmid_dic(args.input_dir, suffix)

    write_output_tsv(dic_plasmids, lst_samples, args.output_file)


def prep_plasmid_dic(in_dir, suf):

    plasmid_dic = {}
    sample_list = []

    with os.scandir(in_dir) as entries:

        for entry in entries:
        #    print(entry.name)
            if entry.is_file() and entry.name.endswith('.tsv'):
                tsv_handle = open(entry, 'r')

                sample=entry.name.rstrip(suf)
                sample_list.append(sample)
                next(tsv_handle)

                for line in tsv_handle:
                    lst_line = line.split('\t')
                    print(lst_line)

                    plasmid_symbol = lst_line[1]

                    # if the plasmid is already in the dic:
                    # add the sample_id to the dic
                    if plasmid_symbol in plasmid_dic.keys():
                        # here the same sample can be added twice or more
                        # if several contigs belong to the same plasmid
                        plasmid_dic[plasmid_symbol].append(sample)
                    else:
                        # new plasmid, add it to the dic first
                        plasmid_dic[plasmid_symbol]=[]
                        plasmid_dic[plasmid_symbol].append(sample)

                tsv_handle.close()

    return(plasmid_dic,sample_list)


def write_output_tsv(plasmid_dic, sample_lst, tsv_file):

#    save_handle = open(tsv_file, 'w')
    with open(tsv_file, 'w', newline='') as csvfile:
        samplewriter = csv.writer(csvfile, delimiter='\t')
        # sort the list of sample_ids to always keep the same order while writing values
        sorted_sample_lst=sorted(sample_lst)
        line =["Plasmid"]
        for sample_id in sorted_sample_lst:
            line.append(sample_id)
        # write the tsv file header
        samplewriter.writerow(line)

        for plasmid, samples in plasmid_dic.items():
            line=[plasmid]

            for s_id in sorted_sample_lst:
                pres=0
                if(s_id in samples):
                    # the number of times that the sample_id is found for this plasmid
                    pres=samples.count(s_id)
                line.append(str(pres))

            samplewriter.writerow(line)


if __name__ == '__main__':
    main()
