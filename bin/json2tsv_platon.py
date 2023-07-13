#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import csv
import json

def main():
    import argparse

    parser = argparse.ArgumentParser(description="To convert Platon json output into tsv")
    parser.add_argument("-j", "--json_input",
                        help="input json file from platon", required=True)
    parser.add_argument("-i", "--inc_types_file",
                        help="output tabular file containing plasmids inc types",
                        required=True)
    parser.add_argument("-p", "--plasmids_file",
                        help="output tabular file containing plasmids hits from platon",
                        required=True)
    parser.add_argument("-a", "--amr_file",
                        help="output tabular file containing ARGs hits from platon",
                        required=False)

    args = parser.parse_args()

    csvfile_inc = open(args.inc_types_file, 'w', newline='')
    # write columns headers for inc_type output file
    inctype_handle = csv.writer(csvfile_inc, delimiter='\t')
    inc_header = ["Contig_id","Inc_type","Identity","Coverage"]
    inctype_handle.writerow(inc_header)

    csvfile_pla = open(args.plasmids_file, 'w', newline='')
    # write columns headers for plasmid hits output file
    plasmid_handle = csv.writer(csvfile_pla, delimiter='\t')
    plasmid_header = ["Contig_id","Plasmid_accession","Identity","coverage", \
            "plasmid_length","contig_start","contig_end","plasmid_start","plasmid_end"]
    plasmid_handle.writerow(plasmid_header)

    # write columns headers for amr hits output file
    if args.amr_file:
        csvfile_amr = open(args.amr_file, 'w', newline='')

        amr_handle = csv.writer(csvfile_amr, delimiter='\t')
        amr_header = ["Contig_id","AMR_type","Gene","Product","Evalue","Start","End","Strand"]
        amr_handle.writerow(amr_header)

    # read the json data
    with open(args.json_input) as fd:
        data = json.load(fd)

    for contig, dics in data.items():

        for key, vals in dics.items():
            if key == "inc_types" and len(vals) > 0:
                # prepare lines for inc_types_file
                for inc in vals:
                    inc_line = [contig,inc['type'],str(inc['identity']),str(inc['coverage'])]
                    inctype_handle.writerow(inc_line)

            if key == "plasmid_hits":
                for ph in vals:
                    #[{'contig_start': 1, 'contig_end': 7073, 'plasmid_start': 3267,
                    # 'plasmid_end': 10339, 'plasmid': {'id': 'NC_021155.1', 'length': 148711},
                    # 'coverage': 1.0, 'identity': 0.9998586172769688}]
                    ph_line = [contig,ph['plasmid']['id'],str(ph['identity']), \
                        str(ph['coverage']),str(ph['plasmid']['length']), \
                        str(ph['contig_start']),str(ph['contig_end']), \
                        str(ph['plasmid_start']),str(ph['plasmid_end'])]
                    plasmid_handle.writerow(ph_line)


            if key == "amr_hits" and args.amr_file:
                for amr in vals:
                    amr_line = [contig,amr['type'],amr['gene'],amr['product'], \
                        str(amr['evalue']),str(amr['start']), \
                        str(amr['end']),str(amr['strand'])]
                    amr_handle.writerow(amr_line)


    csvfile_inc.close()
    csvfile_pla.close()
    csvfile_amr.close()

if __name__ == '__main__':
    main()
