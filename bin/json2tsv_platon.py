#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
#import csv
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

    # write columns headers for inc_type output file
    inctype_handle = open(args.inc_types_file, 'w')
    inc_header = "Contig_id"+"\t"+"Inc_type"+"\t"+"Identity"+"\t"+"Coverage"
    inctype_handle.write(inc_header + os.linesep)

    # write columns headers for plasmid hits output file
    plasmid_handle = open(args.plasmids_file, 'w')
    plasmid_header = "Contig_id"+"\t"+"Plasmid_accession"+"\t"+"Identity"+"\t"+"coverage" \
        +"\t"+"plasmid_length"+"\t"+"contig_start"+"\t"+"contig_end"+"\t" \
        +"plasmid_start"+"\t"+"plasmid_end"
    plasmid_handle.write(plasmid_header + os.linesep)

    # write columns headers for amr hits output file
    if args.amr_file:
        amr_handle = open(args.amr_file, 'w')
        amr_header = "Contig_id"+"\t"+"AMR_type"+"\t"+"Gene"+"\t"+"Product" \
            +"\t"+"Evalue"+"\t"+"Start"+"\t"+"End"+"\t"+"Strand"
        amr_handle.write(amr_header + os.linesep)

    # read the json data
    with open(args.json_input) as fd:
        data = json.load(fd)

    for contig, dics in data.items():

        for key, vals in dics.items():
            if key == "inc_types" and len(vals) > 0:
                # prepare lines for inc_types_file
                for inc in vals:
                    inc_line = contig+"\t"+inc['type']+"\t"+str(inc['identity']) \
                        +"\t"+str(inc['coverage'])
                    inctype_handle.write(inc_line + os.linesep)

            if key == "plasmid_hits":
                for ph in vals:
                    #[{'contig_start': 1, 'contig_end': 7073, 'plasmid_start': 3267,
                    # 'plasmid_end': 10339, 'plasmid': {'id': 'NC_021155.1', 'length': 148711},
                    # 'coverage': 1.0, 'identity': 0.9998586172769688}]
                    ph_line = contig+"\t"+ph['plasmid']['id']+"\t" \
                        +str(ph['identity'])+"\t"+str(ph['coverage'])+"\t" \
                        +str(ph['plasmid']['length'])+"\t"+str(ph['contig_start'])+"\t" \
                        +str(ph['contig_end'])+"\t"+str(ph['plasmid_start'])+"\t" \
                        +str(ph['plasmid_end'])
                    plasmid_handle.write(ph_line + os.linesep)


            if key == "amr_hits" and args.amr_file:
                for amr in vals:
                    amr_line = contig+"\t"+amr['type']+"\t" \
                        +amr['gene']+"\t"+amr['product']+"\t" \
                        +str(amr['evalue'])+"\t"+str(amr['start'])+"\t" \
                        +str(amr['end'])+"\t"+str(amr['strand'])
                    amr_handle.write(amr_line + os.linesep)

    # open the csv file as a DictWriter using those names
    # with open(args.inc_types_file, "w", newline='') as of:
    #     wr = csv.DictWriter(of, names)
    #     wr.writeheader()
    #     for field, vals in data.items():
    #         d['field'] = field
    #         for inner in vals:
    #             for k,v in inner.items():
    #                 d[k] = v
    #         wr.writerow(d)



if __name__ == '__main__':
    main()
