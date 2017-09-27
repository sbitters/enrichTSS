#!/usr/bin/env python3
#  -*- coding: utf-8 -*-

# 2017-07-30
# STB

import sys
import pysam
import regex as re
import pandas as pd
import numpy as np
from argparse import ArgumentParser


MIT_license = """Copyright 2017 Sven T. Bitters (sven.bitters@gmail.com)

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""


def parse_input():

    parser = ArgumentParser()
    parser.add_argument("-l", action="store_true", default=False, dest="license")
    parser.add_argument("-t", type=str, dest="target_gene")
    parser.add_argument("-e", type=int, dest="target_ebe")
    parser.add_argument("-n", type=int, dest="length_ebe")
    parser.add_argument("-w", type=int, default=55, dest="window_size")
    parser.add_argument("-x", type=str, dest="expression_table")
    parser.add_argument("-c", type=str, dest="condition")
    parser.add_argument("-g", type=str, dest="gff_path")
    parser.add_argument("-m", type=str, dest="main_bams")
    parser.add_argument("-s", type=str, nargs="+", dest="set_bams",
                        help="-s <treatment.bam> <control.bam>")
    parser.add_argument("-p", type=float, default=1, dest="pseudocount")
    parser.add_argument("--libsize", type=float, nargs="+", dest="lib_sizes",
                        help="--libsize <#treatment> <#control>")
    parser.add_argument("--is_pe", action="store_true", default=False, dest="is_pe")
    parser.add_argument("-o", type=str, dest="output")

    parserargs = parser.parse_args()

    try:
        if parserargs.license:
            print(MIT_license)

        else:
            target_gene = parserargs.target_gene

            if target_gene is None:
                parser.print_help()
                raise SystemExit

            else:
                target_ebe = parserargs.target_ebe
                length_ebe = parserargs.length_ebe
                window_size = parserargs.window_size
                expression_path = parserargs.expression_table
                experiment_cond = parserargs.condition
                gff_path = parserargs.gff_path
                main_bams = parserargs.main_bams
                set_bams = parserargs.set_bams
                pseudocount = parserargs.pseudocount
                output_file = parserargs.output
                lib_sizes = parserargs.lib_sizes
                is_pe = parserargs.is_pe

                return target_gene, target_ebe, length_ebe, window_size, \
                       expression_path, experiment_cond, gff_path, main_bams, \
                       set_bams, pseudocount, output_file, lib_sizes, is_pe

    except SystemExit:
        sys.exit()


def handle_targets(target_gene, gff_gene, expression_table, window_size, target_ebe, length_ebe):
    target_gene_df_tmp = expression_table[(expression_table["Gene_IDs"] == target_gene)]

    comparison_group_tmp = expression_table[(expression_table["group_TPM"] == target_gene_df_tmp.iloc[0]["group_TPM"])].copy()

    start_list = list()
    end_list = list()
    strand_list = list()
    length_list = list()

    comp_grp_genes = list(comparison_group_tmp["Gene_IDs"])
    for gene in comp_grp_genes:
        gene_info = gff_gene[(gff_gene["gene_name"] == gene)]
        start_list.append(gene_info.iloc[0]["start"])
        end_list.append(gene_info.iloc[0]["end"])
        strand_list.append(gene_info.iloc[0]["strand"])
        length_list.append(gene_info.iloc[0]["end"] - gene_info.iloc[0]["start"])

    comparison_group_tmp.loc[:, "start"] = start_list
    comparison_group_tmp.loc[:, "end"] = end_list
    comparison_group_tmp.loc[:, "strand"] = strand_list
    comparison_group_tmp.loc[:, "length"] = length_list

    target_gene_df = comparison_group_tmp[(comparison_group_tmp["Gene_IDs"] == target_gene)]

    target_length = target_gene_df.iloc[0]["length"]
    comparison_group = comparison_group_tmp[(comparison_group_tmp["length"] <= target_length * 1.2) & (comparison_group_tmp["length"] >= target_length * 0.8)].copy()
    comparison_group = comparison_group.dropna(subset=["Gene_IDs"])
    comparison_group = comparison_group.reset_index()

    if target_gene_df.iloc[0]["strand"] == "+":
        EBE_loc = target_gene_df.iloc[0]["start"] + target_ebe + length_ebe
        window_start = EBE_loc + 35
        window_end = window_start + window_size

        dist_winstart_geneannot = window_start - target_gene_df.iloc[0]["start"]
    else:
        EBE_loc = target_gene_df.iloc[0]["end"] - target_ebe - length_ebe
        window_end = EBE_loc - 35
        window_start = window_end - window_size

        dist_winstart_geneannot = target_gene_df.iloc[0]["end"] - window_end

    return comparison_group, target_gene_df, \
           window_start, window_end, dist_winstart_geneannot


def get_window(current_gene, target_gene, dist_winstart_geneannot,
               target_window_start, target_window_end, window_size, gff_full):

    current_gene_name = current_gene.iloc[0]["Gene_IDs"]
    current_gene_start = current_gene.iloc[0]["start"]
    current_gene_end = current_gene.iloc[0]["end"]
    current_gene_strand = current_gene.iloc[0]["strand"]

    skip_gene = False

    try:
        if current_gene_name != target_gene:
            exon_df = gff_full[(gff_full["type"] == "exon")].copy()
            exon_df = exon_df[(exon_df["gene_name"] == current_gene_name)]

            min_start = min(exon_df["start"])
            max_end = max(exon_df["end"])

            window_start = -1
            window_end = -1
            found_start = False
            found_end = False
            count = 0

            if current_gene_strand == "+":
                exon_df = exon_df.sort_values("start", axis=0)

                if len(exon_df) == 0:
                    skip_gene = True

                else:
                    if dist_winstart_geneannot <= 0:
                        window_start = current_gene_start + dist_winstart_geneannot
                        window_end = window_start + window_size
                    else:

                        if current_gene_start + dist_winstart_geneannot >= current_gene_end:
                            skip_gene = True
                        else:
                            position = current_gene_start
                            prelim_start = position + dist_winstart_geneannot

                            res = check_if_exon(exon_df, prelim_start, count, "+")
                            is_exon = res[2]
                            skip_gene = res[3]

                            if is_exon:
                                window_start = prelim_start
                                position = window_start
                                found_start = True

                                same_exon = check_if_same_exon(window_start, window_start + window_size, exon_df)
                                if same_exon:
                                    window_end = window_start + window_size

                if not skip_gene:
                    while window_end == -1:

                        if not found_start:
                            same_exon = check_if_same_exon(position, position + window_size, exon_df)
                            if same_exon:
                                window_start = position
                                window_end = position + window_size
                                break


                        if not found_start and count == dist_winstart_geneannot:
                            window_start = position
                            count = 0
                            found_start = True

                        if found_start and count == window_size:
                            window_end = position
                            break

                        if (position > max_end) or ((count > (dist_winstart_geneannot * 1.1)) and (count > (window_size * 1.1))):
                            skip_gene = True
                            break

                        res = check_if_exon(exon_df, position, count, "+")
                        count = res[0]
                        position = res[1]
                        skip_gene = res[3]

                        if skip_gene:
                            break

                        position += 1

            else:
                exon_df = exon_df.sort_values("end", axis=0, ascending=False)

                if len(exon_df) == 0:
                    skip_gene = True

                else:
                    if dist_winstart_geneannot <= 0:
                        window_end = current_gene_end - dist_winstart_geneannot
                        window_start = window_end - window_size
                    else:
                        if current_gene_end - dist_winstart_geneannot <= current_gene_start:
                            skip_gene = True
                        else:
                            position = current_gene_end
                            prelim_end = position - dist_winstart_geneannot

                            res = check_if_exon(exon_df, prelim_end, count, "-")
                            is_exon = res[2]
                            skip_gene = res[3]

                            if is_exon:
                                window_end = prelim_end
                                position = window_end
                                found_end = True

                                same_exon = check_if_same_exon(window_end, window_end - window_size, exon_df)
                                if same_exon:
                                    window_start = window_end - window_size

                if not skip_gene:
                    while window_start == -1:

                        if not found_end:
                            same_exon = check_if_same_exon(position, position - window_size, exon_df)
                            if same_exon:
                                window_end = position
                                window_start = position - window_size
                                break

                        if not found_end and count == dist_winstart_geneannot:
                            window_end = position
                            count = 0
                            found_end = True

                        if found_end and count == window_size:
                            window_start = position

                        if (position < min_start) or ((count > (dist_winstart_geneannot * 1.1)) and (count > (window_size * 1.1))):
                            skip_gene = True
                            break

                        res = check_if_exon(exon_df, position, count, "-")
                        count = res[0]
                        position = res[1]
                        skip_gene = res[3]

                        if skip_gene:
                            break

                        position -= 1

        else:
            window_start = target_window_start
            window_end = target_window_end

    except:
        skip_gene = True

    if skip_gene:
        window_start = np.NaN
        window_end = np.NaN

    return window_start, window_end


def check_if_exon(exon_df, position, count, strandedness):
    is_exon = False
    skip_gene = False
    count += 1

    if strandedness == "+":
        for index in range(1, len(exon_df)):
            start = exon_df.iloc[index]["start"]
            end = exon_df.iloc[index]["end"]

            try:
                if position <= end and position >= start:
                    is_exon = True
                    break
                elif position > end and position < exon_df.iloc[index+1]["start"]:
                    position = exon_df.iloc[index+1]["start"]
                    is_exon = False
                    break
            except:
                skip_gene = True

    else:
        for index in range(1, len(exon_df)):
            start = exon_df.iloc[index]["start"]
            end = exon_df.iloc[index]["end"]

            try:
                if position <= end and position >= start:
                    is_exon = True
                    break
                elif position < start and position > exon_df.iloc[index+1]["end"]:
                    position = exon_df.iloc[index+1]["end"]
                    is_exon = False
                    break
            except:
                skip_gene = True

    return [count, position, is_exon, skip_gene]


def check_if_same_exon(pos1, pos2, exon_df):

    same_exon = False

    try:
        for index in range(1, len(exon_df)):
            start = exon_df.iloc[index]["start"]
            end = exon_df.iloc[index]["end"]

            if (pos1 >= start and pos1 <= end) and (pos2 >= start and pos2 <= end):
                same_exon = True
    except:
        same_exon = False

    return same_exon


def get_read_count(bamfiles, set_bams, window_start, window_end, experiment_cond, current_gene_chrom, count_treat_list, count_contr_list):

    for ii in range(0, len(set_bams)):
        bamfile = bamfiles[ii]

        try:
            count = bamfile.count(current_gene_chrom, window_start+21, window_end-19, read_callback="all")
        except:
            count = np.NaN

        if ii == 0:
            count_treat_list.append(count)
        else:
            count_contr_list.append(count)

    return count_contr_list, count_treat_list



def main():
    target_gene, target_ebe, length_ebe, window_size, expression_path, experiment_cond, gff_path, \
    main_bams, set_bams, pseudocount, output_file, lib_sizes, is_pe = parse_input()

    gff_full = pd.read_table(gff_path, header=0, dtype={"seqid": str, "source": str, "type": str,
                                                        "score": str, "strand": str, "phase":str,
                                                        "id": str, "parent": str, "gene_name": str,
                                                        "start": np.int64, "end": np.int64})
    gff_gene = gff_full[(gff_full["type"] == "gene")]

    expression_table = pd.read_table(expression_path, header=0)

    comparison_group, target_gene_df, target_window_start, target_window_end, dist_winstart_geneannot \
        = handle_targets(target_gene, gff_gene, expression_table, window_size, target_ebe, length_ebe)

    start_list = list()
    end_list = list()
    count_treat_list = list()
    count_contr_list = list()

    bampath_1 = main_bams + "/" + set_bams[0]
    bamfile_1 = pysam.AlignmentFile(bampath_1, "rb")

    bampath_2 = main_bams + "/" + set_bams[1]
    bamfile_2 = pysam.AlignmentFile(bampath_2, "rb")

    bamfiles = [bamfile_1, bamfile_2]

    last = 0
    for index in range(0, len(comparison_group)):

        percent = round(((index+1) / len(comparison_group)) * 100, 1)
        if percent != last:
            last = percent
            print("Progress: {:>5} %".format(percent), end="\r")


        current_gene = comparison_group.iloc[[index]]
        current_gene_info = gff_gene[gff_gene["gene_name"] == current_gene.iloc[0]["Gene_IDs"]]
        current_gene_chrom = current_gene_info.iloc[0]["seqid"]

        window_start, window_end = \
            get_window(current_gene, target_gene, dist_winstart_geneannot,
                       target_window_start, target_window_end, window_size, gff_full)

        start_list.append(window_start)
        end_list.append(window_end)

        count_contr_list, count_treat_list = \
            get_read_count(bamfiles, set_bams, window_start, window_end, experiment_cond,
                           current_gene_chrom, count_treat_list, count_contr_list)

    bamfile_1.close()
    bamfile_2.close()

    if is_pe:
        pe_factor = 2
    else:
        pe_factor = 1

    libsize_contr = lib_sizes[1] * pe_factor
    libsize_treat = lib_sizes[0] * pe_factor

    cpm_contr_list = [((count / libsize_contr) * 1e6) for count in count_contr_list]
    cpm_treat_list = [((count / libsize_treat) * 1e6) for count in count_treat_list]

    comparison_group.loc[:, "window_start"] = start_list
    comparison_group.loc[:, "window_end"] = end_list
    comparison_group.loc[:, "count_control"] = count_contr_list
    comparison_group.loc[:, "count_treatment"] = count_treat_list
    comparison_group.loc[:, "cpm_control"] = cpm_contr_list
    comparison_group.loc[:, "cpm_treatment"] = cpm_treat_list

    for ii in range(0, len(comparison_group)):
        if comparison_group.loc[ii, "cpm_control"] > 0:
            comparison_group.set_value(ii, "cpm_enrichment", round((comparison_group.loc[ii, "cpm_treatment"]) / (comparison_group.loc[ii, "cpm_control"]), 4))
        else:
            comparison_group.set_value(ii, "cpm_enrichment", round((comparison_group.loc[ii, "cpm_treatment"]), 4))

        if comparison_group.loc[ii, "cpm_enrichment"] == 0:
            comparison_group.set_value(ii, "cpm_enrichment", 0.0001)

    comparison_group = comparison_group.sort_values("cpm_enrichment", axis=0, ascending=False)

    comparison_group = comparison_group.dropna()

    comparison_group.to_csv(output_file, sep="\t", index=False, header=True)
    print("")

main()
