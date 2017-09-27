#!/usr/bin/env python3
#  -*- coding: utf-8 -*-

# 2017-07-30
# STB

import sys
import regex as re
import pandas as pd
from argparse import ArgumentParser
from file_read_backwards import FileReadBackwards


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
    parser.add_argument("-i", type=str, dest="gff_path")
    parser.add_argument("-o", type=str, dest="output")

    parserargs = parser.parse_args()

    try:
        if parserargs.license:
            print(MIT_license)

        else:
            gff_path = parserargs.gff_path
            output_path = parserargs.output

            if gff_path is None and output_path is None:
                parser.print_help()
                raise SystemExit

            else:
                return gff_path, output_path

    except SystemExit:
        sys.exit()


def read_and_mod_gff(annot_gff):
    # Read the GFF file and replace the values "attributes" column with just the IDs of the annotated elements
    # i.e. something like "ID=id364475;Parent=gene41724;Dbxref=GeneID:19989172;exon_number=1;gbkey=exon;gene=cox2"
    # becomes "id364475"

    id_regex = re.compile(r"(?<=ID=).+?(?=($|;))")
    parent_regex = re.compile(r"(?<=Parent=).+?(?=($|;))")
    loc_regex = re.compile(r"(?<=gene=).+?(?=($|;))")
    name_regex = re.compile(r"(?<=Name=).+?(?=($|;))")

    id_dict = {}
    last_name = ""

    gff_line_list = list()
    print("Reading GFF...")
    gff_in = FileReadBackwards(annot_gff, encoding="utf-8")
    for line in gff_in:

        if re.match("\w", line) and not line.startswith('#'):
            tab_elements = line.split("\t")
            type = tab_elements[2]
            attributes = tab_elements[-1]

            try:
                element_id = re.search(id_regex, attributes)
                element_id = element_id.group()
            except:
                element_id = "."

            if element_id == ".":
                id_dict[type] = id_dict.get(type, 0) + 1
                element_id = type + str(id_dict[type]-1)

            try:
                loc_name = re.search(loc_regex, attributes)
                loc_name = loc_name.group()
            except:
                loc_name = "."

            try:
                gene_name = re.search(name_regex, attributes)
                gene_name = gene_name.group()
                last_name = gene_name
            except:
                gene_name = "."

            if gene_name == ".":
                gene_name = last_name

            if gene_name != "." and loc_name == ".":
                loc_name = gene_name

            try:
                parent_gene = re.search(parent_regex, tab_elements[-1])
                parent_gene = parent_gene.group()
            except:
                parent_gene = "."

            if parent_gene == "." and gene_name != ".":
                parent_gene = gene_name

            tab_elements[-1] = element_id
            tab_elements += [parent_gene, loc_name]

            gff_line_list.append(tab_elements)

    gff_cols = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "id", "parent", "gene_name"]
    gff_df = pd.DataFrame(gff_line_list, columns=gff_cols)

    return gff_df


def main():
    gff_path, output_path = parse_input()

    gff_df = read_and_mod_gff(gff_path)

    print("Saving GFF3...")
    gff_df.to_csv(output_path, sep="\t", index=False)

main()
