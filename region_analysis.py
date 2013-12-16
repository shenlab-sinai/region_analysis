#! /usr/bin/env python

import os
import sys
import json
from itertools import groupby
from optparse import OptionParser
import pybedtools
import pybedtools.featurefuncs


def main():
    opt_parser = OptionParser()
    opt_parser.add_option('-i', '--input', action='store',
                          dest='input_file_name',
                          help='Input region file must assume the first 3 columns contain (chr, start, end)')
    opt_parser.add_option('-d', '--database', action='store',
                          dest='anno_db', help='Choose database: refseq(default) or ensembl',
                          default='refseq')
    opt_parser.add_option('-r', '--rhead', action='store_true',
                          dest='rhead', help='Whether the input file contains column header', default=False)
    opt_parser.add_option('-g', '--genome', action='store',
                          dest='genome', help='Choose genome: mm10(default)',
                          default='mm10')
    try:    
        (options, args) = opt_parser.parse_args(sys.argv)
        script_dir = os.path.dirname(os.path.realpath(__file__))
        db_path = os.path.join(script_dir, "database/")
        input_file_name = options.input_file_name
        anno_db = options.anno_db
        rhead = options.rhead
        genome = options.genome
        if (input_file_name is None) or (len(input_file_name)==0):
            raise SystemExit
    except SystemExit:
        print("Please assign proper input file!")
        opt_parser.print_help()
        return 1

    # create a tmp bed file with index column
    in_f = file(input_file_name)
    input_filtered = [ line  for line in in_f  if not line.startswith("#") ] # filter the comment lines
    input_indexed = [ '%s\t%d\n' % (line.strip(), i) for i, line in enumerate(input_filtered) ] # add index column to the bed lines
    in_f.close()

    # read all annotations into a dictionary, for the further output.
    anno_bed = os.path.join(
        db_path, genome + "." + anno_db + ".biotype_region_ext.bed")

    # use saveas() to convert the BedTool objects to file-based objects,
    # so they could be used multiple times.
    # When debug, we may use saveas("tss.tmp"), and the output of bedtools
    # could be saved.
    pybedtools.set_tempdir("./")
    anno = pybedtools.BedTool(anno_bed).saveas()
    gd = pybedtools.BedTool(
        os.path.join(db_path, genome + "_geneDesert.bed")).saveas()
    pc = pybedtools.BedTool(
        os.path.join(db_path, genome + "_pericentromere.bed")).saveas()
    st = pybedtools.BedTool(
        os.path.join(db_path, genome + "_subtelomere.bed")).saveas()

    # load the input intervals to be annotated
    if rhead == True:
        headlineL = input_indexed[0].strip().split("\t")[:-1]
        del input_indexed[0]
    input_bed = pybedtools.BedTool(
        "".join(input_indexed), from_string=True).saveas()
    list_input = [x.fields[:] for x in input_bed]
    col_no_input = input_bed.field_count()
    # get the midpoint of the intervals
    input_bed_mid = input_bed.each(pybedtools.featurefuncs.midpoint).saveas()

    # intersectBed with annotations
    input_GB = input_bed_mid.intersect(anno, wao=True).saveas()
    list_GB = [x.fields[:] for x in input_GB]
    input_gd = input_bed_mid.intersect(gd, c=True, f=0.5).saveas()
    list_gd = [x.fields[col_no_input + 0] for x in input_gd]
    input_pc = input_bed_mid.intersect(pc, c=True, f=0.5).saveas()
    list_pc = [x.fields[col_no_input + 0] for x in input_pc]
    input_st = input_bed_mid.intersect(st, c=True, f=0.5).saveas()
    list_st = [x.fields[col_no_input + 0] for x in input_st]

    # groupby the intersectBed results based on the index column
    input_idx = key = lambda s: s[col_no_input - 1]
    GB_dict = {}
    for key, GB_hits in groupby(list_GB, key=input_idx):
        GB_dict[key] = list(v for v in GB_hits)

    def getBestHit(anno_db, col_no_input, GB_entry, gd_entry, st_entry, pc_entry):
        """
        format all hits and get the best hit that is neartest to TSS
        """
        best_hit, cur_output, Dis2TSS, Dis2TES = (None, None, None, None)
        formatted = []

        def getDis2TSS(anno_db, cur_input, col_no_input):
            """
            Calculate the distance to TSS and decide the annotation feature of the entry
            """
            # calculate TSS and TES based on strands
            # if strand is not "-" then it will be treated as "+"
            if cur_input[col_no_input + 5] != "-":
                TSS = int(cur_input[col_no_input + 10])
                TES = int(cur_input[col_no_input + 11])
                cur_input[col_no_input + 1] = str(int(TSS))
                cur_input[col_no_input + 2] = str(int(TES))
                Dis2TSS = int(cur_input[1]) - TSS
                Dis2TES = int(cur_input[1]) - TES
            else:
                TSS = int(cur_input[col_no_input + 11])
                TES = int(cur_input[col_no_input + 10])
                cur_input[col_no_input + 1] = str(int(TES))
                cur_input[col_no_input + 2] = str(int(TSS))
                Dis2TSS = TSS - int(cur_input[1])
                Dis2TES = TES - int(cur_input[1])
            if abs(Dis2TSS) <= 250:
                Pos = "ProximalPromoter"
            elif abs(Dis2TSS) <= 1000:
                Pos = "Promoter1k"
            elif abs(Dis2TSS) <= 3000:
                Pos = "Promoter3k"
            else:
                Pos = "genebody"
            if anno_db == "ensembl":  # output is gid
                cur_output = [
                    cur_input[col_no_input+3], cur_input[col_no_input+6], cur_input[col_no_input+5],
                    cur_input[col_no_input+10], cur_input[col_no_input+11], Pos, str(Dis2TSS),
                    cur_input[col_no_input+9]]
            else:  # output is gene symbol
                cur_output = [cur_input[col_no_input+4], cur_input[col_no_input+6], cur_input[col_no_input+5],
                    cur_input[col_no_input+10], cur_input[col_no_input+11], Pos, str(Dis2TSS),
                    cur_input[col_no_input+9]]
            return (cur_output, Dis2TSS, Dis2TES)

        for i in GB_entry:
            # discard null hit of genebody entry
            if not ((i[col_no_input] == ".") and (i[col_no_input + 1] == "-1")):
                cur_output, Dis2TSS, Dis2TES = getDis2TSS(
                    anno_db, i, col_no_input)
                if (best_hit is None) or (abs(int(best_hit[6])) > abs(Dis2TSS)):
                    best_hit = cur_output
                formatted.append(cur_output)
        if gd_entry != "0":
            cur_output = ["", "", "", "", "", "Genedesert", "NA", "No_anno"]
            if (best_hit is None):
                best_hit = cur_output
            formatted.append(cur_output)
        if st_entry != "0":
            cur_output = ["", "", "", "", "", "Subtelomere", "NA", "No_anno"]
            if (best_hit is None):
                best_hit = cur_output
            formatted.append(cur_output)
        if pc_entry != "0":
            cur_output = ["", "", "", "", "",
                          "Pericentromere", "NA", "No_anno"]
            if (best_hit is None):
                best_hit = cur_output
            formatted.append(cur_output)
        if best_hit is None:
            cur_output = ["", "", "", "", "",
                          "Otherintergenic", "NA", "No_anno"]
            best_hit = cur_output
            formatted.append(cur_output)
        return (formatted, best_hit)

    output_file_best = file(input_file_name + ".annotated", "w")
    output_file = file(input_file_name + ".full.annotated", "w")
    output_file_json = file(input_file_name + ".full.annotated.json", "w")
    # Output the header
    if rhead == True:
        output_file.write("\t".join(
            headlineL + ["GName", "TName", "Strand", "TSS", "TES", "Feature", "D2TSS", "Biotype"]) + "\n")
        output_file_best.write("\t".join(
            headlineL + ["GName", "TName", "Strand", "TSS", "TES", "Feature", "D2TSS", "Biotype"]) + "\n")
        start_idx = 1
    else:
        start_idx = 0
    # write to the output: input.bed.annotated, input.bed.full.annotated
    json_dict = {}
    for i in range(start_idx, len(input_bed)):
        output_lineL = list_input[i][:-1]  # original input line
        json_dict[str(i)] = {}
        json_dict[str(i)]["query_interval"] = output_lineL
        formatted, best_hit = getBestHit(
            anno_db, col_no_input, GB_dict[str(i)], list_gd[i], list_st[i], list_pc[i])
        output_file_best.write("\t".join(output_lineL + best_hit) + "\n")
        json_dict[str(i)]["best_hit"] = best_hit
        for j in formatted:
            output_file.write("\t".join(output_lineL + j) + "\n")
        json_dict[str(i)]["all_hits"] = formatted
    output_file_best.close()
    output_file.close()
    json.dump(json_dict, output_file_json, sort_keys=True, indent=2)
    output_file_json.close()
    pybedtools.cleanup()
    return 0

#-------------------------------------------------------------------------
if __name__ == '__main__':
    main()
#-------------------------------------------------------------------------
# EOF
