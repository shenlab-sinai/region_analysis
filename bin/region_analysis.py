#! /usr/bin/env python

import os
import sys
import json
from itertools import groupby
from optparse import OptionParser
import pybedtools
import pybedtools.featurefuncs
import regionanalysis
    
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
    opt_parser.add_option('-v', '--version', action='store_true',
                          dest='version', help='Version of Region_Analysis package')
    try:    
        (options, args) = opt_parser.parse_args(sys.argv)
        if options.version == True:
            print("Region_Analysis Version: %s\n" %regionanalysis.__version__)
        module_dir = os.path.dirname(os.path.realpath(regionanalysis.__file__))
        db_path = os.path.join(module_dir, "database/")
        input_file_name = options.input_file_name
        anno_db = options.anno_db
        rhead = options.rhead
        genome = options.genome
        if (input_file_name is None) or (len(input_file_name)==0):
            raise SystemExit
    except SystemExit:
        sys.stderr.write("Please assign proper input file!")
        opt_parser.print_help()
        return 1

    # create a tmp bed file with index column.
    in_f = file(input_file_name)
    input_filtered = [ line  for line in in_f  if not line.lstrip().startswith("#") ] # filter the comment lines
    # if there is header, store it and remove it from the query BED.
    if rhead == True:
        headlineL = input_filtered[0].strip().split("\t")
        del input_filtered[0]
    input_indexed = [ '%s\t%d\n' % (line.strip(), i) for i, line in enumerate(input_filtered) ] # add index column to the bed lines
    in_f.close()

    # read all annotations into a dictionary, for the further output.
    anno_bed = os.path.join(
        db_path, genome + "." + anno_db + ".biotype_region_ext.bed")
    try:
        if not os.path.exists(anno_bed):
    #         raise regionanalysis.NoGenomeError
    # except regionanalysis.NoGenomeError:
            raise SystemExit
    except SystemExit:
        sys.stderr.write("%s not in the genome database!"%genome)
        return 1

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

    # load the input intervals to be annotated.
    try:
        input_bed = pybedtools.BedTool(
            "".join(input_indexed), from_string=True).saveas()
    except:
        sys.stderr.write("Error in input file! Please check the format!")
        return 1
    list_input = [x.fields[:] for x in input_bed]
    col_no_input = input_bed.field_count()
    # get the midpoint of the intervals.
    # there is a bug in midpoint function of pybedtools 0.6.3, so here an alternative function was used.
    # input_bed_mid = input_bed.each(pybedtools.featurefuncs.midpoint).saveas()
    input_bed_mid = pybedtools.BedTool("".join([regionanalysis.midpoint(x) for x in input_indexed]), from_string=True).saveas()

    # intersectBed with annotations.
    input_GB = input_bed_mid.intersect(anno, wao=True).saveas()
    list_GB = [x.fields[:] for x in input_GB]
    input_gd = input_bed_mid.intersect(gd, c=True, f=0.5).saveas()
    list_gd = [x.fields[col_no_input + 0] for x in input_gd]
    input_pc = input_bed_mid.intersect(pc, c=True, f=0.5).saveas()
    list_pc = [x.fields[col_no_input + 0] for x in input_pc]
    input_st = input_bed_mid.intersect(st, c=True, f=0.5).saveas()
    list_st = [x.fields[col_no_input + 0] for x in input_st]

    # groupby the intersectBed results based on the index column.
    input_idx = key = lambda s: s[col_no_input - 1]
    GB_dict = {}
    for key, GB_hits in groupby(list_GB, key=input_idx):
        GB_dict[key] = list(v for v in GB_hits)

    output_file_best = file(input_file_name + ".annotated", "w")
    output_file = file(input_file_name + ".full.annotated", "w")
    output_file_json = file(input_file_name + ".full.annotated.json", "w")
    # Output the header.
    if rhead == True:
        output_file.write("\t".join(
            headlineL + ["GName", "TName", "Strand", "TSS", "TES", "Feature", "D2TSS", "Biotype"]) + "\n")
        output_file_best.write("\t".join(
            headlineL + ["GName", "TName", "Strand", "TSS", "TES", "Feature", "D2TSS", "Biotype"]) + "\n")
    # write to the output: input.bed.annotated, input.bed.full.annotated.
    json_dict = {}
    for i in range(0, len(input_bed)):
        output_lineL = list_input[i][:-1]  # original input line
        json_dict[str(i)] = {}
        json_dict[str(i)]["query_interval"] = output_lineL
        formatted, best_hit = regionanalysis.getBestHit(
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
