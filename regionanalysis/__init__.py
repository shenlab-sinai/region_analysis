__version__ = "0.1.1"

def midpoint(cur_input):
    """
    Calculate the midpoint of the genomic interval.
    cur_input: a line of input, the first 3 columns are the intervals in BED format.
    """
    mid_input = cur_input.split("\t")
    midpoint = (int(mid_input[1]) + int(mid_input[2]))/2
    mid_input[1] = str(midpoint)
    mid_input[2] = str(midpoint + 1)
    return "\t".join(mid_input)

def getDis2TSS(anno_db, cur_input, col_no_input):
    """
    Calculate the distance to TSS and decide the annotation feature of the entry.
    """
    # calculate TSS and TES based on strands.
    # if strand is not "-" then it will be treated as "+".
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
            cur_input[col_no_input+3], cur_input[col_no_input+5], cur_input[col_no_input+6],
            cur_input[col_no_input+10], cur_input[col_no_input+11], Pos, str(Dis2TSS),
            cur_input[col_no_input+9]]
    else:  # output is gene symbol
        cur_output = [cur_input[col_no_input+4], cur_input[col_no_input+5], cur_input[col_no_input+6],
            cur_input[col_no_input+10], cur_input[col_no_input+11], Pos, str(Dis2TSS),
            cur_input[col_no_input+9]]
    return (cur_output, Dis2TSS, Dis2TES)

def getBestHit(anno_db, col_no_input, GB_entry, gd_entry, st_entry, pc_entry):
    """
    format all hits and get the best hit that is neartest to TSS.
    """
    best_hit, cur_output, Dis2TSS, Dis2TES = (None, None, None, None)
    formatted = []

    for i in GB_entry:
        # discard null hit of genebody entry.
        if not ((i[col_no_input] == ".") and (i[col_no_input + 1] == "-1")):
            cur_output, Dis2TSS, Dis2TES = getDis2TSS(
                anno_db, i, col_no_input)
            if (best_hit is None) or (abs(int(best_hit[6])) > abs(Dis2TSS)):
                best_hit = cur_output
            formatted.append(cur_output)
    if gd_entry != "0":
        cur_output = ["NA", "NA", ".", "NA", "NA", "Genedesert", "NA", "No_anno"]
        if (best_hit is None):
            best_hit = cur_output
        formatted.append(cur_output)
    if st_entry != "0":
        cur_output = ["NA", "NA", ".", "NA", "NA", "Subtelomere", "NA", "No_anno"]
        if (best_hit is None):
            best_hit = cur_output
        formatted.append(cur_output)
    if pc_entry != "0":
        cur_output = ["NA", "NA", ".", "NA", "NA",
                      "Pericentromere", "NA", "No_anno"]
        if (best_hit is None):
            best_hit = cur_output
        formatted.append(cur_output)
    if best_hit is None:
        cur_output = ["NA", "NA", ".", "NA", "NA",
                      "Otherintergenic", "NA", "No_anno"]
        best_hit = cur_output
        formatted.append(cur_output)
    return (formatted, best_hit)