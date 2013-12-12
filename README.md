region_analysis
===============

Usage: region_analysis.py [options]

Options:

  -h, --help            show this help message and exit

  -i INPUT_FILE_NAME, --input=INPUT_FILE_NAME

                        Input region file must assume the first 3 columns

                        contain (chr, start, end)

  -d ANNO_DB, --database=ANNO_DB

                        Choose database: refseq(default) or ensembl

  -r, --rhead           Whether the input file contains column header

  -g GENOME, --genome=GENOME

                        Choose genome: mm10(default)

Output:

  -.annotated: the one-to-one output list, only the annotation entry whose TSS is nearest to the inquiry interval kept.

  -.full.annotated: all hit entries are kept.

  -.full.annotated.json: the json format output of -.full.annotated.
  
  ProximalPromoter: 	+/- 250bp of TSS

  Promoter1k: 	+/- 1kbp of TSS

  Promoter3k: 	+/- 3kbp of TSS

  Genebody: 	Anywhere between a gene's promoter and up to 1kbp downstream of the TES.

  Genedeserts: 	Genomic regions that are depleted with genes and are at least 1Mbp long.

  Pericentromere: 	Between the boundary of a centromere and the closest gene minus 10kbp of that gene's regulatory region.

  Subtelomere: 	Similary defined as pericentromere.

  OtherIntergenic: 	Any region that does not belong to the above categories.
