README.MD
=========

Files:

	For TSS, coordinations are extended to upstream 3K bp. For TES, coordinations are extended to downstream 1K bp. So ext_5' and ext_3' are the extended results based on strand information.

	mm9.ensembl.biotype_region_ext.bed: ensembl annotation.

chrom|ext_5' |ext_3' |ensembl id        |Symbol |strand|transcript_id     |position|Null |biotype   |5'     |3'

col*0|col*1  |col*2  |col*3             |col*4  |col*5 |col*6             |col*7   |col*8|col*9     |col*10 |col*11 

chr1 |3041314|3045814|ENSMUSG00000090025|Gm16088|+     |ENSMUST00000160944|genebody|NA   |pseudogene|3044314|3044814

chr1 |3089097|3093206|ENSMUSG00000064842|U6     |+     |ENSMUST00000082908|genebody|NA   |misc      |3092097|3092206

	mm9.refseq.biotype_region_ext.bed: RefSeq annotation.

chrom|ext_5'   |ext_3'   |NULL |Symbol|strand|transcript_id|position|Null |biotype       |5'       |3'

col*0|col*1    |col*2    |col*3|col*4 |col*5|col*6         |col*7   |col*8|col*9         |col*10   |col*11   

chrX |104591534|104601699|NA   |Itm2a |-    |NM_008409     |genebody|NA   |protein_coding|104592534|104598699

chr7 |53697529 |53898177 |NA   |Sergef|-    |NM_013789     |genebody|NA   |protein_coding|53698529 |53895177
