deg_clean.ipynb
	This notebook aims to merge the DEG lists into a single list. 
	Note that the file FinalGeneList_1.5_fold_NoTrkA contains 69 genes as it contains NTK1 - which was removed. 
	Notably, merging by external_gene_name produces 163 but by ensemble_gene_id_version produces 176.
	This is because 13 external_genes_names have two ensemble id's, which need to be further investigated in ensemble
	to check that they do corrospond to the same gene:

	FNBP1 ['ENSG00000187239.17', 'ENSG00000187239'] 
	ZNF423 ['ENSG00000102935.11', 'ENSG00000102935'] 
	SIX3 ['ENSG00000138083.5', 'ENSG00000138083'] 
	PCDHB14 ['ENSG00000120327.6', 'ENSG00000120327'] 
	TNNT1 ['ENSG00000105048.17', 'ENSG00000105048'] 
	BCAN ['ENSG00000132692.19', 'ENSG00000132692'] 
	DLG2 ['ENSG00000150672.17', 'ENSG00000150672'] 
	TRIM5 ['ENSG00000132256.19', 'ENSG00000132256'] 
	DENND1C ['ENSG00000205744.10', 'ENSG00000205744'] 
	PCDHB2 ['ENSG00000112852.7', 'ENSG00000112852'] 
	RGS10 ['ENSG00000148908', 'ENSG00000148908.15'] 
	SH2D3C ['ENSG00000095370', 'ENSG00000095370.20'] 
	FRAS1 ['ENSG00000138759', 'ENSG00000138759.19'] 

deg_gene_name.csv
	Contains the 169 external_gene_names as used in Arthurs report (the 4 CoV2 genes were appended). 