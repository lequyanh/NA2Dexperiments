=====================================================
| Benchmarking data derived from Zammit et.al. data |
=====================================================
Original article 
"A database of flavivirus RNA structures with a search algorithm for pseudoknots and triple base interactions"

available at:
https://academic.oup.com/bioinformatics/article/37/7/956/5899721?login=true

=========
| Files |
=========
Flavivi_3UTRs_raw.txt
	* 3'UTR regions of all viruses examined by Zammit et. all
	* retrieved from https://rna.liacs.nl

Flavivi_3UTRs.fa
	* Parsed raw Flavivirus 3'UTR data from the webserver into FASTA

zammit_supplementary.pdf
	* Copy of the supplementary tables from original the article
	* Used to derive exact Xrn1-resistant domains from 3'URS

S1_csv/zammit-cISFV.csv
S1_csv/zammit-MBFV.csv
S1_csv/zammit-MBFV-other.csv
S1_csv/zammit-TBFV.csv
	* Decomposed zammit_supplementary.pdf  (tables with locations of SL/DB elements in each virus)
	* CSVs derived from the S1 supplementary table
	* each file describes locations of a particular group of Xrn1-resistant domains

S1_csv/cISV_MBFV_SL1.fa
S1_csv/MBFV_DB2_PK2.fa
S1_csv/TBFV_Y_SL1.fa
S1_csv/TBFV_3_GC_SL.fa
...
	* concrete sequences of the domains derived from the zammit-S1-CSVs above. 
	* the result of parse_zammit_tables.Rmd script


S2_csv/S2-TBFV-other.csv
S2_csv/S2-MBFV.csv
S2_csv/S2-MBFV-other.csv
...
	* Parsed S2 tables from Zammit supplementary.
	* Used regular expressions to create the csv:
		(?<=\D*) (?=\d{2,9})
		(?<=\D*) (?=-\d{1,3}\.\d)
		(?<=-\d*) (?=\d{1,2})
		(?=[A-Z]{1,2}_?\d{5,8})
		(?<=-\d*\.\d)
		(?<=-\d*\.\d)\n

===========
| results |
===========
resuls directory contains 3 folders:
rna.liacs.nl_results/
	* results from the rna.liacs.nl server
	* File for each query used (domain searched)

S1_seqs_NA2D_results/
	* Results from NA2DSearch when exact domain sequences were searched
	* File for each query used (domain searched)

sliding_NA2D_results/
	* Results of "sliding window mode" searching of 3'UTRs
	* File for each query used (domain searched)

