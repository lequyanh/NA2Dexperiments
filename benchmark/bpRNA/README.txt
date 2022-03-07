bpRNA_structures_all_pk.fst
	- structure FASTA created by assembling dot-bracket and FASTA files of individual records from bpRNA site
	- the records were retrieved via a default search; !!! not from the download page - IDs don't match !!!
	
bpRNA_1m_90_structureTypes.csv
	- mapping between RNA original accession (RFAM/CRW accessions, not bpRNA IDs) and their RNA type
	- retrieved from the local installation if bpRNA database. The installation scripts accessible from:
		https://github.com/hendrixlab/bpRNA

##########################################
# Extract structure types of RNA records #
##########################################
PROBLEM1:	MYSQL DB contains RNA types for each record but does not contain RNA dot-bracet structures 
			-> must map the RNA dot-bracket files to the RNA types from DB
PROBLEM2:	bpRNA ID mismatches between MSQL DB version and WEB version - cannot map based on bpRNA ID
			-> The only key to map structures to their types is via their original (e.g. RFAM/CRW) accession
PROBLEM3:	CRW, tmRNA FASTA accessions contain unnecessary info
			-> parse the FASTAS/Dot-brackets and rewrite the headers
PROBLEM4:	CRW records can't be mapped even through original accessions
			bpRNA_CRW_20648
			- DB accession	 	Mollicutes/M.mycoid10 
			- FASTA accession	U26049
		- but it seem they match in bpRNAID => for mapping CRW records to their types we will use bpRNA id

bpRNA_1m_90_structureTypes.csv
	- contains mapping between RNA original accession (not bpRNA) and their RNA type
		"RF00103_AAWU01006939.1_25735-25811"; "mir-1"
		"RF00103_AAWZ02009612.1_57865-57941"; "mir-1"
		"RF00103_AAWZ02017015.1_31896-31821"; "mir-1"

	- created by:
		sudo service mysql start

        cd ..../Load_bpRNA_1m_Mysql_Part1/
        # Fixing the insert script
        # replace double apostrophes with a single one
		sed "s/''/'/g" 5_RNATable_Load.sql > 5_RNATable_Load_escaped.sql

        # replace non-terminating semicolons with colons
        sed -i 's/;\([^$]\)/,\1/g' 5_RNATable_Load_escaped.sql

        # remove rogue apostrophes (those not preceded by a comma and followed by a comma/closing bracket)
        # matching end of SQL expression or a new term in a INSERT statement
        sed -i "s/\([^,]\)'\([^,)]\)/\1\2/g" 5_RNATable_Load_escaped.sql

        # removing 5 occurrences of the "3')_element"
        sed -i "s/3')_/3prime)/g" 5_RNATable_Load_escaped.sql

        # will result in 101470/102318 inserts successful

		mysql -u root -p
		mysql>create database bpRNA;
		mysql>use bpRNA

		mysql>source Create_bpRNA_1m/bpRNA_1m_CreateTablesScript.sql;
		mysql>source Load_bpRNA_1m_Mysql_Part1/1_MethodTable_Load.sql;
		mysql>source Load_bpRNA_1m_Mysql_Part1/2_TypeTable_Load.sql;
		mysql>source Load_bpRNA_1m_Mysql_Part1/3_RefDatabaseTable_Load.sql
		mysql>source Load_bpRNA_1m_Mysql_Part1/4_RefDatabaseTable_Load.sql
		mysql>source Load_bpRNA_1m_Mysql_Part1/5_RNATable_Load_escaped.sql

		mysql>select rna_Desc, type_Name from RNA natural join RNAType into outfile '/var/lib/mysql-files/bpRNA_1m_90_structureTypes_raw.csv' fields terminated by ';';
		mysql>exit

		bash>sudo mv /var/lib/mysql-files/bpRNA_1m_90_structureTypes_raw.csv ./

	- The accession contain junk info. Alter the CSV accessions (extract the field "Original name")
		ORIGINAL: "New Name: bpRNA_RFAM_3, Original Name: RF00001_AB015590.1_1-119 Desc:0"
		NEW: 	  "RF00001_AB015590.1_1-119"
	- created by:
		>cat ./misc/bpRNA_1m_90_structureTypes_raw.csv | grep -Po "Original Name: \K[^\s;,(]*"  > bpRNA_1m_90_structureTypes_accession.csv
		>cat ./misc/bpRNA_1m_90_structureTypes_raw.csv | grep -Po "New Name: \K[^\s;,(]*"  > bpRNA_1m_90_structureTypes_bpRNAID.csv

	- Join the temporary tables into a CSV that maps each RNA to its type
       		> a <- read.csv('bpRNA_1m_90_structureTypes_raw.csv', sep = ';')
        	> accessions <- read.csv('bpRNA_1m_90_structureTypes_accession.csv', sep = ';')
        	> bpRNAids <- read.csv('bpRNA_1m_90_structureTypes_bpRNAID.csv', sep = ';')
        	> out <- data.frame(accession=accessions$Actinobacteria.AB002635, bprna_id=bpRNAids$bpRNA_CRW_1, rna_type=a$X16S)
        	> write.table(out, file='bpRNA_1m_90_structureTypes.csv', row.names = F, sep=';')

#########################################################
# Alter FASTA acessions to match those in structure CSV #
#########################################################
- Run unify_accessions.sh
	1. reads FASTA header of each record to obtain original accession (from bpRNA_fastas/ folder)
	2. parses the header if it contains junk info (CRW and tmRNA)
	3. pastes the cleaned header to the corresponding dot-bracket file (from bpRNA_dbns_new/ folder)
ADD: add the filename to the header
 	for f in bpRNA_dbns_new/*dbn; do file=$(basename $f); bpid=${file%.*}; sed -i "1s/\(.*\)/\1|$bpid/" $f; done

- read all CRW and tmRNA FASTAS that have junk header (read from bpRNA_fasta directory)
- replace the headers with accession only, output to a single composite FASTA
		>cat bpRNA_fasta/*CRW* |  sed 's/.*Accession ID:\(.*\)|.*/>\1/' > bpRNA_CRW_fasta.fa
			ORIGINAL: ">bpRNA_CRW_10000|Accession ID:AY189586|Original Source:Gutell Lab CRW (http://www.rna.ccbb.utexas.edu/DAT/3C/SBPI/index.php)"
			NEW:      ">AY189586"
		>cat bpRNA_fasta/*tmRNA* | sed 's/>\(.*\).*/\n>\1/' | sed  's/.*\/\(.*\).*/>\1/' > bpRNA_tmRNA_fasta.fa
			ORIGINAL: ">../BPSEQFiles/Chla.trac._AE001276_1-417"
			NEW:      ">Chla.trac._AE001276_1-417"