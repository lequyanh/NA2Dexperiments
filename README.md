# Data, NA2DSearch projects, pre-processing and post-processing scripts used for demonstrating NA2DSearch capabilities.

## Data
Data for published experiments are stored in the *benchmark* folder. 
* Each data-source then has its own sub-directory with its name. 
  * here *Zammit* denotes the 3'UTR Flavivirae database published by Zammit et. al. 2021 (https://pubmed.ncbi.nlm.nih.gov/32866223/)
* Each data-source required its own pre-processing steps, which are described in README file within respective directories

```
├── benchmark
│   ├── bpRNA
│   ├── tRNAdb
│   ├── Zammit
```

## NA2DSearch projects
NA2DSearch application can be downloaded from http://labtools.cz/na2dsearch/

NA2DSearch bundles queries and search parameters in a functional unit called *project*.
Projects ensure transferability as they are ready to use right away, once downloaded.

There are two NA2Dsearch projects - the first holding queries for both *bpRNA* and *tRNAdb* data-sources (they both contain tRNA records),
and the second holding queries for the Zammit dataset. 
```
├── bpRNABenchmarkProject
├── ZammitBenchmarkProject
```

## Replicating the results
### Replicating search in the tRNAdb
To replicate the search in *tRNAdb* records do:
1. within NA2DSearch, load the *bpRNABenchmarkProject*
2. navigate to the appropriate search via 
*bpRNABenchmarkProject > Searches > tRNA - structure*
3. Right-click the *Databases* section and select the file `benchmark/tRNAdb/tRNAdb_eukaryota.fst`
4. Run the search (Play button)

### Replicating search in the bpRNA
To replicate the search in *bpRNA* records, you basically repeat the same steps with a few alterations:
1. within NA2DSearch, load the *bpRNABenchmarkProject* (if not loaded already)
2. navigate to the search of choice (e.g., SCARNA1) via 
*bpRNABenchmarkProject > Searches > SCARNA1 - structure*
3. Right-click the *Databases* section and select the file `benchmark/bpRNA/bpRNA_dbns_all_pk.fst`
    * **Note**: As the file is too large, select the option *"Search on disk only"* while loading the database
4. Run the search (Play button)
    
### Replicating search in Flavivirae 3'UTRs sequences (Zammit dataset)
To replicate the search in *Zammit* records:
1. within NA2DSearch, load the *ZammitBenchmarkProject*
2. navigate to the search of choice (e.g., the Tick-Borne Flavivirae, AU-SL structure) via 
*ZammitBenchmarkProject > Searches > TBFV - AU-SL - Zammit descriptor*
3. Right-click the *Databases* section and select the file `benchmark/Zammit/Flavivi_3UTRs.fa`
    * **Note**: These searches will take considerably longer due to the folding sub-routine. To fasten the search, lower the ΔΔG
    folding parameter
4. Run the search (Play button)

## Additional experiments
* To search/filter records of a particular RNA type (e.g., only SCARNA1 records), chose from FASTA files in the `benchmark/bpRNA/misc` directory
  (e.g., mir_10.fst, SCARNA_1.fst)
* The `misc` folder also contains FASTA files with sequence data only (*.fa). These can be supplied to sequence searches 
  (e.g., the one titled *SCARNA_1 - sequence*) 