# Data, NA2DSearch projects, pre-processing and post-processing scripts used for demonstrating NA2DSearch capabilities.

## Data
Data used in the published experiments are stored in the `benchmark` folder (see structure below). 

```
├── benchmark
│   ├── bpRNA
│   ├── tRNAdb
│   ├── Zammit
```

* The sub-directories carry the name of the datasource. 
  * *Zammit* denotes the 3'UTR Flavivirae database published by *Zammit et. al. 2021* (https://pubmed.ncbi.nlm.nih.gov/32866223/)
* Each datasource has been pre-processed differently. The pre-processing steps are listed in README files within the respective directories

## NA2DSearch projects
NA2DSearch application can be downloaded from http://labtools.cz/na2dsearch/ (Java 7 and above is required for successful installation).

NA2DSearch bundles query descriptors and search parameters in a functional unit called a *project*.
Projects ensure transferability as they are ready to use right away, once downloaded.

There are two NA2Dsearch projects - one holding queries for *bpRNA* and *tRNAdb* data-sources (they both contain tRNA records),
and the second holding queries for the Zammit dataset. 
```
├── bpRNABenchmarkProject
├── ZammitBenchmarkProject
```

For convenience, we include example databases in each query to test the search (e.g., a tRNA query has a FASTA of tRNA records pre-attached).

## Replicating the results
### Replicating search in the tRNAdb
To replicate the search in *tRNAdb* records do:
1. within NA2DSearch, load the *bpRNABenchmarkProject*
2. navigate to the appropriate search via 
*bpRNABenchmarkProject > Searches > tRNA - structure*
3. Right-click the *Databases* section and select the FASTA file `benchmark/tRNAdb/tRNAdb_eukaryota.fst`
4. Run the search (Play button)

### Replicating search in the bpRNA
To replicate the search in *bpRNA* records, you basically repeat the same steps with a few alterations:
1. within NA2DSearch, load the *bpRNABenchmarkProject* (if not loaded already)
2. navigate to the search of choice (e.g., the search for SCARNA1 structures) via 
*bpRNABenchmarkProject > Searches > SCARNA1 - structure*
3. Right-click the *Databases* section and select the FASTA file `benchmark/bpRNA/bpRNA_dbns_all_pk.fst`
    * **Note**: As the file is too large, select the option *"Search on disk only"* while loading the database
4. Run the search (Play button)
    
### Replicating search in Flavivirae 3'UTRs sequences (Zammit dataset)
**Note**: This search relies on folded structures and for that, *RNAsubopt* or *RNAshapes* installations are needed.
Make sure to set the path to their binaries in *Options > Global Preferences*

To replicate the search in *Zammit* records:
1. within NA2DSearch, load the *ZammitBenchmarkProject*
2. navigate to the search of choice (e.g., the Tick-Borne Flavivirae, AU-SL structure) via 
*ZammitBenchmarkProject > Searches > TBFV - AU-SL - Zammit descriptor*
3. Right-click the *Databases* section and select the FASTA file `benchmark/Zammit/Flavivi_3UTRs.fa`
    * **Note**: These searches will take considerably longer due to the folding sub-routine. To fasten the search, lower the ΔΔG
    folding parameter
4. Run the search (Play button)
