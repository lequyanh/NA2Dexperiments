source('script_functions.R')
library(data.table)
library(tidyr)

#####################
# Support functions #
#####################
s2_load <- function(csv_path){
  df <- read.csv(csv_path, sep = ';', na.strings = c("NA", "")) %>% 
    fill(everything(vars=c('virus', 'accession', 'ORF3prime_end')), 
         .direction = 'down')
  
  # Add relative structure starts
  df <- df %>% 
    mutate(s2_rel_start = s_start - ORF3prime_end)
}

round2 <- function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5 + sqrt(.Machine$double.eps)
  z = trunc(z)
  z = z/10^n
  z*posneg
}

load_tables <- function(){
  s2_mbfv <<- s2_load('benchmark/Zammit/S2_csv/S2-MBFV.csv')
  s2_mbfv_nkv_isfv <<- s2_load('benchmark/Zammit/S2_csv/S2-MBFV-other.csv')
  s2_tbfv <<- s2_load('benchmark/Zammit/S2_csv/S2-TBFV.csv')
  s2_tbfv_nkv <<- s2_load('benchmark/Zammit/S2_csv/S2-TBFV-other.csv')
  s2_cisfv <- s2_load('benchmark/Zammit/S2_csv/S2-cISFV.csv')
  
  # Fix St. Louis virus coordinates
  s2_mbfv <- s2_mbfv %>% 
    mutate(s2_rel_start = ifelse(accession == 'NC_007580', s2_rel_start + 187, s2_rel_start))
  
  # Fix Powassan virus coordinates
  s2_tbfv <- s2_tbfv %>%
    mutate(s2_rel_start = ifelse(accession == 'NC_003687', s2_rel_start + 230, s2_rel_start))
  
  # Derive a table for each structure
  mbfv_sl_viruses <- full_join(s2_mbfv, s2_mbfv_nkv_isfv)
  mbfv_sl_viruses <<- full_join(mbfv_sl_viruses, s2_cisfv) %>%
    filter(structure_type == 'MBFV SL')
  
  mbfv_db1_nopk_viruses <<- full_join(s2_mbfv, s2_tbfv_nkv) %>%
    filter(structure_type == 'DB1(without PK)')
  
  mbfv_db2pk_viruses <<- full_join(s2_mbfv, s2_mbfv_nkv_isfv) %>%
    filter(structure_type == 'DB2+PK')
  
  y_sl_viruses <<- full_join(s2_tbfv, s2_tbfv_nkv) %>%
    filter(structure_type == 'TBFV Y-SL')
  
  au_sl_viruses <<- full_join(s2_tbfv, s2_tbfv_nkv) %>%
    filter(structure_type == 'TBFV AU-SL')
  
  tbfv_3gc_sl_viruses <<- s2_tbfv %>%
    filter(structure_type == 'TBFV 3GC-SL')
  
  tbfv_5gc_sl_viruses <<- s2_tbfv %>%
    filter(structure_type == 'TBFV 5GC-SL')
}

combineS2withWebapp <- function(s2_table, webapp_results){
  parsed_leiden <- parseLeidenWebResults(webapp_results)
  
  # Round the deltaG from Leiden webapp to match the deltaG format in the S2 table
  parsed_leiden <- parsed_leiden %>% 
    mutate(deltaG = round2(deltaG, n = 1))
  
  # Merge S2 table and Webapp results on 'accession' and 'hit_location'
  # Prepare for rolling join as the hit locations can be shifted by 1-2nt from reference locations (which prevents basic join)
  setDT(s2_table)
  setDT(parsed_leiden)
  
  s2_table[, loc := s2_rel_start]
  parsed_leiden[, loc := hit_loc]
  
  setkey(s2_table, accession, loc)
  setkey(parsed_leiden, accession, loc)
  
  merged_s2_web <- s2_table[parsed_leiden, roll='nearest', rollends=c(TRUE, TRUE)]
  
  # Invalidate joins, where hit location and reference location are too far from each other
  merged_s2_web <- merged_s2_web %>%
    mutate(s2_rel_start = ifelse(abs(s2_rel_start - loc) > 10, NA, s2_rel_start),
           s2_structure_rank = ifelse(abs(s2_rel_start - loc) > 10, NA, s2_structure_rank))
  
  # Re-arrange and rename columns
  merged_s2_web <- merged_s2_web %>%
    relocate(seq_def, .after = accession) %>%
    relocate(structure_type, .after = structure) %>%
    relocate(hit_count, .after = structure_type) %>%
    relocate(s2_ref_deltaG, .before = deltaG) %>%
    select(-loc, -matches('s_'), -ORF3prime_end)
  
  return(list(merged_s2_web = merged_s2_web,
              parsed_leiden = parsed_leiden))
}

# MBFV SL query
getMBFV_SLresults <- function(){
  merged_s2_web <- combineS2withWebapp(
    s2_table = mbfv_sl_viruses,
    webapp_results = 'benchmark/Zammit/results/rna.liacs.nl_results/Leiden_Flavivi3UTRs_MBFV_SL_query.txt'
  )
  
  return(c(merged_s2_web, 
           list(s2_table = mbfv_sl_viruses)))  # return the reference table as well
}

# MBFV DB1 query
getMBFV_DB1results <- function(){
  merged_s2_web <- combineS2withWebapp(
    s2_table = mbfv_db1_nopk_viruses,
    webapp_results = 'benchmark/Zammit/results/rna.liacs.nl_results/Leiden_Flavivi3UTRs_MBFV_DB1_noPK_query.txt'
  )
  
  return(c(merged_s2_web, 
           list(s2_table = mbfv_db1_nopk_viruses)))  # return the reference table as well
}

# MBFV DB2 (YFG) query
getMBFV_DB2_YFGresults <- function(){
  merged_s2_web <- combineS2withWebapp(
    s2_table = mbfv_db2pk_viruses,
    webapp_results = 'benchmark/Zammit/results/rna.liacs.nl_results/Leiden_Flavivi3UTRs_MBFV_DB2_PK_YFG_query.txt'
  )
  
  return(c(merged_s2_web, 
           list(s2_table = mbfv_db2pk_viruses)))  # return the reference table as well
}

# MBFV DB2 (non-YFG) query
getMBFV_DB2_nonYFGresults <- function(){
  merged_s2_web <- combineS2withWebapp(
    s2_table = mbfv_db2pk_viruses,
    webapp_results = 'benchmark/Zammit/results/rna.liacs.nl_results/Leiden_Flavivi3UTRs_MBFV_DB2_PK_non-YFG_query.txt'
  )
  
  return(c(merged_s2_web, 
           s2_table = list(mbfv_db2pk_viruses)))  # return the reference table as well
}

# TBFV Y-SL query
getTBFV_Y_SLresults <- function(){
  merged_s2_web <- combineS2withWebapp(
    s2_table = y_sl_viruses,
    webapp_results = 'benchmark/Zammit/results/rna.liacs.nl_results/Leiden_Flavivi3UTRs_TBFV_Y_SL_query.txt'
  )
  
  return(c(merged_s2_web, 
           s2_table = list(y_sl_viruses)))  # return the reference table as well
}

# TBFV AU-SL query
getTBFV_AU_SLresults <- function(){
  merged_s2_web <- combineS2withWebapp(
    s2_table = au_sl_viruses,
    webapp_results = 'benchmark/Zammit/results/rna.liacs.nl_results/Leiden_Flavivi3UTRs_TBFV_AU_SL_query.txt'
  )
  
  return(c(merged_s2_web, 
           s2_table = list(au_sl_viruses)))  # return the reference table as well
}

# TBFV 3'GC-SL query
getTBFV_3GC_SLresults <- function(){
  merged_s2_web <- combineS2withWebapp(
    s2_table = tbfv_3gc_sl_viruses,
    webapp_results = 'benchmark/Zammit/results/rna.liacs.nl_results/Leiden_Flavivi3UTRs_TBFV_3GC_SL_query.txt'
  )
  
  return(c(merged_s2_web, 
           s2_table = list(tbfv_3gc_sl_viruses)))  # return the reference table as well
}

# TBFV 5'GC-SL query
getTBFV_5GC_SLresults <- function(){
  merged_s2_web <- combineS2withWebapp(
    s2_table = tbfv_5gc_sl_viruses,
    webapp_results = 'benchmark/Zammit/results/rna.liacs.nl_results/Leiden_Flavivi3UTRs_TBFV_5GC_SL_query.txt'
  )
  
  return(c(merged_s2_web, 
           s2_table = list(tbfv_5gc_sl_viruses)))  # return the reference table as well
}

getCombined_S2_webbapp <- function(structure_type){
  if (structure_type == 'TBFV-Y-SL'){
    results <- getTBFV_Y_SLresults()
  } else if (structure_type == 'TBFV-AU-SL'){
    results <- getTBFV_AU_SLresults()
  } else if (structure_type == 'TBFV-3GC-SL'){
    results <- getTBFV_3GC_SLresults()
  } else if (structure_type == 'TBFV-5GC-SL'){
    results <- getTBFV_5GC_SLresults()
  } else if (structure_type == 'MBFV-SL'){
    results <- getMBFV_SLresults()  
  } else if (structure_type == 'MBFV-DB1'){
    results <- getMBFV_DB1results()
  } else if (structure_type == 'MBFV-DB2-YFG'){
    results <- getMBFV_DB2_YFGresults()
  } else if (structure_type == 'MBFV-DB2-non-YFG'){
    results <- getMBFV_DB2_nonYFGresults()
  } else {
    print("Not supported structure")
  }
  
  return(results)
}

simplifyWussStructures <- function(wuss_structures){
  simpl_structs <- sapply(wuss_structures, 
            function(x){
              str_replace_all(x, '[\\]\\-_,:\\}\\{\\[]', '.') %>%
                str_replace_all('<', '(' ) %>%
                str_replace_all('>', ')') %>%
                # Replace the ~ pair indicating tertiary interactions by an ordinary BP
                str_replace('~(?=.*~)', '(') %>%
                str_replace('~', ')')
            }, USE.NAMES = F)
  
  return(simpl_structs)
}

simplifyWussStructures2 <- function(wuss_structures){
  simpl_structs <- sapply(wuss_structures, 
                          function(x){
                            str_replace_all(x, '[\\]\\-_,:\\[]', '.') %>%
                              str_replace_all('[\\(\\)]', '.') %>% # in TBFV case () brackes are used for pseudoknots
                              str_replace_all('[<\\{]', '(' ) %>%  # in TBFV case {} brackes are used base pairs
                              str_replace_all('[>\\}]', ')') %>%
                              # Replace the ~ pair indicating tertiary interactions by an ordinary BP
                              str_replace('~(?=.*~)', '(') %>%
                              str_replace('~', ')')
                          }, USE.NAMES = F)
  
  return(simpl_structs)
}

combineS2WebappNA2D <- function(jointWebappS2Tables, na2d_file, ignore_pk=F){
  # Parse the NA2Dresults to access structures and energies
  parsedNA2Dresults <- parseN2dZammitResults(na2d_file)
  
  if(jointWebappS2Tables[1, 'structure_type'] == 'TBFV Y-SL'){
    # Parse a different format of WUSS
    # {{{{{:.:(((<<<<___>>>>,<~<<_____[[[[___>>~>,}}}}},,,,,,,,,,,,,,,,,,)))::]]]]
    jointWebappS2Tables$simpl_struct <- simplifyWussStructures2(jointWebappS2Tables$structure)
    
  } else {
    # Simplify the structures from webapp so they are comparable to NA2DSearch simple dot-bracket format
    jointWebappS2Tables$simpl_struct <- simplifyWussStructures(jointWebappS2Tables$structure)
  }
  
  if(!ignore_pk){
    # Match NA2D hits with full reference structures, including PK (which has been reduced to a simple single strand)
    # Join by accession and structure to find which NA2D hits match the reference structures
    leidenNA2Djoin <- full_join(jointWebappS2Tables, 
                                parsedNA2Dresults, 
                                by = c("accession", 
                                       "seq_def" = "description",
                                       "simpl_struct" = "na2d_structure"))
  } else {
    # Remove pseudoknoted ends of the reference leiden structures
    jointWebappS2Tables <- jointWebappS2Tables %>% 
      mutate(simpl_struct = str_replace(simpl_struct, '[\\.]{4,20}$', ''))
    
    # Perform simple join by accession (can't by structure, as the S2 structures are longer)
    leidenNA2Djoin <- full_join(jointWebappS2Tables, 
                                parsedNA2Dresults, 
                                by = c("accession", 
                                       "seq_def" = "description"))
    
    # Find in which join records the (simplified) S2 structure is a substring/substructure of the NA2D structure
    mask <- mapply(grepl, 
                   leidenNA2Djoin$simpl_struct, 
                   leidenNA2Djoin$na2d_structure, 
                   fixed=TRUE, USE.NAMES = F)
    
    # When NA2D did not fnd anything, the grep will return NA. Replace NA with TRUE (keep these records)
    mask <- replace(mask, is.na(mask), TRUE)
    
    # Invalidate rows, where the simplified Leiden structure doesn't match the NA2D structure
    leidenNA2Djoin[!mask,] <- leidenNA2Djoin[!mask,] %>% 
      mutate(across(c(s2_structure_rank, hit_loc, simpl_struct, deltaG), ~ NA))
  }
  
  leidenNA2Djoin <- leidenNA2Djoin %>%
    select(-seq_def, -Type) %>%
    relocate(structure, .before = simpl_struct) %>%
    relocate(na2d_energy, .after = s2_ref_deltaG)
  
  return(leidenNA2Djoin)
} # ^combineS2WebappNA2D()

combineS2WebappNA2D_sliding <- function(jointWebappS2Tables, na2d_file){
  if(jointWebappS2Tables[1, 'structure_type'] == 'TBFV Y-SL'){
    # Parse different format of WUSS
    # {{{{{:.:(((<<<<___>>>>,<~<<_____[[[[___>>~>,}}}}},,,,,,,,,,,,,,,,,,)))::]]]]
    jointWebappS2Tables$simpl_struct <- simplifyWussStructures2(jointWebappS2Tables$structure)
  } else {
    # Simplify the structures from webapp so they are comparable to NA2DSearch simple dot-bracket format
    jointWebappS2Tables$simpl_struct <- simplifyWussStructures(jointWebappS2Tables$structure)
  }
  
  # Remove pseudoknoted ends of the leiden structure
  jointWebappS2Tables <- jointWebappS2Tables %>% 
    mutate(simpl_struct = str_replace(simpl_struct, '[\\.]{4,30}$', '')) %>%
    mutate(simpl_struct = str_replace(simpl_struct, '^[\\.]{1,3}', ''))
  
  # Parse the NA2Dresults to access structures and energies
  parsedNA2Dresults <- parseN2dZammitResults(na2d_file)
  
  # Rolling join
  setDT(jointWebappS2Tables)
  setDT(parsedNA2Dresults)
  
  jointWebappS2Tables[, loc := s2_rel_start]
  parsedNA2Dresults[, loc := bnStart]
  
  setkey(jointWebappS2Tables, accession, loc)
  setkey(parsedNA2Dresults, accession, loc)
  
  leidenNA2Djoin <- jointWebappS2Tables[parsedNA2Dresults, roll='nearest', rollends=c(TRUE, TRUE)]
  
  # Invalidate joins, where hit location and reference location are too far from each other
  # Invalidate by setting simpl_struct to NA, hence disabling further joins by structure
  leidenNA2Djoin <- leidenNA2Djoin %>%
    mutate(simpl_struct = ifelse(abs(bnStart - s2_rel_start) > 19, NA, simpl_struct))
  
  # Find in which join records the (simplified) S2 structure is a substring/substructure of the NA2D structure
  mask <- mapply(grepl, 
                 leidenNA2Djoin$simpl_struct, 
                 leidenNA2Djoin$na2d_structure, 
                 fixed=TRUE, USE.NAMES = F)
  
  # When NA2D did not fnd anything, the grep will return NA. Replace NA with TRUE (keep these records)
  mask <- replace(mask, is.na(mask), TRUE)
  
  # Invalidate rows, where the simplified Leiden structure doesn't match the NA2D structure
  leidenNA2Djoin[!mask,] <- leidenNA2Djoin[!mask,] %>% 
    mutate(across(c(s2_structure_rank, hit_loc, structure, simpl_struct, deltaG), ~ NA))
  
  
  leidenNA2Djoin <- leidenNA2Djoin %>%
    relocate(structure, .before = simpl_struct) %>%
    relocate(na2d_energy, .after = s2_ref_deltaG)
  
  return(leidenNA2Djoin)
} # ^combineS2WebappNA2D_sliding()

compareNa2d_S2 <- function(structure_type, ignore_leiden_pk = F, na2d_file_override=NA){
  if (structure_type == 'TBFV-Y-SL'){
    results <- getTBFV_Y_SLresults()
    na2d_file <- 'benchmark/Zammit/results/S1_seqs_NA2D_results/results_TBFV_Y_SL.fa'
  } else if (structure_type == 'TBFV-AU-SL'){
    results <- getTBFV_AU_SLresults()
    na2d_file <- 'benchmark/Zammit/results/S1_seqs_NA2D_results/results_TBFV_NKV_AU_SL.fa'
  } else if (structure_type == 'TBFV-3GC-SL'){
    results <- getTBFV_3GC_SLresults()
    na2d_file <- 'benchmark/Zammit/results/S1_seqs_NA2D_results/results_TBFV_3_GC_SL.fa'
  } else if (structure_type == 'TBFV-5GC-SL'){
    print("not comparable with NA2D")
  } else if (structure_type == 'MBFV-SL'){
    results <- getMBFV_SLresults()  
    na2d_file <- 'benchmark/Zammit/results/S1_seqs_NA2D_results/results_MBFV_SL.fa'
  } else if (structure_type == 'MBFV-DB1'){
    results <- getMBFV_DB1results()
    na2d_file <- 'benchmark/Zammit/results/S1_seqs_NA2D_results/results_MBFV_DB1_noPK.fa'
  } else if (structure_type == 'MBFV-DB2-YFG'){
    results <- getMBFV_DB2_YFGresults()
    na2d_file <- 'benchmark/Zammit/results/S1_seqs_NA2D_results/results_MBFV_DB2_YFG.fa'
  } else if (structure_type == 'MBFV-DB2-non-YFG'){
    results <- getMBFV_DB2_nonYFGresults()
    na2d_file <- 'benchmark/Zammit/results/S1_seqs_NA2D_results/results_MBFV_DB2_nonYFG.fa'
  }
  
  # DEBUG: Can override the NA2D results file (with temporary resuts.fa file)
  if(!is.na(na2d_file_override)){
    na2d_file <- na2d_file_override 
  }
  
  # results have 3 slots
  merged_table <- results$merged_s2_web   # inner join of the S2 table and webapp results (for DEBUG purposes)
  webapp_df <- results$parsed_leiden      # parsed webapp results
  s2_table <- results$s2_table            # reference S2 table (reference structure ranks and locations)
  
  # S2 table enriched with correct structures
  s2_with_structures <- merged_table %>% 
    group_by(accession, s2_rel_start) %>%
    group_modify(~ {.x %>% 
        arrange(deltaG) %>% 
        dplyr::slice(dplyr::first(.x$s2_structure_rank))}
    )
  
  # Join S2 table enriched with structures with NA2D results (structural matches to the descriptors)
  comparisonDf <- combineS2WebappNA2D(s2_with_structures, na2d_file, ignore_leiden_pk) %>%
    select(-matches('bn'), -s2_ref_deltaG, -hit_count)

  return(comparisonDf)
} # ^compareNa2d_S2()

compareSlidingNa2d_S2 <- function(structure_type, na2d_file_override=NA){
  if (structure_type == 'TBFV-Y-SL'){
    results <- getTBFV_Y_SLresults()
    na2d_file <- 'benchmark/Zammit/results/sliding_NA2D_results/results_Y_SL_Flavivi3utrs_sliding_12step_20over_50win_6.3ddg.fa'
  } else if (structure_type == 'TBFV-AU-SL'){
    results <- getTBFV_AU_SLresults()
    na2d_file <- 'benchmark/Zammit/results/sliding_NA2D_results/results_AU_SL_Flavivi3utrs_sliding_8step_15over_57win_1.4ddg.fa'
  } else if (structure_type == 'TBFV-3GC-SL'){
    results <- getTBFV_3GC_SLresults()
    na2d_file <- 'benchmark/Zammit/results/sliding_NA2D_results/results_3GC_SL_Flavivi3utrs_sliding_12step_25over_63win_5ddg.fa'
  } else if (structure_type == 'TBFV-5GC-SL'){
    print("not comparable with NA2D")
  } else if (structure_type == 'MBFV-SL'){
    results <- getMBFV_SLresults()  
    na2d_file <- 'benchmark/Zammit/results/sliding_NA2D_results/results_MBFV_SL_Flavivi3utrs_sliding_12step_20over_55win_8.0ddg.fa'
  } else if (structure_type == 'MBFV-DB1'){
    results <- getMBFV_DB1results()
    na2d_file <- 'benchmark/Zammit/results/sliding_NA2D_results/results_DB1_SL_Flavivi3utrs_sliding_15step_20over_70win_2.6ddg.fa'
  } else if (structure_type == 'MBFV-DB2-YFG'){
    results <- getMBFV_DB2_YFGresults()
    na2d_file <- 'benchmark/Zammit/S1_csv/na2d_results/results_MBFV_DB2_YFG.fa'
  } else if (structure_type == 'MBFV-DB2-non-YFG'){
    results <- getMBFV_DB2_nonYFGresults()
    na2d_file <- 'benchmark/Zammit/results/sliding_NA2D_results/results_nYFG_DB2_SL_Flavivi3utrs_sliding_12step_20over_80win_4.7ddg.fa'
  }
  
  # DEBUG: Can override the NA2D results file (with temporary resuts.fa file)
  if(!is.na(na2d_file_override)){
    na2d_file <- na2d_file_override 
  }
  
  # results have 3 slots
  merged_table <- results$merged_s2_web   # inner join of the S2 table and webapp results (for DEBUG purposes)
  webapp_df <- results$parsed_leiden      # parsed webapp results
  s2_table <- results$s2_table            # reference S2 table (reference structure ranks and locations)
  
  # S2 table enriched with correct structures
  s2_with_structures <- merged_table %>% 
    group_by(accession, s2_rel_start) %>%
    group_modify(~ {.x %>% 
        arrange(deltaG) %>% 
        dplyr::slice(dplyr::first(.x$s2_structure_rank))}
    )
  
  # Join S2 table enriched with structures with NA2D results (structural matches to the descriptors)
  comparisonDf <- combineS2WebappNA2D_sliding(s2_with_structures, na2d_file) %>%
    select(-Type, -structure, -hit_count, -s2_ref_deltaG, -s2_delta_deltaG, -s2_structure_rank, -deltaG, -loc, -description) %>%
    relocate(hit_loc, .before = bnStart)
  
  # Add missing reference record that got dropped out during rolling join
  # When we roll the reference table on NA2D results, the base table is NA2D results DF. If it misses a viruses, so does to join
  missing <- s2_with_structures[!(s2_with_structures$accession %in% comparisonDf$accession), ]
  comparisonDf <- rbind(comparisonDf, missing, fill=T)
  
  return(comparisonDf)
} # ^compareSlidingNa2d_S2()
