# library(genbankr)
library(Biostrings)
library(dplyr)
library(stringr)
library(readtext)


parseLeidenWebResults <- function(web_results, s1_accessions=NULL){
  # Parse Zammit file
  res <- readtext(web_results)
  res <- gsub('\n', ' ', res$text)  # Replace newlines with spaces as str_extract has troubles with them
  
  records <- str_split(res, 'Sequence accession')[[1]]
  records <- records[2:length(records)]
  
  accessions <- sapply(records, 
                       function(x){str_extract(x, '[A-Z]{1,2}_?\\d{5,6}')},
                       USE.NAMES = F)
  
  hit_counts <- sapply(records, 
                       function(x){str_extract(x, '(?<=Hit count\t)\\d+')},
                       USE.NAMES = F) %>% as.numeric()
  
  seq_defs <- sapply(records, 
                     function(x){str_extract(x, '(?<=Sequence definition\t).*(?= \\#\tC)')},
                     USE.NAMES = F)
  
  locations <- sapply(records, 
                      function(x){
                        locs <- str_extract_all(x, '(?<=[0-9]{1,3}\t)[0-9]{1,3}(?=[^\\d.]+)')
                        # replace character(0) fields with NaNs. 
                        # We need equal list lengths and R ignores character(0) fields resulting in shorter lists
                        locs <- ifelse(identical(locs[[1]], character(0)), NA, locs)  
                        return(locs)
                      }, USE.NAMES = F)
  
  energies <- sapply(records, 
                     function(x){
                       energs <- str_extract_all(x, '(?<=\t)-?[0-9]{1,2}[.][0-9]+')
                       energs <- ifelse(identical(energs[[1]], character(0)), NA, energs)
                       return(energs)
                     }, USE.NAMES = F)
  
  structures <- sapply(records, 
                       function(x){
                         structs <- str_extract_all(x, '(?<=\t)[~.:,_\\-<>\\{\\}\\(\\)\\[\\]]{10,200}')
                         structs <- ifelse(identical(structs[[1]], character(0)), NA, structs)
                         return(structs)
                       }, USE.NAMES = F)
  
  times <- sapply(structures, 
                  function(x) ifelse(length(x) == 0, 1, length(x)))
  
  parsed_leiden <- data.frame(accession = rep(accessions, times=times),
                              seq_def = rep(seq_defs, times=times),
                              hit_count = rep(hit_counts, times=times),
                              hit_loc = as.numeric(unlist(locations)),
                              deltaG = as.numeric(unlist(energies)),
                              structure = unlist(structures))
  
  # # Add a flag of whether the virus is in any S1 table
  # parsed_leiden <- parsed_leiden %>% 
  #   mutate(inS1Table = sapply(parsed_leiden$accession, function(x){x %in% s1_accessions}))
  
  return(parsed_leiden)
}

#' Turn NA2D result file into a DF. Parse the description of hits to obtain
#' Zammit data (if the virus is TBFV/cISFV/MBFV, etc.)
#' 
#' @param na2d_file location of NA2DSearch result file
parseN2dZammitResults <- function(na2d_file, view_table=F){
  results_df <- parseNa2dResult(na2d_file)
  
  # Split the accession column into Refseq/Genbank ID column and description column
  # accession   Type     description                            bnStart  hits max_score
  # 1 NC_005039 NKV-MBFV Yokose virus, complete genome          15     5       -14
  # 2 NC_002640 MBFV     Dengue virus 4, complete genome        15     4        -2
  results_df <- results_df %>% 
    tidyr::separate(accession, 
                    sep=' ', 
                    c('accession', 'description'), 
                    extra="merge") %>%
    tidyr::separate(description,
                    sep='\\|',
                    c( 'Type', 'description')) %>%
    mutate(accession = str_replace(accession, '\\.\\d', ''))
  
  if(view_table){
    View(results_df)
  }
  
  # Fix Powasan and St. Louis coordinates
  results_df <- results_df %>%
    mutate(bnStart = ifelse(accession == 'NC_003687', bnStart + 230, bnStart),
           bnEnd = ifelse(accession == 'NC_003687', bnEnd + 230, bnEnd)) %>%
    mutate(bnStart = ifelse(accession == 'NC_007580', bnStart + 187, bnStart),
           bnEnd = ifelse(accession == 'NC_007580', bnEnd + 187, bnEnd)) 
  
  return(results_df)
}

parseNa2dResult <- function(na2d_file, keep_sequences=F){
  na2d_res <- readLines(na2d_file)
  
  # Access information on hits (every third line, start line 1)
  hit_info <- na2d_res[seq(1,length(na2d_res), 3)]
  # Extract structures (every third line; start line 3)
  sequences <- na2d_res[seq(2,length(na2d_res), 3)]
  # Extract structures (every third line; start line 3)
  structures <- na2d_res[seq(3,length(na2d_res), 3)]
  # Get the length of hits (second line; all of them should be equally long)
  match_len <-  na2d_res[2] %>% str_length()
  
  accessions <- sapply(hit_info, 
                       function(x){
                         acc <- str_extract(x, ".*\\|(?=|)")
                         acc <- substr(acc, 2, nchar(acc) - 1)
                         return(acc)
                       }, USE.NAMES = F)
  
  na2d_scores <- sapply(hit_info, 
                        function(x){
                          score <- str_extract(x, '(?<=\\"score\\":\\")-\\d{1,3}.\\d+')
                          return(as.numeric(score))
                        }, USE.NAMES = F)
  
  hit_starts <- sapply(hit_info, 
                       function(x){
                         bnStart <- str_extract(x, '(?<=\\"bnStart\\":\\")\\d*')
                         return(as.numeric(bnStart))
                       }, USE.NAMES = F)
  
  energies <- sapply(hit_info, 
                     function(x){
                       energs <- str_extract(x, '(?<=\\"enAbs\\":\\")-\\d*.\\d*')
                       return(as.numeric(energs))
                     }, USE.NAMES = F)
  
  result_df <- data.frame(
    accession = accessions,
    bnStart = hit_starts,
    bnEnd = hit_starts + match_len,
    na2d_structure = structures,
    na2d_energy = energies,
    na2d_score = na2d_scores
  )
  
  if(keep_sequences){
    result_df$sequence <- sequences
  }
  
  return(result_df)
}

#' Adjusts the default MBFV Zammit et. al table derived from their S1 table.
#' Split potential SL repeats into separate columns, adjust SL locations to be
#' relative to the 3UTR start
#' 
#' @param mbfv DataFrame Zammit S1 table parsed into CSV (MBFV-like viruses)
parseMBFV_SL <- function(mbfv){
  mbfv <- mbfv %>% 
    # Not interested in DB structures - filter these columns out
    select(-matches('DB')) %>%
    # MBFV viruses can have 2 SL elements. Separate the MBFV.SL column into two
    tidyr::separate(MBFV.SL, 
                    sep=' ', 
                    c('MBFV.SL1', 'MBFV.SL2'))
  
  # Each SL element has a start and end -> separate the ranges into own columns
  mbfv <- mbfv %>% 
    tidyr::separate(MBFV.SL1,
                    sep='-',
                    c('MBFV.SL1.start', 'MBFV.SL1.end')) %>%
    tidyr::separate(MBFV.SL2,
                    sep='-',
                    c('MBFV.SL2.start', 'MBFV.SL2.end')) %>%
    mutate(across(matches('MBFV.SL\\d+.'), as.numeric)) %>% 
    # Then subtract 3'UTR start to get the relative coordinates matching ones of NA2D
    # New columns have a suffix '_rel'
    mutate(across(matches('MBFV.SL\\d+.'), list(rel = ~ . - ORF.3.end)))
  
  return(mbfv)
}

#' Adjusts the default MBFV Zammit et. al table (originally S1 table).
#' Split potential Y-SL repeats into separate columns, adjust Y-SL locations to be
#' relative to the 3UTR start
#' 
#' @param mbfv DataFrame Zammit S1 table parsed into CSV (TBFV-like viruses)
parseTBFV_Y_SL <- function(tbfv){
  tbfv <- tbfv %>% 
    # TBFV viruses can have 2 SL elements. Separate the SL column into two
    tidyr::separate(TBFV.Y.SL, 
                    sep=' ', 
                    c('TBFV.Y.SL1', 'TBFV.Y.SL2'))
  
  # Each SL element has a start and end -> separate the ranges into own columns
  tbfv <- tbfv %>% 
    tidyr::separate(TBFV.Y.SL1,
                    sep='-',
                    c('TBFV.Y.SL1.start', 'TBFV.Y.SL1.end')) %>%
    tidyr::separate(TBFV.Y.SL2,
                    sep='-',
                    c('TBFV.Y.SL2.start', 'TBFV.Y.SL2.end')) %>%
    mutate(across(matches('TBFV.Y.SL\\d+.'), as.numeric)) %>% 
    # Then subtract 3'UTR start to get the relative coordinates matching ones of NA2D
    # New columns have a suffix '_rel'
    mutate(across(matches('TBFV.Y.SL\\d+.'), list(rel = ~ . - ORF.3.end)))
  
  return(tbfv)
}

#' Adjusts the default MBFV Zammit et. al table (originally S1 table).
#' Separate the AU-SL range locations to start and end columns relative to the 3UTR start
#' 
#' @param mbfv DataFrame Zammit S1 table parsed into CSV (TBFV-like viruses)
parseTBFV_AU_SL <- function(tbfv){
  tbfv <- tbfv %>% 
    tidyr::separate(TBFV.AU.SL,
                    sep='-',
                    c('TBFV.AU.SL.start', 'TBFV.AU.SL.end')) %>%
    mutate(across(matches('TBFV.AU.SL.'), as.numeric)) %>% 
    # Then subtract 3'UTR start to get the relative coordinates matching ones of NA2D
    # New columns have a suffix '_rel'
    mutate(across(matches('TBFV.AU.SL.'), list(rel = ~ . - ORF.3.end)))
  
  return(tbfv)
}

#' Adjusts the default MBFV Zammit et. al table (originally S1 table).
#' Separate the 3' GC-SL range locations to start and end columns relative to the 3UTR start
#' 
#' @param mbfv DataFrame Zammit S1 table parsed into CSV (TBFV-like viruses)
parseTBFV_3GC_SL <- function(tbfv){
  tbfv <- tbfv %>% 
    tidyr::separate(TBFV.3.GC.SL,
                    sep='-',
                    c('TBFV.3.GC.SL.start', 'TBFV.3.GC.SL.end')) %>%
    mutate(across(matches('TBFV.3.GC.SL.'), as.numeric)) %>% 
    # Then subtract 3'UTR start to get the relative coordinates matching ones of NA2D
    # New columns have a suffix '_rel'
    mutate(across(matches('TBFV.3.GC.SL.'), list(rel = ~ . - ORF.3.end)))
  
  return(tbfv)
}

#' Adjusts the default MBFV Zammit et. al table (originally S1 table).
#' Separate the 5' GC-SL range locations to start and end columns relative to the 3UTR start
#' 
#' @param mbfv DataFrame Zammit S1 table parsed into CSV (TBFV-like viruses)
parseTBFV_5GC_SL <- function(tbfv){
  tbfv <- tbfv %>% 
    tidyr::separate(TBFV.5.GC.SL,
                    sep='-',
                    c('TBFV.5.GC.SL.start', 'TBFV.5.GC.SL.end')) %>%
    mutate(across(matches('TBFV.5.GC.SL.'), as.numeric)) %>% 
    # Then subtract 3'UTR start to get the relative coordinates matching ones of NA2D
    # New columns have a suffix '_rel'
    mutate(across(matches('TBFV.5.GC.SL.'), list(rel = ~ . - ORF.3.end)))
  
  return(tbfv)
}



#' Adjusts the default Zammit et. al tables (originally S1 table).
#' Split potential DB+PK repeats into separate columns, adjust DB+PK locations to be
#' relative to the 3UTR start
#' 
#' @param mbfv DataFrame Zammit S1 table parsed into CSV (MBFV-like viruses)
parseDB_PK <- function(mbfv){
  mbfv <- mbfv %>% 
    select(-MBFV.SL) %>%
    # MBFV viruses can have 2 DB+PK elements. Separate the DB2...PK column into two
    tidyr::separate(DB2...PK, 
                    sep=' ', 
                    c('DB_PK_1', 'DB_PK_2'))
  
  # Each DB element has start and end -> parse the ranges into 2 columns
  mbfv <- mbfv %>% 
    tidyr::separate(DB_PK_1,
                    sep='-',
                    c('DB_PK_1.start', 'DB_PK_1.end')) %>%
    tidyr::separate(DB_PK_2,
                    sep='-',
                    c('DB_PK_2.start', 'DB_PK_2.end')) %>%
    mutate(across(matches('DB_PK_\\d+.'), as.numeric)) %>% 
    # Then subtract 3'UTR start to get the relative coordinates matching ones of NA2D
    # New columns have a suffix '_rel'
    mutate(across(matches('DB_PK_\\d+.'), list(rel = ~ . - ORF.3.end)))
  
  return(mbfv)
}


#' Adjusts the default Zammit et. al tables (originally S1 table).
#' DB_withoutPK elements are not repeated. Adjust location coordinates to be
#' relative to the 3'UTR start
#' 
#' @param mbfv DataFrame Zammit S1 table parsed into CSV (MBFV-like viruses)
parseDB_withoutPK <- function(mbfv){
  # Each DB element has start and end -> parse the ranges into 2 columns
  mbfv <- mbfv %>% 
    tidyr::separate(DB1..without.PK.,
                    sep='-',
                    c('DB_without.PK.start', 'DB_without.PK.end')) %>%
    mutate(across(matches('DB_without.PK'), as.numeric)) %>% 
    # Then subtract 3'UTR start to get the relative coordinates matching ones of NA2D
    # New columns have a suffix '_rel'
    mutate(across(matches('DB_without.PK'), list(rel = ~ . - ORF.3.end)))
  
  return(mbfv)
}

#' Merge NA2DSearch results with the corresponding Zammit et. al table to determine
#' Which hits have been recovered by NA2Dsearch and which viruses were completely missed
#' (columns match_SL1, match_SL2)
#' 
#' @param mbfv Raw table a al Zammit S1 table (for MBFV)
#' @param na2d_results DF with parsed NA2DSearch results
#' 
#' @return Joined DF with columns indicating which structures were hit, which weren't
#'         and what viruses they belong to
combineZammitAndN2D_MBFV_SL <- function(mbfv, na2d_results){
  mbfv_parsed <- parseMBFV_SL(mbfv)
  
  # Join NA2Dhits with the reference table to find out which hits were obtained
  mbfv_join_na2d <- full_join(mbfv_parsed, 
                              na2d_results, 
                              by=c("accession")) %>%
    select(-description, -matches('Type'))
  
  # Add 2 columns telling what SL element (if any) was NA2Dsearch able to match
  mbfv_join_na2d <- mbfv_join_na2d %>% 
    mutate(match_SL1 = bnStart -4 <= MBFV.SL1.start_rel 
           & bnEnd + 4 >= MBFV.SL1.end_rel, .after = virus) %>%
    mutate(match_SL2 = bnStart -4 <= MBFV.SL2.start_rel 
           & bnEnd + 4 >= MBFV.SL2.end_rel, .after = virus) %>%
    # Rearrange columns for better readability
    select(accession, virus, matches('match'), matches('bn'), matches('rel'), matches('MBFV'))
  
  return(mbfv_join_na2d)
}

combineZammitAndN2D_TBFV_Y_SL <- function(tbfv, na2d_results){
  tbfv_parsed <- parseTBFV_Y_SL(tbfv)
  
  # Join NA2Dhits with the reference table to find out which hits were obtained
  tbfv_join_na2d <- full_join(na2d_results, 
                              tbfv_parsed, 
                              by=c("accession")) %>%
    select(-description, -matches('Type'))
  
  # Add 2 columns telling what Y SL elements (if any) NA2Dsearch was able to match
  # Subtract 24 from Y.SL coordinates (long unstructured region (with PK) not captured by query)
  tbfv_join_na2d <- tbfv_join_na2d %>% 
    mutate(match_Y_SL_1 = bnStart - 4 <= TBFV.Y.SL1.start_rel 
           & bnEnd + 4 >= TBFV.Y.SL1.end_rel - 24, .after = virus) %>%
    mutate(match_Y_SL_2 = bnStart - 4 <= TBFV.Y.SL2.start_rel 
           & bnEnd + 4 >= TBFV.Y.SL2.end_rel - 24, .after = virus) %>%
    # Rearrange columns for better readability
    select(accession, virus, matches('match'), matches('bn'), matches('rel'), matches('TBFV'))
  
  return(tbfv_join_na2d)
}

combineZammitAndN2D_TBFV_AU_SL <- function(tbfv, na2d_results){
  tbfv <- parseTBFV_AU_SL(tbfv)
  
  # Keep only records containing AU SL
  tbfv <- tbfv %>% 
    filter(!is.na(TBFV.AU.SL.start))
  
  # Join NA2Dhits with the reference table to find out which hits were obtained
  tbfv_join_na2d <- full_join(na2d_results, 
                              tbfv, 
                              by=c("accession")) %>%
    select(-description, -matches('Type'))
  
  # Add 2 columns telling what Y SL elements (if any) NA2Dsearch was able to match
  tbfv_join_na2d <- tbfv_join_na2d %>% 
    mutate(match_AU_SL = bnStart - 4 <= TBFV.AU.SL.start_rel 
           & bnEnd + 4 >= TBFV.AU.SL.end_rel, .after = virus) %>%
    # Rearrange columns for better readability
    select(accession, virus, matches('match'), matches('bn'), matches('rel'), matches('TBFV'))
  
  return(tbfv_join_na2d)
}

combineZammitAndN2D_MBFV_DB_PK <- function(mbfv, na2d_results){
  mbfv <- parseDB_PK(mbfv)
  
  # Keep only records containing at least one repeat of DB+PK
  mbfv <- mbfv %>% 
    filter(if_any(.cols = matches("DB_PK_\\d.start_rel"),
                  .fns = ~ !is.na(.))) 

  # Join NA2Dhits with the reference table to find out which hits were obtained
  mbfv_join_na2d <- full_join(na2d_results, 
                              mbfv, 
                              by=c("accession")) %>%
    select(-description, -matches('Type'))
  
  # Add 2 columns telling what DB+PK element (if any) was NA2Dsearch able to match
  # We subtract 10-20 as the dangling unstructured (pseudokontted) part is not captured by the query
  mbfv_join_na2d <- mbfv_join_na2d %>% 
    mutate(match_DB2_1 = bnStart -4 <= DB_PK_1.start_rel
           & bnEnd + 4 >= DB_PK_1.end_rel - 20, .after = virus) %>%
    mutate(match_DB2_2 = bnStart -4 <= DB_PK_2.start_rel 
           & bnEnd + 4 >= DB_PK_2.end_rel - 20, .after = virus) %>%
    # Rearrange columns for better readability
    select(accession, virus, matches('match'), matches('bn'), matches('rel'), matches('DB'))
  
  return(mbfv_join_na2d)
}


combineZammitAndN2D_MBFV_DB1_noPK <- function(mbfv, na2d_results){
  mbfv <- parseDB_withoutPK(mbfv)
  
  # Keep only records containing at least one repeat of DB1 noPK
  mbfv <- mbfv %>% 
    filter(!is.na(DB_without.PK.start)) 
  
  # Join NA2Dhits with the reference table to find out which hits were obtained
  mbfv_join_na2d <- full_join(na2d_results, 
                              mbfv, 
                              by=c("accession")) %>%
    select(-description, -matches('Type'))
  
  # Add 2 columns telling what DB+PK element (if any) was NA2Dsearch able to match
  mbfv_join_na2d <- mbfv_join_na2d %>% 
    mutate(match_DB1 = bnStart -4 <= DB_without.PK.start_rel
           & bnEnd + 4 >= DB_without.PK.end_rel, .after = virus) %>%
    # Rearrange columns for better readability
    select(accession, virus, matches('match'), matches('bn'), matches('rel'), matches('DB'))
  
  return(mbfv_join_na2d)
}

combineZammitAndN2D_TBFV_3GC_SL <- function(tbfv, na2d_results){
  tbfv <- parseTBFV_3GC_SL(tbfv)
  
  # Keep only records containing at least one repeat of DB+PK
  tbfv <- tbfv %>% 
    filter(!is.na(TBFV.3.GC.SL.start_rel)) 
  
  # Join NA2Dhits with the reference table to find out which hits were obtained
  tbfv_join_na2d <- full_join(na2d_results, 
                              tbfv, 
                              by=c("accession", "Type")) %>%
    select(-description, -Type)
  
  # Add 2 columns telling what DB+PK element (if any) was NA2Dsearch able to match
  tbfv_join_na2d <- tbfv_join_na2d %>% 
    mutate(match_3GC_SL = bnStart -4 <= TBFV.3.GC.SL.start_rel
           & bnEnd + 4 >= TBFV.3.GC.SL.end_rel, .after = virus) %>%
    # Rearrange columns for better readability
    select(accession, virus, matches('match'), matches('bn'), matches('rel'), matches('TBFV.3.GC'))
  
  return(tbfv_join_na2d)
}

getDetectedElements <- function(zammit_join_na2d){
  # Group the joined table by virus. 
  # Sum the match flags - if there is any correct 3'UTR hit, the summarizing value > 0
  #                     - if the virus doesn't have such structure, the match flags will be NULL
  # Sum bnStart values for each virus - if there are no NA2D hits, the sum will be NULL
  grouped_zammit_join_na2d <- zammit_join_na2d %>% 
    group_by(virus, accession) %>% 
    summarise_all(mean) %>%
    select(virus, accession, matches('match'), matches('bn'), matches('_rel')) %>%
    mutate(across(matches('match'), as.logical))
  
  # Store names of match columns (saying True/False whether the element was matched)
  match_cols <- grep('match', 
                     colnames(grouped_zammit_join_na2d), 
                     value=T)
  
  if(length(match_cols) == 1){
    # Single 'match' column - must be true
    all_found_signature <- c(TRUE)
  } else {
    # Double 'match' column. Due to how the 'match' columns are constructed, if one column is TRUE, the other (grouped) column can be NA 
    # only if the repeat is not present. Care for the NA+NA case, which in contrast signals complete miss (see below)
    all_found_signature <- c(TRUE, NA)
  }
  # The element is lost if match is FALSE (indicating wrong match region) 
  # or NA (in no-repeat structures, the joint table holds only names of viruses containing the element -> hence NA indicates missed element)
  #       (in duplicated structures there can be up to 2 repeats. The NA can mean the second/first duplication is actually not present. 
  #            in any case, any combination of NA+FALSE indicates one or both structures missed)
  all_lost_signature <- c(FALSE, NA)
  
  #-----------------------------------------------------------------------------
  # Find viruses with all elements matched
  detected_all <- grouped_zammit_join_na2d %>% 
    filter(if_all(starts_with('match_'), ~ . %in% all_found_signature)) %>%
    # Remove cases, where both match columns are NA (i.e. no match and thus bnStart is NA too)
    filter(!is.na(bnStart)) %>% 
    # Remove cases, where there is only accession and no virus name (i.e. NA2D matched a record not present in Zammit table)
    filter(!is.na(virus))
  #-----------------------------------------------------------------------------
  # Find viruses with none of the elements matched
  detected_none <- grouped_zammit_join_na2d %>% 
    filter(if_all(starts_with('match_'), ~ . %in% all_lost_signature)) %>%
    # Remove cases, where there is only accession and no virus name (i.e. NA2D matched a record not present in Zammit table)
    filter(!is.na(virus))
  
  # False positives are NA2D accessions that after merge have no virus name (which is supplied by the Zammit table with true positives )
  false_positives <- grouped_zammit_join_na2d %>%
    filter(is.na(virus))
  
  return(list(detected_none = detected_none,
              detected_all = detected_all,
              false_positives = false_positives,
              grouped_results = grouped_zammit_join_na2d
              ))
  
}

#' Process a row of Zammit table S1. If a virus under given accession has a structure
#' (start_ is not NULL), then extract the sub-sequence. Return NA otherwise. 
zammit_subseq <- function(accession, start_, end_){
  if(!exists("zammit_fasta")){
    zammit_fasta <- readRNAStringSet('benchmark/Zammit/Flavivi_3UTRs.fa')
  }
  
  if (is.na(start_)) {
    return(NA)
  }
  
  to_retrieve <- grep(accession, names(zammit_fasta))
  extracted <- zammit_fasta[to_retrieve] %>% 
    as.character() %>% 
    substring(start_, end_)
  
  return(extracted)
}

extract_undetected <- function(parsed_fv_df, 
                               elementColumnName, 
                               na2d_detected_accs, 
                               base_fasta){
  # Get accessions of missed records (selecting only repeats searched by NA2D query)
  missing_elements_acc <- parsed_fv_df %>% 
    filter(!is.na(across(all_of(elementColumnName))) &
             !accession %in% na2d_detected_accs) %>%
    pull(accession)
  
  if(length(missing_elements_acc) == 0){
    return(c())
  }
  
  # Get full FASTA IDs of missed records (they contain accession)
  to_extract <- grep(paste(missing_elements_acc, collapse = '|'), 
                     names(base_fasta), 
                     value=T)
  
  missed_elements_seqs <- base_fasta[names(base_fasta) %in% to_extract]
  
  return(missed_elements_seqs)
}

# Save the undetected sequences for later processing
saveUndetected <- function(detected_none, filepath){
  zammit_fasta <- readRNAStringSet('benchmark/Zammit/Flavivi_3UTRs.fa')
  
  to_retrieve <- sapply(detected_none$accession, function(u) grep(u, names(zammit_fasta))) # get indexes of FASTAs to retrieve
  writeXStringSet(zammit_fasta[to_retrieve], filepath)
}