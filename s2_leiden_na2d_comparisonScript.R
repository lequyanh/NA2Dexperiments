setwd('C:\\Users\\vuleq\\PycharmProjects\\na2dsearch')
source('s2_leiden_na2d_integrationFunctions.R')

load_tables()

# This call joins the appropriate S2 table (MBFV, cISFV,...; holding reference structure locations) 
# with the appropriate webapp result file targeting a specific structure (AU SL, DB+PK, ...)
results <- getCombined_S2_webbapp('MBFV-SL')

# results have 3 slots
merged_table <- results$merged_s2_web   # inner join of the S2 table and webapp results
webapp_df <- results$parsed_leiden      # parsed webapp results
s2_table <- results$s2_table            # reference S2 table (reference structure ranks and locations)

# Assign structures to elements referenced in the S2 table
# Group by accession and structure location
# Each such group has the correct structure encoded in the 'structure rank' column
s2_with_structures <- merged_table %>% 
  group_by(accession, s2_rel_start) %>%
  group_modify(~ {.x %>% 
      arrange(deltaG) %>% 
      dplyr::slice(dplyr::first(.x$s2_structure_rank))}
  )

# Assert that the extracted structure has the same energy as proclaimed in the S2 table
assertthat::are_equal(s2_with_structures$deltaG,
                      s2_with_structures$s2_ref_deltaG)

assertthat::are_equal(s2_table$accession %>% sort,
                      s2_with_structures$accession %>% sort)

# DEBUG: Check for inconsistencies if asserts above fail
s2_table[!s2_table$accession %in%  s2_with_structures$accession, ] %>% View() # Missed          S2 viruses
s2_with_structures %>% filter(deltaG != s2_ref_deltaG) %>% View()             # Wrong-structure S2 viruses

###################################
# LEIDEN FALSE POSITIVES ANALYSIS #
###################################
# False positive locations of S2 viruses
fp_s2 <- merged_table %>% 
  # Filter records with S2 virus but without S2_location 
  #(i.e. the rolling join could not find an S2 location for a given Leiden hit)
  filter(is.na(s2_rel_start) & !is.na(s2_virus)) %>% 
  # Create an artificial grouping column as the locations can be scattered +-2nt
  mutate(approx.loc = 10* round(hit_loc/10)) %>% 
  group_by(s2_virus, approx.loc) %>% 
  summarise(no_false_hits = n(), 
            hit_loc = round(mean(hit_loc)),
            min_deltaG = round(min(deltaG), 1))  %>%
  arrange(min_deltaG)

View(fp_s2)

#-------------------------------------------------------------------------------
# False positive locations of non-S2 viruses
fp_non_s2 <- merged_table %>% 
  # Filter records with S2 virus but without S2_location 
  #(i.e. the rolling join could not find an S2 location for a given Leiden hit)
  filter(!is.na(hit_loc) & is.na(s2_virus)) %>% 
  # Create an artificial grouping column as the locations can be scattered +-2nt
  mutate(approx.loc = 10* round(hit_loc/10)) %>% 
  group_by(accession, seq_def, approx.loc) %>% 
  summarise(no_false_hits = n(), 
            hit_loc = round(mean(hit_loc)),
            min_deltaG = round(min(deltaG), 1))  %>%
  arrange(min_deltaG)

View(fp_non_s2)

fp_non_s2_summary <- fp_non_s2 %>%
  group_by(accession, seq_def) %>%
  summarise(unique_loc_counts = n()) %>%
  arrange(desc(unique_loc_counts))

View(fp_non_s2_summary)

# DEBUG
assertthat::assert_that(intersect(
  fp_non_s2 %>% pull(accession) %>% unique(),
  s2_table %>% pull(accession)
  ) %>% length() == 0)

#######################################################
# COMPARING NA2D RESULTS WITH REFERENCE S2+STRUCUTRES #
#######################################################
# !! Here we assume, that S2 tables will be enriched with correct reference structures !!
comparisonDf <- compareNa2d_S2(structure_type = 'TBFV-Y-SL', 
                               ignore_leiden_pk = T)
# DEBUG: use different NA2D results for comparison with reference
comparisonDf <- compareNa2d_S2(structure_type = 'MBFV-SL', 
                               na2d_file_override = 'NA2DsearchProject/searches/search_13/results.fa',
                               ignore_leiden_pk = T)

# View S2 structures that NA2D missed. Should be empty if NA2D found all relevant S2 structures
missed <- comparisonDf %>% filter(!is.na(s2_virus) & is.na(na2d_score)) %>%
  mutate(ref_simpl_struct = simplifyWussStructures(structure), .after=s2_virus)
View(missed)

# View S2 structures that matched with a NA2D hit (simpl_struct is not NA)
View(comparisonDf %>% filter(!is.na(simpl_struct)))

# Show which viruses are extra. Group by accession and summarize by mean (ignoring NAs)
# Keep only those with non-zero hit_counts. 
# Since structure_rank (from S2 table) has value only if the virus contains the examined structure, NaN indicates false positives
comparisonDf %>%
  group_by(accession, s2_virus) %>%
  summarise(across(where(is.numeric), ~round(mean(.x, na.rm = T), 1), .names="m_{.col}")) %>%
  filter(!is.na(m_na2d_score)) %>% 
  arrange(desc(m_na2d_score)) %>%
  View()
