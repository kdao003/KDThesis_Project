library(tidyverse)

df <- read.delim("2022_BeatAML_Meta.txt", header = TRUE, sep = "\t")

df_subset <- df %>%
  select(dbgap_rnaseq_sample, consensusAMLFusions, otherCytogenetics, NPM1, TP53, variantSummary) %>%
  rename(sample = dbgap_rnaseq_sample)

sAML_mut_subset <- df_subset %>%
  select(everything()) %>%
  filter(str_detect(variantSummary, "SRSF2|U2AF1|SF3B1|ZRSR2|ASXL1|EZH2|BCOR|STAG2"))

write.table(sAML_mut_subset, file = "BEATAML2_sAML_mutation_subset_KD_6_10_24.txt", sep = "\t")

  
neg_mut_subset <- df_subset %>%
  select(everything()) %>%
  filter(str_detect(consensusAMLFusions, "RUNX1-RUNX1T1|CBFB/MYH11") | str_detect(otherCytogenetics, "11q23") 
         | str_detect(NPM1, "positive") | str_detect(TP53, "TP53")) 

write.table(neg_mut_subset, file = "BEATAML2_neg_mut_subset_KD_6_10_24.txt", sep = "\t")























