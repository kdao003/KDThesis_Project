library(tidyverse)

#sAML mutations from 2015 Blood Paper
df <- read.delim("2022_BeatAML_Meta.txt", header = TRUE, sep = "\t")

sAML_mut_subset <- df %>%
  select(everything()) %>%
  filter(str_detect(variantSummary, "SRSF2|U2AF1|SF3B1|ZRSR2|ASXL1|EZH2|BCOR|STAG2")) %>%
  rename(sample = dbgap_rnaseq_sample)


#Inserting myeloid-scores made from table in ssGSEA boxplot
df1 <- read.delim("/Volumes/T7/IDH_4970_I_and_II/B_Cell_Analysis/BeatAML2022_CelltypeScores/2022BeatAML_Inflam_CelltypeScores.txt", header = TRUE, sep = "\t")

sAML_mut_subset <- left_join(sAML_mut_subset, df1, by = "sample")

#Inserting MDS enriched scores
df2 <- read.delim("/Volumes/T7/IDH_4970_I_and_II/B_Cell_Analysis/ssGSEAscore_UpInMDS_StrictCutoff.tsv", header = TRUE, sep = "\t")

sAML_mut_subset <- left_join(sAML_mut_subset, df2, by = "sample")

#writing the table
write.table(sAML_mut_subset, file = "BEATAML2_sAML_mutation_subset_KD_6_10_24.txt", sep = "\t")
  
#Mutations that point to confirmed de novo AML 
neg_mut_subset <- df %>%
  select(everything()) %>%
  filter(str_detect(consensusAMLFusions, "RUNX1-RUNX1T1|CBFB/MYH11") | str_detect(otherCytogenetics, "11q23") 
         | str_detect(NPM1, "positive") | str_detect(TP53, "TP53")) 

write.table(neg_mut_subset, file = "BEATAML2_neg_mut_subset_KD_6_10_24.txt", sep = "\t")























