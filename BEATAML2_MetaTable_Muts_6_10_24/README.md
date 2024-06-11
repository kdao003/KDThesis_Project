This script generates 2 tables. 

The first table are samples that contain a mutation in either SRSF2|U2AF1|SF3B1|ZRSR2|ASXL1|EZH2|BCOR|STAG2. This table is considered patients who 
are likely to have secondary AML.

The second table are samples that contain a RUNX1/RUNX1T1 or CBFB/MYH11 or 11q23 or NPM1 or TP53. This table is considered patients who are likely to have
de novo AML.

Generating the first table required myeloid-lineage scores to be generated in EASY EXPRESSION APP ssGSEA boxplot tab table generating feature using gene sets
that corresponded to HSC, HSC/Prog-like, etc. MDS-enriched (Use UpinMDS_Strict_Cutoff) were generated this same way. Inflammations scores were generated the same way
but hte input was the MSigsDB Hallmark Inflammatory Response gene set. 

