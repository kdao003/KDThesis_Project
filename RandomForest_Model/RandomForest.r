library(randomForest)
library(randomForestExplainer)


data <- read.csv("/Volumes/T7/IDH_4970_I_and_II/B_Cell_Analysis/XgBoost_Model/BEATAML_PriorMDS_MetaData.tsv", sep='\t')

blacklist = c("priorMalignancyNonMyeloid", "priorMalignancyRadiationTx", "cumulativeChemo", 
              "priorMalignancyNonMyeloid", "nonAML_MDSMPN_specificDxAtAcquisition", "specificDxAtAcquisition_MDSMPN", 
              "isDenovo", "isTransformed")
data <- data[, ! names(data) %in% blacklist, drop = F]
# Provided list of columns to convert to factors

#categorical_columns <- c("CBFB_MYH11_Present", "RUNX1_RUNX1T1_Present", "X11q23_Present", "STAG2_Mutation",
#                         "BCOR_Mutation", "EZH2_Mutation", "SF3B1_Mutation", "SRSF2_Mutation", "U2AF1_Mutation",
#                         "ZRSR2_Mutation", "ASXL1_Mutated", "TP53_Mutated", "RUNX1_Mutated", "FLT3_ITD_Mutated",
#                         "isDenovo", "isTransformed", "specificDxAtAcquisition_MDSMPN", 
#                         "nonAML_MDSMPN_specificDxAtAcquisition", "priorMalignancyNonMyeloid", 
#                         "cumulativeChemo", "priorMalignancyRadiationTx", "priorMDS")


categorical_columns <- c("CBFB_MYH11_Present", "RUNX1_RUNX1T1_Present", "X11q23_Present", "STAG2_Mutation",
                         "BCOR_Mutation", "EZH2_Mutation", "SF3B1_Mutation", "SRSF2_Mutation", "U2AF1_Mutation",
                         "ZRSR2_Mutation", "ASXL1_Mutated", "TP53_Mutated", "RUNX1_Mutated", "FLT3_ITD_Mutated",
                         "priorMDS")

# Convert the columns to factors
data[categorical_columns] <- lapply(data[categorical_columns], as.factor)

# Ensure other columns are numeric if they are supposed to be (excluding categorical columns)
numeric_columns <- setdiff(colnames(data), c("dbgap_rnaseq_sample", categorical_columns))
data[numeric_columns] <- lapply(data[numeric_columns], as.numeric)

# Remove columns with only one unique value
data <- data[, sapply(data, function(col) length(unique(col)) > 1)]

# Convert the target variable to a factor with 'y' as the positive class
data$priorMDS <- factor(data$priorMDS, levels = c('n', 'y'))

# Check the structure of the dataset
print(str(data))

set.seed(123) #For reproducibility

#Creaing RandomForest Model
forest <- randomForest(priorMDS ~., data = data, ntree = 100, localIMP = TRUE)
save(forest, file = "MDS_randomForest.rda")
load("MDS_randomForest.rda")

#Plot that examines error associated with number of decision trees
plot(forest, main = "Learning Curve of the Forest")
legend("topright", c("error for 'no priorMDS'", "misclassification error", "error for 'yes prior MDS'"), lty = c(1,1,1),
       col = c("green","black", "red"))

























