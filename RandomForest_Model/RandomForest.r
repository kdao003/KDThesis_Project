library(randomForest)
library(randomForestExplainer)
library(caret)
library(dplyr)
library(ROSE)
library(pROC)
library(DiagrammeR)# For plotting the decision tree
library(iml) #To calculate SHAP for random forest
library(ggplot2) #To visualize SHAP values


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

# Separate features and target label
features <- data[, !names(data) %in% c("dbgap_rnaseq_sample", "priorMDS")]
labels <- data$priorMDS

# Split the data into training and testing sets
set.seed(123) # For reproducibility
trainIndex <- createDataPartition(labels, p = .75, list = FALSE, times = 1)
trainData <- features[trainIndex,]
trainLabels <- labels[trainIndex]
testData <- features[-trainIndex,]
testLabels <- labels[-trainIndex]

# Combine training data and labels for balancing
trainData_combined <- data.frame(trainData, isDenovo = trainLabels)

# Balance the training data
balanced_train_data <- ROSE(isDenovo ~ ., data = trainData_combined, seed = 123)$data

# Separate balanced features and labels
trainData <- balanced_train_data[, !names(balanced_train_data) %in% "isDenovo"]
trainLabels <- balanced_train_data$isDenovo

forest <- randomForest(x = trainData, y = trainLabels, ntree = 1000, localIMP = TRUE)

# Save the model
save(forest, file = "MDS_randomForest.rda")

# Load the model (if needed)
load("MDS_randomForest.rda")

#Learning Curve: Plot the error associated with the number of decision trees
plot(forest, main = "Learning Curve of the Forest")
legend("topright", c("error for 'no priorMDS'", "misclassification error", "error for 'yes prior MDS'"), 
       lty = c(1, 1, 1), col = c("green", "black", "red"))

#Distribution of Minimal Depth
min_depth_frame <- min_depth_distribution(forest)
head(min_depth_frame, n = 10)

#Plotting Minimal Depth
plot_min_depth_distribution(min_depth_frame)

#Plotting Minimal Depth with Tree Cutoff
plot_min_depth_distribution(min_depth_frame, min_no_of_trees = 60, mean_sample = "relevant_trees")

#Exploring Variable Importance Measures
importance_frame <- measure_importance(forest)
save(importance_frame, file = "MDS_randomForest_importance_frame.rda")
load("MDS_randomForest_importance_frame.rda")
head(importance_frame, n = 10)


#SHAP
#*Uses iml package*
# Convert train labels and data to numeric
trainData <- data.frame(sapply(trainData, as.numeric))
trainLabels <- as.numeric(trainLabels) - 1

#New Forest for SHAP that uses only numeric values with no factor values
forest1 <- randomForest(x = trainData, y = trainLabels, ntree = 1000, localIMP = TRUE)

# Create an iml Predictor object
predictor <- Predictor$new(forest1, data = trainData, y = trainLabels)

# Calculate SHAP values
shapley <- Shapley$new(predictor, x.interest = trainData[1, ])

# For all instances in the training set
shapley_all <- FeatureImp$new(predictor, loss = "mse")

# Plot SHAP values using plot function
plot(shapley_all, type = "bar", main = "SHAP Values for Random Forest Model")












