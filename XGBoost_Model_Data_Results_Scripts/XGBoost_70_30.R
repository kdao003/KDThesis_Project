# Load necessary libraries
#install.packages("pROC")
#install.packages("SHAPforxgboost")

library(xgboost)
library(tidyverse)
library(caret)
library(dplyr)
library(SHAPforxgboost)
library(ROSE)
library(pROC)
library(DiagrammeR) 

# Load the data
data <- read.csv("/Volumes/T7/IDH_4970_I_and_II/B_Cell_Analysis/May_to_July_ProjectWork/JuneAndJuly_Work/XG_Boost_Model/MutationData_with_QNE_ssGSEA_Scores.csv", header = TRUE)


#Training and Validating with only BeatAML2 data
#Split is 70% training and 30% Testing
BeatAML2 <- data %>%
  filter(startsWith(sample, "BA")) %>%
  column_to_rownames(var = "sample") 

# Determining Categorical Columns
categorical_columns <- c("CBFB_MYH11_Present", "RUNX1_RUNX1T1_Present", "X11q23_Present", "STAG2_Mutation",
                         "BCOR_Mutation", "EZH2_Mutation", "SF3B1_Mutation", "SRSF2_Mutation", "U2AF1_Mutation",
                         "ZRSR2_Mutation", "ASXL1_Mutated", "TP53_Mutated", "RUNX1_Mutated", "FLT3_ITD_Mutated",
                         "priorMDS")

# Convert the columns to factors
BeatAML2[categorical_columns] <- lapply(BeatAML2[categorical_columns], as.factor)

# Check structure of data
str(BeatAML2)

# Dropping columns not of interest
blacklist <- c("priorMalignancyNonMyeloid", "priorMalignancyRadiationTx", "cumulativeChemo", 
               "nonAML_MDSMPN_specificDxAtAcquisition", "specificDxAtAcquisition_MDSMPN", 
               "isDenovo", "isTransformed", "dataset", "prior_hematologic_disorder_diagnosis_indicator", "leukemia_french_american_british_morphology_code",
               "priorMaligncyNonMyeloid", "priorMaligncyRadiationTx", "cumulativeChemo", 
               "nonAML_MDSMPN_specificDxAtAcquisition", "specificDxAtAcquisition_MDSMPN")

sub_data <- BeatAML2[, !names(BeatAML2) %in% blacklist, drop = FALSE]

# Ensure other columns are numeric if they are supposed to be (excluding categorical columns)
numeric_columns <- setdiff(colnames(sub_data), c("sample", categorical_columns))
sub_data[numeric_columns] <- lapply(sub_data[numeric_columns], as.numeric)

# Remove columns with only one unique value
sub_data <- sub_data[, sapply(sub_data, function(col) length(unique(col)) > 1)]

# Convert the target variable to a factor with 'y' as the positive class
sub_data$priorMDS <- factor(sub_data$priorMDS, levels = c('no', 'yes'))

# Separate features and target label
features <- sub_data[, !names(sub_data) %in% c("dbgap_rnaseq_sample", "priorMDS")]
labels <- sub_data$priorMDS

##Training with 70% data and 30% Testing

# Split the data into training and testing sets, 70% training, 30% testing)
set.seed(123) # For reproducibility
trainIndex <- createDataPartition(labels, p = 0.7, list = FALSE, times = 1)
trainData <- features[trainIndex,]
trainLabels <- labels[trainIndex]
testData <- features[-trainIndex,]
testLabels <- labels[-trainIndex]

# Combine training data and labels for balancing
trainData_combined <- data.frame(trainData, Target_PriorMDS_yes_or_no = trainLabels)

# Balance the training data
balanced_train_data <- ROSE(Target_PriorMDS_yes_or_no ~ ., data = trainData_combined, seed = 123)$data

# Separate balanced features and labels
trainData <- balanced_train_data[, !names(balanced_train_data) %in% "Target_PriorMDS_yes_or_no"]
trainLabels <- balanced_train_data$Target_PriorMDS_yes_or_no

# Convert factors to numeric for xgboost
trainData <- data.frame(sapply(trainData, as.numeric))
testData <- data.frame(sapply(testData, as.numeric))

# Convert labels to numeric
trainLabels <- as.numeric(trainLabels) - 1
testLabels <- as.numeric(testLabels) - 1

# Check for NA values and impute or remove if necessary
trainData[is.na(trainData)] <- 0
testData[is.na(testData)] <- 0

# Prepare data for xgboost
dtrain <- xgb.DMatrix(data = as.matrix(trainData), label = trainLabels)
dtest <- xgb.DMatrix(data = as.matrix(testData), label = testLabels)

# Set parameters for xgboost
params <- list(
  objective = "binary:logistic",
  eval_metric = "logloss",
  max_depth = 6,
  eta = 0.1,
  nthread = 2,
  booster = "gbtree"
)

#Obtaining CV Error


# Train the xgboost model
set.seed(123)
xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 100,
  watchlist = list(train = dtrain, eval = dtest),
  early_stopping_rounds = 10,
  verbose = 1
)

#saveRDS(xgb_model, "xgboost_model_70_30.rds")

#readRDS("xgboost_model_70_30.rds")

# Make predictions
preds <- predict(xgb_model, newdata = dtest)
predictions <- ifelse(preds > 0.5, 1, 0)

# Convert predictions back to factor
predictions <- factor(predictions, levels = c(0, 1), labels = c('no', 'yes'))

# Evaluate the model with 'yes' as the positive class
conf_matrix1 <- confusionMatrix(predictions, factor(testLabels, levels = c(0, 1), labels = c('no', 'yes')), positive = "yes")
print(conf_matrix1)

# Extract precision (Positive Predictive Value)
precision1 <- conf_matrix1$byClass["Pos Pred Value"]

# Extract recall (Sensitivity)
recall1 <- conf_matrix1$byClass["Sensitivity"]

# Calculate F1 score
f1_score1 <- 2 * (precision1 * recall1) / (precision1 + recall1)

# Print precision and F1 score
cat("Precision: ", precision1, "\n")
cat("F1 Score: ", f1_score1, "\n")





# Calculate SHAP values
shap_values <- shap.values(xgb_model = xgb_model, X_train = as.matrix(trainData))

# Plot summary of SHAP values
shap_long <- shap.prep(xgb_model = xgb_model, X_train = as.matrix(trainData))
shap.plot.summary(shap_long)

# Plot ROC curve
roc_curve <- roc(testLabels, preds)
plot(roc_curve, col = "blue", main = "ROC Curve for XGBoost Model")
auc(roc_curve)

# Plot an instance of the decision tree
xgb.plot.tree(model = xgb_model, trees = 0)  # Plot the first tree




























#Removing these mutation features: STAG2, X11q23, U2AF1, SRSF2, SF3B1, RUNX1_RUNXT1, BCOR, EZH2, ZRSR2
#Removing these transcriptional features: Promono, GMP


library(xgboost)
library(tidyverse)
library(caret)
library(dplyr)
library(SHAPforxgboost)
library(ROSE)
library(pROC)
library(DiagrammeR) 

# Load the data
data <- read.csv("/Volumes/T7/IDH_4970_I_and_II/B_Cell_Analysis/May_to_July_ProjectWork/JuneAndJuly_Work/XG_Boost_Model/MutationData_with_QNE_ssGSEA_Scores.csv", header = TRUE)


#Training and Validating with only BeatAML2 data
#Split is 70% training and 30% Testing
BeatAML2 <- data %>%
  filter(startsWith(sample, "BA")) %>%
  column_to_rownames(var = "sample") 

# Determining Categorical Columns
categorical_columns <- c("CBFB_MYH11_Present", "ASXL1_Mutated", "TP53_Mutated", "RUNX1_Mutated", "FLT3_ITD_Mutated",
                         "priorMDS")

# Convert the columns to factors
BeatAML2[categorical_columns] <- lapply(BeatAML2[categorical_columns], as.factor)

# Check structure of data
str(BeatAML2)

# Dropping columns not of interest
blacklist <- c("priorMalignancyNonMyeloid", "priorMalignancyRadiationTx", "cumulativeChemo", 
               "nonAML_MDSMPN_specificDxAtAcquisition", "specificDxAtAcquisition_MDSMPN", 
               "isDenovo", "isTransformed", "dataset", "prior_hematologic_disorder_diagnosis_indicator", "leukemia_french_american_british_morphology_code",
               "priorMaligncyNonMyeloid", "priorMaligncyRadiationTx", "cumulativeChemo", 
               "nonAML_MDSMPN_specificDxAtAcquisition", "specificDxAtAcquisition_MDSMPN", "X11q23_Present", "STAG2_Mutation",
               "SRSF2_Mutation", "U2AF1_Mutation", "SF3B1_Mutation", "RUNX1_RUNX1T1_Present","BCOR_Mutation", "EZH2_Mutation",
               "ZRSR2_Mutation", "GMP.like", "Promono.like")

sub_data <- BeatAML2[, !names(BeatAML2) %in% blacklist, drop = FALSE]

# Ensure other columns are numeric if they are supposed to be (excluding categorical columns)
numeric_columns <- setdiff(colnames(sub_data), c("sample", categorical_columns))
sub_data[numeric_columns] <- lapply(sub_data[numeric_columns], as.numeric)

# Remove columns with only one unique value
sub_data <- sub_data[, sapply(sub_data, function(col) length(unique(col)) > 1)]

# Convert the target variable to a factor with 'y' as the positive class
sub_data$priorMDS <- factor(sub_data$priorMDS, levels = c('no', 'yes'))

# Separate features and target label
features <- sub_data[, !names(sub_data) %in% c("dbgap_rnaseq_sample", "priorMDS")]
labels <- sub_data$priorMDS

##Training with 70% data and 30% Testing

# Split the data into training and testing sets, 70% training, 30% testing)
set.seed(123) # For reproducibility
trainIndex <- createDataPartition(labels, p = 0.7, list = FALSE, times = 1)
trainData <- features[trainIndex,]
trainLabels <- labels[trainIndex]
testData <- features[-trainIndex,]
testLabels <- labels[-trainIndex]

# Combine training data and labels for balancing
trainData_combined <- data.frame(trainData, Target_PriorMDS_yes_or_no = trainLabels)

# Balance the training data
balanced_train_data <- ROSE(Target_PriorMDS_yes_or_no ~ ., data = trainData_combined, seed = 123)$data

# Separate balanced features and labels
trainData <- balanced_train_data[, !names(balanced_train_data) %in% "Target_PriorMDS_yes_or_no"]
trainLabels <- balanced_train_data$Target_PriorMDS_yes_or_no

# Convert factors to numeric for xgboost
trainData <- data.frame(sapply(trainData, as.numeric))
testData <- data.frame(sapply(testData, as.numeric))

# Convert labels to numeric
trainLabels <- as.numeric(trainLabels) - 1
testLabels <- as.numeric(testLabels) - 1

# Check for NA values and impute or remove if necessary
trainData[is.na(trainData)] <- 0
testData[is.na(testData)] <- 0

# Prepare data for xgboost
dtrain <- xgb.DMatrix(data = as.matrix(trainData), label = trainLabels)
dtest <- xgb.DMatrix(data = as.matrix(testData), label = testLabels)

# Set parameters for xgboost
params <- list(
  objective = "binary:logistic",
  eval_metric = "logloss",
  max_depth = 6,
  eta = 0.1,
  nthread = 2,
  booster = "gbtree"
)

#Obtaining CV Error
#set.seed(123)
#xgbcv <- xgb.cv(
#  params = params,
#  data = dtrain,
#  nrounds = 100,
#  nfold = 5,
#  showsd = TRUE,
#  stratified = TRUE,
#  print_every_n = 10,
#  early_stopping_rounds = 20,
#  maximize = FALSE
#)


# Train the xgboost model
set.seed(123)
xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 100,
  watchlist = list(train = dtrain, eval = dtest),
  early_stopping_rounds = 10,
  verbose = 1
)

#saveRDS(xgb_model, "xgboost_model_70_30_RemoveFeatures.rds")

#readRDS("xgboost_model_70_30_RemoveFeatures.rds")

# Make predictions
preds <- predict(xgb_model, newdata = dtest)
predictions <- ifelse(preds > 0.5, 1, 0)

# Convert predictions back to factor
predictions <- factor(predictions, levels = c(0, 1), labels = c('no', 'yes'))

# Evaluate the model with 'yes' as the positive class
conf_matrix1 <- confusionMatrix(predictions, factor(testLabels, levels = c(0, 1), labels = c('no', 'yes')), positive = "yes")
print(conf_matrix1)

# Extract precision (Positive Predictive Value)
precision1 <- conf_matrix1$byClass["Pos Pred Value"]

# Extract recall (Sensitivity)
recall1 <- conf_matrix1$byClass["Sensitivity"]

# Calculate F1 score
f1_score1 <- 2 * (precision1 * recall1) / (precision1 + recall1)

# Print precision and F1 score
cat("Precision: ", precision1, "\n")
cat("F1 Score: ", f1_score1, "\n")





# Calculate SHAP values
shap_values <- shap.values(xgb_model = xgb_model, X_train = as.matrix(trainData))

# Plot summary of SHAP values
shap_long <- shap.prep(xgb_model = xgb_model, X_train = as.matrix(trainData))
shap.plot.summary(shap_long)

# Plot ROC curve
roc_curve <- roc(testLabels, preds)
plot(roc_curve, col = "blue", main = "ROC Curve for XGBoost Model")
auc(roc_curve)

# Plot an instance of the decision tree
xgb.plot.tree(model = xgb_model, trees = 0)  # Plot the first tree






















#Removing these mutation features: STAG2, X11q23, U2AF1, SRSF2, SF3B1, RUNX1_RUNXT1, BCOR, EZH2, ZRSR2
#Removing all transcriptional features


library(xgboost)
library(tidyverse)
library(caret)
library(dplyr)
library(SHAPforxgboost)
library(ROSE)
library(pROC)
library(DiagrammeR) 

# Load the data
data <- read.csv("/Volumes/T7/IDH_4970_I_and_II/B_Cell_Analysis/May_to_July_ProjectWork/JuneAndJuly_Work/XG_Boost_Model/MutationData_with_QNE_ssGSEA_Scores.csv", header = TRUE)


#Training and Validating with only BeatAML2 data
#Split is 70% training and 30% Testing
BeatAML2 <- data %>%
  filter(startsWith(sample, "BA")) %>%
  column_to_rownames(var = "sample") %>%
  select(!cDC.like:UpInMDS_StrictCutoff)

# Determining Categorical Columns
categorical_columns <- c("CBFB_MYH11_Present", "ASXL1_Mutated", "TP53_Mutated", "RUNX1_Mutated", "FLT3_ITD_Mutated",
                         "priorMDS")

# Convert the columns to factors
BeatAML2[categorical_columns] <- lapply(BeatAML2[categorical_columns], as.factor)

# Check structure of data
str(BeatAML2)

# Dropping columns not of interest
blacklist <- c("priorMalignancyNonMyeloid", "priorMalignancyRadiationTx", "cumulativeChemo", 
               "nonAML_MDSMPN_specificDxAtAcquisition", "specificDxAtAcquisition_MDSMPN", 
               "isDenovo", "isTransformed", "dataset", "prior_hematologic_disorder_diagnosis_indicator", "leukemia_french_american_british_morphology_code",
               "priorMaligncyNonMyeloid", "priorMaligncyRadiationTx", "cumulativeChemo", 
               "nonAML_MDSMPN_specificDxAtAcquisition", "specificDxAtAcquisition_MDSMPN", "X11q23_Present", "STAG2_Mutation",
               "SRSF2_Mutation", "U2AF1_Mutation", "SF3B1_Mutation", "RUNX1_RUNX1T1_Present","BCOR_Mutation", "EZH2_Mutation",
               "ZRSR2_Mutation")

sub_data <- BeatAML2[, !names(BeatAML2) %in% blacklist, drop = FALSE]

# Ensure other columns are numeric if they are supposed to be (excluding categorical columns)
numeric_columns <- setdiff(colnames(sub_data), c("sample", categorical_columns))
sub_data[numeric_columns] <- lapply(sub_data[numeric_columns], as.numeric)

# Remove columns with only one unique value
sub_data <- sub_data[, sapply(sub_data, function(col) length(unique(col)) > 1)]

# Convert the target variable to a factor with 'y' as the positive class
sub_data$priorMDS <- factor(sub_data$priorMDS, levels = c('no', 'yes'))

# Separate features and target label
features <- sub_data[, !names(sub_data) %in% c("dbgap_rnaseq_sample", "priorMDS")]
labels <- sub_data$priorMDS

##Training with 70% data and 30% Testing

# Split the data into training and testing sets, 70% training, 30% testing)
set.seed(123) # For reproducibility
trainIndex <- createDataPartition(labels, p = 0.7, list = FALSE, times = 1)
trainData <- features[trainIndex,]
trainLabels <- labels[trainIndex]
testData <- features[-trainIndex,]
testLabels <- labels[-trainIndex]

# Combine training data and labels for balancing
trainData_combined <- data.frame(trainData, Target_PriorMDS_yes_or_no = trainLabels)

# Balance the training data
balanced_train_data <- ROSE(Target_PriorMDS_yes_or_no ~ ., data = trainData_combined, seed = 123)$data

# Separate balanced features and labels
trainData <- balanced_train_data[, !names(balanced_train_data) %in% "Target_PriorMDS_yes_or_no"]
trainLabels <- balanced_train_data$Target_PriorMDS_yes_or_no

# Convert factors to numeric for xgboost
trainData <- data.frame(sapply(trainData, as.numeric))
testData <- data.frame(sapply(testData, as.numeric))

# Convert labels to numeric
trainLabels <- as.numeric(trainLabels) - 1
testLabels <- as.numeric(testLabels) - 1

# Check for NA values and impute or remove if necessary
trainData[is.na(trainData)] <- 0
testData[is.na(testData)] <- 0

# Prepare data for xgboost
dtrain <- xgb.DMatrix(data = as.matrix(trainData), label = trainLabels)
dtest <- xgb.DMatrix(data = as.matrix(testData), label = testLabels)

# Set parameters for xgboost
params <- list(
  objective = "binary:logistic",
  eval_metric = "logloss",
  max_depth = 6,
  eta = 0.1,
  nthread = 2,
  booster = "gbtree"
)

#Obtaining CV Error
#set.seed(123)
#xgbcv <- xgb.cv(
#  params = params,
#  data = dtrain,
#  nrounds = 100,
#  nfold = 5,
#  showsd = TRUE,
#  stratified = TRUE,
#  print_every_n = 10,
#  early_stopping_rounds = 20,
#  maximize = FALSE
#)


# Train the xgboost model
set.seed(123)
xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 100,
  watchlist = list(train = dtrain, eval = dtest),
  early_stopping_rounds = 10,
  verbose = 1
)

#saveRDS(xgb_model, "xgboost_model_70_30_RemoveAllScores_RemoveSomeMut.rds")

#readRDS("xgboost_model_70_30_RemoveAllScores_RemoveSomeMut.rds")

# Make predictions
preds <- predict(xgb_model, newdata = dtest)
predictions <- ifelse(preds > 0.5, 1, 0)

# Convert predictions back to factor
predictions <- factor(predictions, levels = c(0, 1), labels = c('no', 'yes'))

# Evaluate the model with 'yes' as the positive class
conf_matrix1 <- confusionMatrix(predictions, factor(testLabels, levels = c(0, 1), labels = c('no', 'yes')), positive = "yes")
print(conf_matrix1)

# Extract precision (Positive Predictive Value)
precision1 <- conf_matrix1$byClass["Pos Pred Value"]

# Extract recall (Sensitivity)
recall1 <- conf_matrix1$byClass["Sensitivity"]

# Calculate F1 score
f1_score1 <- 2 * (precision1 * recall1) / (precision1 + recall1)

# Print precision and F1 score
cat("Precision: ", precision1, "\n")
cat("F1 Score: ", f1_score1, "\n")





# Calculate SHAP values
shap_values <- shap.values(xgb_model = xgb_model, X_train = as.matrix(trainData))

# Plot summary of SHAP values
shap_long <- shap.prep(xgb_model = xgb_model, X_train = as.matrix(trainData))
shap.plot.summary(shap_long)

# Plot ROC curve
roc_curve <- roc(testLabels, preds)
plot(roc_curve, col = "blue", main = "ROC Curve for XGBoost Model")
auc(roc_curve)

# Plot an instance of the decision tree
xgb.plot.tree(model = xgb_model, trees = 0)  # Plot the first tree





























