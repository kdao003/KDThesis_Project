library(randomForest)
library(randomForestExplainer)
library(caret)
library(dplyr)
library(ROSE)
library(pROC)
#library(DiagrammeR)# For plotting the decision tree
#library(iml) #To calculate SHAP for random forest
#library(ggplot2) #To visualize SHAP values
library(tidyverse)
#library(kernelshap) #To calculate SHAP for random forest
library(shapviz) #To visualize SHAP values
#library(rpart)
#library(rpart.plot)
library(reprtree)
#library(performanceEstimation) #For Balancing Training Dataset
library(treeshap)



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

#Convert the columns yes or y to 1, no or n to 0, Convert the columns to numeric
#yes or no values
data <- data %>%
  mutate(across(c("CBFB_MYH11_Present", "RUNX1_RUNX1T1_Present", "X11q23_Present", "STAG2_Mutation",
                  "BCOR_Mutation", "EZH2_Mutation", "SF3B1_Mutation", "SRSF2_Mutation", "U2AF1_Mutation",
                  "ZRSR2_Mutation"), function(x) ifelse(x == "yes", 1, 0)))
#TRUE or FALSE values
data <- data %>%
  mutate(across(c("ASXL1_Mutated", "TP53_Mutated", "RUNX1_Mutated", "FLT3_ITD_Mutated"), function(x) ifelse(x == TRUE,1,0)))
#y or n values
data <- data %>%
  mutate(priorMDS = ifelse(priorMDS == "y", 1, 0))

data[categorical_columns] <- lapply(data[categorical_columns], as.numeric)

# Ensure other columns are numeric if they are supposed to be (excluding categorical columns)
numeric_columns <- setdiff(colnames(data), c("dbgap_rnaseq_sample", categorical_columns))
data[numeric_columns] <- lapply(data[numeric_columns], as.numeric)

# Remove columns with only one unique value
data <- data[, sapply(data, function(col) length(unique(col)) > 1)]

# Convert the target variable to a factor with '1' as the positive class
data$priorMDS <- factor(data$priorMDS, levels = c('0', '1'))

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

#Combine testing data and labels
testData_combined <- data.frame(testData, isDenovo = testLabels)

# Balance the training data using ROSE
balanced_train_data <- ROSE(isDenovo ~ ., data = trainData_combined, seed = 123)$data

#Balance the training data using SMOTE
#balanced_train_data <- smote(isDenovo ~., data = trainData_combined, perc.over = 4)

# Separate balanced features and labels
#trainData <- balanced_train_data[, !names(balanced_train_data) %in% "isDenovo"]
#trainLabels <- balanced_train_data$isDenovo

set.seed(123) # For reproducibility
forest <- randomForest(isDenovo ~., data = balanced_train_data, ntree = 1000, importance = TRUE, do.trace = TRUE, localIMP = TRUE)
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

#Errorring, need to fix#
#Exploring Variable Importance Measures with mean_measures = "relevant trees"
#importance_frame_reltrees <- measure_importance(forest, mean_sample = "relevant trees")
#save(importance_frame, file = "MDS_randomForest_importance_frame.rda")
##load("MDS_randomForest_importance_frame.rda")
head(importance_frame, n = 10)

#Plot Multi-Way Importance Plot for features that have 
# minimum number of trees being 30 or more, size of points reflects number of nodes split on variable
plot_multi_way_importance(importance_frame, size_measure = "no_of_nodes", min_no_of_trees = 30)

#Multi-Way Importance plot, x is accuracy decrease, y is ginidecrease
plot_multi_way_importance(importance_frame, x_measure = "accuracy_decrease", y_measure = "gini_decrease", size_measure = "p_value")

#Compare measures using ggpairs, method is look for 3 comparisons that least agree with each other 
#then use these 3 measures in multi-way importance plot to select top variables
plot_importance_ggpairs(importance_frame)

#using ggpairs function to plot rankings instead of raw measures
plot_importance_rankings(importance_frame)

#Conditional Minimal Depth
(vars <- important_variables(importance_frame, k = 20, measures = c("mean_min_depth", "no_of_trees")))

interactions_frame <- min_depth_interactions(forest, vars)

head(interactions_frame[order(interactions_frame$occurrences, decreasing = TRUE), ])

save(interactions_frame, file = "RF_Interactions_frame.rda")
load("RF_Interactions_frame.rda")
head(interactions_frame[order(interactions_frame$occurrences, decreasing = TRUE), ])

#Explaining the forest
#explain_forest(forest, interactions = TRUE, data, pred_grid = 80)

plot_min_depth_interactions(interactions_frame)
#Prediction
preds <- predict(forest, testData_combined, type = "prob")[,2]

#SHAP
#*Uses kernelshap package*
#explainer <- kernelshap::kernelshap(forest, X = testData, bg_X = trainData, nsamples = 100)

#Save the explainer object
#saveRDS(explainer, file = "SHAP.rds")

#Load the explainer object back into R 
#explainer <- readRDS("SHAP.rds")

# Extract SHAP values for both classes
#shap_values <- explainer$S 

# Combine SHAP values with the true labels
#combined_shap <- cbind(as.data.frame(shap_values), priorMDS = as.numeric(testLabels) - 1)

#melted_shap <- pivot_longer(combined_shap, cols = -priorMDS, names_to = "Feature", values_to = "SHAP_Value")

#Plotting SHAP
#library(ggplot2)
# Create a function to generate colors based on SHAP values
#generate_color <- function(value) {
#  if (value < 0) {
#    return("#800080")  # Purple for negative SHAP values
#  } else if (value > 0) {
#    return("#FFD700")  # Gold for positive SHAP values
#  } else {
#    return("#FF0000")  # Red for zero SHAP values
#  }
#}

# Modify the generate_color function to handle vectors of SHAP values
#generate_color <- function(values) {
#  ifelse(values < 0, "#800080", ifelse(values > 0, "#FFD700", "#FF0000"))
#}

# Apply the color function to create a new column for colors
#melted_shap <- melted_shap %>%
#  mutate(Color = generate_color(SHAP_Value))

# Plotting bee swarm plot with adjusted colors and flipped vertically
#ggplot(melted_shap, aes(y = SHAP_Value, x = Feature, color = Color)) +
#  geom_point(position = position_jitter(height = 0.2), alpha = 0.5) +
#  scale_color_identity(guide = "none" ) +  # Hide legend for color
#  coord_flip() +  # Flip the plot vertically
#  theme_minimal() +
#  ggtitle("Bee Swarm Plot of SHAP Values for Both Classes") +
#  labs(color = "SHAP Value Color") +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#  scale_y_continuous(limits = c(-2, 2))  # Adjust y-axis limits as needed

#SHAP
#*Uses treeshap package
#testLabels <- as.factor(testLabels)

#Unify the Model
unified <- unify(forest, testData)

#Making shap values
shap_values <- treeshap(unified, testData)

#Making to Shapviz object
shp <- shapviz(shap_values, x = testData_combined)

#Plotting Shap
sv_importance(shp)
sv_importance(shp, kind = "bee")

# Plot ROC curve
roc_curve <- roc(testLabels, preds)
roc_plot <- plot(roc_curve, col = "blue", main = "ROC Curve for Random Forest Model")
par(roc_plot = "s")

#Calculating the Area
auc(roc_curve)

#Plotting Decision Tree
# Extract a single tree from the Random Forest model
tree <- getTree(forest, k = 1, labelVar = TRUE)

#Errorring Here#
reprtree:::plot.getTree(forest)

#SHAP
#*Uses iml package*
#*Not working, did not use*
# Convert train labels and data to numeric
#trainData <- data.frame(sapply(trainData, as.numeric))
#trainLabels <- as.numeric(trainLabels) - 1

# Create an iml Predictor object
#predictor <- Predictor$new(forest, data = trainData, y = trainLabels)

# Calculate SHAP values
#shapley <- Shapley$new(predictor, x.interest = trainData[1, ])

# For all instances in the training set
#shapley_all <- FeatureImp$new(predictor, loss = "mse")

# Plot SHAP values using plot function
#plot(shapley_all, type = "bar", main = "SHAP Values for Random Forest Model")




#Predicting with new dataset
library(tidyverse)
library(caret)

#Load in dataset
prediction <- read.delim("/Volumes/T7/IDH_4970_I_and_II/B_Cell_Analysis/RandomForest_Model/Meta_forMLModel.txt", header = TRUE, sep = "\t")

load("MDS_randomForest.rda")
#Separate test features and target feature
test_features <- prediction[, -which(names(prediction) == "priorMDS")]

#Labeling target feature (prior MDS) as 1 or 0 for y or n
target <- prediction %>%
  select(priorMDS) %>%
  mutate(priorMDS = ifelse(priorMDS == "y", 1, 0))

#Converting target feature to matrix and then factor
target <- as.factor(as.matrix(target))

# Make predictions using the random forest model
pred_new_meta <- predict(forest, newdata = test_features)

#Creating a confusion matrix of the data
confusionMatrix(pred_new_meta, target)







