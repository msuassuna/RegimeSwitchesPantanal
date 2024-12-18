
########################################################################
########################################################################
##
## 7: Forecasting
##
########################################################################
########################################################################

predictionTable <- merge(climaClasse, Lu_Data$PAST, by.x = "Year", by.y = "Time", all.x = TRUE)
colnames(predictionTable)[8] <- "pasture"

predictionTable <- merge(predictionTable, Lu_Data$CULT, by.x = "Year", by.y = "Time", all.x = TRUE)
names(predictionTable)[9] <- "culture"
names(predictionTable)[7] <- "Class_Before"
predictionTable$Class <- as.factor(predictionTable$Class)
predictionTable$Class_Before <- as.factor(predictionTable$Class_Before)

# Prepara data without land use data
data <- predictionTable %>% filter(complete.cases(Class_Before))

# Define the input variables and the target variable
inputVars <- names(data)[c(4:7)]
targetVar <- "Class"

# Select and filter the data table
data <- data %>% dplyr::select(c(inputVars, targetVar)) %>% filter(complete.cases(.))
data[,inputVars[-which(inputVars == "Class_Before")]] <- apply(data[,inputVars[-which(inputVars == "Class_Before")]], 2, scale)

# Cheack the data
head(data)
str(data)

# Split the data into training and testing sets
set.seed(123)
trainIndex <- createDataPartition(data$Class, p = 0.7, list = FALSE)
trainData <- data[trainIndex, ]
testData <- data[-trainIndex, ]


set.seed(123) # For reproducibility
control <- rfeControl(method = "cv", number = 10, verbose = FALSE)
model <- rfe(x = trainData[,inputVars],
             y = trainData$Class, sizes = c(1:6), rfeControl = control, method = "rf")

# Print the selected features
predictors(model)

set.seed(123) # For reproducibility
model_rf <- train(x = trainData[,inputVars], y = trainData$Class,
                  method = "rf",
                  trControl = trainControl(method = "cv", number = 10),
                  tuneGrid = expand.grid(mtry = 2:4))

# Model Evaluation
# Evaluate the performance of the model
predictions <- predict(model_rf, testData[,inputVars])
confusionMatrix(predictions, testData$Class)

# Partial dependence functions
par.culture <- partial(model_rf, pred.var = "culture", chull = TRUE)
plot.culture  <- autoplot(par.culture , contour = TRUE)

par.amo <- partial(model_rf, pred.var = "amo", chull = TRUE)
plot.amo  <- autoplot(par.amo , contour = TRUE)

par.pdo <- partial(model_rf, pred.var = "pdo", chull = TRUE)
plot.pdo  <- autoplot(par.pdo , contour = TRUE)

par.pasture <- partial(model_rf, pred.var = "pasture", chull = TRUE)
plot.pasture  <- autoplot(par.pasture , contour = TRUE)

par.class <- partial(model_rf, pred.var = "Class_Before", chull = TRUE)
plot.class  <- autoplot(par.class , contour = TRUE)

par.nino <- partial(model_rf, pred.var = "nino", chull = TRUE)
plot.nino  <- autoplot(par.nino , contour = TRUE)

par.pdo.enso <- partial(model_rf, pred.var = c("pdo","nino"), chull = TRUE)
plot.pdo.enso  <- autoplot(par.pdo.enso , contour = TRUE)

par.pdo.class <- partial(model_rf, pred.var = c("pdo","Class_Before"), chull = TRUE)
plot.pdo.class  <- autoplot(par.pdo.class , contour = TRUE)

par.nino.class <- partial(model_rf, pred.var = c("nino","Class_Before"), chull = TRUE)
plot.nino.class  <- autoplot(par.nino.class , contour = TRUE)

grid.arrange(plot.nino,
             plot.pdo,
             plot.amo,
             plot.pasture,
             plot.culture,
             plot.class,
             plot.pdo.enso,
             plot.pdo.class,
             plot.pasture.class)

grid.arrange(plot.nino,
             plot.pdo,
             plot.amo,
             plot.class,
             plot.pdo.enso,
             plot.pdo.class)














# Visualize the performance of the model
ggplot(testData, aes(x = Class, fill = predictions)) +
  geom_bar(position = "dodge") +
  theme_classic() + 
  abs(title = "Performance of Random Forest Model",
      x = "True Class",
      y = "Count", fill = "Predicted Class")

















# Multinomial logistic regression with feature selection
stepAICModel <- stepAIC(multinom(as.formula(paste(targetVar, "~", paste(inputVars, collapse = "+"))), data = trainData), direction = "both")
stepAICPredictions <- predict(stepAICModel, testData)
stepAICAccuracy <- mean(stepAICPredictions == testData$Class)
multinom_Variables <- stepAICModel$coefnames[-1]

# Random forest with feature selection
rfModel <- randomForest(as.formula(paste(targetVar, "~", paste(inputVars, collapse = "+"))), data = trainData, ntree = 500)
rfPredictions <- predict(rfModel, testData)
rfAccuracy <- mean(rfPredictions == testData$Class)
rf_Variables <- rownames(importance(rfModel))

# Support vector machine with feature selection
model <- svm(as.formula(paste(targetVar, "~", paste(inputVars, collapse = "+"))), 
             data = trainData, 
             kernel = "linear", # or "radial" for radial kernel
             cost = 1, 
             gamma = 1) 

# Evaluate the importance of different variables
importance <- numeric(length(inputVars))
for (i in seq_along(inputVars)) {
  vars <- setdiff(inputVars, inputVars[i])
  model_i <- svm(as.formula(paste(targetVar, "~", paste(vars, collapse = "+"))), 
                 data = testData, 
                 kernel = "linear", # or "radial" for radial kernel
                 cost = 1, 
                 gamma = 1) 
  importance[i] <- 1 - model_i$tot.nSV / model$tot.nSV
}

# Print the variable importance measures
names(importance) <- inputVars
print(importance)
svmModel <- svm(as.formula(paste(targetVar, "~", paste(inputVars[c(1,2,3,4)], collapse = "+"))), 
                data = testData, 
                kernel = "linear", # or "radial" for radial kernel
                cost = 1, 
                gamma = 1) 

svmPredictions <- predict(svmModel, testData)
svmAccuracy <- mean(svmPredictions == testData$Class)
svm_Variables <- inputVars[c(1,2,3,4)]

# Naive Bayes model
nb_model <- naiveBayes(as.formula(paste(targetVar, "~.")), data = trainData)
nb_varimp <- varSelRF(trainData[, -which(names(trainData) == targetVar)], trainData[[targetVar]], ntree = 200)
nb_varimp_df <- data.frame(varimp = nb_varimp$initialImportances,
                           variables = row.names(nb_varimp$initialImportances))
nb_Variables <- nb_varimp$selected.vars

# KNN model
knn_model <- train(as.formula(paste(targetVar, "~.")), method = "knn", data = trainData,
                   trControl = trainControl(method = "cv", number = 5))

rfe <- rfe(trainData[, -ncol(trainData)], trainData[, "Class"],
           sizes = c(1:ncol(trainData)-1),
           rfeControl = rfeControl(functions = caretFuncs,
                                   method = "cv",
                                   number = 10))

# Fit the k-NN model on the selected variables
knn_Variables <- rfe$optVariables

# Decision Tree model
dt_model <- rpart(as.formula(paste(targetVar, "~.")), data = trainData, method = "class")
dt_varimp <- varImp(dt_model)
dt_varimp_df <- data.frame(varimp = dt_varimp, variables = row.names(dt_varimp))

# Feature selection using recursive feature elimination
ctrl <- rfeControl(functions = caretFuncs, method = "cv", number = 10)
model_rfe <- rfe(trainData[, -ncol(trainData)], trainData$Class,
                 sizes = 1:6, rfeControl = ctrl)
print(model_rfe)
model_rfe$perfNames
dt_Variables <- inputVars[4]

# Gradient Boosting model
set.seed(123)
gbm_model <- gbm(as.formula(paste(targetVar, "~.")), data = trainData,
                 n.trees = 100, interaction.depth = 3, distribution = "multinomial")
summary(gbm_model)
gbm_Variables <- inputVars[c(1,2,3,4)]

# Neural Network model
nn_model <- nnet(as.formula(paste(targetVar, "~.")), data = trainData, size = 5, decay = 0.1)
nn_varimp <- varImp(nn_model)
nn_Variables <- inputVars[c(1,2,3,4)]

###
# Summaries
###
variables <- data.frame(Name = c("Multinomial", "Random Forest",
                                 "Support Vector Machine", "Naive Bayes",
                                 "KNN", "Decision Tree",
                                 "Gradient boosting","Neural Network"),
                        Variables = c(paste(targetVar, "~", paste(multinom_Variables, collapse = "+")),
                                      paste(targetVar, "~", paste(rf_Variables, collapse = "+")),
                                      paste(targetVar, "~", paste(svm_Variables, collapse = "+")),
                                      paste(targetVar, "~", paste(nb_Variables, collapse = "+")),
                                      paste(targetVar, "~", paste(knn_Variables, collapse = "+")),
                                      paste(targetVar, "~", paste(dt_Variables, collapse = "+")),
                                      paste(targetVar, "~", paste(gbm_Variables, collapse = "+")),
                                      paste(targetVar, "~", paste(nn_Variables, collapse = "+")))) %>%
  arrange(Name)
variables

# Evaluate the performance of the model on the test set
# Train the models on the selected variables
multinomialModel <- multinom(as.formula(paste(targetVar, "~", paste(inputVars, collapse = "+"))), data = trainData)
rfModel <- randomForest(as.formula(paste(targetVar, "~", paste(inputVars, collapse = "+"))), data = trainData, ntree = 500)
svmModel <- svm(as.formula(paste(targetVar, "~", paste(inputVars[c(1,2,3,4)], collapse = "+"))), 
                data = trainData, kernel = "linear", cost = 1, gamma = 1) 
nb_model <- naiveBayes(as.formula(paste(targetVar, "~.")), data = trainData)
knn_model <- train(as.formula(paste(targetVar, "~.")), method = "knn", data = trainData,
                   trControl = trainControl(method = "cv", number = 5))
dt_model <- rpart(as.formula(paste(targetVar, "~.")), data = trainData, method = "class")
gbm_model <- gbm(as.formula(paste(targetVar, "~.")), data = trainData,
                 n.trees = 100, interaction.depth = 3, distribution = "multinomial")
nn_model <- nnet(as.formula(paste(targetVar, "~.")), data = trainData, size = 5, decay = 0.1)


# Evaluate the performance of the models on the test set
multinom_pred <- predict(multinomialModel, newdata = testData)
rf_pred <- predict(rfModel, newdata = testData)
svm_pred <- predict(svmModel, newdata = testData)
nb_pred <- predict(nb_model, newdata = testData)
knn_pred <- predict(knn_model, newdata = testData)
dt_pred <- apply(predict(dt_model, newdata = testData), 1, which.max)
gbm_pred <- apply(predict(gbm_model, newdata = testData), 1, which.max) 
nn_pred <- apply(predict(nn_model, newdata = testData), 1, which.max) 

# Estimate accuracies
multinom_Accuracy <- mean(multinom_pred == testData$Class)
rf_Accuracy <- mean(rf_pred == testData$Class)
svm_Accuracy <- mean(svm_pred == testData$Class)
nb_Accuracy <- mean(nb_pred == testData$Class)
knn_Accuracy <- mean(knn_pred == testData$Class)
dt_Accuracy <- mean(dt_pred == testData$Class)
gbm_Accuracy <- mean(gbm_pred == testData$Class)
nn_Accuracy <- mean(nn_pred == testData$Class)

# Evaluate the performance of the model on the train set
# Train the models on the selected variables
multinomialModel <- multinom(as.formula(paste(targetVar, "~", paste(inputVars, collapse = "+"))), data = trainData)
rfModel <- randomForest(as.formula(paste(targetVar, "~", paste(inputVars, collapse = "+"))), data = trainData, ntree = 500)
svmModel <- svm(as.formula(paste(targetVar, "~", paste(inputVars[c(1,2,3,4)], collapse = "+"))), 
                data = trainData, kernel = "linear", cost = 1, gamma = 1) 
nb_model <- naiveBayes(as.formula(paste(targetVar, "~.")), data = trainData)
knn_model <- train(as.formula(paste(targetVar, "~.")), method = "knn", data = trainData,
                   trControl = trainControl(method = "cv", number = 5))
dt_model <- rpart(as.formula(paste(targetVar, "~.")), data = trainData, method = "class")
gbm_model <- gbm(as.formula(paste(targetVar, "~.")), data = trainData,
                 n.trees = 100, interaction.depth = 3, distribution = "multinomial")
nn_model <- nnet(as.formula(paste(targetVar, "~.")), data = trainData, size = 5, decay = 0.1)

# Evaluate the performance of the models on the training set
multinom_pred_train <- predict(multinomialModel, newdata = trainData)
rf_pred_train <- predict(rfModel, newdata = trainData)
svm_pred_train <- predict(svmModel, newdata = trainData)
nb_pred_train <- predict(nb_model, newdata = trainData)
knn_pred_train <- predict(knn_model, newdata = trainData)
dt_pred_train <- apply(predict(dt_model, newdata = trainData), 1, which.max)
gbm_pred_train <- apply(predict(gbm_model, newdata = trainData), 1, which.max) 
nn_pred_train <- apply(predict(nn_model, newdata = trainData), 1, which.max) 

# Estimate accuracies
multinom_Accuracy_train <- mean(multinom_pred_train == trainData$Class)
rf_Accuracy_train <- mean(rf_pred_train == trainData$Class)
svm_Accuracy_train <- mean(svm_pred_train == trainData$Class)
nb_Accuracy_train <- mean(nb_pred_train == trainData$Class)
knn_Accuracy_train <- mean(knn_pred_train == trainData$Class)
dt_Accuracy_train <- mean(dt_pred_train == trainData$Class)
gbm_Accuracy_train <- mean(gbm_pred_train == trainData$Class)
nn_Accuracy_train <- mean(nn_pred_train == trainData$Class)

Accuarcy_Table <- data.frame(Name = c("Multinomial", "Random Forest",
                                      "Support Vector Machine", "Naive Bayes",
                                      "KNN", "Decision Tree",
                                      "Gradient boosting","Neural Network"),
                             Accuracy = c(multinom_Accuracy,rf_Accuracy,
                                          svm_Accuracy,nb_Accuracy,
                                          knn_Accuracy,dt_Accuracy,gbm_Accuracy,
                                          nn_Accuracy),
                             Accuracy_train = c(multinom_Accuracy_train,rf_Accuracy_train,
                                                svm_Accuracy_train,nb_Accuracy_train,
                                                knn_Accuracy_train,dt_Accuracy_train,gbm_Accuracy_train,
                                                nn_Accuracy_train)) %>%
  arrange(Accuracy)

Accuarcy_Table$OF <- Accuarcy_Table$Accuracy_train - Accuarcy_Table$Accuracy

Accuarcy_Table %>% arrange(OF)

########################################################################
########################################################################
##
## 7: Forecasting with Land Use
##
########################################################################
########################################################################

# Prepara data without land use data
data <- predictionTable %>% filter(complete.cases(Class_Before))

# Define the input variables and the target variable
inputVars <- names(data)[c(4:9)]
targetVar <- "Class"

# Select and filter the data table
data <- data %>% dplyr::select(c(inputVars, targetVar)) %>% filter(complete.cases(.))

# Cheack the data
head(data)
str(data)

# Split the data into training and testing sets
set.seed(123)
trainIndex <- createDataPartition(data$Class, p = 0.7, list = FALSE)
trainData <- data[trainIndex, ]
testData <- data[-trainIndex, ]

# Multinomial logistic regression with feature selection
stepAICModel <- stepAIC(multinom(as.formula(paste(targetVar, "~", paste(inputVars, collapse = "+"))), data = trainData), direction = "both")
stepAICPredictions <- predict(stepAICModel, testData)
stepAICAccuracy <- mean(stepAICPredictions == testData$Class)
multinom_Variables <- stepAICModel$coefnames[-1]

# Random forest with feature selection
rfModel <- randomForest(as.formula(paste(targetVar, "~", paste(inputVars, collapse = "+"))), data = trainData, ntree = 500)
rfPredictions <- predict(rfModel, testData)
rfAccuracy <- mean(rfPredictions == testData$Class)
rf_Variables <- rownames(importance(rfModel))

# Support vector machine with feature selection
model <- svm(as.formula(paste(targetVar, "~", paste(inputVars, collapse = "+"))), 
             data = trainData, 
             kernel = "linear", # or "radial" for radial kernel
             cost = 1, 
             gamma = 1) 

# Evaluate the importance of different variables
importance <- numeric(length(inputVars))
for (i in seq_along(inputVars)) {
  vars <- setdiff(inputVars, inputVars[i])
  model_i <- svm(as.formula(paste(targetVar, "~", paste(vars, collapse = "+"))), 
                 data = testData, 
                 kernel = "linear", # or "radial" for radial kernel
                 cost = 1, 
                 gamma = 1) 
  importance[i] <- 1 - model_i$tot.nSV / model$tot.nSV
}

# Print the variable importance measures
names(importance) <- inputVars
print(importance)
svmModel <- svm(as.formula(paste(targetVar, "~", paste(inputVars[c(1,2,3,4)], collapse = "+"))), 
                data = testData, 
                kernel = "linear", # or "radial" for radial kernel
                cost = 1, 
                gamma = 1) 

svmPredictions <- predict(svmModel, testData)
svmAccuracy <- mean(svmPredictions == testData$Class)
svm_Variables <- inputVars

# Naive Bayes model
nb_model <- naiveBayes(as.formula(paste(targetVar, "~.")), data = trainData)
nb_varimp <- varSelRF(trainData[, -which(names(trainData) == targetVar)], trainData[[targetVar]], ntree = 200)
nb_varimp_df <- data.frame(varimp = nb_varimp$initialImportances,
                           variables = row.names(nb_varimp$initialImportances))
nb_Variables <- nb_varimp$selected.vars

# KNN model
knn_model <- train(as.formula(paste(targetVar, "~.")), method = "knn", data = trainData,
                   trControl = trainControl(method = "cv", number = 5))

rfe <- rfe(trainData[, -ncol(trainData)], trainData[, "Class"],
           sizes = c(1:ncol(trainData)-1),
           rfeControl = rfeControl(functions = caretFuncs,
                                   method = "cv",
                                   number = 10))

# Fit the k-NN model on the selected variables
knn_Variables <- rfe$optVariables

# Decision Tree model
dt_model <- rpart(as.formula(paste(targetVar, "~.")), data = trainData, method = "class")
dt_varimp <- varImp(dt_model)
dt_varimp_df <- data.frame(varimp = dt_varimp, variables = row.names(dt_varimp))

# Feature selection using recursive feature elimination
ctrl <- rfeControl(functions = caretFuncs, method = "cv", number = 10)
model_rfe <- rfe(trainData[, -ncol(trainData)], trainData$Class,
                 sizes = 1:6, rfeControl = ctrl)
print(model_rfe)
model_rfe$perfNames
dt_Variables <- inputVars[4]

# Gradient Boosting model
set.seed(123)
gbm_model <- gbm(as.formula(paste(targetVar, "~.")), data = trainData,
                 n.trees = 100, interaction.depth = 3, distribution = "multinomial")
summary(gbm_model)
gbm_Variables <- inputVars[c(1,2,3,4)]

# Neural Network model
nn_model <- nnet(as.formula(paste(targetVar, "~.")), data = trainData, size = 5, decay = 0.1)
nn_varimp <- varImp(nn_model)
nn_Variables <- inputVars[c(1,2,3,4)]

###
# Summaries
###
variables <- data.frame(Name = c("Multinomial", "Random Forest",
                                 "Support Vector Machine", "Naive Bayes",
                                 "KNN", "Decision Tree",
                                 "Gradient boosting","Neural Network"),
                        Variables = c(paste(targetVar, "~", paste(multinom_Variables, collapse = "+")),
                                      paste(targetVar, "~", paste(rf_Variables, collapse = "+")),
                                      paste(targetVar, "~", paste(svm_Variables, collapse = "+")),
                                      paste(targetVar, "~", paste(nb_Variables, collapse = "+")),
                                      paste(targetVar, "~", paste(knn_Variables, collapse = "+")),
                                      paste(targetVar, "~", paste(dt_Variables, collapse = "+")),
                                      paste(targetVar, "~", paste(gbm_Variables, collapse = "+")),
                                      paste(targetVar, "~", paste(nn_Variables, collapse = "+")))) %>%
  arrange(Name)
variables

# Evaluate the performance of the model on the test set
# Train the models on the selected variables
multinomialModel <- multinom(as.formula(paste(targetVar, "~", paste(inputVars, collapse = "+"))), data = trainData)
rfModel <- randomForest(as.formula(paste(targetVar, "~", paste(inputVars, collapse = "+"))), data = trainData, ntree = 500)
svmModel <- svm(as.formula(paste(targetVar, "~", paste(inputVars, collapse = "+"))), 
                data = trainData, kernel = "linear", cost = 1, gamma = 1) 
nb_model <- naiveBayes(as.formula(paste(targetVar, "~.")), data = trainData)
knn_model <- train(as.formula(paste(targetVar, "~.")), method = "knn", data = trainData,
                   trControl = trainControl(method = "cv", number = 5))
dt_model <- rpart(as.formula(paste(targetVar, "~.")), data = trainData, method = "class")
gbm_model <- gbm(as.formula(paste(targetVar, "~.")), data = trainData,
                 n.trees = 100, interaction.depth = 3, distribution = "multinomial")
nn_model <- nnet(as.formula(paste(targetVar, "~.")), data = trainData, size = 5, decay = 0.1)


# Evaluate the performance of the models on the test set
multinom_pred <- predict(multinomialModel, newdata = testData)
rf_pred <- predict(rfModel, newdata = testData)
svm_pred <- predict(svmModel, newdata = testData)
nb_pred <- predict(nb_model, newdata = testData)
knn_pred <- predict(knn_model, newdata = testData)
dt_pred <- apply(predict(dt_model, newdata = testData), 1, which.max)
gbm_pred <- apply(predict(gbm_model, newdata = testData), 1, which.max) 
nn_pred <- apply(predict(nn_model, newdata = testData), 1, which.max) 

# Estimate accuracies
multinom_Accuracy <- mean(multinom_pred == testData$Class)
rf_Accuracy <- mean(rf_pred == testData$Class)
svm_Accuracy <- mean(svm_pred == testData$Class)
nb_Accuracy <- mean(nb_pred == testData$Class)
knn_Accuracy <- mean(knn_pred == testData$Class)
dt_Accuracy <- mean(dt_pred == testData$Class)
gbm_Accuracy <- mean(gbm_pred == testData$Class)
nn_Accuracy <- mean(nn_pred == testData$Class)

# Evaluate the performance of the model on the test set
# Train the models on the selected variables
multinomialModel <- multinom(as.formula(paste(targetVar, "~", paste(inputVars, collapse = "+"))), data = trainData)
rfModel <- randomForest(as.formula(paste(targetVar, "~", paste(inputVars, collapse = "+"))), data = trainData, ntree = 500)
svmModel <- svm(as.formula(paste(targetVar, "~", paste(inputVars[c(1,2,3,4)], collapse = "+"))), 
                data = trainData, kernel = "linear", cost = 1, gamma = 1) 
nb_model <- naiveBayes(as.formula(paste(targetVar, "~.")), data = trainData)
knn_model <- train(as.formula(paste(targetVar, "~.")), method = "knn", data = trainData,
                   trControl = trainControl(method = "cv", number = 5))
dt_model <- rpart(as.formula(paste(targetVar, "~.")), data = trainData, method = "class")
gbm_model <- gbm(as.formula(paste(targetVar, "~.")), data = trainData,
                 n.trees = 100, interaction.depth = 3, distribution = "multinomial")
nn_model <- nnet(as.formula(paste(targetVar, "~.")), data = trainData, size = 5, decay = 0.1)


# Evaluate the performance of the models on the test set
multinom_pred_train <- predict(multinomialModel, newdata = trainData)
rf_pred_train <- predict(rfModel, newdata = trainData)
svm_pred_train <- predict(svmModel, newdata = trainData)
nb_pred_train <- predict(nb_model, newdata = trainData)
knn_pred_train <- predict(knn_model, newdata = trainData)
dt_pred_train <- apply(predict(dt_model, newdata = trainData), 1, which.max)
gbm_pred_train <- apply(predict(gbm_model, newdata = trainData), 1, which.max) 
nn_pred_train <- apply(predict(nn_model, newdata = trainData), 1, which.max) 

# Estimate accuracies
multinom_Accuracy_train <- mean(multinom_pred_train == trainData$Class)
rf_Accuracy_train <- mean(rf_pred_train == trainData$Class)
svm_Accuracy_train <- mean(svm_pred_train == trainData$Class)
nb_Accuracy_train <- mean(nb_pred_train == trainData$Class)
knn_Accuracy_train <- mean(knn_pred_train == trainData$Class)
dt_Accuracy_train <- mean(dt_pred_train == trainData$Class)
gbm_Accuracy_train <- mean(gbm_pred_train == trainData$Class)
nn_Accuracy_train <- mean(nn_pred_train == trainData$Class)

Accuarcy_Table_LU <- data.frame(Name = c("Multinomial", "Random Forest",
                                         "Support Vector Machine", "Naive Bayes",
                                         "KNN", "Decision Tree",
                                         "Gradient boosting","Neural Network"),
                                Accuracy = c(multinom_Accuracy,rf_Accuracy,
                                             svm_Accuracy,nb_Accuracy,
                                             knn_Accuracy,dt_Accuracy,gbm_Accuracy,
                                             nn_Accuracy),
                                Accuracy_train = c(multinom_Accuracy_train,rf_Accuracy_train,
                                                   svm_Accuracy_train,nb_Accuracy_train,
                                                   knn_Accuracy_train,dt_Accuracy_train,gbm_Accuracy_train,
                                                   nn_Accuracy_train)) %>%
  arrange(Accuracy)

Accuarcy_Table_LU$OF <- Accuarcy_Table_LU$Accuracy_train - Accuarcy_Table_LU$Accuracy

write.csv2(
  cbind(
    Accuarcy_Table_LU %>% arrange(Accuracy),
    Accuarcy_Table %>% arrange(Accuracy)
  ),
  "Accuracies.csv", row.names = FALSE
)


########################################################################
########################################################################
##
## 8: Partial Dependence Plots
##
########################################################################
########################################################################

rfModel <- randomForest(as.formula(paste(targetVar, "~", paste(inputVars, collapse = "+"))), data = trainData, ntree = 500)
svmModel <- svm(as.formula(paste(targetVar, "~", paste(inputVars[c(1,2,3,4)], collapse = "+"))), 
                data = trainData, kernel = "linear", cost = 1, gamma = 1) 
nb_model <- naiveBayes(as.formula(paste(targetVar, "~.")), data = trainData)
knn_model <- train(as.formula(paste(targetVar, "~.")), method = "knn", data = trainData,
                   trControl = trainControl(method = "cv", number = 5))
dt_model <- rpart(as.formula(paste(targetVar, "~.")), data = trainData, method = "class")
gbm_model <- gbm(as.formula(paste(targetVar, "~.")), data = trainData,
                 n.trees = 100, interaction.depth = 3, distribution = "multinomial")
nn_model <- nnet(as.formula(paste(targetVar, "~.")), data = trainData, size = 5, decay = 0.1)



par.nino <- partial(multinomialModel, pred.var = "nino", chull = TRUE)
plot.nino  <- autoplot(par.nino , contour = TRUE)

par.pdo <- partial(multinomialModel, pred.var = "pdo", chull = TRUE)
plot.pdo  <- autoplot(par.pdo , contour = TRUE)

par.amo <- partial(multinomialModel, pred.var = "amo", chull = TRUE)
plot.amo  <- autoplot(par.amo , contour = TRUE)

par.pasture <- partial(multinomialModel, pred.var = "pasture", chull = TRUE)
plot.pasture  <- autoplot(par.pasture , contour = TRUE)

par.culture <- partial(multinomialModel, pred.var = "culture", chull = TRUE)
plot.culture  <- autoplot(par.culture , contour = TRUE)

par.class <- partial(multinomialModel, pred.var = "Class_Before", chull = TRUE)
plot.class  <- autoplot(par.class , contour = TRUE)

par.pdo.enso <- partial(multinomialModel, pred.var = c("pdo","nino"), chull = TRUE)
plot.pdo.enso  <- autoplot(par.pdo.enso , contour = TRUE)

par.pdo.class <- partial(multinomialModel, pred.var = c("pdo","Class_Before"), chull = TRUE)
plot.pdo.class  <- autoplot(par.pdo.class , contour = TRUE)

par.nino.class <- partial(multinomialModel, pred.var = c("nino","Class_Before"), chull = TRUE)
plot.nino.class  <- autoplot(par.nino.class , contour = TRUE)

par.pasture.class <- partial(multinomialModel, pred.var = c("pasture","Class_Before"), chull = TRUE)
plot.pasture.class  <- autoplot(par.pasture.class , contour = TRUE)

grid.arrange(plot.nino,
             plot.pdo,
             plot.amo,
             plot.pasture,
             plot.culture,
             plot.class,
             plot.pdo.enso,
             plot.pdo.class,
             plot.pasture.class)


par.nino <- partial(rfModel, pred.var = "nino", chull = TRUE)
plot.nino  <- autoplot(par.nino , contour = TRUE)

par.pdo <- partial(rfModel, pred.var = "pdo", chull = TRUE)
plot.pdo  <- autoplot(par.pdo , contour = TRUE)

par.amo <- partial(rfModel, pred.var = "amo", chull = TRUE)
plot.amo  <- autoplot(par.amo , contour = TRUE)

par.pasture <- partial(rfModel, pred.var = "pasture", chull = TRUE)
plot.pasture  <- autoplot(par.pasture , contour = TRUE)

par.culture <- partial(rfModel, pred.var = "culture", chull = TRUE)
plot.culture  <- autoplot(par.culture , contour = TRUE)

par.class <- partial(rfModel, pred.var = "Class_Before", chull = TRUE)
plot.class  <- autoplot(par.class , contour = TRUE)

par.pdo.enso <- partial(rfModel, pred.var = c("pdo","nino"), chull = TRUE)
plot.pdo.enso  <- autoplot(par.pdo.enso , contour = TRUE)

par.pdo.class <- partial(rfModel, pred.var = c("pdo","Class_Before"), chull = TRUE)
plot.pdo.class  <- autoplot(par.pdo.class , contour = TRUE)

par.nino.class <- partial(rfModel, pred.var = c("nino","Class_Before"), chull = TRUE)
plot.nino.class  <- autoplot(par.nino.class , contour = TRUE)

par.pasture.class <- partial(knn_model, pred.var = c("pasture","Class_Before"), chull = TRUE)
plot.pasture.class  <- autoplot(par.pasture.class , contour = TRUE)

grid.arrange(plot.nino,
             plot.pdo,
             plot.amo,
             plot.pasture,
             plot.culture,
             plot.class,
             plot.pdo.enso,
             plot.pdo.class,
             plot.pasture.class)





########################################################################
########################################################################
##
## 7: Forecasting the Amf
##
########################################################################
########################################################################

predictionTable <- merge(climaClasse, Lu_Data$PAST, by.x = "Year", by.y = "Time", all.x = TRUE)
colnames(predictionTable)[8] <- "pasture"

predictionTable <- merge(predictionTable, Lu_Data$CULT, by.x = "Year", by.y = "Time", all.x = TRUE)
names(predictionTable)[9] <- "culture"
names(predictionTable)[7] <- "Class_Before"
predictionTable$Class <- as.factor(predictionTable$Class)
predictionTable$Class_Before <- as.factor(predictionTable$Class_Before)
predictionTable$Amf_2 <- lag(predictionTable$Amf)

tail(predictionTable,10)

# Prepara data without land use data
data <- predictionTable %>% filter(complete.cases(Class_Before))
tail(data,10)

# Define the input variables and the target variable
inputVars <- names(data)[c(4:6,10)]
targetVar <- "Amf"

# Select and filter the data table
data <- data %>% dplyr::select(c(inputVars, targetVar)) %>% filter(complete.cases(.))
data[,inputVars] <- apply(data[,inputVars], 2, scale)

# Cheack the data
head(data)
str(data)

# Split the data into training and testing sets
set.seed(123)
trainIndex <- createDataPartition(data$Amf, p = 0.7, list = FALSE)
trainData <- data[trainIndex, ]
testData <- data[-trainIndex, ]


# Multinomial logistic regression with feature selection
lm_Model <- lm(as.formula(paste(targetVar, "~", paste(inputVars, collapse = "+"))), data = trainData)
lm_Predictions <- predict(lm_Model, testData)
lm_rmse <- RMSE(lm_Predictions, testData$Amf)
s_result <- summary(lm_Model, test = "F")
lm_Variables <- row.names(s_result$coefficients)[-1][s_result$coefficients[-1,4] < 0.05]

# Random forest with feature selection
rf_Model <- randomForest(as.formula(paste(targetVar, "~", paste(inputVars, collapse = "+"))), data = trainData, ntree = 500)
rf_Predictions <- predict(rf_Model, testData)
rf_rmse <- RMSE(rf_Predictions, testData$Amf)
plot(varImp(rf_Model))


# Support vector machine with feature selection
svm_model <- svm(as.formula(paste(targetVar, "~", paste(inputVars, collapse = "+"))),
                 data = trainData, 
                 kernel = "linear", # or "radial" for radial kernel
                 cost = 1, 
                 gamma = 1) 

# Evaluate the importance of different variables
importance <- numeric(length(inputVars))
for (i in seq_along(inputVars)) {
  vars <- setdiff(inputVars, inputVars[i])
  model_i <- svm(as.formula(paste(targetVar, "~", paste(vars, collapse = "+"))), 
                 data = testData, 
                 kernel = "linear", # or "radial" for radial kernel
                 cost = 1, 
                 gamma = 1) 
  importance[i] <- 1 - model_i$tot.nSV / model$tot.nSV
}

# Print the variable importance measures
names(importance) <- inputVars
print(importance)
svm_Model <- svm(as.formula(paste(targetVar, "~", paste(inputVars[c(1,2,3,4)], collapse = "+"))), 
                 data = testData, 
                 kernel = "linear", # or "radial" for radial kernel
                 cost = 1, 
                 gamma = 1) 

svm_Predictions <- predict(svm_Model, testData)
svm_rmse <- RMSE(svm_Predictions, testData$Amf)

# KNN model
knn_model <- train(as.formula(paste(targetVar, "~.")), method = "knn", data = trainData,
                   trControl = trainControl(method = "cv", number = 5))

rfe <- rfe(trainData[, -ncol(trainData)], trainData[, "Amf"],
           sizes = c(1:ncol(trainData)-1),
           rfeControl = rfeControl(functions = caretFuncs,
                                   method = "cv",
                                   number = 10))

# Fit the k-NN model on the selected variables
knn_Variables <- rfe$optVariables

knn_Predictions <- predict(knn_model, testData)
knn_rmse <- RMSE(knn_Predictions, testData$Amf)


# Decision Tree model
dt_model <- rpart(as.formula(paste(targetVar, "~.")), data = trainData, method = "class")
dt_varimp <- varImp(dt_model)
dt_varimp_df <- data.frame(varimp = dt_varimp, variables = row.names(dt_varimp))

dt_Predictions <- apply(predict(dt_model, testData), 1, which.max)
knn_rmse <- RMSE(knn_Predictions, testData$Amf)


# Feature selection using recursive feature elimination
ctrl <- rfeControl(functions = caretFuncs, method = "cv", number = 10)
model_rfe <- rfe(trainData[, -ncol(trainData)], trainData$Amf,
                 sizes = 1:6, rfeControl = ctrl)

print(model_rfe)
model_rfe$perfNames
dt_Variables <- inputVars[4]

# Gradient Boosting model
set.seed(123)
gbm_model <- gbm(as.formula(paste(targetVar, "~.")), data = trainData,
                 n.trees = 100, interaction.depth = 3, distribution = "multinomial")
summary(gbm_model)
gbm_Variables <- inputVars[c(1,2,3,4)]

# Neural Network model
nn_model <- nnet(as.formula(paste(targetVar, "~.")), data = trainData, size = 5, decay = 0.1)
nn_varimp <- varImp(nn_model)
nn_Variables <- inputVars[c(1,2,3,4)]

###
# Summaries
###
variables <- data.frame(Name = c("Multinomial", "Random Forest",
                                 "Support Vector Machine", "Naive Bayes",
                                 "KNN", "Decision Tree",
                                 "Gradient boosting","Neural Network"),
                        Variables = c(paste(targetVar, "~", paste(multinom_Variables, collapse = "+")),
                                      paste(targetVar, "~", paste(rf_Variables, collapse = "+")),
                                      paste(targetVar, "~", paste(svm_Variables, collapse = "+")),
                                      paste(targetVar, "~", paste(nb_Variables, collapse = "+")),
                                      paste(targetVar, "~", paste(knn_Variables, collapse = "+")),
                                      paste(targetVar, "~", paste(dt_Variables, collapse = "+")),
                                      paste(targetVar, "~", paste(gbm_Variables, collapse = "+")),
                                      paste(targetVar, "~", paste(nn_Variables, collapse = "+")))) %>%
  arrange(Name)
variables

# Evaluate the performance of the model on the test set
# Train the models on the selected variables
multinomialModel <- multinom(as.formula(paste(targetVar, "~", paste(inputVars, collapse = "+"))), data = trainData)
rfModel <- randomForest(as.formula(paste(targetVar, "~", paste(inputVars, collapse = "+"))), data = trainData, ntree = 500)
svmModel <- svm(as.formula(paste(targetVar, "~", paste(inputVars[c(1,2,3,4)], collapse = "+"))), 
                data = trainData, kernel = "linear", cost = 1, gamma = 1) 
nb_model <- naiveBayes(as.formula(paste(targetVar, "~.")), data = trainData)
knn_model <- train(as.formula(paste(targetVar, "~.")), method = "knn", data = trainData,
                   trControl = trainControl(method = "cv", number = 5))
dt_model <- rpart(as.formula(paste(targetVar, "~.")), data = trainData, method = "class")
gbm_model <- gbm(as.formula(paste(targetVar, "~.")), data = trainData,
                 n.trees = 100, interaction.depth = 3, distribution = "gaussian")
nn_model <- nnet(as.formula(paste(targetVar, "~.")), data = trainData, size = 5, decay = 0.1)


# Evaluate the performance of the models on the test set
multinom_pred <- predict(multinomialModel, newdata = testData)
rf_pred <- predict(rfModel, newdata = testData)
svm_pred <- predict(svmModel, newdata = testData)
nb_pred <- predict(nb_model, newdata = testData)
knn_pred <- predict(knn_model, newdata = testData)
dt_pred <- apply(predict(dt_model, newdata = testData), 1, which.max)
gbm_pred <- predict(gbm_model, newdata = testData)
nn_pred <- predict(nn_model, newdata = testData)

# Estimate accuracies
multinom_Accuracy <- RMSE(as.numeric(as.character(multinom_pred)),testData$Amf)
rf_Accuracy <- RMSE(rf_pred,testData$Amf)
svm_Accuracy <- RMSE(svm_pred,testData$Amf)
nb_Accuracy <- RMSE(nb_pred,testData$Amf)
knn_Accuracy <- RMSE(knn_pred,testData$Amf)
dt_Accuracy <- RMSE(dt_pred,testData$Amf)
gbm_Accuracy <- RMSE(gbm_pred,testData$Amf)
nn_Accuracy <- RMSE(nn_pred,testData$Amf)

# Evaluate the performance of the model on the train set
# Train the models on the selected variables
multinomialModel <- multinom(as.formula(paste(targetVar, "~", paste(inputVars, collapse = "+"))), data = trainData)
rfModel <- randomForest(as.formula(paste(targetVar, "~", paste(inputVars, collapse = "+"))), data = trainData, ntree = 500)
svmModel <- svm(as.formula(paste(targetVar, "~", paste(inputVars[c(1,2,3,4)], collapse = "+"))), 
                data = trainData, kernel = "linear", cost = 1, gamma = 1) 
nb_model <- naiveBayes(as.formula(paste(targetVar, "~.")), data = trainData)
knn_model <- train(as.formula(paste(targetVar, "~.")), method = "knn", data = trainData,
                   trControl = trainControl(method = "cv", number = 5))
dt_model <- rpart(as.formula(paste(targetVar, "~.")), data = trainData, method = "class")
gbm_model <- gbm(as.formula(paste(targetVar, "~.")), data = trainData,
                 n.trees = 100, interaction.depth = 3, distribution = "multinomial")
nn_model <- nnet(as.formula(paste(targetVar, "~.")), data = trainData, size = 5, decay = 0.1)

# Evaluate the performance of the models on the training set
multinom_pred_train <- predict(multinomialModel, newdata = trainData)
rf_pred_train <- predict(rfModel, newdata = trainData)
svm_pred_train <- predict(svmModel, newdata = trainData)
nb_pred_train <- predict(nb_model, newdata = trainData)
knn_pred_train <- predict(knn_model, newdata = trainData)
dt_pred_train <- apply(predict(dt_model, newdata = trainData), 1, which.max)
gbm_pred_train <- apply(predict(gbm_model, newdata = trainData), 1, which.max) 
nn_pred_train <- apply(predict(nn_model, newdata = trainData), 1, which.max) 

# Estimate accuracies
multinom_Accuracy_train <- mean(multinom_pred_train == trainData$Class)
rf_Accuracy_train <- mean(rf_pred_train == trainData$Class)
svm_Accuracy_train <- mean(svm_pred_train == trainData$Class)
nb_Accuracy_train <- mean(nb_pred_train == trainData$Class)
knn_Accuracy_train <- mean(knn_pred_train == trainData$Class)
dt_Accuracy_train <- mean(dt_pred_train == trainData$Class)
gbm_Accuracy_train <- mean(gbm_pred_train == trainData$Class)
nn_Accuracy_train <- mean(nn_pred_train == trainData$Class)

Accuarcy_Table <- data.frame(Name = c("Multinomial", "Random Forest",
                                      "Support Vector Machine", "Naive Bayes",
                                      "KNN", "Decision Tree",
                                      "Gradient boosting","Neural Network"),
                             Accuracy = c(multinom_Accuracy,rf_Accuracy,
                                          svm_Accuracy,nb_Accuracy,
                                          knn_Accuracy,dt_Accuracy,gbm_Accuracy,
                                          nn_Accuracy),
                             Accuracy_train = c(multinom_Accuracy_train,rf_Accuracy_train,
                                                svm_Accuracy_train,nb_Accuracy_train,
                                                knn_Accuracy_train,dt_Accuracy_train,gbm_Accuracy_train,
                                                nn_Accuracy_train)) %>%
  arrange(Accuracy)

Accuarcy_Table$OF <- Accuarcy_Table$Accuracy_train - Accuarcy_Table$Accuracy

Accuarcy_Table %>% arrange(OF)