library(MethylPipeR)
library(doMC)
library(doParallel)
registerDoMC(cores = 3)
registerDoParallel(3)

runInformation <- beginPipelineRun(note = 'Random Forest CV with lasso features 10-year survival incremental model. Random Forest timestamp: 2021_05_03_17_51_26. Variables: CpGs, age, sex, bmi, self-reported parent/sibling diabetes, self-reported hypertension. NA covariate rows removed. No class weights', 
                                   logFolderPath = '/path/to/log/folder/',
                                   log = TRUE)

trainingData <- loadData('/path/to/training/data')
trainingTarget <- loadData('/path/to/training/target')

trainingRowsToKeep <- !is.na(trainingTarget$family_diabetes)
trainingData <- trainingData[trainingRowsToKeep, ]
trainingTarget <- trainingTarget[trainingRowsToKeep, ]
row.names(trainingTarget) <- NULL

print(paste0('Training set - before removal of deaths in controls: n = ', nrow(trainingTarget)))
removeDeathsFromControlsResult <- removeDeathsFromControls(trainingTarget, list(trainingData))
trainingData <- removeDeathsFromControlsResult$objectsFiltered[[1]]
trainingTarget <- removeDeathsFromControlsResult$targetFiltered
removeDeathsFromControlsResult <- NULL
gc()
print(paste0('Training set - after removal of deaths in controls: n = ', nrow(trainingTarget)))

threshold <- 10

thresholdTTEResult <- thresholdTTE(trainingTarget,
                                   list(trainingData = trainingData),
                                   threshold)

trainingData <- thresholdTTEResult$objectsFiltered[[1]]
trainingTarget <- thresholdTTEResult$targetFiltered
row.names(trainingTarget) <- NULL
thresholdTTEResult <- NULL
gc()

trainingFolds <- assignTrainingFolds(trainingTarget, 'Event', 3, runInformation)

colsToKeep <- setdiff(colnames(trainingData), c('age', 'sex', 'bmi'))

trainingTarget$age <- trainingData[, 'age']
trainingTarget$sex <- trainingData[, 'sex']
trainingTarget$bmi <- trainingData[, 'bmi']

trainingData <- trainingData[, colsToKeep]


gc()

testData <- loadData('/path/to/test/data')
testTarget <- loadData('/path/to/test/target')

rfResponse <- readRDS(paste0(runInformation[['log folder path']], 'random_forest_lasso_test_response_2021_05_03_17_51_26.rds'))

testRowsToKeep <- !is.na(testTarget$family_diabetes)
testData <- testData[testRowsToKeep, ]
testTarget <- testTarget[testRowsToKeep, ]
row.names(testTarget) <- NULL

print(paste0('Test set - before removal of deaths in controls: n = ', nrow(testTarget)))
removeDeathsFromControlsResult <- removeDeathsFromControls(testTarget, list(testData))
testData <- removeDeathsFromControlsResult$objectsFiltered[[1]]
testTarget <- removeDeathsFromControlsResult$targetFiltered
removeDeathsFromControlsResult <- NULL
gc()
print(paste0('Test set - after removal of deaths in controls: n = ', nrow(testTarget)))

thresholdTTEResult <- thresholdTTE(testTarget,
                                   list(testData = testData, rfResponse = rfResponse),
                                   threshold)

testData <- thresholdTTEResult$objectsFiltered[[1]]
testTarget <- thresholdTTEResult$targetFiltered
row.names(testTarget) <- NULL
thresholdTTEResult <- NULL
gc()

testTarget$age <- testData[, 'age']
testTarget$sex <- testData[, 'sex']
testTarget$bmi <- testData[, 'bmi']

testData <- testData[, colsToKeep]

testTarget$methyl <- rfResponse

nullModel <- glm(Event ~ age + as.factor(sex) + bmi + as.factor(family_diabetes) + as.factor(high_BP), family = 'binomial', data = testTarget)
fullModel <- glm(Event ~ age + as.factor(sex) + bmi + as.factor(family_diabetes) + as.factor(high_BP) + methyl, family = 'binomial', data = testTarget)

nullModelResponse <- predict(nullModel, type = 'response')
fullModelResponse <- predict(fullModel, type = 'response')

if (runInformation[['log']]) {
  saveRDS(nullModel, paste0(runInformation[['log folder path']], 'incremental_rf_lasso_cv_null_model_', runInformation[['start timestamp']], '.rds'))
  saveRDS(fullModel, paste0(runInformation[['log folder path']], 'incremental_rf_lasso_cv_full_model_', runInformation[['start timestamp']], '.rds'))
  saveRDS(nullModelResponse, paste0(runInformation[['log folder path']], 'incremental_rf_lasso_cv_null_model_response_', runInformation[['start timestamp']], '.rds'))
  saveRDS(fullModelResponse, paste0(runInformation[['log folder path']], 'incremental_rf_lasso_cv_full_model_response_', runInformation[['start timestamp']], '.rds'))
}


nullAUC <- pROC::roc(testTarget$Event, nullModelResponse)$auc
fullAUC <- pROC::roc(testTarget$Event, fullModelResponse)$auc
nullPRAUC <- MLmetrics::PRAUC(nullModelResponse, testTarget$Event)
fullPRAUC <- MLmetrics::PRAUC(fullModelResponse, testTarget$Event)

tictoc::tic(paste0('Null model AUC: ', nullAUC, ', full model AUC: ', fullAUC, ', null model PRAUC: ', nullPRAUC, ', full model PRAUC: ', fullPRAUC))
tictoc::toc(log = TRUE)

endPipelineRun(runInformation)
gc()

