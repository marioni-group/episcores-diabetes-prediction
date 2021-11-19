options(java.parameters = "-Xmx50g")
library(bartMachine)
library(MethylPipeR)

runInformation <- beginPipelineRun(note = 'Probit BART (binary classification) with features selected with logistic lasso methylation only (experiment timestamp: 2021_03_29_11_54_07).', 
                                   logFolderPath = '/path/to/log/folder/', 
                                   log = TRUE)

# Here we select rows -1 to exclude the intercept.
selectedFeatures <- read.csv(paste0(runInformation[['log folder path']], 'logistic_lasso_model_non_zero_model_coefficients_2021_03_29_11_54_07.csv'))[-1, 'X']

trainingData <- loadData('/path/to/training/data')
trainingTarget <- loadData('/path/to/training/target')

trainingRowsToKeep <- !is.na(trainingTarget$family_diabetes)
trainingData <- trainingData[trainingRowsToKeep, ]
trainingTarget <- trainingTarget[trainingRowsToKeep, ]
row.names(trainingTarget) <- NULL

trainingData <- cbind(trainingData, as.matrix(trainingTarget[, c('family_diabetes', 'high_BP')]))
trainingData <- trainingData[, selectedFeatures]

removeDeathsFromControlsResult <- removeDeathsFromControls(trainingTarget, list(trainingData))
trainingData <- removeDeathsFromControlsResult$objectsFiltered[[1]]
trainingTarget <- removeDeathsFromControlsResult$targetFiltered
removeDeathsFromControlsResult <- NULL
gc()

threshold <- 10

thresholdTTEResult <- thresholdTTE(trainingTarget,
                                   list(trainingData = trainingData),
                                   threshold)
trainingData <- thresholdTTEResult$objectsFiltered[[1]]
trainingTarget <- thresholdTTEResult$targetFiltered
row.names(trainingTarget) <- NULL
thresholdTTEResult <- NULL
gc()

testData <- loadData('/path/to/test/data')
testTarget <- loadData('/path/to/test/target')

testRowsToKeep <- !is.na(testTarget$family_diabetes)
testData <- testData[testRowsToKeep, ]
testTarget <- testTarget[testRowsToKeep, ]
row.names(testTarget) <- NULL

testData <- cbind(testData, as.matrix(testTarget[, c('family_diabetes', 'high_BP')]))

testData <- testData[, selectedFeatures]
gc()

removeDeathsFromControlsResult <- removeDeathsFromControls(testTarget, list(testData))
testData <- removeDeathsFromControlsResult$objectsFiltered[[1]]
testTarget <- removeDeathsFromControlsResult$targetFiltered
removeDeathsFromControlsResult <- NULL
gc()

thresholdTTEResult <- thresholdTTE(testTarget,
                                   list(testData = testData),
                                   threshold)

testData <- thresholdTTEResult$objectsFiltered[[1]]
testTarget <- thresholdTTEResult$targetFiltered
row.names(testTarget) <- NULL
thresholdTTEResult <- NULL
gc()

set_bart_machine_num_cores(8)

trainingData <- data.frame(trainingData)

bartModel <- bartMachine(X = trainingData,
                           y = factor(trainingTarget$Event, levels = c(1, 0)),
                           num_trees = 100,
                           verbose = TRUE,
                           num_burn_in = 10000,
                           num_iterations_after_burn_in = 20000,
                           seed = runInformation[['random seed']])

colnames(testData) <- selectedFeatures
trainingPredictions <- predict(bartModel, trainingData)
testPredictions <- predict(bartModel, data.frame(testData))

aucTrain <- pROC::roc(trainingTarget$Event, trainingPredictions)$auc
aucTest <- pROC::roc(testTarget$Event, testPredictions)$auc

praucTrain <- MLmetrics::PRAUC(trainingPredictions, trainingTarget$Event)
praucTest <- MLmetrics::PRAUC(testPredictions, testTarget$Event)

if (runInformation[['log']]) {
  saveRDS(bartModel, paste0(runInformation[['log folder path']], 'class_bart_model_', runInformation[['start timestamp']], '.rds'))
  saveRDS(trainingPredictions, paste0(runInformation[['log folder path']], 'class_bart_training_response_', runInformation[['start timestamp']], '.rds'))
  saveRDS(testPredictions, paste0(runInformation[['log folder path']], 'class_bart_test_response_', runInformation[['start timestamp']], '.rds'))
}

endPipelineRun(runInformation)
