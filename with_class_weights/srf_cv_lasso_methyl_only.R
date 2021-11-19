library(MethylPipeR)
library(doMC)
library(doParallel)
registerDoMC(cores = 3)
registerDoParallel(3)
library(randomForestSRC)

runInformation <- beginPipelineRun(note = 'Survival random forest with lasso-selected variables methylation only, with weights (lasso experiment timestamp: 2021_04_01_13_03_03). Train on wave 3, test on wave 1. NA covariates removed', 
                                   logFolderPath = '/path/to/log/folder/',
                                   log = TRUE)

selectedFeatures <- read.csv(paste0(runInformation[['log folder path']], 'logistic_lasso_model_non_zero_model_coefficients_2021_04_01_13_03_03.csv'))[-1, 'X']

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

# thresholdTTEResult <- thresholdTTE(trainingTarget,
#                                    list(trainingData = trainingData),
#                                    threshold)

# trainingData <- thresholdTTEResult$objectsFiltered[[1]]
# trainingTarget <- thresholdTTEResult$targetFiltered
# row.names(trainingTarget) <- NULL
# thresholdTTEResult <- NULL
# gc()

trainingFolds <- assignTrainingFolds(trainingTarget, 'Event', 3, runInformation)
trainingWeights <- createClassWeights(trainingTarget$Event)

trainingData <- cbind(trainingData, as.matrix(trainingTarget[, c('family_diabetes', 'high_BP')]))

trainingData <- trainingData[, selectedFeatures]

testData <- loadData('/path/to/test/data')
testTarget <- loadData('/path/to/test/target')

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
                                   list(testData = testData),
                                   threshold)

testData <- thresholdTTEResult$objectsFiltered[[1]]
testTarget <- thresholdTTEResult$targetFiltered
row.names(testTarget) <- NULL
thresholdTTEResult <- NULL
gc()

testData <- cbind(testData, as.matrix(testTarget[, c('family_diabetes', 'high_BP')]))

testData <- testData[, selectedFeatures]

trainingData <- data.frame(trainingData)

trainingData$time_to_event <- trainingTarget$time_to_event
trainingData$Event <- trainingTarget$Event


srfTrainAndTest <- function(trainDF, testDF, metric = 'AUC', ntree, mtry, nodesize, pipelineRunInformation) {
  srfModel <- rfsrc(Surv(time_to_event, Event) ~ ., trainDF, seed = pipelineRunInformation[['random seed']], ntree = ntree, mtry = mtry, nodesize = nodesize)

  # Threshold for binary prediction
  thresholdTTEResult <- thresholdTTE(testDF,
                                     NULL,
                                     threshold)

  testDF <- thresholdTTEResult$targetFiltered
  row.names(testDF) <- NULL
  thresholdTTEResult <- NULL
  gc()

  testPredictions <- predict(srfModel, testDF)

  testResultIndex <- match(threshold, testPredictions$time.interest)

  testOnsetPredictions <- 1 - testPredictions$survival[, testResultIndex]
  if (metric == 'AUC') {
    metricResult <- pROC::roc(testDF$Event, testOnsetPredictions)$auc
  } else if (metric == 'PRAUC') {
    metricResult <- MLmetrics::PRAUC(testOnsetPredictions, testDF$Event)
  }
  print(paste0('ntree = ', ntree, ', mtry = ', mtry, ', nodesize = ', nodesize, ', metricResult = ', metricResult))
  metricResult
}

srfCVIter <- function(dataDF, foldIDs, metric, ntree, mtry, nodesize, pipelineRunInformation) {
  foldIDSet <- unique(foldIDs)
  nFolds <- length(foldIDSet)
  metricResults <- sapply(foldIDSet, function(foldID) {
    testIndex <- foldIDs == foldID
    testDF <- dataDF[testIndex, ]
    trainDF <- dataDF[!testIndex, ]
    srfTrainAndTest(trainDF, testDF, metric, ntree, mtry, nodesize, pipelineRunInformation)
  })
  meanMetricResult <- mean(metricResults)
  print(paste0('ntree = ', ntree, ', mtry = ', mtry, ', nodesize = ', nodesize, ', mean_metric = ', meanMetricResult))
  meanMetricResult
}

srfCV <- function(dataDF, foldIDs, metric = 'AUC', ntrees, mtrys, nodesizes, pipelineRunInformation) {
  meanMetricResults <- data.frame(ntree = integer(0), mtry = integer(0), nodesize = integer(0), mean_metric = double(0))
  for (ntree in ntrees) {
    for (mtry in mtrys) {
      for (nodesize in nodesizes) {
        cvIterResult <- srfCVIter(dataDF, foldIDs, metric, ntree, mtry, nodesize, pipelineRunInformation)
        meanMetricResults <- rbind(meanMetricResults, list(ntree = ntree, mtry = mtry, nodesize = nodesize, mean_metric = cvIterResult))
      }
    }
  }
  meanMetricResults
}

srfCVResult <- srfCV(dataDF = trainingData, foldIDs = trainingFolds, metric = 'AUC',
                     ntrees = c(1600, 3200, 6400),
                     mtrys = c(4, 8, 16),
                     nodesizes = c(32, 64, 96, 128),
                     pipelineRunInformation = runInformation)
# sort
srfCVResult <- srfCVResult[order(srfCVResult$mean_metric, decreasing = TRUE),]

ntreeBest <- srfCVResult[1, 'ntree']
mtryBest <- srfCVResult[1, 'mtry']
nodesizeBest <- srfCVResult[1, 'nodesize']

survivalRFModel <- rfsrc(Surv(time_to_event, Event) ~ ., trainingData, seed = runInformation[['random seed']], ntree = ntreeBest, mtry = mtryBest, nodesize = nodesizeBest)

testPredictions <- predict(survivalRFModel, data.frame(testData))

# Find index of first time point that is greater than or equal to 10
testResultIndex <- match(TRUE, testPredictions$time.interest >= threshold)

testOnsetPredictions <- 1 - testPredictions$survival[, testResultIndex]

# Threshold TTE for training set for binary prediction metrics
thresholdTTEResult <- thresholdTTE(trainingTarget,
                                   list(trainingData = trainingData),
                                   threshold)

trainingData <- thresholdTTEResult$objectsFiltered[[1]]
trainingTarget <- thresholdTTEResult$targetFiltered
row.names(trainingTarget) <- NULL
thresholdTTEResult <- NULL
gc()

trainingData$time_to_event <- NULL
trainingData$Event <- NULL
trainPredictions <- predict(survivalRFModel, trainingData)

# Find index of first time point that is greater than or equal to 10
trainResultIndex <- match(TRUE, trainPredictions$time.interest >= threshold)

trainOnsetPredictions <- 1 - trainPredictions$survival[, trainResultIndex]

trainROCResult <- pROC::roc(trainingTarget$Event, trainOnsetPredictions)$auc
trainPRAUCResult <- MLmetrics::PRAUC(trainOnsetPredictions, trainingTarget$Event)
testROCResult <- pROC::roc(testTarget$Event, testOnsetPredictions)$auc
testPRAUCResult <- MLmetrics::PRAUC(testOnsetPredictions, testTarget$Event)

if (runInformation[['log']]) {
  tictoc::tic(paste0('Training AUC: ', trainROCResult, ', training PRAUC: ', trainPRAUCResult, ', test AUC: ', testROCResult, ', test PRAUC: ', testPRAUCResult))
  tictoc::toc(log = TRUE)
  tictoc::tic('Save training and test responses.')
  saveRDS(trainOnsetPredictions, paste0(runInformation[['log folder path']], 'survival_random_forest_train_response_', runInformation[['start timestamp']], '.rds'))
  saveRDS(testOnsetPredictions, paste0(runInformation[['log folder path']], 'survival_random_forest_test_response_', runInformation[['start timestamp']], '.rds'))
  tictoc::toc(log = TRUE)
}

endPipelineRun(runInformation)
gc()

