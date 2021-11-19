library(MethylPipeR)
library(doMC)
library(doParallel)
registerDoMC(cores = 3)
registerDoParallel(3)
library(randomForest)

runInformation <- beginPipelineRun(note = 'Random forest with lasso-selected variables, methylation only (Experiment timestamp: 2021_03_31_10_31_42. No class weights). Train on wave 3, test on wave 1.', 
                                   logFolderPath = '/path/to/log/folder/',
                                   log = TRUE)

selectedFeatures <- read.csv(paste0(runInformation[['log folder path']], 'logistic_lasso_model_non_zero_model_coefficients_2021_03_31_10_31_42.csv'))[-1, 'X']

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

rfTrainAndTest <- function(trainXs, trainY, testXs, testY, metric = 'AUC', ntree = 500, mtry = NULL, nodesize = 1, pipelineRunInformation) {
  # If mtry is not set, set default mtry for randomForest
  if (is.null(mtry)) {
    mtry <- if (!is.null(trainY) && !is.factor(trainY)) max(floor(ncol(trainXs)/3), 1) else floor(sqrt(ncol(trainXs)))
  }
  set.seed(pipelineRunInformation[['random seed']])
  model <- randomForest(trainXs, y = trainY, xtest = testXs, ytest = testY, ntree = ntree, mtry = mtry, nodesize = nodesize)
  if (metric == 'AUC') {
    metricResult <- pROC::roc(testY, model$test$votes[, 2])$auc
  } else if (metric == 'PRAUC') {
    metricResult <- MLmetrics::PRAUC(rfModel$test$votes[, 2], testTarget$Event)
  }
  print(paste0('ntree = ', ntree, ', mtry = ', mtry, ', nodesize = ', nodesize, ', metricResult = ', metricResult))
  metricResult
}

rfCVIter <- function(dataset, labels, foldIDs, metric, ntree, mtry, nodesize, pipelineRunInformation) {
  foldIDSet <- unique(foldIDs)
  nFolds <- length(foldIDSet)
  metricResults <- sapply(foldIDSet, function(foldID) {
    testIndex <- foldIDs == foldID
    testXs <- dataset[testIndex, ]
    testY <- labels[testIndex]
    trainXs <- dataset[!testIndex, ]
    trainY <- labels[!testIndex]
    rfTrainAndTest(trainXs, trainY, testXs, testY, metric, ntree, mtry, nodesize, pipelineRunInformation)
  })
  meanMetricResult <- mean(metricResults)
  print(paste0('ntree = ', ntree, ', mtry = ', mtry, ', nodesize = ', nodesize, ', mean_metric = ', meanMetricResult))
  meanMetricResult
}

rfCV <- function(dataset, labels, foldIDs, metric = 'AUC', ntrees, mtrys, nodesizes, pipelineRunInformation) {
  meanMetricResults <- data.frame(ntree = integer(0), mtry = integer(0), nodesize = integer(0), mean_metric = double(0))
  for (ntree in ntrees) {
    for (mtry in mtrys) {
      for (nodesize in nodesizes) {
        cvIterResult <- rfCVIter(dataset, labels, foldIDs, metric, ntree, mtry, nodesize, pipelineRunInformation)
        meanMetricResults <- rbind(meanMetricResults, list(ntree = ntree, mtry = mtry, nodesize = nodesize, mean_metric = cvIterResult))
      }
    }
  }
  meanMetricResults
}



rfCVResult <- rfCV(dataset = trainingData, labels = as.factor(trainingTarget$Event), foldIDs = trainingFolds, metric = 'AUC', 
                   ntrees = c(400, 800, 1600, 3200, 6400),
                   mtrys = c(2, 4, 8, 16, 32),
                   nodesizes = c(2, 4, 8, 16, 32),
                   pipelineRunInformation = runInformation)

# Sort
rfCVResult <- rfCVResult[order(rfCVResult$mean_metric, decreasing = TRUE),]

ntreeBest <- rfCVResult[1, 'ntree']
mtryBest <- rfCVResult[1, 'mtry']
nodesizeBest <- rfCVResult[1, 'nodesize']

rfModel <- randomForest(trainingData, y = as.factor(trainingTarget$Event), xtest = testData, ytest = as.factor(testTarget$Event), ntree = ntreeBest, mtry = mtryBest)

rfTrainResponse <- rfModel$votes[, 2]
rfTestResponse <- rfModel$test$votes[, 2]

if (runInformation[['log']]) {
  saveRDS(rfTrainResponse, paste0(runInformation[['log folder path']], 'random_forest_lasso_train_response_', runInformation[['start timestamp']], '.rds'))
  saveRDS(rfTestResponse, paste0(runInformation[['log folder path']], 'random_forest_lasso_test_response_', runInformation[['start timestamp']], '.rds'))
}

aucTrain <- pROC::roc(trainingTarget$Event, rfModel$votes[, 2])$auc
aucTest <- pROC::roc(testTarget$Event, rfModel$test$votes[, 2])$auc
praucTrain <- MLmetrics::PRAUC(rfModel$votes[, 2], trainingTarget$Event)
praucTest <- MLmetrics::PRAUC(rfModel$test$votes[, 2], testTarget$Event)

endPipelineRun(runInformation)
gc()

