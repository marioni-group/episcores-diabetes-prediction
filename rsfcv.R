library(randomForestSRC)

srfTrainAndTest <- function(trainDF, testDF, metric = 'AUC', ntree, mtry, nodesize, seed, threshold = 10) {
  srfModel <- rfsrc(Surv(time_to_event, Event) ~ ., trainDF, seed = seed, ntree = ntree, mtry = mtry, nodesize = nodesize)

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

srfCVIter <- function(dataDF, foldIDs, metric, ntree, mtry, nodesize, seed, threshold = 10) {
  foldIDSet <- unique(foldIDs)
  nFolds <- length(foldIDSet)
  metricResults <- sapply(foldIDSet, function(foldID) {
    testIndex <- foldIDs == foldID
    testDF <- dataDF[testIndex, ]
    trainDF <- dataDF[!testIndex, ]
    srfTrainAndTest(trainDF, testDF, metric, ntree, mtry, nodesize, seed, threshold)
  })
  meanMetricResult <- mean(metricResults)
  print(paste0('ntree = ', ntree, ', mtry = ', mtry, ', nodesize = ', nodesize, ', mean_metric = ', meanMetricResult))
  meanMetricResult
}

srfCV <- function(dataDF, foldIDs, metric = 'AUC', ntrees, mtrys, nodesizes, seed, threshold = 10) {
  meanMetricResults <- data.frame(ntree = integer(0), mtry = integer(0), nodesize = integer(0), mean_metric = double(0))
  for (ntree in ntrees) {
    for (mtry in mtrys) {
      for (nodesize in nodesizes) {
        cvIterResult <- srfCVIter(dataDF, foldIDs, metric, ntree, mtry, nodesize, seed, threshold = threshold)
        meanMetricResults <- rbind(meanMetricResults, list(ntree = ntree, mtry = mtry, nodesize = nodesize, mean_metric = cvIterResult))
      }
    }
  }
  meanMetricResults
}


