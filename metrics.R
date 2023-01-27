library(formattable)

metricsPlot <- function(responses, classes, nThresholds = 101, metricsToPlot = c('Sensitivity', 'Specificity', 'Pos Pred Value', 'Neg Pred Value'), positive = 'TRUE') {
  thresholds <- seq(from = 0, to = 1, length.out = nThresholds)
  predictionsAtThresholds <- lapply(thresholds, function(threshold) {
    factor(responses > threshold)
  })

  predictionCounts <- lapply(predictionsAtThresholds, function(predictions) {
    positives <- sum(predictions == 'TRUE')
    negatives <- sum(predictions == 'FALSE')
    list(positives = positives, negatives = negatives)
  })

  predictedPositives <- sapply(predictionCounts, function(predictionCount) {
    predictionCount$positives
  })

  predictedNegatives <- sapply(predictionCounts, function(predictionCount) {
    predictionCount$negatives
  })

  predictionCountDF <- data.frame(list(Threshold = thresholds, Positives = predictedPositives, Negatives = predictedNegatives))

  # Threshold index (keep 11 rows)
  thresholdIndex <- seq(from = 1, to = nThresholds, by = floor(nThresholds/10))
  
  predictionCountDFToDisplay <- predictionCountDF[thresholdIndex, ]

  predictionCountDFToDisplay <- data.frame(t(predictionCountDFToDisplay))

  colnames(predictionCountDFToDisplay) <- predictionCountDFToDisplay['Threshold',]

  # Remove Threshold row
  predictionCountDFToDisplay <- predictionCountDFToDisplay[row.names(predictionCountDFToDisplay) != 'Threshold',]

  # Make a new column named 'Threshold' with elements 'Positives' and 'Negatives' for display purposes
  predictionCountDFToDisplay$Threshold <- row.names(predictionCountDFToDisplay)
  row.names(predictionCountDFToDisplay) <- NULL
  
  # Move Threshold column to beginning
  colnamesReordered <- c('Threshold', colnames(predictionCountDFToDisplay)[1:length(predictionCountDFToDisplay)-1])
  predictionCountDFToDisplay <- predictionCountDFToDisplay[, colnamesReordered]

  predictionCountFormattedTable <- formattable(predictionCountDFToDisplay)

  confusionMatrixResults <- lapply(predictionsAtThresholds, function(predictions) {
    confusionMatrix(predictions, classes, positive = positive)$byClass
  })
  
  plotDF <- data.frame(lapply(metricsToPlot, function(metricToPlot) {
    sapply(confusionMatrixResults, function(confusionMatrixResult) {
      confusionMatrixResult[[metricToPlot]]
    })
  }))

  colnames(plotDF) <- metricsToPlot

  plotDF <- melt(plotDF)

  plotDF$Threshold <- rep(thresholds, length(metricsToPlot))

  p <- ggplot(data = plotDF, aes(x = Threshold, y = value, colour = variable)) +
         geom_line() + theme_minimal()

  list(plotDF = plotDF, confusionMatrixResults = confusionMatrixResults, predictionsAtThresholds = predictionsAtThresholds, p = p, predictionCountDF = predictionCountDF, predictionCountFormattedTable = predictionCountFormattedTable)
}
