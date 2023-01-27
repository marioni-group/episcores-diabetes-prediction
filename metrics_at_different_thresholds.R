source('metrics.R')
library(MethylPipeR)
library(caret)
library(reshape2)
library(patchwork)

# Load Wave 1 target
w1Target <- readRDS('/path/to/file')

# Event should be 0 if tte is > 10
w1Target$Event <- sapply(1:nrow(w1Target), function(i) {
  if (w1Target$time_to_event[[i]] > 10) {
    0
  } else {
    w1Target$Event[[i]]
  }
})

coxTestResults <- readRDS('/path/to/file')

nullResponse <- coxTestResults$r$onsetPredictions
fullResponse <- coxTestResults$d$onsetPredictions

incrementalCoxLassoResponse <- fullResponse

nullMetrics0.5 <- MLmetrics::ConfusionMatrix(nullResponse > 0.5, w1Target$Event)
incrementalCoxLassoMetrics0.5 <- MLmetrics::ConfusionMatrix(incrementalCoxLassoResponse > 0.5, w1Target$Event)

caretNullMetrics0.5 <- caret::confusionMatrix(as.factor(nullResponse > 0.5), as.factor(w1Target$Event == 1), positive = 'TRUE')
caretIncrementalCoxLassoMetrics0.5 <- caret::confusionMatrix(as.factor(incrementalCoxLassoResponse > 0.5), as.factor(w1Target$Event == 1), positive = 'TRUE')

calculateMetricsAtThreshold <- function(predicted, actual, threshold) {
  caret::confusionMatrix(as.factor(predicted > threshold), as.factor(actual == 1), positive = 'TRUE')
}

nullMetrics0.5 <- calculateMetricsAtThreshold(nullResponse, w1Target$Event, 0.5)
incrementalCoxLassoMetrics0.5 <- calculateMetricsAtThreshold(incrementalCoxLassoResponse, w1Target$Event, 0.5)

thresholds <- seq(0.1, 1, 0.1)

incrementalCoxLassoMetrics <- lapply(thresholds, function(threshold) {list(
  nullModel = calculateMetricsAtThreshold(nullResponse, w1Target$Event, threshold),
  fullModel = calculateMetricsAtThreshold(incrementalCoxLassoResponse, w1Target$Event, threshold)
)})

names(incrementalCoxLassoMetrics) <- thresholds

nullTruePositives <- sapply(thresholds, function(threshold) {
  incrementalCoxLassoMetrics[[toString(threshold)]]$nullModel$table[2,2]
})

fullTruePositives <- sapply(thresholds, function(threshold) {
  incrementalCoxLassoMetrics[[toString(threshold)]]$fullModel$table[2,2]
})

nullTrueNegatives <- sapply(thresholds, function(threshold) {
  incrementalCoxLassoMetrics[[toString(threshold)]]$nullModel$table[1,1]
})

fullTrueNegatives <- sapply(thresholds, function(threshold) {
  incrementalCoxLassoMetrics[[toString(threshold)]]$fullModel$table[1,1]
})

nullFalsePositives <- sapply(thresholds, function(threshold) {
  incrementalCoxLassoMetrics[[toString(threshold)]]$nullModel$table[2,1]
})

fullFalsePositives <- sapply(thresholds, function(threshold) {
  incrementalCoxLassoMetrics[[toString(threshold)]]$fullModel$table[2,1]
})

nullFalseNegatives <- sapply(thresholds, function(threshold) {
  incrementalCoxLassoMetrics[[toString(threshold)]]$nullModel$table[1,2]
})

fullFalseNegatives <- sapply(thresholds, function(threshold) {
  incrementalCoxLassoMetrics[[toString(threshold)]]$fullModel$table[1,2]
})

tables <- lapply(thresholds, function(threshold) {
  incrementalCoxLassoMetrics[[toString(threshold)]]$fullModel$table
})

par(mfrow = c(2,2))

plot(thresholds, fullTruePositives, type = 'o', col = 'blue', pch='o',
     main = 'True positives', 
     xlab = 'Threshold', ylab = 'N',
     ylim = c(0, 150))
points(thresholds, nullTruePositives, col = 'red', pch='*')
lines(thresholds, nullTruePositives, col = 'red')
legend('topright', legend = c('Risk factors only', 'Full'), col = c('red', 'blue'), lty = 1)

# plot(thresholds, fullTrueNegatives, type = 'o', col = 'blue', pch='o',
#      main = 'True negatives', 
#      xlab = 'Threshold', ylab = 'Number of true negatives',
#      ylim = c(4100, 4600))
# points(thresholds, nullTrueNegatives, col = 'red', pch='*')
# lines(thresholds, nullTrueNegatives, col = 'red')
# legend('topleft', legend = c('null', 'full'), col = c('red', 'blue'), lty = 1)

plot(thresholds, fullFalseNegatives, type = 'o', col = 'blue', pch='o',
     main = 'False negatives',
     xlab = 'Threshold', ylab = 'N',
     ylim = c(50, 250))
points(thresholds, nullFalseNegatives, col = 'red', pch='*')
lines(thresholds, nullFalseNegatives, col = 'red')

plot(thresholds, fullFalsePositives, type = 'o', col = 'blue', pch='o',
     main = 'False positives', 
     xlab = 'Threshold', ylab = 'N',
     ylim = c(0, 500))
points(thresholds, nullFalsePositives, col = 'red', pch='*')
lines(thresholds, nullFalsePositives, col = 'red')
# legend('topright', legend = c('null', 'full'), col = c('red', 'blue'), lty = 1)

plot(thresholds, fullTrueNegatives, type = 'o', col = 'blue', pch='o',
     main = 'True negatives',
     xlab = 'Threshold', ylab = 'N',
     ylim = c(4100, 4600))
points(thresholds, nullTrueNegatives, col = 'red', pch='*')
lines(thresholds, nullTrueNegatives, col = 'red')

# plot(thresholds, fullFalseNegatives, type = 'o', col = 'blue', pch='o',
#      main = 'False negatives', 
#      xlab = 'Threshold', ylab = 'Number of false negatives',
#      ylim = c(50, 250))
# points(thresholds, nullFalseNegatives, col = 'red', pch='*')
# lines(thresholds, nullFalseNegatives, col = 'red')
# legend('topleft', legend = c('null', 'full'), col = c('red', 'blue'), lty = 1)

# incrementalClassBartLassoResponse <- readRDS(paste0(logFolderPath, '/path/to/file'))

# incrementalCoxMetricsPlot <- metricsPlot(incrementalCoxLassoResponse, factor(w1Target$Event == 1))
# incrementalBartMetricsPlot <- metricsPlot(incrementalClassBartLassoResponse, factor(w1Target$Event == 1))

# incrementalCombinedPlot <- incrementalCoxMetricsPlot$p / incrementalBartMetricsPlot$p

# bartResponse <- readRDS(paste0(logFolderPath, '/path/to/file'))
# logisticElnetResponse <- readRDS(paste0(logFolderPath, '/path/to/file'))

# logisticElnetMetricsPlot <- metricsPlot(logisticElnetResponse, factor(w1Target$Event == 1))
# bartMetricsPlot <- metricsPlot(bartResponse, factor(w1Target$Event == 1))

# combinedPlot <- logisticElnetMetricsPlot$p / bartMetricsPlot$p

fullTPRs <- fullTruePositives / (fullTruePositives + fullFalseNegatives)
nullTPRs <- nullTruePositives / (nullTruePositives + nullFalseNegatives)

fullTNRs <- fullTrueNegatives / (fullTrueNegatives + fullFalsePositives)
nullTNRs <- nullTrueNegatives / (nullTrueNegatives + nullFalsePositives)

fullPPVs <- fullTruePositives / (fullTruePositives + fullFalsePositives)
nullPPVs <- nullTruePositives / (nullTruePositives + nullFalsePositives)

fullNPVs <- fullTrueNegatives / (fullTrueNegatives + fullFalseNegatives)
nullNPVs <- nullTrueNegatives / (nullTrueNegatives + nullFalseNegatives)


