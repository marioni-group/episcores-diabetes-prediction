library(tidyverse)
library(MLmetrics)

targetW1 <- readRDS('/path/to/file')

testResults <- readRDS('/path/to/file')

# Event should be 0 if tte is > 10
targetW1$Event <- sapply(1:nrow(targetW1), function(i) {
  if (targetW1$time_to_event[[i]] > 10) {
    0
  } else {
    targetW1$Event[[i]]
  }
})

onsetPredictions <- lapply(testResults, function(x) {x$onsetPredictions})

onlyCasesAfter2YearsIndex <- (targetW1$tte > 2) | (targetW1$Event == 0)

targetFiltered <- targetW1[onlyCasesAfter2YearsIndex, ]
onsetPredictionsFiltered <- lapply(testResults, function(testResult) {
  testResult$onsetPredictions[onlyCasesAfter2YearsIndex]
})


metrics <- lapply(onsetPredictionsFiltered, function(x) {
  auc <- MLmetrics::AUC(y_pred = x, y_true = targetFiltered$Event)
  prauc <- MLmetrics::PRAUC(y_pred = x, y_true = targetFiltered$Event)
  list(auc = auc, prauc = prauc)
})

metricsOriginal <- lapply(onsetPredictions, function(x) {
  auc <- MLmetrics::AUC(y_pred = x, y_true = targetW1$Event)
  prauc <- MLmetrics::PRAUC(y_pred = x, y_true = targetW1$Event)
  list(auc = auc, prauc = prauc)
})

aucs <- sapply(metrics, function(x) {
  x$auc
})

praucs <- sapply(metrics, function(x) {
  x$prauc
})

aucsPrevious <- sapply(metricsOriginal, function(x) {
  x$auc
})

praucsPrevious <- sapply(metricsOriginal, function(x) {
  x$prauc
})

aucDF <- data.frame(original_auc = aucsPrevious, two_year_threshold_auc = aucs)
row.names(aucDF) <- c('risk_factors', 'protein_episcore', 'direct_episcore', 'full_model')

praucDF <- data.frame(original_prauc = praucsPrevious, two_year_threshold_prauc = praucs)
row.names(praucDF) <- c('risk_factors', 'protein_episcore', 'direct_episcore', 'full_model')

write.csv(aucDF, 'cases_after_2_years_aucs.csv')
write.csv(praucDF, 'cases_after_2_years_praucs.csv')

