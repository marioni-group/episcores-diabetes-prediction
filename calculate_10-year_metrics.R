library(survival)

predictCoxPHOnset <- function(dataDF, coxPHModel, threshold = 10) {
  uniqueTimes <- sort(unique(dataDF$time_to_event))
  thresholdIndex <- match(threshold, uniqueTimes)
  cumulativeBaseHaz <- gbm::basehaz.gbm(dataDF$time_to_event, dataDF$Event, predict(coxPHModel), uniqueTimes)
  survivalPredictions <- exp(-cumulativeBaseHaz[[thresholdIndex]]) ^ exp(predict(coxPHModel))
  onsetPredictions <- 1 - survivalPredictions

  # Event should be 0 if tte is > 10
  dataDF$Event <- sapply(1:nrow(dataDF), function(i) {
    if (dataDF$time_to_event[[i]] > 10) {
      0
    } else {
      dataDF$Event[[i]]
    }
  })

  auc <- MLmetrics::AUC(y_pred = onsetPredictions, y_true = dataDF$Event)
  prauc <- MLmetrics::PRAUC(y_pred = onsetPredictions, y_true = dataDF$Event)
  roc <- pROC::roc(response = dataDF$Event, predictor = onsetPredictions)
  list(cumulativeBaseHaz = cumulativeBaseHaz, onsetPredictions = onsetPredictions, auc = auc, prauc = prauc, roc = roc)
}

calculateModels <- function(targetW1) {
  riskFactorsOnlyCoxPH <- coxph(Surv(time_to_event, Event) ~ Age + sex + bmi + high_BP + family_diabetes, targetW1)
  proteinCoxPH <- coxph(Surv(time_to_event, Event) ~ Age + sex + bmi + high_BP + family_diabetes + pScore, targetW1)
  directCoxPH <- coxph(Surv(time_to_event, Event) ~ Age + sex + bmi + high_BP + family_diabetes + dScore, targetW1)
  fullCoxPH <- coxph(Surv(time_to_event, Event) ~ Age + sex + bmi + high_BP + family_diabetes + pScore + dScore, targetW1)

  list(r = riskFactorsOnlyCoxPH, p = proteinCoxPH, d = directCoxPH, f = fullCoxPH)
}


calculateResults <- function(targetW1) {
  models <- calculateModels(targetW1)
  testResults <- lapply(models, function(m) {predictCoxPHOnset(targetW1, m)})

  aucs <- sapply(testResults, function(r) {r$auc})
  praucs <- sapply(testResults, function(r) {r$prauc})

  # The 5th column corresponds to the p-value in the coefficients table
  pModelPValue <- summary(models$p)$coefficients['pScore', 5]
  dModelPValue <- summary(models$d)$coefficients['dScore', 5]
  fModelPScorePValue <- summary(models$f)$coefficients['pScore', 5]
  fModelDScorePValue <- summary(models$f)$coefficients['dScore', 5]

  pValues <- c(pModel = pModelPValue, dModel = dModelPValue, fModelPScore = fModelPScorePValue, fModelDScore = fModelDScorePValue)

  predictions <- lapply(testResults, function(r) {r$onsetPredictions})

  rocTestResults <- sapply(2:length(predictions), function(i) {
    pROC::roc.test(response = targetW1$Event, predictor1 = predictions[[1]], predictor2 = predictions[[i]])$p.value
  })
  names(rocTestResults) <- c('p', 'd', 'f')

  metricsTable <- data.frame(AUC = aucs, PRAUC = praucs)

  row.names(metricsTable) <- c('Risk factors only', 'Risk factors + protein EpiScore', 'Risk factors + direct EpiScore', 'Full model')
  list(testResults = testResults, metricsTable = metricsTable, models = models, pValues = pValues, rocTestResults = rocTestResults)
}


targetW1Cox <- readRDS('/path/to/file')
targetW1BART <- readRDS('/path/to/file')
targetW1RSF <- readRDS('/path/to/file')
targets <- list(cox = targetW1Cox, rsf = targetW1RSF, bart = targetW1BART)

results <- lapply(targets, calculateResults)

tables <- lapply(results, function(result) {result$metricsTable})

aucs <- lapply(tables, function(t) {t$AUC})
praucs <- lapply(tables, function(t) {t$PRAUC})

aucDF <- data.frame(aucs)
praucDF <- data.frame(praucs)

pValues <- lapply(results, function(result) {result$pValues})
pValuesDF <- data.frame(pValues)

rocTestResults <- lapply(results, function(result) {result$rocTestResults})
rocTestResultsDF <- data.frame(rocTestResults)

row.names(aucDF) <- row.names(tables$cox)
row.names(praucDF) <- row.names(tables$cox)



# write.csv(aucDF, '/path/to/file')
# write.csv(praucDF, '/path/to/file')
# write.csv(pValuesDF, '/path/to/file')
# write.csv(rocTestResultsDF, '/path/to/file')

