library(MethylPipeR)
library(survival)
library(BART)

startTimestamp <- format(Sys.time(), '%Y_%m_%d_%H_%M_%S')

initLogs('/path/')

set.seed(42)

cpgs <- readRDS('/path/to/file')

targetW4 <- readRDS('/path/to/file')
methylW4 <- readRDS('/path/to/file')
methylW4 <- methylW4[, cpgs]
gc()

targetW3 <- readRDS('/path/to/file')
methylW3 <- readRDS('/path/to/file')
methylW3 <- methylW3[, cpgs]
gc()

targetW1 <- readRDS('/path/to/file')
methylW1 <- readRDS('/path/to/file')
methylW1 <- methylW1[, cpgs]

# Filter out NA BMI rows from targetW1
targetW1NonNABMIIndex <- !is.na(targetW1$bmi)
methylW1 <- methylW1[targetW1NonNABMIIndex, ]
targetW1 <- targetW1[targetW1NonNABMIIndex, ]
row.names(targetW1) <- NULL

getFilterByVarianceIDs <- function(data, numberOfFeatures) {
  vars <- apply(data, 2, var)
  sortedVars <- sort(vars, index.return = TRUE, decreasing = TRUE)
  cpgIDs <- colnames(data)[sortedVars$ix]
  cpgIDs[1:numberOfFeatures]
}

p <- 200000
topCpGsByVariance <- getFilterByVarianceIDs(methylW3, p)
gc()

methylW4 <- methylW4[, topCpGsByVariance]
methylW3 <- methylW3[, topCpGsByVariance]
methylW1 <- methylW1[, topCpGsByVariance]
gc()

pedigree <- read.csv('/path/to/file')

targetW4 <- merge(targetW4, pedigree[, c('volid', 'famid')], by.x = 'Sample_Name', by.y = 'volid', all.x = TRUE)

targetW3 <- merge(targetW3, pedigree[, c('volid', 'famid')], by.x = 'Sample_Name', by.y = 'volid', all.x = TRUE)

targetW1 <- merge(targetW1, pedigree[, c('volid', 'famid')], by.x = 'Sample_Name', by.y = 'volid', all.x = TRUE)

writeLines('Loaded data')

set.seed(42)

w4ShuffleIndex <- sample(1:nrow(targetW4))
methylW4 <- methylW4[w4ShuffleIndex, ]
targetW4 <- targetW4[w4ShuffleIndex, ]
row.names(targetW4) <- NULL

w1ShuffleIndex <- sample(1:nrow(targetW1))
methylW1 <- methylW1[w1ShuffleIndex, ]
targetW1 <- targetW1[w1ShuffleIndex, ]
row.names(targetW1) <- NULL

set.seed(42)

methylW4W3 <- rbind(methylW4, methylW3)
targetW4W3 <- rbind(targetW4, targetW3)

w4w3ShuffleIndex <- sample(1:nrow(targetW4W3))
methylW4W3 <- methylW4W3[w4w3ShuffleIndex, ]
targetW4W3 <- targetW4W3[w4w3ShuffleIndex, ]
row.names(targetW4W3) <- NULL

set.seed(42)
foldIDsW4W3 <- getGroupCVFoldIDs(targetW4W3$famid, 9)$foldIDs

writeLines('Fitting proteinLasso')
proteinEpiScoresW4 <- read.csv('/path/to/file', row.names = 1)
proteinEpiScoresW3W1 <- read.csv('/path/to/file', row.names = 1)
proteinEpiScores <- rbind(proteinEpiScoresW4, proteinEpiScoresW3W1)

proteinEpiScores <- as.matrix(proteinEpiScores)

proteinEpiScoresW4W3 <- proteinEpiScores[match(row.names(methylW4W3), row.names(proteinEpiScores)), ]
proteinEpiScoresW1 <- proteinEpiScores[match(row.names(methylW1), row.names(proteinEpiScores)), ]

proteinLasso <- fitMPRModelCV('survival', 'glmnet', proteinEpiScoresW4W3, targetW4W3[, c('Event', 'time_to_event')], seed = 42, nFolds = 9, foldID = foldIDsW4W3)

proteinPredictionW1 <- predictMPRModel(proteinLasso, proteinEpiScoresW1, s = 'lambda.min')

targetW1$pScore <- proteinPredictionW1



writeLines('Fitting direct T2D Lasso')


directLasso <- readRDS('/path/to/file')
directLassoCoefs <- coef(directLasso$model, s = 'lambda.min')
directLassoCpGs <- names(directLassoCoefs[directLassoCoefs[, 1] != 0, ])

methylW4W3 <- methylW4W3[, directLassoCpGs]
methylW1 <- methylW1[, directLassoCpGs]

gc()

modelDirectT2DW4W3 <- mc.surv.bart(x.train = methylW4W3, times = targetW4W3$time_to_event, delta = targetW4W3$Event, x.test = methylW1, mc.cores = 8, seed = 42)

gc()

predictionTimePoint <- 10
numberOfTimePoints <- modelDirectT2DW4W3$K
predictionTimePointIndex <- match(predictionTimePoint, modelDirectT2DW4W3$times)

# Given an index for an individual (a row in the dataset), a time point and the total number of time points, returns the corresponding column index in bart.survival.model$surv.test (or any result with the same structure).
getBartResultColumn <- function(individual, timePointIndex, nTimePoints) {
  (individual - 1) * nTimePoints + timePointIndex
}

survivalMeanPredictions <- sapply(1:nrow(targetW1), function(individual) {
  modelDirectT2DW4W3$surv.test.mean[[getBartResultColumn(individual, predictionTimePointIndex, numberOfTimePoints)]]
})

onsetPredictions <- 1 - survivalMeanPredictions

targetW1$dScore <- onsetPredictions

riskFactorsOnlyCoxPH <- coxph(Surv(time_to_event, Event) ~ Age + sex + bmi + high_BP + family_diabetes, targetW1)
proteinCoxPH <- coxph(Surv(time_to_event, Event) ~ Age + sex + bmi + high_BP + family_diabetes + pScore, targetW1)
directCoxPH <- coxph(Surv(time_to_event, Event) ~ Age + sex + bmi + high_BP + family_diabetes + dScore, targetW1)
fullCoxPH <- coxph(Surv(time_to_event, Event) ~ Age + sex + bmi + high_BP + family_diabetes + pScore + dScore, targetW1)

models <- list(r = riskFactorsOnlyCoxPH, p = proteinCoxPH, d = directCoxPH, f = fullCoxPH)

predictCoxPHOnset <- function(dataDF, coxPHModel, threshold = 10) {
  uniqueTimes <- sort(unique(targetW1$time_to_event))
  thresholdIndex <- match(threshold, uniqueTimes)
  cumulativeBaseHaz <- gbm::basehaz.gbm(dataDF$time_to_event, dataDF$Event, predict(coxPHModel), uniqueTimes)
  survivalPredictions <- exp(-cumulativeBaseHaz[[thresholdIndex]]) ^ exp(predict(coxPHModel))
  onsetPredictions <- 1 - survivalPredictions
  auc <- MLmetrics::AUC(y_pred = onsetPredictions, y_true = dataDF$Event)
  prauc <- MLmetrics::PRAUC(y_pred = onsetPredictions, y_true = dataDF$Event)
  roc <- pROC::roc(response = dataDF$Event, predictor = onsetPredictions)
  list(cumulativeBaseHaz = cumulativeBaseHaz, onsetPredictions = onsetPredictions, auc = auc, prauc = prauc, roc = roc)
}

testResults <- lapply(models, function(m) {predictCoxPHOnset(targetW1, m)})

aucs <- sapply(testResults, function(r) {r$auc})
praucs <- sapply(testResults, function(r) {r$prauc})
metricsTable <- data.frame(AUC = aucs, PRAUC = praucs)

row.names(metricsTable) <- c('Risk factors only', 'Risk factors + protein EpiScore', 'Risk factors + direct EpiScore', 'Full model')

sessionStartTimestamp <- getOption('mprSessionStartTimestamp')
sessionLogFolder <- getOption('mprSessionLogFolder')
folderPath <- paste0(sessionLogFolder, 'output_', sessionStartTimestamp, '/')

dir.create(paste0(folderPath))

saveRDS(testResults, paste0(folderPath, 'testResults.rds'))
saveRDS(metricsTable, paste0(folderPath, 'metricsTable.rds'))
saveRDS(targetW1, paste0(folderPath, 'targetW1.rds'))
saveRDS(methylW1, paste0(folderPath, 'methylW1.rds'))
saveRDS(modelDirectT2DW4W3, paste0(folderPath, 'modelDirectT2DW4W3.rds'))


methylW4W3 <- NULL
methylW4 <- NULL
methylW3 <- NULL
methylW1 <- NULL

gc()

