# Required packages: pROC, MLmetrics, DescTools, caret

# Loading in model coefficients. Make sure these files are present in the current working directory or change the paths to the correct directory.
directEpiScoreCoefficients <- readRDS('directEpiScoreCoefficients.rds')
dModelCoefficients <- readRDS('dModelCoefficients.rds')
nullModelCoefficients <- readRDS('nullModelCoefficients.rds')

cumulativeBaselineHazardNull <- 0.02596089 # Baseline cumulative hazard at t = 10.
cumulativeBaselineHazardDirect <- 0.02033703

# Insert column names for each of the variables corresponding to the names in the data table.
ageColname <- '' # Age in years
sexColname <- '' # Binary 1=male, 0=female
bmiColname <- '' # BMI
familyDiabetesColname <- '' # Parent/sibling diabetes (binary). 1=yes, 0=no
hypertensionColname <- '' # History of hypertension (binary). 1=yes, 0=no

diabetesColname <- '' # Incident diabetes onset (binary). 1=yes, 0=no

# Read in table with CpG data. Each column corresponds to a CpG. Each row corresponds to an individual.
methylationTable <-  # read in CpG table

# Read in table with phenotype data where each column corresponds to a variable listed above. Each row corresponds to an individuals.
phenotypeTable <- # read in phenotype table


## Edit above
#########################################################################################################

# Calculate direct EpiScore
directCpGColnames <- names(directEpiScoreCoefficients)

directExcludedCpGs <- setdiff(directCpGColnames, colnames(methylationTable))

directCpGColnames <- intersect(directCpGColnames, colnames(methylationTable))

directCpGs <- methylationTable[, directCpGColnames]

# Result of matrix multiplication is a nrow(data) x 1 matrix. Select first column for result as a vector.
directEpiScore <- (as.matrix(directCpGs) %*% as.numeric(directEpiScoreCoefficients[directCpGColnames]))[, 1]

age <- phenotypeTable[, ageColname]
sex <- phenotypeTable[, sexColname]
bmi <- phenotypeTable[, bmiColname]
familyDiabetes <- phenotypeTable[, familyDiabetesColname]
hypertension <- phenotypeTable[, hypertensionColname]
diabetes <- phenotypeTable[, diabetesColname]

meanAge <- 48.28064
meanSex <- 0.3796568
meanBMI <- 26.83208
meanFamilyDiabetes <- 0.2015488
meanHypertension <- 0.1391796
meanDirectEpiScore <- 12.43837

dLinearPredictor <- (age - meanAge) * dModelCoefficients[['Age']] +
                (sex - meanSex) * dModelCoefficients[['sexM']] +
                (bmi - meanBMI) * dModelCoefficients[['bmi']] +
                (familyDiabetes - meanFamilyDiabetes) * dModelCoefficients[['family_diabetes']] +
                (hypertension - meanHypertension) * dModelCoefficients[['high_BP']] +
                (directEpiScore - meanDirectEpiScore) * dModelCoefficients[['dScore']]

dSurvival <- exp(-cumulativeBaselineHazardDirect) ^ exp(dLinearPredictor)
dPrediction <- 1 - dSurvival

nullLinearPredictor <- (age - meanAge) * dModelCoefficients[['Age']] +
                (sex - meanSex) * dModelCoefficients[['sexM']] +
                (bmi - meanBMI) * dModelCoefficients[['bmi']] +
                (familyDiabetes - meanFamilyDiabetes) * dModelCoefficients[['family_diabetes']] +
                (hypertension - meanHypertension) * dModelCoefficients[['high_BP']]

nullSurvival <- exp(-cumulativeBaselineHazardNull) ^ exp(nullLinearPredictor)
nullPrediction <- 1 - nullSurvival


nullAUC <- pROC::roc(diabetes, nullPrediction)$auc
dAUC <- pROC::roc(diabetes, dPrediction)$auc

nullPRAUC <- MLmetrics::PRAUC(nullPrediction, diabetes)
dPRAUC <- MLmetrics::PRAUC(dPrediction, diabetes)


results <- list(nullAUC = nullAUC,
                dAUC = dAUC,
                nullPRAUC = nullPRAUC,
                dPRAUC = dPRAUC
           )

resultsDF <- data.frame(results)

write.csv(resultsDF, 'direct_kora_validation_results.csv')
write.csv(data.frame(excluded_cpgs = directExcludedCpGs), 'direct_kora_logistic_excluded_cpgs.csv')


dataDF <- data.frame(age = age, sex = sex, bmi = bmi, familyDiabetes = familyDiabetes, hypertension = hypertension, directEpiScore = directEpiScore, diabetes = diabetes)

logisticNullModel <- glm(diabetes ~ age + as.factor(sex) + bmi + as.factor(familyDiabetes) + as.factor(hypertension), family = 'binomial', data = dataDF)
logisticFullModel <- glm(diabetes ~ age + as.factor(sex) + bmi + as.factor(familyDiabetes) + as.factor(hypertension) + directEpiScore, family = 'binomial', data = dataDF)

logisticNullPredictions <- predict(logisticNullModel, type = 'response')
logisticFullPredictions <- predict(logisticFullModel, type = 'response')

nullLogisticAUC <- pROC::roc(diabetes, logisticNullPredictions)$auc
fullLogisticAUC <- pROC::roc(diabetes, logisticFullPredictions)$auc

nullLogisticPRAUC <- MLmetrics::PRAUC(logisticNullPredictions, diabetes)
fullLogisticPRAUC <- MLmetrics::PRAUC(logisticFullPredictions, diabetes)

logisticResults <- list(nullAUC = nullLogisticAUC,
                        dAUC = fullLogisticAUC,
                        nullPRAUC = nullLogisticPRAUC,
                        dPRAUC = fullLogisticPRAUC)

logisticResultsDF <- data.frame(logisticResults)
write.csv(logisticResultsDF, 'direct_kora_validation_logistic_results.csv')

# Calculate confusion matrix metrics

calculatePredictions <- function(response1, response2, actual) {
  thresholds <- seq(from = 0, to = 1, by = 0.1)
  
  binaryPredictions1 <- lapply(thresholds, function(threshold) {
    as.factor(as.numeric(response1 >= threshold))
  })

  names(binaryPredictions1) <- thresholds

  binaryPredictions2 <- lapply(thresholds, function(threshold) {
    as.factor(as.numeric(response2 >= threshold))
  })

  names(binaryPredictions2) <- thresholds

  confusionMatrices1 <- lapply(binaryPredictions1, function(x) {
    caret::confusionMatrix(x, actual, positive = '1')
  })

  confusionMatrices2 <- lapply(binaryPredictions2, function(x) {
    caret::confusionMatrix(x, actual, positive = '1')
  })

  metrics1 <- do.call(rbind, lapply(confusionMatrices1, function(confusionMatrix) {confusionMatrix[[4]][1:4]}))

  colnames(metrics1) <- c('Full Sensitivity', 'Full Specificity', 'Full PPV', 'Full NPV')

  metrics2 <- do.call(rbind, lapply(confusionMatrices2, function(confusionMatrix) {confusionMatrix[[4]][1:4]}))
  colnames(metrics2) <- c('Null Sensitivity', 'Null Specificity', 'Null PPV', 'Null NPV')

  metricsDF <- as.data.frame(cbind(metrics1, metrics2))

  rocTest <- pROC::roc.test(response = actual, predictor1 = response1, predictor2 = response2)$p.value

  list(metricsDF = metricsDF, rocTest = rocTest)
}

calculateAtPrevalence <- function(metricsTable, prevalence = 0.33, n = 10000) {

  metrics <- metricsTable

  nCases <- floor(n * prevalence)
  nControls <- n - nCases

  metrics$diffTP <- metrics[, 'Full Sensitivity'] * nCases - metrics[, 'Null Sensitivity'] * nCases
  metrics$diffFN <- (1 - metrics[, 'Full Sensitivity']) * nCases - (1 - metrics[, 'Null Sensitivity']) * nCases
  metrics$diffTN <- metrics[, 'Full Specificity'] * nControls - metrics[, 'Null Specificity'] * nControls
  metrics$diffFP <- (1 - metrics[, 'Full Specificity']) * nControls - (1 - metrics[, 'Null Specificity']) * nControls
  metrics$diffCorrect <- metrics$diffTP + metrics$diffTN

  metrics
}

calculateCalibration <- function(response1, response2, actual) {
  calibrationDF <- data.frame(nullResponse = response1, fullResponse = response2, actual = actual)

  calibrationResult <- caret::calibration(actual ~ fullResponse + nullResponse, data = calibrationDF)

  calibrationResult
}

metricsAndROCTest <- calculatePredictions(nullPrediction, dPrediction, as.factor(diabetes))
confusionMatrixMetrics <- calculateAtPrevalence(metricsAndROCTest$metricsDF)
calibrationTest <- calculateCalibration(nullPrediction, dPrediction, factor(diabetes, labels = c('1', '0'), levels = c(1, 0)))

metricsAndROCTestLogistic <- calculatePredictions(logisticNullPredictions, logisticFullPredictions, as.factor(diabetes))
confusionMatrixMetricsLogistic <- calculateAtPrevalence(metricsAndROCTestLogistic$metricsDF)
calibrationTestLogistic <- calculateCalibration(logisticNullPredictions, logisticFullPredictions, factor(diabetes, labels = c('1', '0'), levels = c(1, 0)))

resultsObject <- list(metricsAndROCTest = metricsAndROCTest,
                      confusionMatrixMetrics = confusionMatrixMetrics,
                      calibrationTest = calibrationTest,
                      metricsAndROCTestLogistic = metricsAndROCTestLogistic,
                      confusionMatrixMetricsLogistic = confusionMatrixMetricsLogistic,
                      calibrationTestLogistic = calibrationTestLogistic)

saveRDS(resultsObject, 'resultsObject.rds')

