library(MethylPipeR)
library(formattable)
runInformation <- beginPipelineRun(note = 'Generate test results tables. All models', log = TRUE,
       logFolderPath = '/path/to/log/folder/')

logFolderPath <- '/path/to/log/folder/'

logisticNullTestResponse <- readRDS(paste0(logFolderPath, 'test_response_logistic_model_predictions_2021_04_16_12_17_04.rds'))
logisticLassoTestResponse <- readRDS(paste0(logFolderPath, 'test_response_logistic_model_predictions_2021_04_13_19_18_37.rds'))
logisticElnetTestResponse <- readRDS(paste0(logFolderPath, 'test_response_logistic_model_predictions_2021_04_13_19_09_15.rds'))

classBartLassoTestResponse <- readRDS(paste0(logFolderPath, 'class_bart_test_response_2021_04_27_11_35_53.rds'))
classBartElnetTestResponse <- readRDS(paste0(logFolderPath, 'class_bart_test_response_2021_04_27_11_40_47.rds'))

survivalBartLassoTestResponse <- readRDS(paste0(logFolderPath, 'bart_survival_model_from_2021_04_14_10_42_19_event_predictions_2021_04_14_17_06_06.rds'))
survivalBartElnetTestResponse <- readRDS(paste0(logFolderPath, 'bart_survival_model_from_2021_04_15_15_27_59_event_predictions_2021_04_15_23_00_05.rds'))

randomForestCVLassoTestResponse <- readRDS(paste0(logFolderPath, 'random_forest_lasso_test_response_2021_04_23_16_37_18.rds'))
randomForestCVElnetTestResponse <- readRDS(paste0(logFolderPath, 'random_forest_elnet_test_response_2021_04_23_16_34_24.rds'))

survivalRandomForestLassoTestResponse <- readRDS(paste0(logFolderPath, 'survival_random_forest_test_response_2021_04_15_18_10_08.rds'))
survivalRandomForestElnetTestResponse <- readRDS(paste0(logFolderPath, 'survival_random_forest_test_response_2021_04_15_19_08_02.rds'))

coxLassoTestResponse <- readRDS(paste0(logFolderPath, 'cox_lasso_onset_predictions_test_2021_04_14_12_19_11.rds'))
coxElnetTestResponse <- readRDS(paste0(logFolderPath, 'cox_elnet_onset_predictions_test_2021_04_14_11_19_15.rds'))

testResponses <- list(logistic_null = logisticNullTestResponse,
                      logistic_lasso = logisticLassoTestResponse,
                      logistic_elnet = logisticElnetTestResponse,
                      class_bart_lasso = classBartLassoTestResponse,
                      class_bart_elnet = classBartElnetTestResponse,
                      survival_bart_lasso = survivalBartLassoTestResponse,
                      survival_bart_elnet = survivalBartElnetTestResponse,
                      random_forest_cv_lasso = randomForestCVLassoTestResponse,
                      random_forest_cv_elnet = randomForestCVElnetTestResponse,
                      survival_random_forest_lasso = survivalRandomForestLassoTestResponse,
                      survival_random_forest_elnet = survivalRandomForestElnetTestResponse,
                      cox_lasso = coxLassoTestResponse,
                      cox_elnet = coxElnetTestResponse)

# Load Wave 1 target
w1Target <- loadData('path/to/w1/target')
w1Target <- w1Target[!is.na(w1Target$family_diabetes),]
row.names(w1Target) <- NULL

removeDeathsFromControlsResult <- removeDeathsFromControls(w1Target, NULL)
w1Target <- removeDeathsFromControlsResult$targetFiltered

threshold <- 10

thresholdTTEResult <- thresholdTTE(w1Target,
                                      NULL,
                                      threshold)

w1Target <- thresholdTTEResult$targetFiltered
row.names(w1Target) <- NULL

calculateAUC <- function(response, target) {
  pROC::roc(target, response)$auc
}

calculatePRAUC <- function(response, target) {
  MLmetrics::PRAUC(response, target)
}

testAUCs <- lapply(testResponses, function(response) {
  signif(unlist(calculateAUC(response, w1Target$Event)), digits = 3)
})

testPRAUCs <- lapply(testResponses, function(response) {
  signif(unlist(calculatePRAUC(response, w1Target$Event)), digits = 3)
})

testModelNames <- c('Logistic Null (Lasso)', 'Logistic Lasso', 'Logistic Elastic-net',
                    'Classification BART Lasso', 'Classification BART Elastic-Net', 'Survival BART Lasso', 'Survival BART Elastic-Net', 'Random Forest CV Lasso', 'Random Forest CV Elastic-net', 'Survival Random Forest Lasso',
                    'Survival Random Forest Elastic-net', 'Cox Lasso', 'Cox Elastic-net')

testTable <- data.frame(Model = testModelNames, AUC = unlist(testAUCs), PRAUC = unlist(testPRAUCs))

formattedTestTable <- formattable(testTable, align = c("l",rep("r", NCOL(testTable) - 1)), list(
  `Model` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")),
  area(col = c(2)) ~ color_tile("#DeF7E9", "#71CA97"), area(col = c(3)) ~ color_tile("#DeF7E9", "#71CA97")))

# Incremental models:
incrementalSurvivalRandomForestLassoResponse <- readRDS(paste0(logFolderPath, 'incremental_srf_lasso_cv_full_model_response_2021_05_05_01_47_49.rds'))
incrementalSurvivalRandomForestElnetResponse <- readRDS(paste0(logFolderPath, 'incremental_srf_elnet_cv_full_model_response_2021_05_05_01_02_03.rds'))
incrementalRandomForestLassoResponse <- readRDS(paste0(logFolderPath, 'incremental_rf_lasso_cv_full_model_response_2021_05_04_14_25_51.rds'))
incrementalRandomForestElnetResponse <- readRDS(paste0(logFolderPath, 'incremental_rf_elnet_cv_full_model_response_2021_05_04_15_17_49.rds'))
incrementalCoxLassoResponse <- readRDS(paste0(logFolderPath, 'incremental_cox_lasso_full_model_response_2021_05_04_10_44_41.rds'))
incrementalCoxElnetResponse <- readRDS(paste0(logFolderPath, 'incremental_cox_elnet_full_model_response_2021_05_04_14_49_10.rds'))
incrementalClassBartLassoResponse <- readRDS(paste0(logFolderPath, 'incremental_class_bart_lasso_full_model_response_2021_05_04_12_51_50.rds'))
incrementalClassBartElnetResponse <- readRDS(paste0(logFolderPath, 'incremental_class_bart_elnet_full_model_response_2021_05_04_13_17_20.rds'))
incrementalBartElnetResponse <- readRDS(paste0(logFolderPath, 'incremental_bart_elnet_full_model_response_2021_05_04_14_07_15.rds'))
incrementalBartLassoResponse <- readRDS(paste0(logFolderPath, 'incremental_bart_lasso_full_model_response_2021_04_26_11_48_56.rds'))
incrementalLogisticLassoResponse <- readRDS(paste0(logFolderPath, 'incremental_logistic_lasso_full_model_response_2021_04_01_13_03_03.rds'))
incrementalLogisticElnetResponse <- readRDS(paste0(logFolderPath, 'incremental_logistic_elnet_full_model_response_2021_04_01_13_22_48.rds'))
incrementalLogisticNullResponse <- readRDS(paste0(logFolderPath, 'incremental_logistic_lasso_null_model_response_2021_04_01_13_03_03.rds'))

incrementalSurvivalRandomForestLassoModel <- readRDS(paste0(logFolderPath, 'incremental_srf_lasso_cv_full_model_2021_05_05_01_47_49.rds'))
incrementalSurvivalRandomForestElnetModel <- readRDS(paste0(logFolderPath, 'incremental_srf_elnet_cv_full_model_2021_05_05_01_02_03.rds'))
incrementalRandomForestLassoModel <- readRDS(paste0(logFolderPath, 'incremental_rf_lasso_cv_full_model_2021_05_04_14_25_51.rds'))
incrementalRandomForestElnetModel <- readRDS(paste0(logFolderPath, 'incremental_rf_elnet_cv_full_model_2021_05_04_15_17_49.rds'))
incrementalCoxLassoModel <- readRDS(paste0(logFolderPath, 'incremental_cox_lasso_full_model_2021_05_04_10_44_41.rds'))
incrementalCoxElnetModel <- readRDS(paste0(logFolderPath, 'incremental_cox_elnet_full_model_2021_05_04_14_49_10.rds'))
incrementalClassBartLassoModel <- readRDS(paste0(logFolderPath, 'incremental_class_bart_lasso_full_model_2021_05_04_12_51_50.rds'))
incrementalClassBartElnetModel <- readRDS(paste0(logFolderPath, 'incremental_class_bart_elnet_full_model_2021_05_04_13_17_20.rds'))
incrementalBartElnetModel <- readRDS(paste0(logFolderPath, 'incremental_bart_elnet_full_model_2021_05_04_14_07_15.rds'))
incrementalBartLassoModel <- readRDS(paste0(logFolderPath, 'incremental_bart_lasso_full_model_2021_04_26_11_48_56.rds'))
incrementalLogisticLassoModel <- readRDS(paste0(logFolderPath, 'incremental_logistic_lasso_full_model_2021_04_01_13_03_03.rds'))
incrementalLogisticElnetModel <- readRDS(paste0(logFolderPath, 'incremental_logistic_elnet_full_model_2021_04_01_13_22_48.rds'))
incrementalLogisticNullModel <- readRDS(paste0(logFolderPath, 'incremental_logistic_lasso_null_model_2021_04_01_13_03_03.rds'))

incrementalResponses <- list(srf_lasso = incrementalSurvivalRandomForestLassoResponse,
                             srf_elnet = incrementalSurvivalRandomForestElnetResponse,
                             rf_lasso = incrementalRandomForestLassoResponse,
                             rf_elnet = incrementalRandomForestElnetResponse,
                             cox_lasso = incrementalCoxLassoResponse,
                             cox_elnet = incrementalCoxElnetResponse,
                             class_bart_lasso = incrementalClassBartLassoResponse,
                             class_bart_elnet = incrementalClassBartElnetResponse,
                             bart_lasso = incrementalBartLassoResponse,
                             bart_elnet = incrementalBartElnetResponse,
                             logistic_lasso = incrementalLogisticLassoResponse,
                             logistic_elnet = incrementalLogisticElnetResponse,
                             logistic_null = incrementalLogisticNullResponse)

incrementalModels <- list(srf_lasso = incrementalSurvivalRandomForestLassoModel,
                          srf_elnet = incrementalSurvivalRandomForestElnetModel,
                          rf_lasso = incrementalRandomForestLassoModel,
                          rf_elnet = incrementalRandomForestElnetModel,
                          cox_lasso = incrementalCoxLassoModel,
                          cox_elnet = incrementalCoxElnetModel,
                          class_bart_lasso = incrementalClassBartLassoModel,
                          class_bart_elnet = incrementalClassBartElnetModel,
                          bart_lasso = incrementalBartLassoModel,
                          bart_elnet = incrementalBartElnetModel,
                          logistic_lasso = incrementalLogisticLassoModel,
                          logistic_elnet = incrementalLogisticElnetModel,
                          logistic_null = incrementalLogisticNullModel)

incrementalAUCs <- sapply(incrementalResponses, function(response) {
  signif(unlist(calculateAUC(response, w1Target$Event)), digits = 3)
})

incrementalPRAUCs <- sapply(incrementalResponses, function(response) {
  signif(unlist(calculatePRAUC(response, w1Target$Event)), digits = 3)
})

incrementalMcFaddenAdjs <- sapply(incrementalModels, function(model) {
  signif(DescTools::PseudoR2(model, 'McFaddenAdj'), digits = 3)
})

incrementalMethylPValues <- sapply(incrementalModels, function(model) {
  coefData <- coef(summary(model))
  if ('methyl' %in% row.names(coefData)) {
    signif(coefData['methyl', 4], digits = 3)
  } else {
    NA
  }
})

incrementalNames <- c('Survival Random Forest Lasso', 'Survival Random Forest Elastic-net', 'Random Forest Lasso', 'Random Forest Elastic-net', 'Cox Lasso', 'Cox Elastic-net', 'Class BART Lasso', 'Class BART Elastic-net', 'Survival BART Lasso', 'Survival BART Elastic-net', 'Logistic Lasso', 'Logistic Elastic-net', # 'Logistic Ridge', 
                      'Logistic Null')

incrementalTable <- data.frame(Model = incrementalNames, AUC = incrementalAUCs, PRAUC = incrementalPRAUCs, "McFaddenAdjR2" = incrementalMcFaddenAdjs, 'MethylPValue' = incrementalMethylPValues)

# Sort by Incremental R2
sortIndex <- order(incrementalTable$McFaddenAdjR2, decreasing = TRUE)
incrementalTable <- incrementalTable[order(incrementalTable$McFaddenAdjR2, decreasing = TRUE),]

incrementalResponses <- incrementalResponses[sortIndex]
incrementalModels <- incrementalModels[sortIndex]
incrementalNames <- incrementalNames[sortIndex]

row.names(incrementalTable) <- NULL


formattedIncrementalTable <- formattable(incrementalTable, align = c("l",rep("r", NCOL(incrementalTable) - 1)), list(
  `Model` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")), 
  area(col = c(2)) ~ color_tile("#DeF7E9", "#71CA97"), area(col = c(3)) ~ color_tile("#DeF7E9", "#71CA97"), area(col = c(4)) ~ color_tile("#DeF7E9", "#71CA97"), area(col = c(5)) ~ color_tile("#DeF7E9", "#71CA97")))

if (runInformation[['log']]) {
  saveRDS(testTable, paste0(runInformation[['log folder path']], 'ml_test_table_', runInformation[['start timestamp']], '.rds'))
  saveRDS(incrementalTable, paste0(runInformation[['log folder path']], 'incremental_table_', runInformation[['start timestamp']], '.rds'))
}

numberOfMethods <- length(incrementalModels)
labels <- rep(list(w1Target[, 'Event']), numberOfMethods)
names(labels) <- names(incrementalResponses)

print('Plotting ROCs')

rocPlot <- multipleROC(labels, incrementalResponses, 'Incremental Models Prediction ROC', runInformation)
prPlot <- multipleROC(labels, incrementalResponses, 'Incremental Models Prediction Precision-Recall', runInformation, 'prec', 'rec', auc = 'PR')

endPipelineRun(runInformation)
