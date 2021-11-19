library(MethylPipeR)
runInformation <- beginPipelineRun(note = 'Plot test ROC curves. All models', log = FALSE,
		   logFolderPath = '/path/to/log/folder/')

w3Target <- loadData('/path/to/w3/target')

removeDeathsFromControlsResult <- removeDeathsFromControls(w3Target, NULL)
w3Target <- removeDeathsFromControlsResult$targetFiltered

threshold <- 10

thresholdTTEResult <- thresholdTTE(w3Target,
                                      NULL,
                                      threshold)

w3Target <- thresholdTTEResult$targetFiltered
row.names(w3Target) <- NULL

# Load Wave 1 target
w1Target <- loadData('/path/to/w1/target')
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

lassoResponse <- readRDS(paste0(runInformation[['log folder path']], 'test_response_logistic_model_predictions_2021_03_16_14_34_18.rds'))
elnetResponse <- readRDS(paste0(runInformation[['log folder path']], 'test_response_logistic_model_predictions_2021_03_16_14_17_33.rds'))
ridgeResponse <- readRDS(paste0(runInformation[['log folder path']], 'test_response_logistic_model_predictions_2021_03_16_14_46_51.rds'))
lassoNullResponse <- readRDS(paste0(runInformation[['log folder path']], 'test_response_logistic_model_predictions_2021_03_16_16_55_07.rds'))
bartSurvivalResponse <- readRDS(paste0(runInformation[['log folder path']], 'bart_survival_model_from_2021_03_17_10_51_34_event_predictions_2021_03_18_09_05_37.rds'))
srfResponse <- readRDS(paste0(runInformation[['log folder path']], 'survival_random_forest_test_response_2021_03_18_15_33_15.rds'))
coxElnetResponse <- readRDS(paste0(runInformation[['log folder path']], 'cox_elnet_onset_predictions_test_2021_03_18_09_10_52.rds'))

# Ordered by AUC
methods <- c('logistic_lasso_null', 'logistic_lasso', 'cox_elnet', 'logistic_elnet', 'survival_bart', 'survival_random_forest', 'logistic_ridge')
responses <- list(lassoNullResponse, lassoResponse, coxElnetResponse, elnetResponse, bartSurvivalResponse, srfResponse, ridgeResponse)

numberOfMethods <- length(methods)
labels <- rep(list(w1Target[, 'Event']), numberOfMethods)

names(responses) <- methods
names(labels) <- methods

rocPlot <- multipleROC(labels, responses, 'Test Set Prediction ROC', runInformation)
prPlot <- multipleROC(labels, responses, 'Test Set Prediction Precision-Recall', runInformation, 'prec', 'rec', auc = 'PR')

calibrationDF <- data.frame(responses)
calibrationDF$actual <- factor(w1Target[, 'Event'], levels = c(1, 0))

calibrationResult <- caret::calibration(actual ~ logistic_lasso + logistic_elnet + logistic_ridge + logistic_lasso_null + survival_bart + survival_random_forest + cox_elnet,
                                        data = calibrationDF)

logisticCalibrationResult <- caret::calibration(actual ~ logistic_lasso + logistic_elnet + logistic_ridge + logistic_lasso_null, data = calibrationDF)
otherCalibrationResult <- caret::calibration(actual ~ survival_bart + survival_random_forest + cox_elnet, data = calibrationDF)

endPipelineRun(runInformation)
