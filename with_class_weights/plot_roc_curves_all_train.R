library(MethylPipeR)
runInformation <- beginPipelineRun(note = 'Plot train ROC curves. All models (with class weights applied)', log = FALSE,
		   logFolderPath = 'path/to/log/folder/')

w3Target <- loadData('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/incident_diabetes_pipeline/data/test_target.csv')

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

lassoResponse <- readRDS(paste0(runInformation[['log folder path']], 'train_response_logistic_model_predictions_2021_03_09_12_25_07.rds'))
elnetResponse <- readRDS(paste0(runInformation[['log folder path']], 'train_response_logistic_model_predictions_2021_03_09_11_50_27.rds'))
ridgeResponse <- readRDS(paste0(runInformation[['log folder path']], 'train_response_logistic_model_predictions_2021_03_09_12_38_54.rds'))
lassoNullResponse <- readRDS(paste0(runInformation[['log folder path']], 'train_response_logistic_model_predictions_2021_03_09_12_49_53.rds'))
srfResponse <- readRDS(paste0(runInformation[['log folder path']], 'survival_random_forest_train_response_2021_03_18_11_49_31.rds'))
coxElnetResponse <- readRDS(paste0(runInformation[['log folder path']], 'cox_elnet_onset_predictions_train_2021_03_16_10_13_40.rds'))


methods <- c('survival_random_forest', 'cox_elnet', 'logistic_elnet', 'logistic_lasso', 'logistic_lasso_null')
responses <- list(srfResponse, coxElnetResponse, elnetResponse, lassoResponse, lassoNullResponse)

numberOfMethods <- length(methods)
labels <- rep(list(w3Target[, 'Event']), numberOfMethods)

names(responses) <- methods
names(labels) <- methods

rocPlot <- multipleROC(labels, responses, 'Training Set Prediction ROC', runInformation)
prPlot <- multipleROC(labels, responses, 'Training Set Prediction Precision-Recall', runInformation, 'prec', 'rec', auc = 'PR')

calibrationDF <- data.frame(responses)
calibrationDF$actual <- factor(w3Target[, 'Event'], levels = c(1, 0))

calibrationResult <- caret::calibration(actual ~ logistic_lasso + logistic_elnet + logistic_lasso_null + survival_random_forest + cox_elnet,
                                        data = calibrationDF)

logisticCalibrationResult <- caret::calibration(actual ~ logistic_lasso + logistic_elnet + logistic_lasso_null, data = calibrationDF)
otherCalibrationResult <- caret::calibration(actual ~ survival_random_forest + cox_elnet, data = calibrationDF)

endPipelineRun(runInformation)
