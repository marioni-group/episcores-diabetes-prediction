library(MethylPipeR)
library(DescTools)
runInformation <- beginPipelineRun(note = 'Plot incremental ROC curves.', log = FALSE,
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

lassoResponse <- readRDS(paste0(runInformation[['log folder path']], 'incremental_logistic_lasso_full_model_response_2021_03_31_10_31_42.rds'))
elnetResponse <- readRDS(paste0(runInformation[['log folder path']], 'incremental_logistic_elnet_full_model_response_2021_03_31_10_45_31.rds'))
bartLassoResponse <- readRDS(paste0(runInformation[['log folder path']], 'incremental_bart_lasso_full_model_response_2021_03_31_10_46_34.rds'))
coxLassoResponse <- readRDS(paste0(runInformation[['log folder path']], 'incremental_cox_lasso_full_model_response_2021_03_31_12_31_29.rds'))
nullResponse <- readRDS(paste0(runInformation[['log folder path']], 'incremental_logistic_lasso_null_model_response_2021_03_31_10_31_42.rds'))

methods <- c('lasso', 'elnet', 'bart_lasso', 'cox_lasso', 'null')
responses <- list(lassoResponse, elnetResponse, bartLassoResponse, coxLassoResponse, nullResponse)


numberOfMethods <- length(methods)
labels <- rep(list(w1Target[, 'Event']), numberOfMethods)

names(responses) <- methods
names(labels) <- methods

rocPlot <- multipleROC(labels, responses, 'Incremental Models Prediction ROC', runInformation)
prPlot <- multipleROC(labels, responses, 'Incremental Models Prediction Precision-Recall', runInformation, 'prec', 'rec', auc = 'PR')

calibrationDF <- data.frame(responses)
calibrationDF$actual <- factor(w1Target[, 'Event'], levels = c(1, 0))

calibrationResult <- caret::calibration(actual ~ lasso + elnet + bart_lasso + cox_lasso + null,
                                        data = calibrationDF)

# Model analysis
lassoModel <- readRDS(paste0(runInformation[['log folder path']], 'incremental_logistic_lasso_full_model_2021_03_31_10_31_42.rds'))
elnetModel <- readRDS(paste0(runInformation[['log folder path']], 'incremental_logistic_elnet_full_model_2021_03_31_10_45_31.rds'))
bartLassoModel <- readRDS(paste0(runInformation[['log folder path']], 'incremental_bart_lasso_full_model_2021_03_31_10_46_34.rds'))
coxLassoModel <- readRDS(paste0(runInformation[['log folder path']], 'incremental_cox_lasso_full_model_2021_03_31_12_31_29.rds'))
nullModel <- readRDS(paste0(runInformation[['log folder path']], 'incremental_logistic_lasso_null_model_2021_03_31_10_31_42.rds'))

models <- list(lasso = lassoModel, elnet = elnetModel, bart_lasso = bartLassoModel, cox_lasso = coxLassoModel, null = nullModel)
summaries <- lapply(models, summary)

adjR2McFadden <- lapply(models, function(model) {PseudoR2(model, 'McFaddenAdj')})
r2Nagelkerke <- lapply(models, function(model) {PseudoR2(model, 'Nagelkerke')})


endPipelineRun(runInformation)
