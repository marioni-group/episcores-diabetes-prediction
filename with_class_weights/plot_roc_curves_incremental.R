library(MethylPipeR)
library(DescTools)
runInformation <- beginPipelineRun(note = 'Plot incremental ROC curves (models with class weights).', log = FALSE,
		   logFolderPath = '/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/incident_diabetes_pipeline/using_methylpiper/logs/with_class_weights_final_pipeline_instances/')

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
w1Target <- loadData('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/incident_diabetes_pipeline/data/training_target.csv')
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


lassoResponse <- readRDS(paste0(runInformation[['log folder path']], 'incremental_logistic_lasso_full_model_response_2021_04_01_13_03_03.rds'))
elnetResponse <- readRDS(paste0(runInformation[['log folder path']], 'incremental_logistic_elnet_full_model_response_2021_04_01_13_22_48.rds'))
nullResponse <- readRDS(paste0(runInformation[['log folder path']], 'incremental_logistic_elnet_null_model_response_2021_04_01_13_22_48.rds'))

methods <- c('lasso', 'elnet', 'null')
responses <- list(lassoResponse, elnetResponse, nullResponse)


numberOfMethods <- length(methods)
labels <- rep(list(w1Target[, 'Event']), numberOfMethods)

names(responses) <- methods
names(labels) <- methods

rocPlot <- multipleROC(labels, responses, 'Incremental Models Prediction ROC', runInformation)
prPlot <- multipleROC(labels, responses, 'Incremental Models Prediction Precision-Recall', runInformation, 'prec', 'rec', auc = 'PR')

calibrationDF <- data.frame(responses)
calibrationDF$actual <- factor(w1Target[, 'Event'], levels = c(1, 0))

calibrationResult <- caret::calibration(actual ~ lasso + elnet + null,
                                        data = calibrationDF)


# Model analysis
lassoModel <- readRDS(paste0(runInformation[['log folder path']], 'incremental_logistic_lasso_full_model_2021_04_01_13_03_03.rds'))
elnetModel <- readRDS(paste0(runInformation[['log folder path']], 'incremental_logistic_elnet_full_model_2021_04_01_13_22_48.rds'))
nullModel <- readRDS(paste0(runInformation[['log folder path']], 'incremental_logistic_elnet_null_model_2021_04_01_13_22_48.rds'))

models <- list(lasso = lassoModel, elnet = elnetModel, null = nullModel)
summaries <- lapply(models, summary)

adjR2McFadden <- lapply(models, function(model) {PseudoR2(model, 'McFaddenAdj')})
r2Nagelkerke <- lapply(models, function(model) {PseudoR2(model, 'Nagelkerke')})


endPipelineRun(runInformation)
