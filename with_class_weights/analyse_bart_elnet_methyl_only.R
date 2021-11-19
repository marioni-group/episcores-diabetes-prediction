library(pROC)

library(MethylPipeR)

runInformation <- beginPipelineRun(note = 'Analysis of BART survival model from 2021_05_04_10_13_42 with features selected with logistic elnet - with class weights, methylation only. NA covariate rows removed. Train on wave 3, test on wave 1. Test data is filtered by TTE threshold. Deaths removed from controls. event.predictions object saved.', 
                                   logFolderPath = '/path/to/log/folder/',
                                   log = TRUE)

selectedFeatures <- read.csv(paste0(runInformation[['log folder path']], 'logistic_elnet_model_non_zero_model_coefficients_2021_04_01_13_22_48.csv'))[-1, 'X']

testData <- loadData('/path/to/test/data')
testTarget <- loadData('/path/to/test/target')

testRowsToKeep <- !is.na(testTarget$family_diabetes)
testData <- testData[testRowsToKeep, ]
testTarget <- testTarget[testRowsToKeep, ]
row.names(testTarget) <- NULL

testData <- cbind(testData, as.matrix(testTarget[, c('family_diabetes', 'high_BP')]))

testData <- testData[, selectedFeatures]
gc()

removeDeathsFromControlsResult <- removeDeathsFromControls(testTarget, list(testData))
testData <- removeDeathsFromControlsResult$objectsFiltered[[1]]
testTarget <- removeDeathsFromControlsResult$targetFiltered
removeDeathsFromControlsResult <- NULL
gc()

# Resetting row names so they match the row number after removal of deaths from controls.
rownames(testTarget) <- NULL

bartSurvivalModel <- readRDS(paste0(runInformation[['log folder path']], 'bart_survival_model_2021_05_04_10_13_42.rds'))

predictionTimePoint <- 10
numberOfTimePoints <- bartSurvivalModel$K
predictionTimePointIndex <- match(predictionTimePoint, bartSurvivalModel$times)

# Given an index for an individual (a row in the dataset), a time point and the total number of time points, returns the corresponding column index in bart.survival.model$surv.test (or any result with the same structure).
getBartResultColumn <- function(individual, timePointIndex, nTimePoints) {
  (individual - 1) * nTimePoints + timePointIndex
} 

survivalMeanPredictions <- sapply(1:nrow(testData), function(individual) {
  bartSurvivalModel$surv.test.mean[[getBartResultColumn(individual, predictionTimePointIndex, numberOfTimePoints)]]
})

thresholdTTEResult <- thresholdTTE(testTarget,
				      list(testData,
					   survivalMeanPredictions),
				      predictionTimePoint)
testTarget <- thresholdTTEResult$targetFiltered
row.names(testTarget) <- NULL
testData <- thresholdTTEResult$objectsFiltered[[1]]
survivalMeanPredictions <- thresholdTTEResult$objectsFiltered[[2]]
thresholdTTECounts <- thresholdTTEResult$counts
thresholdTTEResult <- NULL
gc()

# event.probability is calculated as 1 - survival probability
eventPredictions <- 1 - survivalMeanPredictions

if (runInformation[['log']]) {
  saveRDS(eventPredictions, paste0(runInformation[['log folder path']], 'bart_survival_model_from_2021_05_04_10_13_42_event_predictions_', runInformation[['start timestamp']], '.rds'))
}

rocResult <- roc(testTarget[, 'Event'], eventPredictions)
# plot(rocResult)

aucResult <- auc(rocResult)

tictoc::tic(paste0('Test AUC: ', aucResult))
tictoc::toc(log = TRUE)

endPipelineRun(runInformation)
