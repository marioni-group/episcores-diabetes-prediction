library(MethylPipeR)

runInformation <- beginPipelineRun(note = 'BART survival model - features selected with logistic elnet (methylation only, with class weights), rows with na covariates removed (timestamp: 2021_04_01_13_22_48). Train on wave 3', 
                                   logFolderPath = '/path/to/log/folder/')

trainingData <- loadData('/path/to/training/data')
trainingTarget <- loadData('/path/to/training/target')

# Here we select rows -1 to exclude the intercept.
selectedFeatures <- read.csv(paste0(runInformation[['log folder path']], 'logistic_elnet_model_non_zero_model_coefficients_2021_04_01_13_22_48.csv'))[-1, 'X']

# trainingTarget$family_diabetes[is.na(trainingTarget$family_diabetes)] <- 0
# trainingTarget$high_BP[is.na(trainingTarget$high_BP)] <- 0
trainingData <- cbind(trainingData, as.matrix(trainingTarget[, c('family_diabetes', 'high_BP')]))
trainingData <- trainingData[, selectedFeatures]

trainingRowsToKeep <- !is.na(trainingTarget$family_diabetes)
trainingData <- trainingData[trainingRowsToKeep, ]
trainingTarget <- trainingTarget[trainingRowsToKeep, ]
row.names(trainingTarget) <- NULL

gc()

removeDeathsFromControlsResult <- removeDeathsFromControls(trainingTarget, list(trainingData))
trainingData <- removeDeathsFromControlsResult$objectsFiltered[[1]]
trainingTarget <- removeDeathsFromControlsResult$targetFiltered
removeDeathsFromControlsResult <- NULL
gc()

trainingFolds <- assignTrainingFolds(trainingTarget, 'Event', 3, runInformation)

testData <- loadData('/path/to/test/data')
testTarget <- loadData('/path/to/test/target')

# testTarget$family_diabetes[is.na(testTarget$family_diabetes)] <- 0
# testTarget$high_BP[is.na(testTarget$high_BP)] <- 0
testData <- cbind(testData, as.matrix(testTarget[, c('family_diabetes', 'high_BP')]))
testData <- testData[, selectedFeatures]

testRowsToKeep <- !is.na(testTarget$family_diabetes)
testData <- testData[testRowsToKeep, ]
testTarget <- testTarget[testRowsToKeep, ]
row.names(testTarget) <- NULL
gc()

removeDeathsFromControlsResult <- removeDeathsFromControls(testTarget, list(testData))
testData <- removeDeathsFromControlsResult$objectsFiltered[[1]]
testTarget <- removeDeathsFromControlsResult$targetFiltered
removeDeathsFromControlsResult <- NULL
gc()

bartSurvivalModel <- fitBARTSurvival(trainXs = trainingData, trainTarget = trainingTarget, testXs = testData, pipelineRunInformation = runInformation)

endPipelineRun(runInformation)
