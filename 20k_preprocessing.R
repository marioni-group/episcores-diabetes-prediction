library(tictoc)
tic()

writeLines('Obtaining IDs for Wave 4+3 and Wave 1 individuals\n')
# Filter out wave 4 individuals with family members in wave 1

w4Target <- readRDS('/path/to/file')
gs10kTarget <- readRDS('/path/to/file')

w1Target <- gs10kTarget[gs10kTarget$Set == 'W1',]
w3Target <- gs10kTarget[gs10kTarget$Set == 'W3',]

pedigree <- read.csv('/path/to/file')

w4IDs <- w4Target$Sample_Name
w1IDs <- w1Target$Sample_Name
w3IDs <- w3Target$Sample_Name

writeLines(paste0('Initial n individuals in Wave 4 = ', length(w4IDs)))
writeLines(paste0('Initial n individuals in Wave 3 = ', length(w3IDs)))
writeLines(paste0('Initial n individuals in Wave 1 = ', length(w1IDs), '\n'))

coxTable <- read.csv('/path/to/file')
coxTableEvents <- coxTable[coxTable$Event == 1, ]

writeLines(paste0('Initial n cases individuals in Wave 4 = ', sum(sapply(w4IDs, function(x) {x %in% coxTableEvents$Sample_Name}))))
writeLines(paste0('Initial n cases individuals in Wave 3 = ', sum(sapply(w3IDs, function(x) {x %in% coxTableEvents$Sample_Name}))))
writeLines(paste0('Initial n cases individuals in Wave 1 = ', sum(sapply(w1IDs, function(x) {x %in% coxTableEvents$Sample_Name})), '\n'))

w1Pedigree <- pedigree[pedigree$volid %in% w1IDs, ]
w4Pedigree <- pedigree[pedigree$volid %in% w4IDs, ]

w4Pedigree$inW1 <- sapply(w4Pedigree$famid, function(fid) {fid %in% w1Pedigree$famid})

w4PedigreeFiltered <- w4Pedigree[!w4Pedigree$inW1, ]

w4IDs <- w4PedigreeFiltered$volid

writeLines(paste0('Individuals in Wave 4 after filtering out those with family in Wave 1 = ', length(w4IDs)))
writeLines(paste0('Cases in Wave 4 after filtering out those with family in Wave 1 = ', sum(sapply(w4IDs, function(x) {x %in% coxTableEvents$Sample_Name})), '\n'))

trainingIDs <- unique(c(w4IDs, w3IDs))
testIDs <- w1IDs

writeLines(paste0('Number of training IDs (Waves 4+3) = ', length(trainingIDs)))
writeLines(paste0('Number of test IDs (Wave 1) = ', length(testIDs), '\n'))

writeLines('Preprocessing Cox tables\n')

trainingCoxTable <- coxTable[sapply(coxTable$Sample_Name, function(sid) {sid %in% trainingIDs}), ]
testCoxTable <- coxTable[sapply(coxTable$Sample_Name, function(sid) {sid %in% testIDs}), ]

trainingCoxTable$Set <- sapply(trainingCoxTable$Sample_Name, function(sid) {
  if (sid %in% w4IDs) {
    'W4'
  } else if (sid %in% w3IDs) {
    'W3'
  } else {
    NA
  }
})

if (any(is.na(trainingCoxTable$Set))) {
  stop('Training Cox table contains individuals with ID not matching Wave 3 or 4')
}

testCoxTable$Set <- sapply(testCoxTable$Sample_Name, function(sid) {
  if (sid %in% w1IDs) {
    'W1'
  } else {
    NA
  }
})

if (any(is.na(testCoxTable$Set))) {
  stop('Test Cox table contains individuals with ID not matching Wave 1')
}

writeLines(paste0('Initial n individuals in training Cox table = ', nrow(trainingCoxTable)))
writeLines(paste0('Initial n cases in training Cox table = ', nrow(trainingCoxTable[trainingCoxTable$Event == 1,])))
writeLines(paste0('Initial n individuals in test Cox table = ', nrow(testCoxTable)))
writeLines(paste0('Initial n cases in test Cox table = ', nrow(testCoxTable[testCoxTable$Event == 1,]), '\n'))

bodyTable <- read.csv('/path/to/file')

# For the training table, all.x = TRUE as we are not removing individuals from the training set with NA covariates
trainingCoxTable <- merge(trainingCoxTable, bodyTable[, c('id', 'bmi')], all.x = TRUE, by.x = 'Sample_Name', by.y = 'id')
testCoxTable <- merge(testCoxTable, bodyTable[, c('id', 'bmi')], by.x = 'Sample_Name', by.y = 'id')

writeLines(paste0('Number of individuals in training Cox table after body table (BMI) merge = ', nrow(trainingCoxTable)))
writeLines(paste0('Number of cases in training Cox table after body table (BMI) merge = ', nrow(trainingCoxTable[trainingCoxTable$Event == 1, ])))
writeLines(paste0('Number of individuals in test Cox table after body table (BMI) merge = ', nrow(testCoxTable)))
writeLines(paste0('Number of cases in test Cox table after body table (BMI) merge = ', nrow(testCoxTable[testCoxTable$Event == 1, ]), '\n'))

diseaseTable <- read.csv('/path/to/file')

# For the training table, all.x = TRUE as we are not removing individuals from the training set with NA covariates
# trainingCoxTable <- merge(trainingCoxTable, diseaseTable[, c('ID', 'high_BP_Y')], all.x = TRUE, by.x = 'Sample_Name', by.y = 'ID')
# testCoxTable <- merge(testCoxTable, diseaseTable[, c('ID', 'high_BP_Y')], by.x = 'Sample_Name', by.y = 'ID')
trainingCoxTable$high_BP <- trainingCoxTable$high_BP_Y
# trainingCoxTable$high_BP_Y <- NULL

testCoxTable <- testCoxTable[!is.na(testCoxTable$high_BP_Y), ]
testCoxTable$high_BP <- testCoxTable$high_BP_Y
# testCoxTable$high_BP_Y <- NULL

writeLines(paste0('Number of individuals in training Cox table after disease table (high blood pressure) merge = ', nrow(trainingCoxTable)))
writeLines(paste0('Number of cases in training Cox table after disease table (high blood pressure) merge = ', nrow(trainingCoxTable[trainingCoxTable$Event == 1, ])))
writeLines(paste0('Number of individuals in test Cox table after disease table (high blood pressure) merge = ', nrow(testCoxTable)))
writeLines(paste0('Number of cases in test Cox table after disease table (high blood pressure) merge = ', nrow(testCoxTable[testCoxTable$Event == 1, ]), '\n'))

trainingCoxTable$family_diabetes <- apply(trainingCoxTable[, c('diabetes_M', 'diabetes_F', 'diabetes_BS', 'diabetes_B', 'diabetes_S')], 1, max)
testCoxTable$family_diabetes <- apply(testCoxTable[, c('diabetes_M', 'diabetes_F', 'diabetes_BS', 'diabetes_B', 'diabetes_S')], 1, max)

trainingCoxTable <- trainingCoxTable[!is.na(trainingCoxTable$time_to_event), ]
testCoxTable <- testCoxTable[!is.na(testCoxTable$time_to_event), ]

writeLines(paste0('Number of individuals in training Cox table after removing NA time_to_event = ', nrow(trainingCoxTable)))
writeLines(paste0('Number of cases in training Cox table after removing NA time_to_event = ', nrow(trainingCoxTable[trainingCoxTable$Event == 1, ])))
writeLines(paste0('Number of individuals in test Cox table after removing NA time_to_event = ', nrow(testCoxTable)))
writeLines(paste0('Number of cases in test Cox table after removing NA time_to_event = ', nrow(testCoxTable[testCoxTable$Event == 1, ]), '\n'))

trainingCoxTable <- trainingCoxTable[trainingCoxTable$time_to_event > 0, ]
testCoxTable <- testCoxTable[testCoxTable$time_to_event > 0, ]

writeLines(paste0('Number of individuals in training Cox table after removing time_to_event <= 0 = ', nrow(trainingCoxTable)))
writeLines(paste0('Number of cases in training Cox table after removing time_to_event <= 0 = ', nrow(trainingCoxTable[trainingCoxTable$Event == 1, ])))
writeLines(paste0('Number of individuals in test Cox table after removing NA time_to_event <= 0 = ', nrow(testCoxTable)))
writeLines(paste0('Number of cases in test Cox table after removing time_to_event <= 0 = ', nrow(testCoxTable[testCoxTable$Event == 1, ]), '\n'))

# Incorporate DNAm data

writeLines('Loading w4Methyl')
w4Methyl <- readRDS('/path/to/file')
writeLines('Transpose w4Methyl')
w4Methyl <- t(w4Methyl)

gc()

writeLines('Loading w3Methyl')
w3Methyl <- readRDS('/path/to/file')
writeLines('Transpose w3Methyl')
w3Methyl <- t(w3Methyl)

gc()

writeLines('Loading w1Methyl')
w1Methyl <- readRDS('/path/to/file')
writeLines('Transpose w1Methyl')
w1Methyl <- t(w1Methyl)

w4CpGs <- colnames(w4Methyl)
w3CpGs <- colnames(w3Methyl)
w1CpGs <- colnames(w1Methyl)


cpgIntersection <- intersect(intersect(w1CpGs, w3CpGs), w4CpGs)

writeLines('Filter w4Methyl to CpG intersection')
w4Methyl <- w4Methyl[, cpgIntersection]
writeLines('Filter w3Methyl to CpG intersection')
w3Methyl <- w3Methyl[, cpgIntersection]
writeLines('Filter w1Methyl to CpG intersection')
w1Methyl <- w1Methyl[, cpgIntersection]

gc()

# Function for converting M-values to beta values
mToBeta <- function(m) {
  2^m/(2^m + 1)
}

writeLines('Convert w4Methyl M to Beta')
w4Methyl <- mToBeta(w4Methyl)
writeLines('Convert w3Methyl M to Beta')
w3Methyl <- mToBeta(w3Methyl)
writeLines('Convert w1Methyl M to Beta')
w1Methyl <- mToBeta(w1Methyl)

gc()

# Mean imputation function for methylation data
imputeMean <- function(data) {
  for (i in 1:ncol(data)) {
    data[is.na(data[, i]), i] <- mean(data[, i], na.rm = TRUE)
  }
  data
}

writeLines('w4Methyl mean imputation')
w4Methyl <- imputeMean(w4Methyl)
writeLines('w3Methyl mean imputation')
w3Methyl <- imputeMean(w3Methyl)
writeLines('w1Methyl mean imputation')
w1Methyl <- imputeMean(w1Methyl)

gc()

w4SampleSentrixIDs <- data.frame(Sample_Name = w4Target$Sample_Name, Sample_Sentrix_ID = row.names(w4Target))
w3SampleSentrixIDs <- w3Target[, c('Sample_Name', 'Sample_Sentrix_ID')]
w1SampleSentrixIDs <- w1Target[, c('Sample_Name', 'Sample_Sentrix_ID')]

trainingSampleSentrixIDs <- rbind(w4SampleSentrixIDs, w3SampleSentrixIDs)
testSampleSentrixIDs <- w1SampleSentrixIDs

trainingCoxTable <- merge(trainingCoxTable, trainingSampleSentrixIDs, by = 'Sample_Name')
testCoxTable <- merge(testCoxTable, testSampleSentrixIDs, by = 'Sample_Name')

w4SentrixIntersection <- intersect(trainingCoxTable$Sample_Sentrix_ID, row.names(w4Methyl))
w3SentrixIntersection <- intersect(trainingCoxTable$Sample_Sentrix_ID, row.names(w3Methyl))
w1SentrixIntersection <- intersect(testCoxTable$Sample_Sentrix_ID, row.names(w1Methyl))

w4CoxTableIndex <- match(w4SentrixIntersection, trainingCoxTable$Sample_Sentrix_ID)
w3CoxTableIndex <- match(w3SentrixIntersection, trainingCoxTable$Sample_Sentrix_ID)
w1CoxTableIndex <- match(w1SentrixIntersection, testCoxTable$Sample_Sentrix_ID)

w4MethylIndex <- match(w4SentrixIntersection, row.names(w4Methyl))
w3MethylIndex <- match(w3SentrixIntersection, row.names(w3Methyl))
w1MethylIndex <- match(w1SentrixIntersection, row.names(w1Methyl))

writeLines('Creating w4CoxTable and filtering w4Methyl to match w4CoxTable by row')
w4CoxTable <- trainingCoxTable[w4CoxTableIndex, ]
w4Methyl <- w4Methyl[w4MethylIndex, ]

writeLines('Creating w3CoxTable and filtering w3Methyl to match w3CoxTable by row')
w3CoxTable <- trainingCoxTable[w3CoxTableIndex, ]
w3Methyl <- w3Methyl[w3MethylIndex, ]

writeLines('Creating w1CoxTable and filtering w1Methyl to match w1CoxTable by row')
w1CoxTable <- testCoxTable[w1CoxTableIndex, ]
w1Methyl <- w1Methyl[w1MethylIndex, ]

gc()

w4FinalNCases <- nrow(w4CoxTable[w4CoxTable$Event == 1, ])
w4FinalNControls <- nrow(w4CoxTable[w4CoxTable$Event == 0, ])
writeLines(paste0('Final n cases wave 4: ', w4FinalNCases, ', n controls: ', w4FinalNControls))

w3FinalNCases <- nrow(w3CoxTable[w3CoxTable$Event == 1, ])
w3FinalNControls <- nrow(w3CoxTable[w3CoxTable$Event == 0, ])
writeLines(paste0('Final n cases wave 3: ', w3FinalNCases, ', n controls: ', w3FinalNControls))


w1FinalNCases <- nrow(w1CoxTable[w1CoxTable$Event == 1, ])
w1FinalNControls <- nrow(w1CoxTable[w1CoxTable$Event == 0, ])
writeLines(paste0('Final n cases wave 1: ', w1FinalNCases, ', n controls: ', w1FinalNControls))


# writeLines('Save locally to reduce load times')

# writeLines('Wave 4')
# saveRDS(w4CoxTable, '/path/to/file')
# saveRDS(w4Methyl, '/path/to/file')

# writeLines('Wave 3')
# saveRDS(w3CoxTable, '/path/to/file')
# saveRDS(w3Methyl, '/path/to/file')

# writeLines('Wave 1')
# saveRDS(w1CoxTable, '/path/to/file')
# saveRDS(w1Methyl, '/path/to/file')

gc()

toc()
