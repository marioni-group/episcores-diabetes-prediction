library(stringr)
library(glmnet)

set.seed(42)

methylW1 <- readRDS('/path/to/files')

proteinEpiScoresW3W1 <- read.csv('/path/to/file', row.names = 1)
proteinEpiScoresW3W1 <- as.matrix(proteinEpiScoresW3W1)
proteinEpiScoresW1 <- proteinEpiScoresW3W1[match(row.names(methylW1), row.names(proteinEpiScoresW3W1)), ]


targetW1 <- readRDS('/path/to/file')

hospTable <- read.csv('/path/to/file')
hospTable$dtString <- as.character(hospTable$dt)
hospTable$year <- sapply(hospTable$dtString, function(s) {str_sub(s, 1, 4)})
hospTable$month <- sapply(hospTable$dtString, function(s) {str_sub(s, 5,6)})

hospTable$dateDecimal <- as.integer(hospTable$year) + as.integer(hospTable$month)/12

testTable <- read.csv('/path/to/file')
testTable$dtString <- as.character(testTable$dt)
testTable$year <- sapply(testTable$dtString, function(s) {str_sub(s, 1, 4)})
testTable$month <- sapply(testTable$dtString, function(s) {str_sub(s, 5,6)})
testTable$dateDecimal <- as.integer(testTable$year) + as.integer(testTable$month)/12

positiveTestIndex <- sapply(targetW1$Sample_Name, function(id) {id %in% testTable$id})
targetW1 <- targetW1[positiveTestIndex, ]
methylW1 <- methylW1[positiveTestIndex, ]
proteinEpiScoresW1 <- proteinEpiScoresW1[positiveTestIndex, ]
row.names(targetW1) <- NULL

targetW1$date_of_event <- targetW1$yob + targetW1$mob/12 + targetW1$age_at_event

targetW1$covid_hosp <- sapply(targetW1$Sample_Name, function(id) {ifelse(id %in% hospTable$id, 1, 0)})

# Calculate closest test date to each individual
targetW1$closestTestDate <- sapply(1:nrow(targetW1), function(i) {
  id <- targetW1$Sample_Name[[i]]
  if (id %in% hospTable$id) {
    hospDecimal <- hospTable[hospTable$id == id, 'dateDecimal']
  } else {
    hospDecimal <- NULL
  }
  # For subsetting testTable to only rows that match id
  testTableIndex <- testTable$id == id
  testTableIndex[is.na(testTableIndex)] <- FALSE
  dateDecimals <- testTable[testTableIndex, 'dateDecimal']
  # print(dateDecimals)
  if (!is.null(hospDecimal)) {
    distances <- sapply(dateDecimals, function(d) {abs(d - hospDecimal)})
    # print(distances)
    minIndex <- which.min(distances)
    closestDate <- dateDecimals[minIndex]
  } else {
    closestDate <- min(dateDecimals)
  }
  closestDate
})

diabetesAfterCovid <- (targetW1$Event == 1) & (targetW1$date_of_event > targetW1$closestTestDate)
targetW1[diabetesAfterCovid, 'Event'] == 0

targetW1$baselineDecimal <- targetW1$date_of_event - targetW1$time_to_event
targetW1$timeToTest <- targetW1$closestTestDate - targetW1$baselineDecimal

# Calculate scores at t = timeToTest
directEpiScoreModel <- readRDS('/path/to/file')
proteinEpiScoreModel <- readRDS('/path/to/file')

uniqueTimes <- sort(unique(targetW1$timeToTest))

dPredictions <- predict(directEpiScoreModel$model, s = 'lambda.min', newx = methylW1)
pPredictions <- predict(proteinEpiScoreModel$model, s = 'lambda.min', newx = proteinEpiScoresW1)

dCumulativeBaseHaz <- gbm::basehaz.gbm(targetW1$time_to_event, targetW1$Event, dPredictions, uniqueTimes)
pCumulativeBaseHaz <- gbm::basehaz.gbm(targetW1$time_to_event, targetW1$Event, pPredictions, uniqueTimes)

directSurvivalPredictions <- sapply(1:nrow(targetW1), function(i) {
  thresholdIndex <- match(targetW1$timeToTest[[i]], uniqueTimes)
  cumulativeBaseHaz <- dCumulativeBaseHaz[[thresholdIndex]]
  dPrediction <- dPredictions[[i]]
  survivalPrediction <- exp(-cumulativeBaseHaz) ^ exp(dPrediction)
  onsetPrediction <- 1 - survivalPrediction
  onsetPrediction
})

targetW1$D <- directSurvivalPredictions

proteinSurvivalPredictions <- sapply(1:nrow(targetW1), function(i) {
  thresholdIndex <- match(targetW1$timeToTest[[i]], uniqueTimes)
  cumulativeBaseHaz <- pCumulativeBaseHaz[[thresholdIndex]]
  pPrediction <- pPredictions[[i]]
  survivalPrediction <- exp(-cumulativeBaseHaz) ^ exp(pPrediction)
  onsetPrediction <- 1 - survivalPrediction
  onsetPrediction
})

targetW1$P <- proteinSurvivalPredictions

targetW1 <- targetW1[!is.na(targetW1$D), ]

nullModel <- glm(covid_hosp ~ Age + sex + bmi + family_diabetes + high_BP + Event, family = 'binomial', data = targetW1)
pModel <- glm(covid_hosp ~ Age + sex + bmi + family_diabetes + high_BP + Event + P, family = 'binomial', data = targetW1)
dModel <- glm(covid_hosp ~ Age + sex + bmi + family_diabetes + high_BP + Event + D, family = 'binomial', data = targetW1)
fullModel <- glm(covid_hosp ~ Age + sex + bmi + family_diabetes + high_BP + Event + P + D, family = 'binomial', data = targetW1)

nullModelCoefSummary <- summary(nullModel)$coefficients
pModelCoefSummary <- summary(pModel)$coefficients
dModelCoefSummary <- summary(dModel)$coefficients
fullModelCoefSummary <- summary(fullModel)$coefficients

write.csv(nullModelCoefSummary, '/path/to/file')
write.csv(pModelCoefSummary, '/path/to/file')
write.csv(dModelCoefSummary, '/path/to/file')
write.csv(fullModelCoefSummary, '/path/to/file')

targetW1Cases <- targetW1[targetW1$covid_hosp == 1, ]
targetW1Controls <- targetW1[targetW1$covid_hosp == 0, ]

tables <- list(cases = targetW1Cases, controls = targetW1Controls)

n <- sapply(tables, nrow)

meanAgeAtBaseline <- sapply(tables, function(t) {mean(t$Age)})
sdAgeAtBaseline <- sapply(tables, function(t) {sd(t$Age)})

meanTimeToCovid <- sapply(tables, function(t) {mean(t$timeToTest)})
sdTimeToCovid <- sapply(tables, function(t) {sd(t$timeToTest)})

nMale <- sapply(tables, function(t) {sum(t$sex == 'M')})
percentMale <- sapply(tables, function(t) {sum(t$sex == 'M') / nrow(t)})

meanBMI <- sapply(tables, function(t) {mean(t$bmi, na.rm = TRUE)})
sdBMI <- sapply(tables, function(t) {sd(t$bmi, na.rm = TRUE)})

nFamilyDiabetes <- sapply(tables, function(t) {sum(t$family_diabetes)})
percentFamilyDiabetes <- sapply(tables, function(t) {sum(t$family_diabetes) / nrow(t)})

nHypertension <- sapply(tables, function(t) {sum(t$high_BP)})
percentHypertension <- sapply(tables, function(t) {sum(t$high_BP) / nrow(t)})

nT2D <- sapply(tables, function(t) {sum(t$Event)})
percentT2D <- sapply(tables, function(t) {sum(t$Event) / nrow(t)})

df <- data.frame(n = n, meanAgeAtBaseline = meanAgeAtBaseline, sdAgeAtBaseline = sdAgeAtBaseline, meanTimeToCovid = meanTimeToCovid, sdTimeToCovid = sdTimeToCovid, nMale = nMale, percentMale = percentMale, meanBMI = meanBMI, sdBMI = sdBMI, nFamilyDiabetes = nFamilyDiabetes, percentFamilyDiabetes = percentFamilyDiabetes, nHypertension = nHypertension, percentHypertension = percentHypertension, nT2D = nT2D, percentT2D = percentT2D)

df <- as.data.frame(t(df))

write.csv(df, '/path/to/file')
