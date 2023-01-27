targetW4 <- readRDS('/path/to/file')
targetW3 <- readRDS('/path/to/file')
targetW1 <- readRDS('/path/to/file')

# Filter out NA BMI rows from targetW1
targetW1NonNABMIIndex <- !is.na(targetW1$bmi)
targetW1 <- targetW1[targetW1NonNABMIIndex, ]
row.names(targetW1) <- NULL

targetW4W3 <- rbind(targetW4, targetW3)

trainingTarget <- targetW4W3
testTarget <- targetW1

trainingTargetControls <- trainingTarget[trainingTarget$Event == 0, ]
trainingTargetCases <- trainingTarget[trainingTarget$Event == 1, ]

testTargetControls <- testTarget[testTarget$Event == 0, ]
testTargetCases <- testTarget[testTarget$Event == 1, ]

targetW4Cases <- targetW4[targetW4$Event == 1, ]
targetW4Controls <- targetW4[targetW4$Event == 0, ]

targetW3Cases <- targetW3[targetW3$Event == 1, ]
targetW3Controls <- targetW3[targetW3$Event == 0, ]

targetW1Cases <- targetW1[targetW1$Event == 1, ]
targetW1Controls <- targetW1[targetW1$Event == 0, ]

tables <- list(trainingControls = trainingTargetControls, trainingCases = trainingTargetCases, testControls = testTargetControls, testCases = testTargetCases,
               targetW4Controls = targetW4Controls, targetW4Cases = targetW4Cases,
               targetW3Controls = targetW3Controls, targetW3Cases = targetW3Cases,
               targetW1Controls = targetW1Controls, targetW1Cases = targetW1Cases)

n <- sapply(tables, nrow)

meanTTE <- sapply(tables, function(t) {mean(t$time_to_event)})
sdTTE <- sapply(tables, function(t) {sd(t$time_to_event)})

meanAge <- sapply(tables, function(t) {mean(t$age_at_event)})
sdAge <- sapply(tables, function(t) {sd(t$age_at_event)})

nMale <- sapply(tables, function(t) {sum(t$sex == 'M')})
percentMale <- sapply(tables, function(t) {sum(t$sex == 'M') / nrow(t)})

meanBMI <- sapply(tables, function(t) {mean(t$bmi, na.rm = TRUE)})
sdBMI <- sapply(tables, function(t) {sd(t$bmi, na.rm = TRUE)})

nFamilyDiabetes <- sapply(tables, function(t) {sum(t$family_diabetes == 1, na.rm = TRUE)})
percentFamilyDiabetes <- sapply(tables, function(t) {sum(t$family_diabetes == 1, na.rm = TRUE) / nrow(t)})

nHypertension <- sapply(tables, function(t) {sum(t$high_BP == 1, na.rm = TRUE)})
percentHypertension <- sapply(tables, function(t) {sum(t$high_BP == 1, na.rm = TRUE) / nrow(t)})

df <- data.frame(n = n, meanTTE = meanTTE, sdTTE = sdTTE, meanAge = meanAge, sdAge = sdAge, nMale = nMale, percentMale = percentMale, meanBMI = meanBMI, sdBMI = sdBMI, nFamilyDiabetes = nFamilyDiabetes, percentFamilyDiabetes = percentFamilyDiabetes, nHypertension = nHypertension, percentHypertension = percentHypertension)

df <- as.data.frame(t(df))

write.csv(df, '/path/to/file')

