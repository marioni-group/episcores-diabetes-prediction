library(MethylPipeR)
library(ggplot2)
library(precrec)
library(ROCR)
library(RColorBrewer)
library(MLmetrics)

w1Target <- readRDS('/path/to/file')
# w1Target$Event <- as.factor(w1Target$Event)

coxTestResults <- readRDS('/path/to/file')
bartTestResults <- readRDS('/path/to/file')
rsfTestResults <- readRDS('/path/to/file')

nullResponse <- coxTestResults$r$onsetPredictions
fullResponse <- coxTestResults$f$onsetPredictions
bartResponse <- bartTestResults$f$onsetPredictions
rsfResponse <- rsfTestResults$f$onsetPredictions

calibrationDFLasso <- data.frame(list('Full_model' = fullResponse, 'Risk_factors_only' = nullResponse))

w1Target$Event <- sapply(1:nrow(w1Target), function(i) {
  if (w1Target$time_to_event[[i]] > 10) {
    0
  } else {
    w1Target$Event[[i]]
  }
})

calibrationDFLasso$actual <- as.factor(w1Target$Event) # factor(w1Target$Event, levels = c(1, 0))

calibrationResultLasso <- caret::calibration(actual ~ Full_model + Risk_factors_only,
                                             data = calibrationDFLasso, class = '1')


pCox <- ggplot(calibrationResultLasso)


scores <- join_scores(nullResponse, bartResponse, rsfResponse, fullResponse)

curves <- evalmod(scores = scores, labels = w1Target$Event, calc_avg = FALSE, modnames = c('Risk factors only', 'Survival BART full model', 'Random survival forest full model', 'Cox full model'))
plot(curves)


multipleROC <- function(labelsList, predictionsList, title, yMetric = 'tpr', xMetric = 'fpr', auc = 'ROC', legendPosition = 'bottomright') {
  # Create ROCR prediction and performance objects
  methodsList <- names(labelsList)
  predictionObjects <- lapply(methodsList, function(method) {
    prediction(predictionsList[[method]], labelsList[[method]])
  })
  performanceObjects <- lapply(predictionObjects, function(predictionObject) {
    performance(predictionObject, yMetric, xMetric)
  })

  # aucValues <- lapply(predictionObjects, function(predictionObject) {
  #   performance(predictionObject, 'auc')@y.values[[1]]
  # })

  if (auc == 'ROC') {
    aucValues <- lapply(predictionObjects, function(predictionObject) {
      performance(predictionObject, 'auc')@y.values[[1]]
    })
  } else if (auc == 'PR') {
    aucValues <- lapply(methodsList, function(method) {
      PRAUC(predictionsList[[method]], labelsList[[method]])
    })
  }

  # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  # colourPalette = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  # Remove first 4 colours
  # colourPalette <- colourPalette[5:length(colourPalette)]

  colourPalette <- brewer.pal(n = 8, 'Dark2')

  plot(performanceObjects[[1]], col = colourPalette[[1]], main = title)
  if (length(methodsList) > 1) {
    lapply(2:length(methodsList), function(methodIndex) {
      plot(performanceObjects[[methodIndex]], col = colourPalette[[methodIndex]], add = TRUE)
    })
  }

  if (!is.null(auc)) {
    legendText <- lapply(1:length(methodsList), function(i) {
      paste0(methodsList[[i]], '; AUC = ', signif(aucValues[[i]], 3))
    })
  } else {
    legendText <- methodsList
  }


  # legendText <- lapply(1:length(methodsList), function(i) {
  #   paste0(methodsList[[i]], '; AUC = ', signif(aucValues[[i]], 3))
  # })

  legend(legendPosition, legend = legendText, bty = 'n', col = colourPalette[1:length(methodsList)], lty = rep(1, length(methodsList)),
         lwd = rep(5, length(methodsList))
         )
  # if (plotPRBaseline & auc == 'PR') {
  #   baselinePrecision <- sum(labelsList[[1]] == 1) / length(labelsList[[1]])
  #   abline(h = baselinePrecision)
  # }
  rocPlot <- recordPlot()

  list(rocPlot = rocPlot, performanceObjects = performanceObjects, aucValues = aucValues)
}

predictions <- list('Risk factors only' = nullResponse, 'Survival BART' = bartResponse, 'Random survival forest' = rsfResponse, 'Cox lasso' = fullResponse)
labels <- list('Risk factors only' = w1Target$Event, 'Survival BART' = w1Target$Event, 'Random survival forest' = w1Target$Event, 'Cox lasso' = w1Target$Event)
rocResult <- multipleROC(labels, predictions, 'Receiver Operating Characteristic')
prResult <- multipleROC(labels, predictions, 'Precision-Recall', 'prec', 'rec', auc = 'PR', legendPosition = 'topright')


gbm::calibrate.plot(y = w1Target$Event, p = fullResponse, distribution = 'bernoulli', shade.col = 'darkgrey')
# set.seed(42)
# pdf(file = 'full_model_calibration_curve.pdf')
# fullModelCalibrationCurve <- CalibrationCurves::val.prob.ci.2(fullResponse, w1Target$Event, statloc = FALSE, CL.BT = TRUE, d0lab = '', d1lab = '', legendloc = FALSE)
# dev.off()
# pdf(file = 'null_model_calibration_curve.pdf')
# nullModelCalibrationCurve <- CalibrationCurves::val.prob.ci.2(nullResponse, w1Target$Event, statloc = FALSE, CL.BT = TRUE, d0lab = '', d1lab = '', legendloc = FALSE)
# dev.off()

# Final
# set.seed(42)
# pdf(file = 'full_model_calibration_curve_zoomed.pdf')
# fullModelCalibrationCurveZoomed <- CalibrationCurves::val.prob.ci.2(fullResponse, w1Target$Event, statloc = FALSE, CL.BT = TRUE, d0lab = '', d1lab = '', legendloc = FALSE, xlim = c(-0.02, 0.2), ylim = c(-0.15, 0.4))
# dev.off()
# pdf(file = 'null_model_calibration_curve_zoomed.pdf')
# nullModelCalibrationCurve <- CalibrationCurves::val.prob.ci.2(nullResponse, w1Target$Event, statloc = FALSE, CL.BT = TRUE, d0lab = '', d1lab = '', legendloc = FALSE, xlim = c(-0.02, 0.2), ylim = c(-0.15, 0.4))
# dev.off()
