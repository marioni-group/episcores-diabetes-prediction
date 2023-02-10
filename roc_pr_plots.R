library(MethylPipeR)
library(ROCR)
library(RColorBrewer)
library(MLmetrics)

w1Target <- readRDS('/path/to/file')

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


multipleROC <- function(labelsList, predictionsList, title, yMetric = 'tpr', xMetric = 'fpr', auc = 'ROC', legendPosition = 'bottomright') {
  # Create ROCR prediction and performance objects
  methodsList <- names(labelsList)
  predictionObjects <- lapply(methodsList, function(method) {
    prediction(predictionsList[[method]], labelsList[[method]])
  })
  performanceObjects <- lapply(predictionObjects, function(predictionObject) {
    performance(predictionObject, yMetric, xMetric)
  })

  if (auc == 'ROC') {
    aucValues <- lapply(predictionObjects, function(predictionObject) {
      performance(predictionObject, 'auc')@y.values[[1]]
    })
  } else if (auc == 'PR') {
    aucValues <- lapply(methodsList, function(method) {
      PRAUC(predictionsList[[method]], labelsList[[method]])
    })
  }

  colourPalette <- brewer.pal(n = 8, 'Dark2')

  plot(performanceObjects[[1]], col = colourPalette[[1]], main = title, axes = FALSE)
  if (length(methodsList) > 1) {
    lapply(2:length(methodsList), function(methodIndex) {
      plot(performanceObjects[[methodIndex]], col = colourPalette[[methodIndex]], add = TRUE, axes = FALSE)
    })
  }

  if (!is.null(auc)) {
    legendText <- lapply(1:length(methodsList), function(i) {
      paste0(methodsList[[i]], '; AUC = ', signif(aucValues[[i]], 3))
    })
  } else {
    legendText <- methodsList
  }


  legend(legendPosition, legend = legendText, bty = 'n', col = colourPalette[1:length(methodsList)], lty = rep(1, length(methodsList)),
         lwd = rep(5, length(methodsList))
         )
  rocPlot <- recordPlot()

  list(rocPlot = rocPlot, performanceObjects = performanceObjects, aucValues = aucValues)
}

predictions <- list('Risk factors only' = nullResponse, 'Survival BART' = bartResponse, 'Random survival forest' = rsfResponse, 'Cox lasso' = fullResponse)
labels <- list('Risk factors only' = w1Target$Event, 'Survival BART' = w1Target$Event, 'Random survival forest' = w1Target$Event, 'Cox lasso' = w1Target$Event)

pdf(file = 't2d_nat_aging_figure2.pdf', width = 6, height = 12)
par(mfrow = c(2, 1))
rocResult <- multipleROC(labels, predictions, 'Receiver Operating Characteristic')
prResult <- multipleROC(labels, predictions, 'Precision-Recall', 'prec', 'rec', auc = 'PR', legendPosition = 'topright')
dev.off()

pdf(file = 'roc_plot.pdf', width = 6, height = 6)
rocResult <- multipleROC(labels, predictions, 'Receiver Operating Characteristic')
dev.off()

