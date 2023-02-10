library(MethylPipeR)
library(rms)
library(CalibrationCurves)


source('val.prob.ci.2.R')
source('auc.nonpara.mw.R')
source('ci.auc.R')
source('BT.samples.R')

w1Target <- readRDS('/path/to/file')

coxTestResults <- readRDS('/path/to/file')

nullResponse <- coxTestResults$r$onsetPredictions
fullResponse <- coxTestResults$f$onsetPredictions

w1Target$Event <- sapply(1:nrow(w1Target), function(i) {
  if (w1Target$time_to_event[[i]] > 10) {
    0
  } else {
    w1Target$Event[[i]]
  }
})

pdf(file = 'calibration_curve_full.pdf')
set.seed(42)
fullModelCalibrationCurve <- CalibrationCurves::val.prob.ci.2(fullResponse, w1Target$Event, statloc = FALSE, CL.BT = TRUE, d0lab = '', d1lab = '', legendloc = FALSE, xlab = '')
dev.off()

pdf(file = 'calibration_curve_full_zoomed.pdf')
set.seed(42)
fullModelCalibrationCurveZoomed <- val.prob.ci.2(fullResponse, w1Target$Event, statloc = FALSE, CL.BT = TRUE, d0lab = '', d1lab = '', legendloc = FALSE, xlim = c(-0.02, 0.2), ylim = c(-0.15, 0.4), xlab = '', ylab = '')
dev.off()

pdf(file = 'calibration_curve_null.pdf')
set.seed(42)
nullModelCalibrationCurve <- CalibrationCurves::val.prob.ci.2(nullResponse, w1Target$Event, statloc = FALSE, CL.BT = TRUE, d0lab = '', d1lab = '', legendloc = FALSE)
dev.off()

pdf(file = 'calibration_curve_null_zoomed.pdf')
set.seed(42)
nullModelCalibrationCurveZoomed <- val.prob.ci.2(nullResponse, w1Target$Event, statloc = FALSE, CL.BT = TRUE, d0lab = '', d1lab = '', legendloc = FALSE, xlim = c(-0.02, 0.2), ylim = c(-0.15, 0.4), ylab = '')
dev.off()
