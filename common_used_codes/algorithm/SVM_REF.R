library(e1071)

source('/Users/stead/Documents/SourceTree/R/common_used_codes/algorithm/SVM_REF/msvmRFE.R')
load('/Users/stead/Documents/SourceTree/R/common_used_codes/algorithm/SVM_REF/demo/input.Rdata')
svmRFE(input, k=10, halve.above = 100)

nfold = 10
nrows = nrow(input)
folds = rep(1:nfold, len = nrows)[sample(nrows)]
results = lapply(folds, svmRFE.wrap, input, k=10, halve.above=100)

###get top features
top.features = WriteFeatures(results, input, save = F)
featsweep = lapply(1 : 20, FeatSweep.wrap, results, input)
head(top.features)

no.info = min(prop.table(table(input[, 1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

dev.new(width=4, height=4, bg='white')
PlotErrors(errors, no.info=no.info)
dev.off()