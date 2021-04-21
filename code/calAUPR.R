# function to calculate the AUCR and AUC
calAUPR <- function(obsLabel, predProb) {
  unique_labels <- unique(obsLabel)
	# cat("unique labels:", unique_labels, "\n")
	# flush.console()
  library(caret)
  if (length(unique_labels) != 2) stop("The first argument 'obsLabel' should be two classes!\n")
	# calculate AUC using pROC
  library(pROC)
  roc1 <- pROC::roc(obsLabel, predProb)
  auc1 = as.numeric(auc(roc1))
  thresholds = coords(roc1,'all')['threshold']
  len = length(thresholds[[1]])
  f1 = 0
  acc = 0
  recall = 0
  specificity = 0
  pre = 0
  for (i in 1:len){
    cm  = table(as.integer(predProb >= thresholds[[1]][i]), obsLabel)
    if (nrow(cm) > 1){
      TP = cm['1','1']
      TN = cm['0','0']
      FP = cm['1','0']
      FN = cm['0','1']

      acc = acc + (TP+TN)/(TP+TN+FN+FP) 
      recall = recall + TP/(TP+FN)
      specificity = specificity + TN/(FP+TN)
      pre = pre + TP/(TP+FP)
      
      if (TP == 0){
        f1 = f1
      }else{
        f1 = f1 + 2*((TP/(TP+FP))*(TP/(TP+FN)))/((TP/(TP+FP))+(TP/(TP+FN)))
      }
      
    }else if(0 %in% rownames(cm)){
      TP = 0
      TN = cm['0','0']
      FP = 0
      FN = cm['0','1']
      acc = acc + (TP+TN)/(TP+TN+FN+FP) 
      recall = recall
      specificity = specificity + TN/(FP+TN)
      f1 = f1
      pre = pre
    }else{
      TP = cm['1','1']
      TN = 0
      FP = cm['1','0']
      FN = 0
      acc = acc + (TP+TN)/(TP+TN+FN+FP) 
      specificity = specificity + TN/(FP+TN)
      pre = pre+ TP/(TP+FP)
      if (TP == 0){
        f1 = f1
        recall = recall
        pre = pre
      }else{
        f1 = f1 + 2*(TP/(TP+FP))/(TP/(TP+FP)+1)
        recall = recall + TP/(TP+FN)
        pre = pre+ TP/(TP+FP)
      }
    }
  }
  
  acc = acc / len
  f1 = f1 / len
  recall = recall / len
  specificity = specificity / len
  pre = pre / len

  # Calculate AUPR using auprc
  library(auprc)
  aupr <- auprc(predProb, obsLabel, 1)

	# Save the result
	statRes <- matrix(0, nrow = 1, ncol = 7)
	colnames(statRes) <- c("auc", "aupr", "acc", "recall", "pre", "specificity", "f1")
  statRes[, "auc"] <- auc1
	statRes[, "aupr"] <- aupr
	statRes[, "acc"] <- acc
	statRes[, "recall"] <- recall
	statRes[, "pre"] <- pre
	statRes[, "specificity"] <- specificity
	statRes[, "f1"] <- f1
	return(statRes)
}
