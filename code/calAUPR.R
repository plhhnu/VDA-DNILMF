# function to calculate the AUCR and AUC
calAUPR <- function(obsLabel, predProb) {
  unique_labels <- unique(obsLabel)
	# cat("unique labels:", unique_labels, "\n")
	# flush.console()
  library(caret)
  if (length(unique_labels) != 2) stop("The first argument 'obsLabel' should be two classes!\n")
	# calculate AUC using ROCR
	source('calculate_measures_by_threshold.R')
  pred <- ROCR::prediction(predProb, obsLabel)
  perf <- ROCR::performance(pred, "auc")
  auc <- as.numeric(perf@y.values)
	
  # Calculate AUPR using auprc
  library(auprc)
  aupr <- auprc(predProb, obsLabel, 1)

	# Save the result
	statRes <- matrix(0, nrow = 1, ncol = 7)
	
	thresholds = calculate_measures_by_threshold(predProb, obsLabel, 1)['thresh']
	f1 = 0
	acc = 0
	recall = 0
	specificity = 0
	pre = 0
	for (i in 1:1000){
	  pred = as.integer(predProb >= thresholds[i,1])
	  table_matrix = table(pred,obsLabel)
	
	  if (nrow(table_matrix) > 1){
	    TP = table_matrix['1','1']
	    TN = table_matrix['0','0']
	    FP = table_matrix['1','0']
	    FN = table_matrix['0','1']
	    acc = acc + (TP+TN)/(TP+TN+FN+FP) 
	    recall = recall + TP/(TP+FN)
	    specificity = specificity + TN/(FP+TN)
	    pre = pre+ TP/(TP+FP)
	    # f1 = f1 + 2*((TP/(TP+FP))*(TP/(TP+FN)))/((TP/(TP+FP))+(TP/(TP+FN)))
	    
	    if (TP == 0){
	      f1 = f1
	    }else{
	      f1 = f1 + 2*((TP/(TP+FP))*(TP/(TP+FN)))/((TP/(TP+FP))+(TP/(TP+FN)))
	    }
	    
	  }else if(0 %in% rownames(table_matrix)){
	    TP = 0
	    TN = table_matrix['0','0']
	    FP = 0
	    FN = table_matrix['0','1']
	    acc = acc + (TP+TN)/(TP+TN+FN+FP) 
	    recall = recall + TP/(TP+FN)
	    specificity = specificity + TN/(FP+TN)
	    f1 = f1
	    pre = pre
	  }else{
	    TP = table_matrix['1','1']
	    TN = 0
	    FP = table_matrix['1','0']
	    FN = 0
	    acc = acc + (TP+TN)/(TP+TN+FN+FP) 
	    specificity = specificity + TN/(FP+TN)
	    pre = pre+ TP/(TP+FP)
	    if (TP == 0){
	      f1 = f1
	      recall = recall
	    }else{
	      f1 = f1 + 2*(TP/(TP+FP))/(TP/(TP+FP)+1)
	      recall = recall + TP/(TP+FN)
	    }
	  }
	}

	f1 = f1/1000
	acc = acc/1000
	recall = recall/1000
	specificity = specificity/1000
	pre = pre/1000
	colnames(statRes) <- c("auc", "aupr","acc","recall","pre","specificity","f1")
  statRes[, "auc"] <- auc
	statRes[, "aupr"] <- aupr
	statRes[, "acc"] <- acc
	statRes[, "recall"] <- recall
	statRes[, "pre"] <- pre
	statRes[, "specificity"] <- specificity
	statRes[, "f1"] <- f1
	return(statRes)
}








