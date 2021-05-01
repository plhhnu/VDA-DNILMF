## The following code is taken from Hao et al. [1].
## [1] Predicting drug-target interactions by dual-network integrated logistica matrix factorization.


## setwd("YourDir\\VDA-KLMF")
library(Rcpp)
# current data set name
db <- 'vd'

switch (db,
  vd = {
    cat("vd data\n")
    flush.console()
    sd <- read.table("../datasets/vd_simmat_dc.txt")
    sd <- as.matrix(sd)
    st <- read.table("../datasets/vd_simmat_dg.txt")
    st <- as.matrix(st)
    Y <- read.table("../datasets/vd_admat_dgc.txt")
    Y <- as.matrix(Y)
    Y <- Y
  },
  vd2 = {
    cat("vd2 data\n")
    flush.console()
    sd <- read.table("../datasets/vd2_simmat_dc.txt")
    sd <- as.matrix(sd)
    st <- read.table("../datasets/vd2_simmat_dg.txt")
    st <- as.matrix(st)
    Y <- read.table("../datasets/vd2_admat_dgc.txt")
    Y <- as.matrix(Y)
    Y <- Y
  },
  vd3 = {
    cat("vd3 data\n")
    flush.console()
    sd <- read.table("../datasets/vd3_simmat_dc.txt")
    sd <- as.matrix(sd)
    st <- read.table("../datasets/vd3_simmat_dg.txt")
    st <- as.matrix(st)
    Y <- read.table("../datasets/vd3_admat_dgc.txt")
    Y <- as.matrix(Y)
    Y <- Y
  },
  stop("db should be one of the follows: {vd, vd2 or vd3}\n")
)

## load required packages
pkgs <- c("matrixcalc", "data.table", "Rcpp", "ROCR", "Bolstad2", "MESS")
rPkgs <- lapply(pkgs, require, character.only = TRUE)

## source required R files
rSourceNames <- c("doCrossValidation.R",
                  "doCrossValidationByRow.R",
                  "doCrossValidationByCol.R",
                  "constrNeig.R", 
                  "inferZeros.R",
                  "calcLogLik.R",
                  "calcDeriv.R",
                  "updateUV.R",
                  "calcPredScore.R",
                  "calAUPR.R",
                  "calcPredScore.R")

rSN <- lapply(rSourceNames, source, verbose = FALSE)

## sourceCPP required C++ files
cppSourceNames <- c("fastKF.cpp", "fastKgipMat.cpp", "log1pexp.cpp", "sigmoid.cpp")
cppSN <- lapply(cppSourceNames, sourceCpp, verbose = FALSE)



## convert to kernel
isKernel <- TRUE
if (isKernel) {
  if (!isSymmetric(sd)) {
    sd <- (sd + t(sd)) / 2
  }
  epsilon <- 0.1
  while (!is.positive.semi.definite(sd)) {
    sd <- sd + epsilon * diag(nrow(sd))
  }
  if (!isSymmetric(st)) {
    st <- (st + t(st)) / 2
  }
  epsilon <- 0.1
  while (!is.positive.semi.definite(st)) {
    st <- st + epsilon * diag(nrow(st))
  }
}

## do cross-validation
kfold <- 5
numSplit <- 1

## split training and test sets
savedFolds <- doCrossValidationByCol(Y, kfold = kfold, numSplit = numSplit)

## hyper-parameters

isDefaultPara <- TRUE
if (isDefaultPara) {
  numLat <- 50
  cc <- 5
  thisAlpha <- 0.5
  lamU <- 5
  lamV <- 1
  K1 <- 5
} else {
  numLat <- 90
  cc <- 6
  thisAlpha <- 0.4
  lamU <- 2
  lamV <- 2
  K1 <- 2
}


## values according to hyper-parameters
thisBeta <- (1 - thisAlpha)/2
thisGamma <- 1 - thisAlpha - thisBeta

## for saving results
AUPRVec <- vector(length = kfold)
AUCVec <- vector(length = kfold)
accVec <- vector(length = kfold)
recallVec <- vector(length = kfold)
preVec <- vector(length = kfold)
f1Vec <- vector(length = kfold)
specVec <- vector(length = kfold)
finalResult <- matrix(NA, nrow = numSplit, ncol = 7)
colnames(finalResult) <- c("AUPR", "AUC","acc", "recall", "pre", "specificity", "f1")

# main loop
for (i in 1:numSplit) {
  for (j in 1:kfold) {
    cat("numSplit:", i, "/", numSplit, ";", "kfold:", j, "/", kfold, "\n")
    flush.console()
    Y <- savedFolds[[i]][[j]][[7]]
    Yr <- inferZeros(Y, sd, K = K1)
    Yc <- inferZeros(t(Y), st, K = K1)
    KgipD <- fastKgipMat(Yr, 1)
    KgipT <- fastKgipMat(Yc, 1)
    # nNeig = 3, nIter = 2
    sd <- fastKF(KgipD, sd, 3, 2)
    st <- fastKF(KgipT, st, 3, 2)
    lap <- constrNeig(sd, st, K = K1)
    lapD <- lap$lapD
    lapT <- lap$lapT
    simD <- lap$simD
    simT <- lap$simT
    ## use AdaGrid to update U and V
    UV <- updateUV(
      cc = cc,
      inMat = Y,
      thisAlpha = thisAlpha,
      thisBeta = thisBeta,
      Sd = simD,
      thisGamma = thisGamma,
      St = simT,
      lamU = lamU,
      lamV = lamV,
      numLat = numLat,
      initMethod = "useNorm",
      thisSeed = 123,
      maxIter = 100)
    
    U <- UV$U
    V <- UV$V
    
    knownDrugIndex <- savedFolds[[i]][[j]][[5]]
    knownTargetIndex <- savedFolds[[i]][[j]][[6]] 
    testIndexRow = savedFolds[[i]][[j]][[3]]      
    testIndexCol = savedFolds[[i]][[j]][[4]]      
    testLabel = savedFolds[[i]][[j]][[1]]         
    # result
    results <- calcPredScore(
      U = U,
      V = V,
      simDrug = simD,
      simTarget = simT,
      knownDrugIndex = knownDrugIndex,
      knownTargetIndex = knownTargetIndex,
      testIndexRow = testIndexRow,
      testIndexCol = testIndexCol,
      K = K1,
      testLabel = testLabel,
      thisAlpha = thisAlpha, 
      thisBeta = thisBeta,   
      thisGamma = thisGamma  
    )
    AUPRVec[j] <- results[1,"aupr"]
    AUCVec[j] <- results[1,"auc"]
    accVec[j] <- results[1,"acc"]
    recallVec[j] <- results[1,"recall"] 
    preVec[j] <- results[1,"pre"]
    specVec[j] <- results[1,"specificity"]
    f1Vec[j] <- results[1,"f1"]
  }
  AUPR <- mean(AUPRVec)
  AUC <- mean(AUCVec)
  acc <- mean(accVec)
  recall <- mean(recallVec)
  pre <- mean(preVec)
  spec <- mean(specVec)
  f1 <- mean(f1Vec)
  
  finalResult[i, "AUPR"] <- AUPR
  finalResult[i, "AUC"] <- AUC
  finalResult[i, "acc"] <- acc
  finalResult[i, "recall"] <- recall
  finalResult[i, "f1"] <- f1
  finalResult[i, "specificity"] <- spec
  finalResult[i, "pre"] <- pre
}


## print the result
cat(
  "\n======================\n\n",
  "db is: ", db, "\n",
  ## hyper-parameters
  "numLat = ", numLat, "\n",
  "cc = ", cc, "\n",
  "thisAlpha = ", thisAlpha, "\n",
  "lamU = ", lamU, "\n",
  "lamV = ", lamV, "\n",
  "K1 = ", K1, "\n",
  "\n=====================\n")
cat(numSplit, "trails 5-fold CV", "\n")
print(summary(finalResult))
cat("\n\n mean values:\n")
print(apply(finalResult, 2, mean))
cat("\n\n sd values:\n")
print(apply(finalResult, 2, sd))

# save to file
curDate <- format(Sys.time(), format = "%Y-%m-%d")
curTime <- format(Sys.time(), format =  "%H.%M.%S")
savedFileName <- paste0(db, "_", curDate, "_", curTime, ".RData")
cat("\n\n")
print(savedFileName)
#save.image(file = savedFileName)
