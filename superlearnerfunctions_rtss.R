sitetzScreen <- function(X,...){
  ind <- rep(FALSE, ncol(X))
  ind[which(colnames(X) %in% c("site1","site2","site3","site4","site5","Z","t"))] <- TRUE
  ind
}

sitezScreen <- function(X,...){
  ind <- rep(FALSE, ncol(X))
  ind[which(colnames(X) %in% c("site1","site2","site3","site4","site5","Z"))] <- TRUE
  ind
}

tzScreen <- function(X,...){
  ind <- rep(FALSE, ncol(X))
  ind[which(colnames(X) %in% c("Z","t"))] <- TRUE
  ind
}

screen.corRank5 <- function(...,rank=5){
  screen.corRank(...,rank=rank)
}

screen.corRank10 <- function(...,rank=10){
  screen.corRank(...,rank=rank)
}

trtOnlyScreen <- function(X,...){
  ind <- rep(FALSE, ncol(X))
  ind[which(colnames(X)=="Z")] <- TRUE
  ind
}

noIntScreen <- function(X,...){
  ind <- rep(TRUE, ncol(X))
  ind[grep("Int",colnames(X))] <- FALSE
  ind
}

predict.SL.randomForest.itMeans <- 
function (object, newdata, ...) 
{
  SuperLearner:::.SL.require("randomForest")
  if (!object$binOutcome) {
    pred <- predict(object$object, newdata = newdata, type = "response")
  }else{
    pred <- predict(object$object, newdata = newdata, type = "vote")[, 
                                                                     2]
  }
  pred
}

SL.gbm.itMeans <- function (Y, X, newX, family, obsWeights, gbm.trees = 1000, 
                            interaction.depth = 2, n.cores=1, ...) 
{
  SuperLearner:::.SL.require("gbm")
  gbm.model <- as.formula(paste("Y~", paste(colnames(X), collapse = "+")))
  binOutcome  <- all(Y %in% c(0,1))
  if (!binOutcome) {
    fit.gbm <- gbm::gbm(formula = gbm.model, data = X, distribution = "gaussian", 
                        n.trees = gbm.trees, interaction.depth = interaction.depth, 
                        cv.folds = 0, keep.data = TRUE, weights = obsWeights, 
                        verbose = FALSE, n.cores=n.cores)
  }else{
    fit.gbm <- gbm::gbm(formula = gbm.model, data = X, distribution = "bernoulli", 
                        n.trees = gbm.trees, interaction.depth = interaction.depth, 
                        cv.folds = 0, keep.data = TRUE, verbose = FALSE, 
                        weights = obsWeights, n.cores=n.cores)
  }
  best.iter <- gbm::gbm.perf(fit.gbm, method = "OOB", plot.it = FALSE)
  pred <- predict(fit.gbm, newdata = newX, best.iter, type = "response")
  fit <- list(object = fit.gbm, n.trees = best.iter, binOutcome=binOutcome)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.gbm")
  return(out)
}


SL.randomForest.itMeans <- function (Y, X, newX, family, mtry = ifelse(family$family == 
                                                                         "gaussian", floor(sqrt(ncol(X))), max(floor(ncol(X)/3), 1)), 
                                     ntree = 1000, nodesize = ifelse(family$family == "gaussian", 
                                                                     5, 1), ...) 
{
  SuperLearner:::.SL.require("randomForest")
  binOutcome  <- all(Y %in% c(0,1))
  if (!binOutcome) {
    fit.rf <- randomForest::randomForest(Y ~ ., data = X, 
                                         ntree = ntree, xtest = newX, 
                                         keep.forest = TRUE, 
                                         mtry = mtry, nodesize = nodesize)
    pred <- fit.rf$test$predicted
    fit <- list(object = fit.rf, binOutcome = binOutcome)
  }else{
    fit.rf <- randomForest::randomForest(y = as.factor(Y), 
                                         x = X, ntree = ntree, 
                                         xtest = newX, keep.forest = TRUE, 
                                         mtry = mtry, nodesize = nodesize)
    pred <- fit.rf$test$votes[, 2]
    fit <- list(object = fit.rf, binOutcome = binOutcome)
  }
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.randomForest.itMeans")
  return(out)
}

