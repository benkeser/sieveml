#---------------------------------------------------------------------
# Super Learner functions defined for use in the RV144 data analysis
#---------------------------------------------------------------------
# this function defines a version of the function SL.randomForest that is 
# compatible with the iterative means implementation of the TMLE
SL.randomForest1 <- function (Y, X, newX, family, mtry = ifelse(family$family == 
                                              "gaussian", floor(sqrt(ncol(X))), max(floor(ncol(X)/3), 1)), 
          ntree = 1000, nodesize = ifelse(family$family == "gaussian", 
                                          5, 1), ...) 
{
  if (length(unique(Y))>2) {
    fit.rf <- randomForest(Y ~ ., data = X, ntree = ntree, 
                           xtest = newX, keep.forest = TRUE, mtry = mtry, nodesize = nodesize)
    pred <- fit.rf$test$predicted
    fit <- list(object = fit.rf)
  }
  if (length(unique(Y))<=2) {
    fit.rf <- randomForest(y = as.factor(Y), x = X, ntree = ntree, 
                           xtest = newX, keep.forest = TRUE, mtry = mtry, nodesize = nodesize)
    pred <- fit.rf$test$votes[, 2]
    fit <- list(object = fit.rf)
  }
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.randomForest")
  return(out)
}

# this function defines a version of the SL.caret function from the SuperLearner package
# that is compatible with the iterative means implementation
SL.caret1 <- function (Y, X, newX, family, obsWeights, method = "rf", tuneLength = 3, 
          trControl = trainControl(method = "cv", number = 10, verboseIter = FALSE), 
          metric, 
          ...) 
{
  # change from family = specification of original function to 
  # check the unique values of Y to decide whether to use MSE or llik loss
  if (length(unique(Y))>2){
    if(is.matrix(Y)) Y <- as.numeric(Y)
    metric <- "RMSE"
    fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                              metric = metric, method = method, tuneLength = tuneLength, 
                              trControl = trControl)
    pred <- predict(fit.train, newdata = newX, type = "raw")
  }
  if (length(unique(Y))<=2) {
    cat("length <=2 \n")
    metric <- "Accuracy"
    Y.f <- as.factor(Y)
    levels(Y.f) <- c("A0", "A1")
    fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights, 
                              metric = metric, method = method, tuneLength = tuneLength, 
                              trControl = trControl)
    pred <- predict(fit.train, newdata = newX, type = "prob")[,2]
  }
  fit <- list(object = fit.train)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.caret")
  return(out)
}

# this function defines a regression tree tuned using caret for use
# with SuperLearner()
SL.rpart.caret1 <- function(...,method="rpart",tuneLength = 8){
  SL.caret1(...,method=method,tuneLength=tuneLength)
}
# this function defines a random forest tuned using caret for use
# with SuperLearner()
SL.rf.caret1 <- function(...,method="rf",tuneLength=8){
  SL.caret1(...,method=method,tuneLength=tuneLength)
}
