

PredictCellType <- function(object, classifier.path="1. Library preprocessing/data/celltype.classifier.rds", column.prefix="cell.type") {
  suppressPackageStartupMessages(require(glmnet))
  classifier <- readRDS(classifier.path)
  name <- classifier$name
  message(paste0("Classifying cell type using fitted ", name, " classifier"))
  classifier <- classifier$model
  
  # Arrange design matrix
  feature.space <- intersect(rownames(classifier$glmnet.fit$beta[[1]]), rownames(object))
  missing.features <- setdiff(rownames(classifier$glmnet.fit$beta[[1]]), rownames(object))
  
  message(paste("Not all features were found in given dataset. Missing features are treated as if contain a zero value.",
                " Missing features are:", paste0(missing.features, collapse = ", ")))
  
  
  temp <- NormalizeData(object, verbose = F) %>% ScaleData(features=feature.space, verbose = F)
  X <- t(GetAssayData(temp, "scale.data"))
  X <- cbind(X, matrix(0, nrow=nrow(X), ncol = length(missing.features), dimnames = list(rownames(X), missing.features)))
  X <- X[,rownames(classifier$glmnet.fit$beta[[1]])]
  rm(temp)
  
  
  # Perform prediction
  pred <- data.frame(predict(classifier, X, s="lambda.min", type="response"))
  colnames(pred) <- gsub(".1", "", colnames(pred))
  rm(X, feature.space, missing.features)
  
  max.val <- apply(pred, 1, max)
  nxt.max.val <- apply(pred, 1, function(v) v[order(v, decreasing = T)[2]])
  
  argmax <- colnames(pred)[apply(pred, 1, which.max)]
  next.argmax <- apply(pred, 1, function(v) colnames(pred)[order(v, decreasing = T)[2]])
  entropy <- apply(pred, 1, function(i) -sum(i*log2(i)))
  
  df <- data.frame(cell.type=argmax, cell.type.prob=max.val, 
                   cell.type.next=next.argmax, cell.type.next.prob=nxt.max.val, 
                   cell.type.entropy=entropy)
  colnames(df) <- gsub("cell.type", column.prefix, colnames(df))
  rm(max.val, nxt.max.val, argmax, next.argmax, entropy)
  
  # Append results to object
  object <- AddMetaData(object, df)
  Tool(object) <- list(classifier.name=name, classifier=classifier.path, prediction=pred)
  return(object)
}
