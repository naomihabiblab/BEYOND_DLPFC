

# -------------------------------------------------- #
# Receive file name from command-line input          #
# -------------------------------------------------- #
suppressPackageStartupMessages(library(optparse))
args <- parse_args(OptionParser(usage="%prog [options]", 
                                option_list=list(
                                  make_option(c("-d","--data.path"), action = "store",  
                                              help="Path to rds file containing train and test sets"),
                                  make_option(c("-n","--name"), action = "store", help="Name of model fitted"),
                                  make_option(c("-a","--alpha"), action = "store",
                                              help="Alpha parameter of fitted model. Passed to glmnet"),
                                  make_option(c("-c", "--ncores"), action="store",
                                              help="Number of cores to use while fitting model"),
                                  make_option(c("-o", "--output.path"), action="store",
                                              help="Output path for files and models"))))
unloadNamespace("optparse")




# -------------------------------------------------- #
# Fit classifier to given data                       #
# -------------------------------------------------- #
library(uuid)
library(doMC)
library(glmnet)

.trainModel <- function(train, alpha, name, lambdas, ncores=2, id=UUIDgenerate(), output.fig.path=NULL) {
  message(paste0(Sys.time(), "\t", id, "\tFitting model with alpha=",alpha, ". Using ", ncores, " cores"))
  
  st <- Sys.time()
  registerDoMC(cores = ncores)
  fit <- cv.glmnet(train$x, train$y, weights = train$w, family = "multinomial", 
                   type.measure="deviance", alpha=alpha, lambda = lambdas, parallel = T)
  message(paste0(Sys.time(), "\t", id, "\tTraining time was: ", Sys.time() - st))
  message(paste0(Sys.time(), "\t", id, "\tSelected labmda ", fit$lambda.min))
  
  if (!is.null(output.fig.path)) {
    plot(fit)
    title(paste0("CV Results - alpha=", alpha, ", lambda=",fit$lambda.min))
    dev.copy(png, paste0(output.fig.path, "model.", alpha, ".cv.res.png"))
    dev.off()
  }
  return(fit)
}

.evaluateModel <- function(test, model, model.name, id=UUIDgenerate()) {
  message(paste0(Sys.time(), "\t", id, "\tEvaluating model ", model.name, " on test test"))
  
  pred.prob <- predict(model, test$x, s="lambda.min", type="response")
  pred.type <- colnames(pred.prob)[apply(pred.prob, 1, which.max)]
  
  lam_id <- match(model$lambda.min, model$lambda)
  coefs <- do.call(cbind, lapply(model$glmnet.fit$beta, function(b) b[, lam_id]))
  rownames(coefs)[1] <- "Intercept"
  
  return(list(
    name=model.name,
    model=model,
    lambda=model$lambda.min,
    n.coef=model$nzero[lam_id],
    coef=coefs,
    class.prob <- pred.prob,
    confusion.matrix=table(test$y, pred.type, dnn = c("Predicted","True Classes")),
    accuracy=round(100 * mean(pred.type == test$y), 4),
    normalized.accuracy=round(100 * sum((pred.type == test$y)/table(test$y)[test$y])/length(unique(test$y)), 4)))
}


id <- UUIDgenerate()
lambdas <- 10^seq(-8, 3, length.out = 200)
datasets <- readRDS(args$data.path)

m <- .trainModel(datasets$train, args$alpha, args$name, lambdas, ncores = args$ncores, id = id, output.fig.path = args$output.path)
m <- .evaluateModel(datasets$test, m, args$name, id=id)

saveRDS(m, file.path(args$output.path, paste0(args$name, ".model.rds")))



# Create Classifier report 
coefs <- m$coef
coefs <- coefs[Matrix::rowSums(coefs) != 0,]

summary <- data.frame(m$lambda, m$n.coef, m$accuracy, m$normalized.accuracy)
colnames(summary) <- c("Selected Lambda", "# Non-Zero Features","Accuracy", "Normalized Accuracy")

library(cowplot)
library(gridExtra)
pdf(file.path(args$output.path, paste0(m$name, ".report.pdf")), width = 10, height = 10)
plot_grid(
  ggdraw() + draw_label(paste(m$name, "Report"), fontface='bold'),
  tableGrob(summary, rows = c()),
  plot_grid(
    pheatmap::pheatmap(m$confusion.matrix/Matrix::rowSums(m$confusion.matrix), 
                       display_numbers = m$confusion.matrix, 
                       cluster_cols = F, cluster_rows = F, silent = T,
                       color = colorRampPalette(c("navy","white","red3"))(40), 
                       main = "Prediction Confusion Map\n(predictions as rows)")[[4]],
    pheatmap::pheatmap(coefs, 
                       color = colorRampPalette(c("navy","white","red3"))(40), 
                       breaks = unique(c(seq(min(coefs), 0,length.out = 20), seq(0, max(coefs),length.out = 20))),
                       cluster_cols = F, show_rownames = F, silent = T,
                       main = "Fitted Coefficients Heatmap\n")[[4]]),
  
  tableGrob(do.call(rbind, lapply(m$model$glmnet.fit$beta,function(d) paste0(names(rev(tail(sort(d[,m$model$index[1]]),n = 10))), collapse = ", "))),
            cols = c("Features With Highest Coefficients"),
            theme = gridExtra::ttheme_default(base_size=8)),
  
  tableGrob(do.call(rbind, lapply(m$model$glmnet.fit$beta, function(d) paste0(names(sort(d[,m$model$index[1]])[1:10]), collapse = ", "))),
            cols = c("Features With Lowest Coefficients"),
            theme = gridExtra::ttheme_default(base_size=8)),
  ncol = 1, rel_heights = c(.1,.2,1.5,1,1))
while (!is.null(dev.list()))  dev.off()


