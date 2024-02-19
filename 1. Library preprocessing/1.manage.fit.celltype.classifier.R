

####################################################################################################################
##                                          #  Create Cell Type Classifier  #                                     ##
####################################################################################################################



## ----------------------------------------- ##
#          Prepare train and test sets       ##
## ----------------------------------------- ##
library(Seurat)
obj <- readRDS("1. Library preprocessing/data/Cain.2023.DLPFC.rds")
obj@assays$RNA@scale.data <- matrix(0)
obj <- subset(obj, features = rownames(obj)[!grepl("^(AL|AC|LINC)\\d+", rownames(obj))])

obj <- NormalizeData %>%
  FindVariableFeatures(nfeatures=700) %>%
  ScaleData()

x <- t(obj@assays$RNA@scale.data)
y <- as.character(obj$cell.type)

x <- x[!is.na(y),]
y <- y[!is.na(y)]


idx <- sample(seq_len(nrow(x)), size=floor(.75 * nrow(x)))
train <- list(x=x[idx,], y=y[idx], w=c(1/(table(y[idx])[y[idx]] * length(unique(y[idx])))))
test <- list(x=x[-idx,], y=y[-idx], w=c(1/(table(y[-idx])[y[-idx]] * length(unique(y[-idx])))))

saveRDS(list(train=train, test=test), 
        file = "1. Library preprocessing/data/celltype.classifier.datasets.rds")
rm(object, x, y, idx, train, test)
gc()



## ----------------------------------------- ##
#          Train & evaluate classifiers      ##
## ----------------------------------------- ##

models <- list("Ridge"=0, "ElasticNet.0.25"=.25, "ElasticNet.0.5"=.5, "ElasticNet.0.75"=.75, "Lasso"=1)
ram    <- c(45, 45, 45, 45, 45)
cores  <- c(10, 10, 10, 10, 10)


for(i in seq_along(models)) {
  exe.name <- file.path("1. Library preprocessing/scripts", paste0("classifier.", names(models)[[i]], ".sh"))
  
  if(file.exists(exe.name)) file.remove(exe.name)
  lapply(c("#!/bin/bash",
           "#SBATCH --partition=elsc.q",
           "#SBATCH --nodes=1",
           "#SBATCH --ntasks=1",
           paste0("#SBATCH --cpus-per-task=", cores[[i]]),
           paste0("#SBATCH --mem=",ram[[i]], "GB"),
           "#SBATCH --output=%x.o%j",
           paste0('Rscript "1. Library preprocessing/1.fit.celltype.classifier.R"',
                  ' -d "1. Library preprocessing/data/celltype.classifier.datasets.rds"', 
                  ' -n "', names(models)[[i]],'"',   
                  ' -a "', models[[i]], '"',
                  ' -c "', cores[[i]], '"', 
                  ' -o "1. Library preprocessing/data/"')),
         write, exe.name, append=T, ncolumns=1000)
}



## ----------------------------------------- ##
#            Select Best Classifier          ##
## ----------------------------------------- ##
models <- lapply(list.files("1. Library preprocessing/data", 
                            pattern = ".*model\\.rds", full.names = T), readRDS)

stats <- do.call(rbind, lapply(models, function(m) data.frame(m$lambda, m$n.coef, m$accuracy, m$normalized.accuracy)))
colnames(stats) <- c("Selected Lambda", "# Non-Zero Features","Accuracy", "Normalized Accuracy")
rownames(stats) <- lapply(models, function(m) m$name)
message("Models Comparison:")
print(stats)
 
# id <- ID OF SELECTED MODEL
# saveRDS(models[[id]], file="1. Library preprocessing/data/celltype.classifier.rds")
