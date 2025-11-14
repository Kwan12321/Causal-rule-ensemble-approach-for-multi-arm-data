

variable.importance <- function(object, dat){
  
  ### Base function importance for Group lasso ###
  beta.first <- object$beta.first
  rules.first <- object$rules.first
  base.importance.first <- matrix(beta.first, nrow = length(rules.first), ncol = object$group.num)
  rownames(base.importance.first) <- rules.first
  colnames(base.importance.first) <- paste0("tau", 0:(object$group.num - 1))
  #base.importance.first <- base.importance.first[,-1] - base.importance.first[,1]
  # Calculate differences between columns
  differences <- list()
  for (i in 1:(ncol(base.importance.first) - 1)) {
    for (j in (i + 1):ncol(base.importance.first)) {
      differences[[paste(i, "-", j, sep = "")]] <- base.importance.first[, i] - base.importance.first[, j]
    }
  }
  base.importance.first <- do.call(cbind, differences)
  
  var.base.importance.first <- base.importance.first[(nrow(base.importance.first) - ncol(dat[,-c(1:2)]) + 1):(nrow(base.importance.first)),]
  rules.base.importance.first <- base.importance.first[1:(nrow(base.importance.first) - ncol(dat[,-c(1:2)])),]
  var.base.importance.first <- abs(var.base.importance.first)*apply(winsorize(dat), 2, sd)
  rules.base.importance.first <- abs(rules.base.importance.first)*apply(tran_rules(dat, rownames(rules.base.importance.first)), 2, sd)
  base.importance.first <- rbind(rules.base.importance.first, var.base.importance.first)
  
  ### Base function importance for Adaptive Group lasso ###
  beta <- object$beta
  rules <- object$rules
  base.importance <- matrix(beta, nrow = length(rules), ncol = object$group.num)
  rownames(base.importance) <- rules
  colnames(base.importance) <- paste0("tau", 0:(object$group.num - 1))
  #base.importance <- base.importance[,-1] - base.importance[,1]
  # Calculate differences between columns
  differences <- list()
  for (i in 1:(ncol(base.importance) - 1)) {
    for (j in (i + 1):ncol(base.importance)) {
      differences[[paste(i, "-", j, sep = "")]] <- base.importance[, i] - base.importance[, j]
    }
  }
  base.importance <- do.call(cbind, differences)
  
  var.base.importance <- base.importance[(nrow(base.importance) - ncol(dat[,-c(1:2)]) + 1):(nrow(base.importance)),]
  rules.base.importance <- base.importance[1:(nrow(base.importance) - ncol(dat[,-c(1:2)])),]
  var.base.importance <- abs(var.base.importance)*apply(winsorize(dat), 2, sd)
  rules.base.importance <- abs(rules.base.importance)*apply(tran_rules(dat, rownames(rules.base.importance)), 2, sd)
  base.importance <- rbind(rules.base.importance, var.base.importance)

  ### Variable importance ###
  base.counts.first <- sapply(rownames(base.importance.first), function(x){
    num <- length(unlist(strsplit(x," & ")))
    num
  })
  base.counts <- sapply(rownames(base.importance), function(x){
    num <- length(unlist(strsplit(x," & ")))
    num
  })
  
  var.imp.first <- var.imp <- NULL
  for(var in colnames(dat[,-c(1:2)])){
    var.base.first <- apply((base.importance.first/base.counts.first)[grep(var, rules.first),,drop = FALSE], 2, sum)
    var.base <- apply((base.importance/base.counts)[grep(var, rules),,drop = FALSE], 2, sum) 
    var.imp.first <- rbind(var.imp.first, var.base.first)
    var.imp <- rbind(var.imp, var.base)
  }
  rownames(var.imp.first) <- rownames(var.imp) <- colnames(dat[,-c(1:2)])
  
  return(list(base.imp.first = base.importance.first,
              base.imp = base.importance,
              var.imp.first = var.imp.first,
              var.imp = var.imp
              ))
}
