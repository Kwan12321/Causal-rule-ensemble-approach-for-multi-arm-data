
### Main function for propsoed method ###

mcre_agenet <- function(X, Y, W, W.hat = NULL,
                        ntrees = 333, depth = 2,
                        gen.tree = "gbm", learnrate = learnrate,
                        init = "glasso", fin = "glasso",
                        tune.init = "CV", tune.fin = "CV",
                        nfolds = 10L, seed = 1) {
  
  
  ### Input:
  ### -------------------------------------
  
  ############
  ### Data ###
  ############
  
  # X: Covariates
  # Y: Outcome
  # W: Treatment indicators
  
  ##################################
  ### Rule generation parameters ###
  ##################################
  
  # gen.tree: Select from "gbm", "ctree", "gbm.b", or "ctree.b"
  #   - "gbm":  PS trimming + GBM
  #   - "ctree": PS trimming + Conditional Inference Tree
  #   - "gbm.b": CBPS + GBM
  #   - "ctree.b": CBPS + Conditional Inference Tree
  
  # ntrees: Number of base learners for the boosting model
  
  # depth: Maximum depth of each base learner
  
  # learnrate: Learning rate for the boosting model
  
  ################################
  ### Rule ensemble parameters ###
  ################################
  
  # init: First-stage group regularization method (corresponds to group lasso)
  
  # fin: Second-stage group regularization method (corresponds to adaptive group lasso)
  
  # init/fin: Choose from:
  #   - "genet": Group elastic-net
  #   - "glasso": Group lasso
  
  # tune.init / tune.fin: Method to select the optimal lambda
  #   Choose from "CV", "EBIC", "BIC", "AIC", or "GCV"
  
  # nfolds: Number of folds used for cross-validation (if "CV" is selected)
  
  ### -------------------------------------------------
  
  start <- proc.time()
  
  ### Stage 1: Rule generation
  set.seed(seed)
  rules <- NULL
  
  if(gen.tree == "gbm") {
    rules <- rule_gen.gbm(X = X, Y = Y, W = W, W.hat = W.hat, ntrees = ntrees, depth = depth, learnrate = learnrate, seed = seed)
    print("gbm")
  } else if(gen.tree == "cf") {
    rules <- rule_gen.cf(X = X, Y = Y, W = W, ntrees = ntrees, seed = seed)
    print("cf")
  } else if(gen.tree == "ctree") {
    rules <- rule_gen.ctree(X = X, Y = Y, W = W, W.hat = W.hat, ntrees = ntrees, depth = depth, learnrate = learnrate, seed = seed)
    print("ctree")
  } else if(gen.tree == "gbm.b") {
    rules <- rule_gen.gbm.b(X = X, Y = Y, W = W, ntrees = ntrees, depth = depth, learnrate = learnrate, seed = seed)
    print("gbm.b")
  } else if(gen.tree == "ctree.b") {
    rules <- rule_gen.ctree.b(X = X, Y = Y, W = W, ntrees = ntrees, depth = depth, learnrate = learnrate, seed = seed)
    print("ctree.b")  
  }
  
  if (is.null(rules)) {
    stop("Error in rule generation: rules are NULL")
  }
  
  rules <- unique(rules)
  print(paste("Number of unique rules:", length(rules)))
  
  # Remove duplicate and complementary rules
  #rules <- remove_complement(X, rules)[[1]]
  if (length(rules) == 0) {
    stop("Error after removing complements: no rules remain.")
  }
  
  ### Stage 2: Rule ensemble 
  
  # Grouping the base function for group lasso
  grp.list <- group.train(rules = rules, data = data.frame(Y = Y, W = W, X))
  if (length(grp.list) == 0) {
    stop("Error: grp.list is empty.")
  }
  
  # Fit group lasso
  fit <- agenet(X = grp.list$base, y = Y, group = grp.list$index, init = init, fin = fin,
                tune.init = tune.init, tune.fin = tune.fin, nfolds = nfolds, seed = seed)
  
  end <- proc.time() - start
  
  intercept.fin <- fit$beta[1,]
  beta.fin <- fit$beta[-1,] / grp.list$scale
  intercept.first <- fit$beta.first[1,]
  beta.first <- fit$beta.first[-1,] / grp.list$scale
  time.first <- fit$time.ini
  
  return(list(intercept = intercept.fin, beta = beta.fin, rules = c(rules, colnames(X)),
              intercept.first = intercept.first, beta.first = beta.first,
              rules.first = c(rules, colnames(X)),
              group.num = as.numeric(table(grp.list$index)[1]),
              time.first = time.first[3], time = end[3],
              model.first = fit$model.first, model = fit$model))
}

### Prediction function for proposed method ###

predict.mcre <- function(object, newdata) {
  
  rules.fin <- object$rules
  intercept.fin <- object$intercept
  beta.fin <- object$beta
  rules.first <- object$rules.first
  intercept.first <- object$intercept.first
  beta.first <- object$beta.first
  group.num <- object$group.num
  
  grp.list.fin <- group.test(rules.fin, newdata, group.num)
  grp.list.first <- group.test(rules.first, newdata, group.num)
  
  if (length(grp.list.fin) == 0 || length(grp.list.first) == 0) {
    stop("Error: Group terms generation failed in predict function.")
  }
  
  tau.can.fin <- NULL
  tau.can.first <- NULL
  for (g in 1:group.num) {
    tau.can.fin <- cbind(tau.can.fin, grp.list.fin[[g]] %*% beta.fin + intercept.fin)
    tau.can.first <- cbind(tau.can.first, grp.list.first[[g]] %*% beta.first + intercept.first)
  }
  tau.fin_hat <- tau.can.fin[, -1] - tau.can.fin[, 1]
  tau.first_hat <- tau.can.first[, -1] - tau.can.first[, 1]
  res <- list(tau = tau.fin_hat, tau.first = tau.first_hat, y = tau.can.fin, y.first = tau.can.first)
  res
}

####################################################
##### 1: Helpful functions for Rule generation #####
####################################################

### 1: Causal Forest 

rule_gen.cf <- function(X, Y, W, ntrees = 500, seed = 1){
  
  set.seed(seed)
  
  fit <- grf::multi_arm_causal_forest(Y = Y, W = as.factor(W), X = X, num.trees = ntrees)
  rules <- unlist(sapply(1:fit$`_num_trees`,function(x){    
    tree <- grf::get_tree(fit, index = x)
    get.rules.cf(tree)}))
  
  #return(list(rules = rules, var.imp = t(var.imp)))
  rules
}

### 2: Conditional Inference Tree 

rule_gen.ctree <- function(X, Y, W, W.hat, ntrees = 500, depth = 2, sampfrac = NULL, learnrate = 0.01, seed = 1){
  
  set.seed(seed)
  
  ## trimming the propensity score () 
  r_upper <- apply(W.hat, 2, quantile, 0.9)
  r_lower <- apply(W.hat, 2, quantile, 0.1)
  
  for(c in 1:length(unique(W))){
    W.hat[,c][W.hat[,c] >= r_upper[c]] <- r_upper[c]
    W.hat[,c][W.hat[,c] <= r_lower[c]] <- r_lower[c]
  }
  
  # Transform outcome
  W.ind <- sapply(sort(unique(W)), function(z){as.numeric(W==z)})
  Y.tilde <- (Y*W.ind/W.hat)[,-1] - (Y*W.ind/W.hat)[,1]
  Y.names <- paste0("Y",1:ncol(Y.tilde))
  colnames(Y.tilde) <- Y.names
  X.names <- colnames(X)
  
  # Train outcome & Formula
  data_with_y_learn <- data.learn <- cbind(Y.tilde, X)
  formula <- as.formula(paste0(paste(Y.names, collapse = "+"),"~",paste(X.names, collapse = "+")))
  
  # Determine the sample fraction size for training the base learner (Default using the setting in Friedman and Popescu (2008))
  n <- nrow(X)
  if (is.null(sampfrac)) {
    size <- min(n / 2, 100 + 6 * sqrt(n))
  } else {
    size <- sampfrac * n
  }
  
  # Create the index for training data for each base learner
  subsample <- list()
  subsample <- mapply(function(i) {
    if (size == n) {
      subsample[[i]] <- sample(1:n, size = size, replace = FALSE)
    } 
    else if (size < n) {
      subsample[[i]] <- sample(1:n, size = round(size), replace = FALSE)
    }
  }, i = 1:ntrees, SIMPLIFY = FALSE)
  
  # Calculate the number of terminal nodes for each tree using an exponential distribution
  maxdepth <- ceiling(log(2 + floor(rexp(ntrees, rate = 1 / (2^depth - 2))), base = 2))
  
  # Initialize estimates
  y <- data.learn[ , Y.names]
  eta_0 <- apply(y, 2, weighted.mean, weights = rep(1, nrow(y)))
  eta <- t(replicate(n = nrow(y), expr = eta_0))
  
  # Calculate the pseudo-response variable using the negative gradient function
  data_with_y_learn[, Y.names] <- y - eta
  
  # Initialize an empty vector to store the rules
  rules <- c()
  
  for(i in 1:ntrees) {
    
    # Set up controls for ctree
    ctree.control <- partykit::ctree_control(maxdepth = maxdepth[i])
    
    # Fit the individual tree
    tree <- partykit::ctree(formula = formula, control = ctree.control, data = data_with_y_learn[subsample[[i]], ])
    
    # Decompose tree into rules
    paths <- list.rules(tree)
    
    # Append new rules to existing rules
    rules <- c(rules, paths)
    
    # Update eta using predictions from the new tree
    eta <- eta + learnrate * predict(tree, newdata = data_with_y_learn)
    
    # Update pseudo-response variable using the new eta
    data_with_y_learn[ ,Y.names] <- y - eta
    
  }
  
  rules
}  

### 3: Gradient Boosting Tree (GBT) 

rule_gen.gbm <- function(X, Y, W, W.hat, ntrees = 500, depth = 2, sampfrac = NULL, learnrate = 0.01, seed = 1){
  
  set.seed(seed)
  
  ## trimming the propensity score () 
  r_upper <- apply(W.hat, 2, quantile, 0.9)
  r_lower <- apply(W.hat, 2, quantile, 0.1)
  
  for(c in 1:length(unique(W))){
    W.hat[,c][W.hat[,c] >= r_upper[c]] <- r_upper[c]
    W.hat[,c][W.hat[,c] <= r_lower[c]] <- r_lower[c]
  }
  
  # Transform outcome
  W.ind <- sapply(sort(unique(W)), function(z){as.numeric(W==z)})
  Y.tilde <- (Y*W.ind/W.hat)[,-1] - (Y*W.ind/W.hat)[,1]
  Y.names <- paste0("Y",1:ncol(Y.tilde))
  colnames(Y.tilde) <- Y.names
  X.names <- colnames(X)
  
  # Train outcome & Formula
  data_with_y_learn <- data.learn <- cbind(Y.tilde, X)
  formula <- as.formula(paste0(paste(Y.names, collapse = "+"),"~",paste(X.names, collapse = "+")))
  
  # Determine the sample fraction size for training the base learner (Default using the setting in Friedman and Popescu (2008))
  n <- nrow(X)
  if (is.null(sampfrac)) {
    size <- min(n / 2, 100 + 6 * sqrt(n))
  } else {
    size <- sampfrac * n
  }
  
  # Create the index for training data for each base learner
  subsample <- list()
  subsample <- mapply(function(i) {
    if (size == n) {
      subsample[[i]] <- sample(1:n, size = size, replace = FALSE)
    } 
    else if (size < n) {
      subsample[[i]] <- sample(1:n, size = round(size), replace = FALSE)
    }
  }, i = 1:ntrees, SIMPLIFY = FALSE)
  
  # Calculate the number of terminal nodes for each tree using an exponential distribution
  maxdepth <- ceiling(log(2 + floor(rexp(ntrees, rate = 1 / (2^depth - 2))), base = 2))
  
  # Initialize estimates
  y <- data.learn[ , Y.names]
  eta_0 <- apply(y, 2, weighted.mean, weights = rep(1, nrow(y)))
  eta <- t(replicate(n = nrow(y), expr = eta_0))
  
  # Calculate the pseudo-response variable using the negative gradient function
  data_with_y_learn[, Y.names] <- y - eta
  
  # Initialize an empty vector to store the rules
  rules <- c()
  #var.imp <- numeric(ncol(X))
  #names(var.imp) <- X.names
  
  for(i in 1:ntrees) {
    
    # Set up controls for rpart tree
    tree.control <- rpart::rpart.control(maxdepth = maxdepth[i])
    
    # Fit the individual tree
    tree <- rpart::rpart(formula, control = tree.control, data = data_with_y_learn[subsample[[i]], ])
    
    # Decompose the tree into rules
    paths <- rpart::path.rpart(tree, nodes = rownames(tree$frame), print.it = FALSE, pretty = 0)
    paths <- unname(sapply(sapply(paths, `[`, index = -1), paste, collapse = " & ")[-1])
    
    # Append new rules to existing rules
    rules <- c(rules, paths)
    
    # Update eta using predictions from the new tree
    eta <- eta + learnrate * predict(tree, newdata = data_with_y_learn)
    
    # Update pseudo-response variable using the new eta
    data_with_y_learn[ ,Y.names] <- y - eta
    
  }
  
  rules
}

### 2-1: Conditional Inference Tree (Balancing)

rule_gen.ctree.b <- function(X, Y, W, ntrees = 500, depth = 2, sampfrac = NULL, learnrate = 0.01, seed = 1){
  
  set.seed(seed)
  
  W.hat <- WeightIt::weightit(W ~., data = data.frame(as.factor(W), X), method = "npcbps", estimand = "ATE")
  W.hat <- W.hat$weight
  
  # Transform outcome
  W.ind <- sapply(sort(unique(W)), function(z){as.numeric(W==z)})
  Y.tilde <- (Y*W.hat*W.ind)[,-1] - (Y*W.hat*W.ind)[,1]
  Y.names <- paste0("Y",1:ncol(Y.tilde))
  colnames(Y.tilde) <- Y.names
  X.names <- colnames(X)
  
  # Train outcome & Formula
  data_with_y_learn <- data.learn <- cbind(Y.tilde, X)
  formula <- as.formula(paste0(paste(Y.names, collapse = "+"),"~",paste(X.names, collapse = "+")))
  
  # Determine the sample fraction size for training the base learner (Default using the setting in Friedman and Popescu (2008))
  n <- nrow(X)
  if (is.null(sampfrac)) {
    size <- min(n / 2, 100 + 6 * sqrt(n))
  } else {
    size <- sampfrac * n
  }
  
  # Create the index for training data for each base learner
  subsample <- list()
  subsample <- mapply(function(i) {
    if (size == n) {
      subsample[[i]] <- sample(1:n, size = size, replace = FALSE)
    } 
    else if (size < n) {
      subsample[[i]] <- sample(1:n, size = round(size), replace = FALSE)
    }
  }, i = 1:ntrees, SIMPLIFY = FALSE)
  
  # Calculate the number of terminal nodes for each tree using an exponential distribution
  maxdepth <- ceiling(log(2 + floor(rexp(ntrees, rate = 1 / (2^depth - 2))), base = 2))
  
  # Initialize estimates
  y <- data.learn[ , Y.names]
  eta_0 <- apply(y, 2, weighted.mean, weights = rep(1, nrow(y)))
  eta <- t(replicate(n = nrow(y), expr = eta_0))
  
  # Calculate the pseudo-response variable using the negative gradient function
  data_with_y_learn[, Y.names] <- y - eta
  
  # Initialize an empty vector to store the rules
  rules <- c()
  
  for(i in 1:ntrees) {
    
    # Set up controls for ctree
    ctree.control <- partykit::ctree_control(maxdepth = maxdepth[i])
    
    # Fit the individual tree
    tree <- partykit::ctree(formula = formula, control = ctree.control, data = data_with_y_learn[subsample[[i]], ])
    
    # Decompose tree into rules
    paths <- list.rules(tree)
    
    # Append new rules to existing rules
    rules <- c(rules, paths)
    
    # Update eta using predictions from the new tree
    eta <- eta + learnrate * predict(tree, newdata = data_with_y_learn)
    
    # Update pseudo-response variable using the new eta
    data_with_y_learn[ ,Y.names] <- y - eta
    
  }
  
  rules
}  

### 3-1: Gradient Boosting Tree (GBT)(Balancing)

rule_gen.gbm.b <- function(X, Y, W, ntrees = 500, depth = 2, sampfrac = NULL, learnrate = 0.01, seed = 1){
  
  set.seed(seed)
  
  W.hat <- WeightIt::weightit(W ~., data = data.frame(as.factor(W), X), method = "npcbps", estimand = "ATE")
  W.hat <- W.hat$weight
  
  # Transform outcome
  W.ind <- sapply(sort(unique(W)), function(z){as.numeric(W==z)})
  Y.tilde <- (Y*W.hat*W.ind)[,-1] - (Y*W.hat*W.ind)[,1]
  Y.names <- paste0("Y",1:ncol(Y.tilde))
  colnames(Y.tilde) <- Y.names
  X.names <- colnames(X)
  
  # Train outcome & Formula
  data_with_y_learn <- data.learn <- cbind(Y.tilde, X)
  formula <- as.formula(paste0(paste(Y.names, collapse = "+"),"~",paste(X.names, collapse = "+")))
  
  # Determine the sample fraction size for training the base learner (Default using the setting in Friedman and Popescu (2008))
  n <- nrow(X)
  if (is.null(sampfrac)) {
    size <- min(n / 2, 100 + 6 * sqrt(n))
  } else {
    size <- sampfrac * n
  }
  
  # Create the index for training data for each base learner
  subsample <- list()
  subsample <- mapply(function(i) {
    if (size == n) {
      subsample[[i]] <- sample(1:n, size = size, replace = FALSE)
    } 
    else if (size < n) {
      subsample[[i]] <- sample(1:n, size = round(size), replace = FALSE)
    }
  }, i = 1:ntrees, SIMPLIFY = FALSE)
  
  # Calculate the number of terminal nodes for each tree using an exponential distribution
  maxdepth <- ceiling(log(2 + floor(rexp(ntrees, rate = 1 / (2^depth - 2))), base = 2))
  
  # Initialize estimates
  y <- data.learn[ , Y.names]
  eta_0 <- apply(y, 2, weighted.mean, weights = rep(1, nrow(y)))
  eta <- t(replicate(n = nrow(y), expr = eta_0))
  
  # Calculate the pseudo-response variable using the negative gradient function
  data_with_y_learn[, Y.names] <- y - eta
  
  # Initialize an empty vector to store the rules
  rules <- c()
  #var.imp <- numeric(ncol(X))
  #names(var.imp) <- X.names
  
  for(i in 1:ntrees) {
    
    # Set up controls for rpart tree
    tree.control <- rpart::rpart.control(maxdepth = maxdepth[i])
    
    # Fit the individual tree
    tree <- rpart::rpart(formula, control = tree.control, data = data_with_y_learn[subsample[[i]], ])
    
    # Decompose the tree into rules
    paths <- rpart::path.rpart(tree, nodes = rownames(tree$frame), print.it = FALSE, pretty = 0)
    paths <- unname(sapply(sapply(paths, `[`, index = -1), paste, collapse = " & ")[-1])
    
    # Append new rules to existing rules
    rules <- c(rules, paths)
    
    # Update eta using predictions from the new tree
    eta <- eta + learnrate * predict(tree, newdata = data_with_y_learn)
    
    # Update pseudo-response variable using the new eta
    data_with_y_learn[ ,Y.names] <- y - eta
    
  }
  
  rules
}

### 4: Extract Rules from Causal Forest

get.rules.cf <- function(tree, edge = NULL, index = 1){
  
  ### Current Node (If not leaf)###
  node <- tree$nodes[[index]]
  if(node$is_leaf){
    res <- NULL
    return(res)
  }
  variable_name <- tree$columns[node$split_variable]
  
  ### Rules of left edge ###
  if(!node$is_leaf){
    edge_info_left <- paste(variable_name, "<=", round(node$split_value,digits = 2))
    if(index > 1){
      edge_info_left <- paste(edge,"&",edge_info_left)
    }
  } else {
    edge_info_left < NULL
  } 
  
  ### Rules of right edge ###
  if(!node$is_leaf){
    edge_info_right <- paste(variable_name, ">", round(node$split_value,digits = 2))
    if(index > 1){
      edge_info_right <- paste(edge,"&",edge_info_right)
    }
  } else {
    edge_info_right < NULL
  } 
  
  ### Rules created up to the current node ###  
  this_lines <- c(edge, edge_info_left, edge_info_right)
  
  ### Get rules recursively ###
  right_child_lines <- get.rules.cf(tree, edge = edge_info_right, index = node$right_child)
  left_child_lines <- get.rules.cf(tree, edge = edge_info_left, index = node$left_child)
  
  ### Summary the result ###
  res <- c(this_lines,left_child_lines,right_child_lines)
  
  return(res)
  
}

### 5: Transform "ctree" into list of rules (Fokkema et al., 2020; R package "pre") 

list.rules <- function (x, i = NULL, removecomplements = TRUE, 
                        singleconditions = FALSE, ...) {
  
  ## singleconditions can take values FALSE (rules code only node membership);
  ##   TRUE (rules code both node membership and all individual splits; 
  ##   "only" (rules code only individual splits)
  
  if (is.null(i)) 
    i <- partykit::nodeids(x, terminal = TRUE)
  if (length(i) > 1) {
    # ret <- sapply(i, list.rules, x = x)
    # TODO: Benjamin Christoffersen changed this part. This can be done smarter
    # than finding all and then removing duplicates. I guess the computational
    # cost is low, though
    if (isTRUE(singleconditions)) {
      ret <- lapply(i, list.rules, x = x, simplify = FALSE)
      ret <- c(ret, lapply(i, list.rules, x = x, simplify = FALSE, 
                           singleconditions = TRUE))
    } else if (singleconditions == "only") {
      ret <- lapply(i, list.rules, x = x, simplify = FALSE, 
                    singleconditions = TRUE)
    } else {
      ret <- lapply(i, list.rules, x = x, simplify = FALSE)
    }
    
    # Find the first rules. We will only keep one of these
    
    ## TODO: If we apply non-negativity constraints,
    ## the rule that is kept should correlate positively
    ## with the outcome, if we apply negativity constraint,
    ## the rule that is kep should correlate negatively with 
    ## the response.
    ## I.e., if  'lower.limits = 0' or 'upper.limits = 0' was used in calling pre()
    ##
    ## Easier solution may be to just not remove first rule here
    ## E..g, employ rm.firstrule argument (which is true by default)
    if (removecomplements) {
      first_rules <- unique(sapply(ret, "[[", 1))
      first_rule_remove <- first_rules[2]
    }
    
    # Make list of final rules
    ret <- unlist(ret)
    ret <- ret[!duplicated(ret)]
    if (removecomplements) {
      ret <- ret[ret != first_rule_remove]
    }
    ## of rules with only single conditions, retain one of each pair
    if (isTRUE(singleconditions) && removecomplements) {
      mcrs <- ret[grep(" & ", ret)] ## multi-condition rules
      scrs <- ret[-grep(" & ", ret)] ## single-condition rules
      ## keep only odd numbered of those 
      ret <- c(mcrs, scrs[1:length(scrs) %% 2 == 1])
    } else if (singleconditions == "only" && removecomplements) {
      ## keep only odd numbered rules 
      ret <- c(ret[1:length(ret) %% 2 == 1])      
    }
    
    # TODO: this still leaves us with complements for non-terminal rules
    # names(ret) <- if (is.character(i)) 
    #   i else names(x)[i]
    return(ret) # Root node returns here
  }
  
  # Non-root nodes starts here
  if (is.character(i) && !is.null(names(x))) 
    i <- which(names(x) %in% i)
  #stopifnot(length(i) == 1 & is.numeric(i))
  #stopifnot(i <= length(x) & i >= 1)
  i <- as.integer(i)
  # dat <- partykit::data_party(x, i)
  # if (!is.null(x$fitted)) {
  #   findx <- which("(fitted)" == names(dat))[1]
  #   fit <- dat[, findx:ncol(dat), drop = FALSE]
  #   dat <- dat[, -(findx:ncol(dat)), drop = FALSE]
  #   if (ncol(dat) == 0) 
  #     dat <- x$data
  # }
  # else {
  #   fit <- NULL
  #   dat <- x$data
  # }
  dat <- x$data
  rule <- c()
  recFun <- function(node) {
    # if (partykit::id_node(node) == i) {
    #   return(NULL)
    # }
    if (node$id == i) {
      return(NULL)
    }
    # kid <- sapply(partykit::kids_node(node), partykit::id_node)
    kid <- sapply(node$kids, function(x) x$id)
    whichkid <- max(which(kid <= i))
    #split <- partykit::split_node(node)
    split <- node$split
    # ivar <- partykit::varid_split(split)
    ivar <- split$varid
    svar <- names(dat)[ivar]
    # index <- partykit::index_split(split)
    index <- split$index
    if (is.factor(dat[, svar])) {
      # if (is.null(index)) 
      #   index <- ((1:nlevels(dat[, svar])) > partykit::breaks_split(split)) + 1
      if (is.null(index))
        index <- ((1:nlevels(dat[, svar])) > split$breaks) + 1
      slevels <- levels(dat[, svar])[index == whichkid]
      # factor levels not occurring in the node will be coded as NA
      # and should be removed from rule description:
      slevels <- slevels[!is.na(slevels)]
      srule <- paste(svar, " %in% c(\"", paste(slevels, 
                                               collapse = "\", \"", sep = ""), "\")", sep = "")
    } else {
      if (is.null(index)) {
        index <- 1:length(kid)
      }
      # breaks <- cbind(c(-Inf, partykit::breaks_split(split)), c(partykit::breaks_split(split), 
      #                                                           Inf))
      breaks <- cbind(c(-Inf, split$breaks), c(split$breaks, Inf))
      sbreak <- breaks[index == whichkid, ]
      # right <- partykit::right_split(split)
      right <- split$right
      srule <- c()
      if (is.finite(sbreak[1])) {
        srule <- c(srule, paste(svar, ifelse(right, ">", 
                                             ">="), sbreak[1]))
      }
      if (is.finite(sbreak[2])) { 
        srule <- c(srule, paste(svar, ifelse(right, "<=", 
                                             "<"), sbreak[2]))
      }
      srule <- paste(srule, collapse = " & ")
    }
    rule <<- c(rule, srule)
    return(recFun(node[[whichkid]]))
  }
  # node <- recFun(partykit::node_party(x))
  node <- recFun(x$node)
  # paste(rule, collapse = " & ")
  
  if(is.null(rule))
    return(character())
  
  if (isTRUE(singleconditions) || singleconditions == "only") {
    ## keep conditions separate
    rule
  } else {
    ## combine conditions into single rule
    sapply(seq_along(rule), function(r) paste(rule[1:r], collapse = " & "))
  }
  
}

### 6: Transforms specified rules in the dataset into binary variables
tran_rules <- function(dat, rules) {
  tryCatch({
    expr <- parse(text = paste0("cbind(", paste0(rules, collapse = ", "), ")"))
    rulevars <- eval(expr, dat) + 0
    colnames(rulevars) <- rules
    return(rulevars)
  }, error = function(e) {
    stop("Error in tran_rules: ", e$message)
  })
}

### 7: Remove the complementary rules
remove_complement <- function(dat, rules) {
  rulevars <- tran_rules(dat, rules)
  vars <- apply(rulevars, 2, sd)
  
  vars_distinct <- lapply(unique(vars), function(x) {
    idx <- which(is_almost_eq(x, vars))
    list(var = x, n = length(idx), idx = idx)
  })
  
  complements <- logical(ncol(rulevars))
  
  for (va in vars_distinct) {
    if (va$n < 2L) next
    idx <- setdiff(va$idx, which(complements))
    if (length(idx) < 2) next
    
    n_idx <- length(idx)
    for (j in 1:(n_idx - 1)) {
      if (complements[idx[j]]) next
      this_val <- rulevars[, idx[j]]
      is_compl <- which(apply(rulevars[, idx[(j + 1):n_idx], drop = FALSE], 2, function(x) all(x != this_val))) + j
      if (length(is_compl) > 0) complements[idx[is_compl]] <- TRUE
    }
  }
  
  rules <- rules[!complements]
  rulevars <- rulevars[, !complements, drop = FALSE]
  return(list(rules = rules, rulevars = rulevars))
}

### 7.1 Check near equality
is_almost_eq <- function(x, y, tolerance = sqrt(.Machine$double.eps)) {
  stopifnot(is.numeric(x), length(x) == 1L)
  x_abs <- abs(x)
  xy <- if (x_abs > tolerance) abs(x - y) / x_abs else abs(x - y)
  return(xy <= tolerance)
}

####################################################
#####  2: Helpful functions for Rule ensemble  #####
####################################################

### 1: Grouping the base function for train data

group.train <- function(rules, data){
  
  # Modified version of linear terms ("Winsorized" linear terms)
  var.trunc <- winsorize(data)
  
  # Transform rules to binary matrix
  rulevars <- tran_rules(dat = data, rules = rules)
  
  # Combine the rule terms and linear terms
  base.fun <- cbind(rulevars, var.trunc)
  base.names <- c(rules, colnames(var.trunc))
  colnames(base.fun) <- base.names
  
  # Indicator variables (Dummy)
  grp.ind <- fastDummies::dummy_cols(data[,2])[,-1]
  grp.class <- unique(sort(data[,2]))
  
  # Grouping the base function
  group.info.list <- list()
  group.var <- NULL
  
  ## Scale of the base function (Freidman & Popescu, 2008)
  var.sd <- apply(data[,-c(1:2)], 2, sd)
  x.scale <- c(rep(1,ncol(rulevars)), var.sd/0.4)
  
  # Normalized the base functions 
  base.fun <- scale(base.fun, scale = x.scale, center = FALSE) 
  
  # Grouping the base function
  group.base <- NULL
  for(g in 1:ncol(grp.ind)){
    base.fun.can <- base.fun*(grp.ind[,g])
    colnames(base.fun.can) <-  paste0(base.names,"_t", grp.class[g])
    group.base <- cbind(group.base, base.fun.can)
  }
  
  # Summarizing the grouping information
  group.info <- list(base = group.base, 
                     index = rep(1:ncol(base.fun),ncol(grp.ind)), 
                     scale = rep(x.scale,ncol(grp.ind)))
  group.info
}

### 2: Grouping the base function for test data

group.test <- function(rules, data, group.num){
  
  # Transform rules to binary matrix
  base.fun <- tran_rules(dat = data, rules = rules)
  
  # Group class
  grp.class <- (1:group.num) - 1 
  
  # Grouping the base function
  group.base.list <- list()
  
  for(gg in 1:group.num){
    
    grp.ind <- data.frame(matrix(0, nrow(data), group.num))
    grp.ind[,gg] <- 1
    
    # Grouping the base function
    group.base <- NULL
    for(g in 1:ncol(grp.ind)){
      base.fun.can <- base.fun*(grp.ind[,g])
      colnames(base.fun.can) <-  paste0(colnames(base.fun),"_t", grp.class[g])
      group.base <- cbind(group.base, base.fun.can)
    }
    group.base.list[[gg]] <- group.base
  }  
  group.base.list
}

### 3: Function for "Winsorization" (truncate outlier)

winsorize <- function(data, win = 0.025) {
  # Apply winsorization to each column in the linear terms data frame
  linear <- data[, -c(1:2)]
  linear.win <- apply(linear, 2, function(x) {
    if (length(unique(x)) > 3) {
      upper <- quantile(x, 1 - win, na.rm = TRUE)
      bottom <- quantile(x, win, na.rm = TRUE)
      x[x >= upper] <- upper
      x[x <= bottom] <- bottom
    }
    return(x)
  })
  return(linear.win)
}

### 4: Function of adaptive group elastic-net (Modified based on R packages "msaenet")###

agenet <- function(
    X, y,
    init = c("genet", "glasso"),
    fin = c("genet", "glasso"),
    alphas = seq(0.05, 0.95, 0.05),
    tune.init = c("CV", "EBIC", "BIC", "AIC", "GCV"),
    tune.fin = c("CV", "EBIC", "BIC", "AIC", "GCV"),
    nfolds = 10L,
    group = 1:ncol(X),
    penalty.factor.init = sqrt(table(group)),
    scale = 1,
    seed = 1,
    verbose = FALSE
){
  
  start <- proc.time()
  
  init <- match.arg(init)
  tune.init <- match.arg(tune.init)
  fin <- match.arg(fin)
  tune.fin <- match.arg(tune.fin)
  call <- match.call()
  
  if (verbose) cat("Starting step 1 ...\n")
  
  if (init == "genet") {
    genet.cv <- msaenet.tune.grpreg(
      X = X, y = y,
      group = group,
      alphas = alphas,
      tune = tune.init,
      nfolds = nfolds,
      penalty.factor = as.numeric(penalty.factor.init),
      seed = seed
    )
  }
  
  if (init == "glasso") {
    genet.cv <- msaenet.tune.grpreg(
      X = X, y = y,
      group = group,
      alphas = 1,
      tune = tune.init,
      nfolds = nfolds,
      penalty.factor = penalty.factor.init,
      seed = seed
    )
  }
  
  best.alpha.genet <- genet.cv$"best.alpha"
  best.lambda.genet <- genet.cv$"best.lambda"
  step.criterion.genet <- genet.cv$"step.criterion"
  
  genet.full <- grpreg::grpreg(
    X = X, y = y,
    group = group,
    alpha = best.alpha.genet,
    lambda = best.lambda.genet,
    group.multiplier = penalty.factor.init,
    returnX = TRUE
  )
  
  end <- proc.time() - start
  
  bhat <- as.matrix(genet.full$"beta")[-1]/(genet.full$XG$scale + .Machine$double.eps*2)
  beta.weight <- abs(bhat)
  
  adpen <- c()
  for(num in 1:length(unique(group))){
    g <- unique(group)[num]
    adpen.each <- (sqrt(sum((beta.weight[which(g == group)])^2)))^(-scale)
    adpen <- c(adpen, adpen.each)
  }
  
  ### Adaptive Shrinkage ####
  
  if (verbose) cat("Starting step 2 ...\n")
  
  if (fin == "genet") {
    
    ### Adaptive group elastic net ####
    
    agenet.cv <- msaenet.tune.grpreg(
      X = X, y = y,
      group = group,
      alphas = alphas,
      tune = tune.fin,
      nfolds = nfolds,
      seed = seed + 1L,
      penalty.factor = adpen
    )
    
  }    
  
  if (fin == "glasso") {
    
    ### Adaptive group lasso ####
    
    agenet.cv <- msaenet.tune.grpreg(
      X = X, y = y,
      group = group,
      alphas = 1,
      tune = tune.fin,
      nfolds = nfolds,
      penalty.factor = penalty.factor.init,
      seed = seed + 1L
    )
    
  } 
  
  best.alpha.agenet <- agenet.cv$"best.alpha"
  best.lambda.agenet <- agenet.cv$"best.lambda"
  step.criterion.agenet <- agenet.cv$"step.criterion"
  
  agenet.full <- grpreg::grpreg(
    X = X, y = y,
    group = group,
    alpha = best.alpha.agenet,
    lambda = best.lambda.agenet,
    group.multiplier = adpen
  )
  
  # final beta stored as sparse matrix
  bhat.full <- agenet.full$"beta"
  
  
  agenet.model <- list(
    "beta" = bhat.full,
    "model" = agenet.full,
    "beta.first" = genet.full$"beta",
    "model.first" = genet.full,
    "best.alpha.genet" = best.alpha.genet,
    "best.alpha.agenet" = best.alpha.agenet,
    "best.lambda.genet" = best.lambda.genet,
    "best.lambda.agenet" = best.lambda.agenet,
    "step.criterion" = c(step.criterion.genet, step.criterion.agenet),
    "adpen" = adpen,
    "seed" = seed,
    "time.ini" = end
  )
  
  agenet.model
}

### 4.1 Function of adaptive weight

msaenet.tune.grpreg <- function(
    X, y,
    group=1:ncol(X),
    alphas,
    tune,
    nfolds,
    penalty.factor,
    seed) {
  
  if (tune == "CV") {
    
    model.list <- vector("list", length(alphas))
    for (i in 1L:length(alphas)) {
      set.seed(seed)
      model.list[[i]] <- grpreg::cv.grpreg(
        X = X, y = y,
        group = group,
        nfolds = nfolds, alpha = alphas[i], nlamba = 200,
        group.multiplier = penalty.factor
      )
    }
    
    errors <- unlist(lapply(model.list, function(x) min(x$"cve")))
    errors.min.idx <- which.min(errors)
    
    best.model <- model.list[[errors.min.idx]]
    best.alpha <- alphas[errors.min.idx]
    best.lambda <- best.model$"lambda.min"
    
    step.criterion <- errors[errors.min.idx]
    
    print(tune)
    
  } else {
    
    model.list <- vector("list", length(alphas))
    for (i in 1L:length(alphas)) {
      set.seed(seed)
      model.list[[i]] <- grpreg::grpreg(
        X = X, y = y,
        group = group,
        alpha = alphas[i],
        nlamba = 200,
        group.multiplier = penalty.factor
      )
    }
    
    crits <- unlist(lapply(model.list, function(x, tune){min(grpreg::select(x, criterion = tune)$IC)}, tune = tune))
    crits.min.idx <- which.min(crits)
    
    best.model <- model.list[[crits.min.idx]]
    best.alpha <- alphas[crits.min.idx]
    best.lambda <- grpreg::select(best.model, criterion = tune)$lambda
    
    step.criterion <- crits[crits.min.idx]
    
    print(tune)
  }
  
  list(
    "best.model" = best.model,
    "best.alpha" = best.alpha,
    "best.lambda" = best.lambda,
    "step.criterion" = step.criterion
  )
}




