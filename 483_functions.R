#' run Function
#'
#' This function loads or runs and saves (if not previously done) a model, table,
#' or any computationally intensive task. It's designed to avoid redundant computations
#' by reusing previously saved results when possible.
#'
#' @param expr An expression representing the computationally intensive task to be executed.
#' @param path A string specifying the file path where the result should be saved or loaded from.
#' @param reuse A logical flag indicating whether to attempt reusing saved results to avoid recomputation.
#' @return The result of evaluating `expr`. If `reuse` is TRUE and a saved result exists, 
#'         that result is returned; otherwise, `expr` is evaluated.
#' @examples
#' # Assuming lm_result.Rds does not exist, this will compute the linear model and save it.
#' run(lm(mpg ~ cyl, data = mtcars), path = "lm_result", reuse = TRUE)

run <- function(expr, path, reuse = TRUE) {
  fit <- NULL
  if (reuse) {
    path <- paste0(path, ".Rds")
    fit <- suppressWarnings(try(readRDS(path), silent = TRUE))
    if (inherits(fit, "try-error")) {
      fit <- NULL
    }
    }
  if (is.null(fit)) {
    fit <- eval(substitute(expr))
    if (reuse && !is.null(path) && nzchar(path)) {
      saveRDS(fit, file = path)
    }
    }
  return(fit)
  }





#' perm_over 
#' 
#' 
perm_over <- function(outcome, 
                      predictors, 
                      nperm = 1000, 
                      seed = 123){
  
  pred_name <- list()
  r2 <- list()
  p_val <- list()
  
  # outcome = outcomes
  # predictors = data_predict
  # nperm = 100
  # seed = 483
  # 
  # i = 2

  for (i in seq_along(predictors)){
    
    set.seed(seed)
    tr <- vegan::adonis2(outcome ~ predictors[,i],
                  by = 'margin', 
                  permutations = nperm,
                  method = 'euclidean')
      
      pred_name[i] <- colnames(predictors[i])
      r2[[i]] <- tr$R2[1]
      p_val[[i]] <- tr$`Pr(>F)`[1]
  }
  
  tre <- data.frame(
    predictor = unlist(pred_name),
    r2 = unlist(r2),
    p_val = unlist(p_val)
    )
  
 return(tre)
}


## Univariate models

glm_iori_poisson <- function(outcomes, data){
  
  # outcomes <- outcomes
  # data <- data
  
  Outcome <- c()
  Smoking_logFC <- c()
  Smoking_P_val <- c()
  cmv_logFC <- c()
  cmv_P_val <- c()
  receiver_age_30y_logFC <- c()
  receiver_age_30y_P_val <- c()
  dialysis_logFC <- c()
  dialysis_P_val <- c()
  ascvd_logFC <- c()
  ascvd_P_val <- c()
  
  for (i in 1:length(outcomes)){
    
    data$outcome <- data[[outcomes[i]]]
    
    model <- glmmTMB::glmmTMB(
      outcome ~ 
        Smoking + 
        CMV + 
        receiver_age_30y + 
        dialysis + 
        ASCVD +
        (1|id),
      data = data,
      family = poisson(link = 'log'),
      control = glmmTMBControl(
        optCtrl = list(iter.max = 10000, eval.max = 10000)
        )
      )
    
    Smoking_logFC[i] <- summary(model)$coefficients$cond[2,1]
    cmv_logFC[i] <- summary(model)$coefficients$cond[3,1]
    receiver_age_30y_logFC[i] <- summary(model)$coefficients$cond[4,1]
    dialysis_logFC[i] <- summary(model)$coefficients$cond[5,1]
    ascvd_logFC[i] <- summary(model)$coefficients$cond[6,1]
    
    Smoking_P_val[i] <- summary(model)$coefficients$cond[2,4]
    cmv_P_val[i] <- summary(model)$coefficients$cond[3,4]
    receiver_age_30y_P_val[i] <- summary(model)$coefficients$cond[4,4]
    dialysis_P_val[i] <- summary(model)$coefficients$cond[5,4]
    ascvd_P_val[i] <- summary(model)$coefficients$cond[6,4]
    
    Outcome[i] <- outcomes[i]
  }
  
  res <- data.frame(
    Outcome, 
    Smoking_logFC, Smoking_P_val,
    cmv_logFC, cmv_P_val, 
    receiver_age_30y_logFC, receiver_age_30y_P_val,
    dialysis_logFC, dialysis_P_val,
    ascvd_logFC, ascvd_P_val
  )
  
  return(res)
}

replace_zeros <- function(x) {
  if (is.numeric(x)) {
    min_nonzero <- min(x[x > 0])
    x[x == 0] <- (2/3) * min_nonzero
  }
  return(x)
}


lm_lognormal <- function(outcomes, data){
  
  # outcomes <- outcomes
  # data <- data
  
  Outcome <- c()
  Smoking_logFC <- c()
  Smoking_P_val <- c()
  cmv_logFC <- c()
  cmv_P_val <- c()
  receiver_age_30y_logFC <- c()
  receiver_age_30y_P_val <- c()
  dialysis_logFC <- c()
  dialysis_P_val <- c()
  ascvd_logFC <- c()
  ascvd_P_val <- c()
  
  for (i in 1:length(outcomes)){
    
    data$outcome <- data[[outcomes[i]]]
    
    model <- lm(
      log(outcome) ~ 
        Smoking + 
        CMV + 
        receiver_age_30y + 
        dialysis + 
        ASCVD,
      data = data)
    
    Smoking_logFC[i] <- summary(model)$coefficients[2,1]
    cmv_logFC[i] <- summary(model)$coefficients[3,1]
    receiver_age_30y_logFC[i] <- summary(model)$coefficients[4,1]
    dialysis_logFC[i] <- summary(model)$coefficients[5,1]
    ascvd_logFC[i] <- summary(model)$coefficients[6,1]
    
    Smoking_P_val[i] <- summary(model)$coefficients[2,4]
    cmv_P_val[i] <- summary(model)$coefficients[3,4]
    receiver_age_30y_P_val[i] <- summary(model)$coefficients[4,4]
    dialysis_P_val[i] <- summary(model)$coefficients[5,4]
    ascvd_P_val[i] <- summary(model)$coefficients[6,4]
    
    Outcome[i] <- outcomes[i]
  }
  
  res <- data.frame(
    Outcome, 
    Smoking_logFC, Smoking_P_val,
    cmv_logFC, cmv_P_val, 
    receiver_age_30y_logFC, receiver_age_30y_P_val,
    dialysis_logFC, dialysis_P_val,
    ascvd_logFC, ascvd_P_val
  )
  
  return(res)
}

replace_zeros <- function(x) {
  if (is.numeric(x)) {
    min_nonzero <- min(x[x > 0])
    x[x == 0] <- (2/3) * min_nonzero
  }
  return(x)
}





