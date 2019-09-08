## treezy functions for GBM

#' Plot partial dependence for BRT models
#'
#' Uses the "partial_dependence" function to plot partial dependence for BRT models. Future work will be into finding a way to generalize these methods to rpart and randomForest models, as an S3 method. This code is bespoke at the moment, and isn't designed as a flexible way to create plots, so I would recommend that people who want to plot their own partial plots just use `partial_dependence` and go from there.
#'
#' @param x The GBM model to be used
#'
#' @param vars The variables used in the GBM model, this is a character vector
#'
#' @return a faceted ggplot plot of the variables
#'
#' @examples
#'
#' \dontrun{
#'
#' # using gbm.step from the dismo package
#'
#' library(gbm)
#' library(dismo)
#'
#' # load data
#'
#' data(Anguilla_train)
#' anguilla_train <- Anguilla_train[1:200,]
#'
#' # fit model
#' angaus_tc_5_lr_01 <- gbm.step(data = anguilla_train,
#'                               gbm.x = 3:14,
#'                               gbm.y = 2,
#'                               family = "bernoulli",
#'                               tree.complexity = 5,
#'                               learning.rate = 0.01,
#'                               bag.fraction = 0.5)
#'
#' gg_partial_plot(angaus.tc5.lr01,
#'                    var = c("SegSumT",
#'                            "SegTSeas"))
#'
#'}
#'
#' @export

gg_partial_plot <- function(x,
                            vars, var.labs){
  
  df_box <- list("vector", length(vars))
  
  for (i in (1:length(vars))){
    
    df_box[[i]] <- partial_dependence(x, vars[[i]])
    
  }
  
  df <- dplyr::bind_rows(df_box)
  
  # make another
  
  df_mean <-
    df %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(mean = mean(fitted_function))
  df_mean$variable <- factor(df_mean$variable, levels = unique(df_mean$variable), labels = var.labs)
  #browser()
  df$variable <- factor(df$variable, levels = unique(df$variable), labels = var.labs)
  ggplot2::ggplot(data = df,
                  ggplot2::aes(x = value,
                               y = fitted_function)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~variable,
                        ncol = 2,
                        scales = "free_x") +
    ggplot2::geom_hline(data = df_mean,
                        ggplot2::aes(yintercept = mean),
                        colour = "red",
                        linetype = "dashed",
                        alpha = 0.75) +
    ggplot2::labs(x = "Parameter Values",
                  y = "Model Predicted Values")
  
}  
  #' partial_dependence
  #'
  #' @description : Some code that returns the partial dependence values for a given set of variables for a gbm.step model. In the future this function will work for other decision trees
  #'
  #' @param x a gbm.step object
  #' @param var a set of variables you want to retrieve partial dependence for
  #' @param ... extra arguments you might want to pass downstream
  
  #'
  #' @note This requires the loading of the `gbm.step` function. Hopefully sometime soom I can just write this in vanilla R myself. Future extensions will allow for this function to work for `rpart`, `gbm`, `gbm.step`, and `randomForest`.
  #'
  #'
  #' @export
partial_dependence <- function(x, var, ...) UseMethod("partial_dependence") 
  
  #' Default method for partial dependence
  #'
  #' @param x a decision tree object
  #' @param var a set of variables you want to retrieve partial dependence for
  #' @param ... extra arguments you might want to pass downstream
  #' @export
  
  partial_dependence.default <- function(x, var, ...){
    
    # grab the name sof the variables in the dataframe used in the model, and give their vector columns position to `i`
    ### x = ###angaus.tc5.lr01
    ### var = ###"SegSumT"
    i <- which(x$var.names == var)
    
    # get the matrix out of the `plot.gbm`
    response_matrix <- gbm::plot.gbm(x,
                                     i.var = i,
                                     n.trees = x$n.trees,
                                     return.grid = TRUE)
    
    # make a dataframe, which contains the observed calues, and the fitted function values, and then adds another column containing the variable name.
    df <- data.frame(value = as.numeric(response_matrix[ , 1]),
                     fitted_function = response_matrix[ , 2])
    df$variable <- x$var.names[i]
    # df$variable <- factor(df$variable, levels = as.character(df$variable))
    
    return(df)
    
  } # end of thingy
  
  #' Partian dependence method for caret::train objects
  #'
  #' @param x an object of class train
  #' @param var a set of variables you want to retrieve partial dependence for
  #' @param ... extra arguments you might want to pass downstream
  #'
  #' @return a dataframe of partial dependence
  #' @export
  #'
  partial_dependence.train <- function(x, var, ...){
    
    # grab the name sof the variables in the dataframe used in the model, and give their vector columns position to `i`
    
    i <- which(x$finalModel$var.names == var)
    
    # get the matrix out of the `plot.gbm`
    response_matrix <- gbm::plot.gbm(x,
                                     i.var = i,
                                     n.trees = x$n.trees,
                                     return.grid = TRUE)
    
    # make a dataframe, which contains the observed calues, and the fitted function values, and then adds another column containing the variable name.
    df <- data.frame(value = as.numeric(response_matrix[ , 1]),
                     fitted_function = response_matrix[ , 2]) %>%
      dplyr::mutate(variable = x$var.names[i])
    
    return(df)
    
  } # end of thingy
  
# end function
  
  #' importance_plot
  #'
  #' importance_plot make a graph of variable importance
  #'
  #' takes an `rpart` or `gbm.step` fitted object and makes a plot of variable importance
  #'
  #' @param x is an rpart or gbm.step object
  #' @param ... extra functions or arguments
  #'
  #' @return a ggplot plot of the variable importance
  #'
  #' @examples
  #'
  #' \dontrun{
  #'  # an rpart object
  #'  library(rpart)
  #'  library(treezy)
  #' fit.rpart <- rpart(Kyphosis ~ Age + Number + Start, data = kyphosis)
  #'
  #' importance_plot(fit.rpart)
  #'
  #' # you can even use piping
  #'
  #' fit.rpart %>% importance_plot
  #'
  #'  # a randomForest object
  #'
  #'  set.seed(131)
  #'   ozone.rf <- randomForest(Ozone ~ ., data=airquality, mtry=3,
  #'                            importance=TRUE, na.action=na.omit)
  #'   print(ozone.rf)
  #'   ## Show "importance" of variables: higher value mean more important:
  #'   importance(ozone.rf)
  #'
  #'   ## use importance_table
  #'
  #'   importance_table(ozone.rf)
  #'
  #'   # now plot it
  #'
  #'   importance_plot(ozone.rf)
  #'
  #'}
  #' @export
  #'
  importance_plot <- function(x, ...) UseMethod("importance_plot")
  
  #' @export
  importance_plot.default <- function(x, ...){
    # x = fit_rpart_kyp
    # library(dplyr)
    importance_table(x) %>%
      ggplot2::ggplot(ggplot2::aes(x = stats::reorder(variable,
                                                      importance),
                                   # make sure the plot is ordered by most important
                                   y = importance))+
      ggplot2::geom_bar(stat="identity",
                        position="dodge",
                        width = 0,
                        colour = "black") +
      ggplot2::geom_point() +
      ggplot2::labs(x = "Variables",
                    y = "Importance Score") +
      ggplot2::coord_flip()
    
  } # end function
  
  #' @export
  importance_plot.randomForest <- function(x, ...){
    
    # get names of columns (which changes according to the type of RF model)
    # new_cols <-
    # importance_table(x) %>%
    #     colnames %>%
    #     grep("variable",
    #          x = .,
    #          invert = TRUE,
    #          value=TRUE)
    
    # importance_table(x) %>%
    #     tidyr::gather_(data = .,
    #                    key_col = "importance_metric",
    #                    value_col = "importance",
    #                    gather_cols = new_cols) %>%
    importance_table(x = x,
                     importance_metric = TRUE) %>%
      ggplot2::ggplot(ggplot2::aes(x = stats::reorder(variable,
                                                      importance),
                                   y = importance)) +
      # ggplot2::geom_point() +
      ggplot2::geom_bar(stat="identity",
                        position="dodge",
                        width = 0,
                        colour = "black") +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(~ importance_metric,
                          scales = "free_x") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45)) +
      ggplot2::labs(x = "Variables",
                    y = "Importance Score") +
      ggplot2::coord_flip()
    
  } # end function

  #' importance_table
  #'
  #' importance_table returns a data_frame of variable importance for decision trees methods rpart randomForest, and gbm.step (from the dismo package).
  #'
  #' @note treezy currently only works for rpart and gbm.step functions. In the
  #' future more features will be added so that it works for many
  #' decision trees
  #'
  #' @param x An rpart, randomForest, or gbm.step, or model
  #' @param ... extra functions or arguments
  #' @return A tibble containing the importance score made with the intention of turning it into a text table with `knitr` or `xtable`
  #'
  #' @examples
  #'
  #' # retrieve a tibble of the variable importance from a decision tree model
  #' \dontrun{
  #'  # rpart object
  #'  library(rpart)
  #' fit_rpart <- rpart(Kyphosis ~ Age + Number + Start, data = kyphosis)
  #'
  #' importance_table(fit_rpart)
  #'
  #' # you can even use piping
  #'
  #' fit_rpart %>% importance_table
  #'
  #' # gbm.step object
  #'
  #' library(dismo)
  #' library(gbm)
  #'
  #' fit_gbm_step <- gbm.step(data = iris,
  #'                          gbm.x = c(1:3),
  #'                          gbm.y = 4,
  #'                          tree.complexity = 1,
  #'                          family = "gaussian",
  #'                          learning.rate = 0.01,
  #'                          bag.fraction = 0.5)
  #'
  #' importance_table(fit_gbm_step)
  #'
  #' # with piping
  #' fit.gbm.step %>% importance_table
  #'
  #' Unfortunately it cannot yet run a gbm object:
  #'
  #'
  #' gbm.fit <- gbm(Sepal.Width ~ .,
  #'                distribution = "gaussian",
  #'                data = iris)
  #'
  #' importance_table(gbm.fit)
  #'
  #'
  #' #A randomForest object
  #'     set.seed(1)
  #'     data(iris)
  #'     fit_rf <- randomForest(Species ~ ., iris,
  #'                             proximity=TRUE,
  #'                             keep.forest=FALSE)
  #'
  #' importance_table(fit_rf)
  #'
  #' }
  #' @note https://github.com/dgrtwo/broom
  #'
  #' @export
  importance_table <- function(x, ...) UseMethod("importance_table")
  
  #' @export
  importance_table.NULL <- function(x, ...) NULL
  
  #' @export
  importance_table.default <- function(x, ...) {
    
    stop("alas, importance_table does not know how to deal with data of class ", class(x), call. = FALSE)
    
  }
  
  # rpart -----------------------------------------------------------------------
  
  #' @export
  importance_table.rpart <- function(x, ...){
    
    # Some trees are stumps, we need to skip those that are NULL (stumps)
    # so here we say, "If variable importance is NOT NULL, do the following"
    # Another option would be to only include those models which are not null.
    
    if (is.null(x$variable.importance) == FALSE) {
      
      x <-
        x$variable.importance %>%
        data.frame(variable = names(x$variable.importance),
                   importance = as.vector(x$variable.importance),
                   row.names = NULL) %>%
        dplyr::select(variable,
                      importance) %>%
        dplyr::as_data_frame()
      
    } else {
      
      # if rpart_frame just contains a decision stump, make NULL datasets.
      
      x <- data.frame(variable = NULL,
                      importance = NULL,
                      row.names = NULL) %>%
        dplyr::as_data_frame()
    } # end else
    
    # no need to modify the class of the object.
    # res <- x
    # class(res) <- c("imp_tbl", class(res))
    # return(res)
    return(x)
    
  }
  
  # gbm -------------------------------------------------------------------------
  
  #' @export
  importance_table.gbm <- function(x, ...){
    
    x <-
      x$contributions %>%
      # make it a dataframe
      dplyr::as_data_frame() %>%
      # rename the variables
      dplyr::rename(variable = var,
                    importance = rel.inf) %>%
      dplyr::as_data_frame()
    
    # no need to return this info
    # res <- x
    # class(res) <- c("imp_tbl", class(res))
    # return(res)
    return(x)
  }
  
  # random forests --------------------------------------------------------------
  
  #' @note you can pass importance_metric = FALSE or TRUE, if you want to show the importance metrics
  #' @export
  importance_table.randomForest <- function(x, importance_metric = FALSE, ...){
    
    
    # get the names of the variables used
    variable <- randomForest::importance(x) %>%
      row.names %>%
      data.frame(variable = .)
    
    pt1 <- randomForest::importance(x) %>%
      as.data.frame(row.names = F) %>%
      dplyr::as_data_frame()
    
    imp_tbl <- dplyr::bind_cols(variable,
                                pt1) %>%
      dplyr::as_data_frame()
    
    # don't rearrange / gather
    if (importance_metric == FALSE) {
      
      return(imp_tbl)
      
      # gather to get different info
    } else if(importance_metric == TRUE){
      
      new_cols <- imp_tbl %>%
        colnames %>%
        grep("variable",
             x = .,
             invert = TRUE,
             value=TRUE)
      
      aug_imp_tbl <-
        tidyr::gather_(data = imp_tbl,
                       key_col = "importance_metric",
                       value_col = "importance",
                       gather_cols = new_cols)
      
      return(aug_imp_tbl)
      
    }
    
    
    
    
  }
  
  
  # caret -----------------------------------------------------------------------
  
  #' @export
  
  importance_table.train <- function(x, ...){
    
    x <- caret::varImp(x)$importance %>%
      data.frame(variable = row.names(.),
                 importance = .,
                 row.names = NULL) %>%
      dplyr::rename(importance = Overall) %>%
      dplyr::as_data_frame()
    
    # no need to return this info
    # res <- x
    # class(res) <- c("imp_tbl", class(res))
    # return(res)
    return(x)
    
  }
  
  
  # future code for gbm objects, not gbm.step...
  ## new method for gbm to get importance value for .gbm methods (not gbm.step)
  #
  # view_importance <- function(x){
  #
  #   relative.influence(x) %>%
  #     as.data.frame %>%
  #     data.frame(importance = .[,1],
  #                variable = row.names(.),
  #                row.names = NULL) %>%
  #     select(variable,
  #            importance) %>%
  #     arrange(-importance)
  #