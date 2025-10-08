#' analyze_time_against_chance()
#'
#' Runs a test on each time-bin of \code{time_sequence_data}. Supports \code{t.test},
#' \code{wilcox.test}, \code{(g)lm}, and \code{(g)lmer}. Also includes support for
#' the "bootstrapped-splines" test (see \code{?make_boot_splines_data} and
#' \href{http://www.eyetracking-r.com/vignettes/divergence}{the divergence vignette} for more info).
#' By default, this function uses 'proportion-looking' (\code{Prop}) as the DV, which can be changed
#' by manually specifying the formula. Results can be plotted to see how test-results or parameters
#' estimates vary over time. P-values can be adjusted for multiple comparisons with \code{p_adjust_method}.
#'
#' @export

analyze_time_against_chance = function(data, ...) {
  UseMethod("analyze_time_against_chance")
}
#' @describeIn analyze_time_against_chance
#'
#' @param data   The output of the 'make_time_sequence_data' function
#' @param predictor_column  The variable whose test statistic you are interested in. If you are not
#'   interested in a predictor, but the intercept, you can enter "intercept" for this argument.
#'   Interaction terms are not currently supported.
#' @param test              What type of test should be performed in each time bin? Supports
#'   \code{t.test}, \code{wilcox.test}, \code{(g)lm}, and \code{(g)lmer}. Also includes support for
#' the "bootstrapped-splines" test (see \code{?make_boot_splines_data} and
#' \href{http://www.eyetracking-r.com/vignettes/divergence}{the divergence vignette} for more info).
#' @param threshold         Value of statistic used in determining significance
#' @param alpha             Alpha value for determining significance, ignored if threshold is given
#' @param aoi               Which AOI should be analyzed? If not specified (and dataframe has multiple AOIs),
#'                          then AOI should be a predictor/covariate in your model (so `formula` needs
#'                          to be specified).
#' @param formula           What formula should be used for the test? Optional for all but
#'   \code{(g)lmer}, if unset will use \code{Prop ~ [predictor_column]}. Change this if you want to use a custom DV.
#' @param treatment_level   If your predictor is a factor, regression functions like `lm` and `lmer` by default will
#'                          treatment-code it. One option is to sum-code this predictor yourself before entering it
#'                          into this function. Another is to use the `treatment_level` argument, which specifies the
#'                          level of the predictor. For example, you are testing a model where `Target` is a predictor,
#'                          which has two levels, 'Animate' and 'Inanimate'. R will code 'Animate' as the reference
#'                          level, and code 'Inanimate' as the treatment level. You'd therefore want to set
#'                          `treatment_level = Inanimate`.
#' @param p_adjust_method   Method to adjust p.values for multiple corrections (default="none").
#'                          See \code{p.adjust.methods}.
#' @param quiet             Should messages and progress bars be suppressed? Default is to show
#' @param ...               Any other arguments to be passed to the selected 'test' function (e.g.,
#'   paired, var.equal, etc.)
#'
#' @examples
#' \dontrun{
#' data(word_recognition)
#' data <- make_eyetrackingr_data(word_recognition,
#'                                participant_column = "ParticipantName",
#'                                trial_column = "Trial",
#'                                time_column = "TimeFromTrialOnset",
#'                                trackloss_column = "TrackLoss",
#'                                aoi_columns = c('Animate','Inanimate'),
#'                                treat_non_aoi_looks_as_missing = TRUE
#' )
#' response_time <- make_time_sequence_data(data, time_bin_size = 250,
#'                                          predictor_columns = c("MCDI_Total"),
#'                                          aois = "Animate", summarize_by = "ParticipantName")
#' tb_analysis <- analyze_time_bins(response_time, predictor_column = "MCDI_Total",
#'                                  test = "lm", threshold = 2)
#' summary(tb_analysis)
#' }
#'
#' @export
#' @return A dataframe indicating the results of the test at each time-bin.

analyze_time_against_chance.time_sequence_data <- function(data,
                                                 test = 't.test',
                                                 threshold = NULL,
                                                 alpha = NULL,
                                                 aoi = NULL,
                                                 formula = NULL,
                                                 mu = 0.5,
                                                 treatment_level = NULL,
                                                 p_adjust_method = "none",
                                                 quiet = FALSE,
                                                 ...)
{

  ## Helper:
  .check_fix_nobs = function(data, data_options, predictor_column, test, dv, paired) {

    # if (g)lmer, don't need to do anything:
    if (test %in% c("lmer","glmer")) return(data)

    # otherwise, need to check if within-subjects:
    data_summarized_by <-  attr(data, "eyetrackingR")$summarized_by
    if (test %in% c("t.test", "wilcox.test") & identical(paired, TRUE)) {
      # if so, need to convert implicit NAs (not even in the data.frame) to explicit NAs:
      df_completed_cases <- expand.grid(Time = unique(data[["Time"]]),
                                        SumBy = unique(data[[attr(data, "eyetrackingR")$summarized_by]]),
                                        Pred = unique(data[[predictor_column]])
      )
      colnames(df_completed_cases) <- c("Time", data_summarized_by, predictor_column )
      save_attr <- attr(data, "eyetrackingR")
      save_class <- class(data)
      data <- left_join(df_completed_cases, data, by = colnames(df_completed_cases))
      attr(data, "eyetrackingR") <- save_attr
      class(data) <- save_class
    }

    # for lm, t.test, wilcox, need to check that there's no more than one observation per
    # participant per predictor-level per time bin.
    # df_grouped <- group_by_(data, .dots = paste0("`",c("Time", data_summarized_by, predictor_column), "`"))
    df_grouped <- group_by(data,
                           Time, !!!syms(data_summarized_by))
    df_summarized <- summarize(df_grouped, N = n()) %>% ungroup()
    if (any(df_summarized$N > 1)) {
      warning(test, " almost always requires no more than one observation per ", data_summarized_by,
              ", per level of ", predictor_column, ". However, your data has more than this for some time-bins. ",
              "Did you pass additional predictors to `make_time_sequence_data` that aren't being used here?")
      if (!quiet) print(arrange(df_summarized, -N))
    }
    return(data)
  }

  # Data info:
  data_options <- attr(data, "eyetrackingR")$data_options
  if (is.null(data_options)) stop("Dataframe has been corrupted.") # <----- TO DO: more informative?
  dots <- lazyeval::lazy_dots(...)

  # Need either alpha or threshold:
  if (!is.null(threshold) & !is.null(alpha)) stop("Please only specify alpha or threshold, not both.")
  if (is.null(threshold) & is.null(alpha)) stop("Please specify either alpha or threshold.")
  if (is.null(alpha)) {
    if (test %in% c('boot_splines', 'wilcox.test')) stop("This test requires `alpha` rather than `threshold`.")
  }

  # Which Test?:
  test <- match.arg(test, c("t.test","wilcox.test","lm","lmer","glm","glmer","boot_splines"))
  if (!requireNamespace("lme4", quietly = TRUE)) {
    if (test %in% c("lmer","glmer")) stop("Please install the 'lme4' package to use this method.")
  }

  # Multiple aois:
  if (!'AOI' %in% colnames(data)) stop("'AOI' column is missing from data.")
  unique_aois <- unique(data[['AOI']])
  if ( length(unique_aois) > 1 ) {           # DF has more than one AOI...
    if (is.null(aoi)) {                      # they haven't specified which is of interest
      if (predictor_column != "AOI") {       # its not the main predictor....
        if (!is.null(formula)) {             # they did specify a formula
          if (!grepl(pattern = 'AOI', x = deparse(formula))) {
            warning("There are multiple AOIs in your data frame, but your model does not use AOI as a predictor/covariate!",
                    immediate. = TRUE)
          }                                  # the formula has AOI-- no problem
        } else {                             # they did NOT specify a formula
          stop("If multiple AOIs in data, and `aoi` argument not specified, then you must manually specify a formula.")
        }
      }                                      # predictor column *was* AOI-- no problem.
    } else {                                 # they *did* specify which AOI is of interest
      data <- filter(data, AOI == aoi)
      if (nrow(data)==0) stop("AOI not found in data.")
    }
  }

  # Check that data is collapsed (by e.g. participants):
  if (!test %in% c("lmer","glmer")) {
    attrs <- attr(data, "eyetrackingR")
    summarized_by <- attrs$summarized_by
    if (is.null(summarized_by)) stop(test, " requires summarized data. ",
                                     "When using the 'make_time_sequence_data' function, please select an argument for 'summarize_by'",
                                     " (e.g., the participant column).")
  }

  # Formula Checks:
  # if (grepl("intercept", predictor_column, ignore.case = TRUE)) {
  #   if (is.null(formula)) stop("If testing intercept, please manually specify formula")
  #   predictor_column <- "(Intercept)"
  # }
  if (is.null(formula)) { # auto-make a formula, unless they specified one
    if (test%in% c("lmer","glmer")) stop("Must specify a formula if using (g)lmer.")
    formula <- as.formula(paste("Prop ~ 1"))
    dv <- "Prop"
  } else {
    dv <- gsub(formula[2], pattern = "()", replacement = "", fixed = TRUE)
  }

  # Run a model for each time-bin
  if (test!="boot_splines") {
    if (!quiet) message("Computing ", test, " for each time bin...")

    # Check for Number of Observations for Each Time-Bin x PPT (or whatever is summarize-by).
    # also, fix incomplete cases if paired = TRUE
    paired <- eval(dots[["paired"]]$expr)
    data <- .check_fix_nobs(data = data,
                            data_options = data_options,
                            test = test,
                            dv = dv,
                            paired = paired)


    # Get Testing Function
    if (test %in% c("lmer","glmer")) {
      the_func <- get(test, envir = getNamespace("lme4"))
    } else {
      the_func <- get(test)
    }

    # Create Test/Summarize function, which also catches errors:
    the_test <- function(formula, data, na.action, ...) {

      warn <- err <- NULL
      res <- withCallingHandlers(
        tryCatch(the_func(formula = formula, data = data, na.action = na.action, ... = ...),
                 error=function(e) {
                   err <<- conditionMessage(e)
                   NULL
                 }),
        warning=function(w) {
          warn <<- append(warn, conditionMessage(w))
          invokeRestart("muffleWarning")
        })

      res_err_warn <- list(res=res, warn=warn, err=err)

      if (length(res_err_warn$err)>0) {
        df_err <- as.data.frame( do.call(cbind, as.list(res_err_warn$err)) )
        colnames(df_err) <- paste0("ErrorMsg", seq_along(res_err_warn$err))
        out <- df_err
      } else {
        if (test %in% c("lmer","glmer")) {
          out <- broom.mixed::tidy(res_err_warn$res, effects = 'fixed')
        } else {
          out <- broom::tidy(res_err_warn$res)
        }
        # need to manually extract degrees of freedom for some tests
        # (not needed if func resulted in error)
        if (test %in% c("lm","glm")) out$parameter <- df.residual(res_err_warn$res)
      }

      if (length(res_err_warn$warn)>0) {
        df_warn <- as.data.frame( do.call(cbind, as.list(res_err_warn$warn)) )
        colnames(df_warn) <- paste0("WarningMsg", seq_along(res_err_warn$warn))
        df_warn <- bind_rows( lapply(X = 1:nrow(out), FUN = function(x) df_warn) )
        out <- bind_cols(out, df_warn)
      }

      out
    }

    # Run models:
    na_action <- ifelse(identical(paired,TRUE), 'na.pass', unlist(options("na.action")))

    if(test == 't.test'){
      # formula <- as.formula(paste("Prop ~", predictor_column))
      # dv <- "Prop"
      # lllev = unique(data[[predictor_column]])
      # formula = as.formula(paste('Pair(',dv,'[',predictor_column, '== \"',lllev[[1]],'\"] ,', dv,'[',predictor_column,' == \"',lllev[[2]],'\"]) ~ 1', sep = ''))

      df_models <- data %>%
        group_by(Time) %>%
        do(the_test(formula = formula, mu = mu, data = ., na.action = na_action)) %>%
        as.data.frame()
    }else{
      #patch for paired = F which is still not allowed in formula syntax
      # dotty <- list(...)
      # if(paired == F){
      #   dotty$paired = NULL
      # }
      # cat('FALSE') #test this here
      # cat('args are ', dotty)
      df_models <- data %>%
        group_by(Time) %>%
        do(the_test(formula = formula, data = ., na.action = na_action, ... = ...)) %>%
        as.data.frame()
    }
    # Warn about warnings
    warn_cols <- grep("WarningMsg", colnames(df_models))
    if ( length(warn_cols)>0 ) {
      if (!quiet) message("At least one time-bin produced warnings--be sure to check 'Warning' col in output.")
      for (col in warn_cols) {
        unique_msgs <- unique(as.character(df_models[,col]))
        if (length(unique_msgs)==1) warning("All time-bins produced the same warning: '", unique_msgs, "'")
      }
    }
    err_cols <- grep("ErrorMsg", colnames(df_models))
    if ( length(err_cols)>0 ) {
      if (!quiet) message("At least one time-bin produced errors--be sure to check 'Error' col in output.")
      for (col in err_cols) {
        unique_msgs <- unique(as.character(df_models[,col]))
        if (length(unique_msgs)==1) stop("All time-bins produced same error: '", unique_msgs, "'")
      }
    }

    # Make model specifications the same for all test-types:
    if (test == "t.test") {
      df_models$std.error <- with(df_models, (conf.high-conf.low)/(1.96*2))
      if (!identical(paired, TRUE)) df_models$estimate <- df_models$estimate1 - df_models$estimate2
    } else if (test %in% c('glmer', 'lmer')) {
      if (!quiet) message("Using the normal approximation for p-value on parameter in ", test,".")
      df_models$p.value <- with(df_models, 2*pnorm(q = abs(statistic), lower.tail = FALSE))
      df_models$parameter <- NA
    } else if (test == "wilcox.test") {
      df_models$statistic <- NA
      df_models$std.error <- NA
      df_models$estimate <- NA
      df_models$parameter <- NA
    }

    # Find Critical Value:
    if (is.null(threshold)) {
      df_models$CritStatisticPos <- .get_threshold(alpha, test, df_models$parameter, quiet)
    } else {
      df_models$CritStatisticPos <- ifelse(sign(threshold)==1,  threshold, -threshold)
    }
    df_models$CritStatisticNeg <- -df_models$CritStatisticPos

    # Filter Param-of-Interest:
    if (test %in% c('t.test', 'wilcox.test')) {
      df_models_this_param <- df_models
    } else {
      terms <- paste0("'", unique(df_models$term), "'", collapse= ", ")

      if (is.null(treatment_level)) {

        if (! predictor_column %in% df_models$term)
          stop("\nThe term '", predictor_column, "' was not found in your model.",
               "\nThis can happen if your predictor is treatment-coded, ",
               "in which case you should specify the `treatment_level` argument.",
               "\nFor example, if 'Target' is treatment coded, then it becomes 'TargetInanimate'",
               "in the model, and you should set `treatment_level = 'Inanimate'`.",
               "\nThe terms in the model were: ", terms)

        df_models_this_param <- filter(df_models, term == predictor_column )

      } else {

        if (! paste0(predictor_column, treatment_level) %in% df_models$term )
          stop("\nThe term '", paste0(predictor_column, treatment_level), "' was not found in your model.",
               "\nThe terms in the model were: ", terms)

        df_models_this_param <- filter(df_models, term == paste0(predictor_column, treatment_level))
      }

    }

    # Generate Output:
    new_cols <- list(Estimate = quote(estimate),
                     StdErr = quote(std.error),
                     Statistic = quote(statistic),
                     Prob = quote(p.value),
                     DF = quote(parameter),
                     quote(Time), quote(CritStatisticPos), quote(CritStatisticNeg)
    )
    new_cols <- append(new_cols, grep("WarningMsg", colnames(df_models), value = TRUE))
    new_cols <- append(new_cols, grep("ErrorMsg", colnames(df_models), value = TRUE))
    #out <- select_(.data = df_models_this_param, .dots = new_cols)
    out <- select(.data = df_models_this_param, !!!(new_cols))

  } else if (test=="boot_splines") {

    # arg-check:
    if (is.null( dots$within_subj$expr )) stop("Method 'boot_splines' requires you specify `within_subj`.")

    # Make Boot Splines:
    bs_samples <- eval(dots$bs_samples$expr)
    if (is.null(bs_samples)) {
      if (!is.null(dots$samples$expr)) {
        bs_samples <- eval(dots$samples$expr)
        warning("The 'samples' argument is deprecated in boot-splines, please use 'bs_samples' in the future.",
                call. = FALSE)
      }
    }
    the_args <- list(data = data, predictor_column = predictor_column, aoi = aoi, alpha = alpha)
    for (this_arg in names(dots)) {
      the_args[[this_arg]] <- dots[[this_arg]]$expr
    }
    the_args[["bs_samples"]] <- ifelse(is.null(bs_samples), 1000, bs_samples)
    bs_dat <- do.call(make_boot_splines_data, the_args)

    # Get Estimates:
    bs_anal <- analyze_boot_splines(bs_dat)

    out <- select(bs_anal,
                  Time, Estimate = MeanDiff, StdErr = SE, Statistic, Significant, ConfLow = CI_low, ConfHigh = CI_high)
    out <- mutate(out,
                  Prob = NA,
                  DF = NA,
                  CritStatisticNeg =  NA,
                  CritStatisticPos =  NA)

  }

  # Adjust P-val.
  if (p_adjust_method != "none") {
    if (test == "boot_splines") stop("The p_adjust_method must be 'none' for boot-splines test.")
    if (is.null(alpha)) stop("If specifying p-value adjustment, must give alpha, not threshold.")
    out$Prob <- p.adjust(p = out$Prob, method = p_adjust_method)
    out$CritStatisticNeg <- out$CritStatisticPos <- NA
  }

  # Compute Information about Runs:
  if (test == "boot_splines") {
    out$PositiveRuns <- .label_consecutive( out$Significant & out$Statistic>0 )
    out$NegativeRuns <- .label_consecutive( out$Significant & out$Statistic<0 )
  } else {
    if (is.null(alpha)) {
      out$PositiveRuns <- .label_consecutive(out$Statistic>out$CritStatisticPos)
      out$NegativeRuns <- .label_consecutive(out$Statistic<out$CritStatisticNeg)
    } else {
      out$PositiveRuns <- .label_consecutive((alpha>out$Prob) & out$Statistic>0)
      out$NegativeRuns <- .label_consecutive((alpha>out$Prob) & out$Statistic<0)
    }
  }

  positive_runs = lapply(unique(na.omit(out$PositiveRuns)), function(run) {
    list(start_time = with(out, min(Time[which(PositiveRuns==run)], na.rm=TRUE)),
         stop_time  = with(out, max(Time[which(PositiveRuns==run)], na.rm=TRUE)) + attr(data, "eyetrackingR")$time_bin_size
    )
  })
  negative_runs = lapply(unique(na.omit(out$NegativeRuns)), function(run) {
    list(start_time = with(out, min(Time[which(NegativeRuns==run)], na.rm=TRUE)),
         stop_time  = with(out, max(Time[which(NegativeRuns==run)], na.rm=TRUE)) + attr(data, "eyetrackingR")$time_bin_size
    )
  })

  # Rename Class, add attributes
  out <- as.data.frame(out)
  class(out) <- unique(c('bin_analysis', 'eyetrackingR_df', class(out)))
  attr(out,"eyetrackingR") <- c(
    attr(data, "eyetrackingR"),
    list(formula= formula,
         p_adjust_method = p_adjust_method,
         alpha = alpha,
         threshold = threshold,
         test = test,
         predictor = predictor_column,
         positive_runs = positive_runs,
         negative_runs = negative_runs)
  )
  out = out[order(out$Time),] # arrange chronologically if not already
  out
}

#' make_time_cluster_against_chance()
#'
#' Takes data that has been summarized into time-bins by \code{make_time_sequence_data()}, finds adjacent time
#' bins that pass some test-statistic threshold, and assigns these adjacent bins into groups (clusters).
#' Output is ready for a cluster permutation-based analyses (Maris & Oostenveld, 2007). Supports \code{t.test},
#' \code{wilcox.test}, \code{(g)lm}, and \code{(g)lmer}. Also includes support for
#' the "bootstrapped-splines" test (see \code{?make_boot_splines_data} and
#' \href{http://www.eyetracking-r.com/vignettes/divergence}{the divergence vignette} for more info).
#' By default, this function uses 'proportion-looking' (\code{Prop}) as the DV, which can be changed
#' by manually specifying the formula.
#'
#' @export

make_time_cluster_against_chance <-function(data, ...) {
  UseMethod("make_time_cluster_against_chance")
}
#' @describeIn make_time_cluster_against_chance Make data for time cluster analysis
#' @param data   The output of the \code{make_time_sequence_data} function
#' @param predictor_column  The column name containing the variable whose test statistic you are interested in.
#' @param aoi               Which AOI should be analyzed? If not specified (and dataframe has multiple AOIs),
#'                          then AOI should be a predictor/covariate in your model (so `formula` needs
#'                          to be specified).
#' @param test              What type of test should be performed in each time bin? Supports
#'                          \code{t.test}, \code{(g)lm}, or \code{(g)lmer}. Also includes experimental support for
#'                          the "bootstrapped-splines" test (see \code{?make_boot_splines_data} and
#'                          \href{http://www.eyetracking-r.com/vignettes/divergence}{the divergence vignette}
#'                          for more info). Does not support \code{wilcox.test}.
#' @param threshold         Time-bins with test-statistics greater than this amount will be grouped into clusters.
#' @param formula           What formula should be used for test? Optional (for all but \code{(g)lmer}), if unset
#'   uses \code{Prop ~ [predictor_column]}
#' @param treatment_level   If your predictor is a factor, regression functions like `lm` and `lmer` by default will
#'                          treatment-code it. One option is to sum-code this predictor yourself before entering it
#'                          into this function. Another is to use the `treatment_level` argument, which specifies the
#'                          level of the predictor. For example, you are testing a model where `Target` is a predictor,
#'                          which has two levels, 'Animate' and 'Inanimate'. R will code 'Animate' as the reference
#'                          level, and code 'Inanimate' as the treatment level. You'd therefore want to set
#'                          `treatment_level = Inanimate`.
#' @param ...               Any other arguments to be passed to the selected 'test' function (e.g., paired,
#'   var.equal, etc.)
#'
#' @examples
#' \dontrun{
#' data(word_recognition)
#' data <- make_eyetrackingr_data(word_recognition,
#'                                participant_column = "ParticipantName",
#'                                trial_column = "Trial",
#'                                time_column = "TimeFromTrialOnset",
#'                                trackloss_column = "TrackLoss",
#'                                aoi_columns = c('Animate','Inanimate'),
#'                                treat_non_aoi_looks_as_missing = TRUE )
#' response_window <- subset_by_window(data, window_start_time = 15500, window_end_time = 21000,
#'                                     rezero = FALSE)
#'
#' # identify clusters in the sequence data using a t-test with
#' # threshold t-value of 2
#'
#' # (note: t-tests require a summarized dataset)
#' response_time <- make_time_sequence_data(response_window, time_bin_size = 500, aois = "Animate",
#'                                          predictor_columns = "Sex",
#'                                          summarize_by = "ParticipantName")
#'
#' time_cluster_data <- make_time_cluster_data(data = response_time,
#'                                             predictor_column = "Sex",
#'                                             aoi = "Animate",
#'                                             test = "t.test",
#'                                             threshold = 2
#' )
#'
#' # identify clusters in the sequence data using an lmer() random-effects
#' # model with a threshold t-value of 1.5.
#'
#' # random-effects models don't require us to summarize
#' response_time <- make_time_sequence_data(response_window, time_bin_size = 500, aois = "Animate",
#'                                          predictor_columns = "Sex")
#'
#' # but they do require a formula to be specified
#' time_cluster_data <- make_time_cluster_data(data = response_time,
#'                            predictor_column = "SexM",
#'                            aoi = "Animate",
#'                            test = "lmer",
#'                            threshold = 1.5,
#'                            formula = LogitAdjusted ~ Sex + (1|Trial) + (1|ParticipantName)
#' )
#' }
#'
#' @export
#' @return The original data, augmented with information about clusters. Calling summary on this data will
#'   describe these clusters. The dataset is ready for the \code{\link{analyze_time_clusters}} method.
#'
make_time_cluster_against_chance.time_sequence_data <- function(data,
                                                      aoi = NULL,
                                                      threshold = NULL,
                                                      test = 't.test',
                                                      mu = 0.5,
                                                      treatment_level = NULL,
                                                      ...) {


  attrs <- attr(data, "eyetrackingR")
  data_options <- attrs$data_options
  dots <- lazyeval::lazy_dots(...)
  if (is.null(data_options)) stop("Dataframe has been corrupted.") # <----- TO DO: fix later
  if (test == "wilcox.test") stop("The wilcox.test is not supported for cluster-analysis.")
  if (test == "boot_splines") {
    if (is.null(dots$alpha$expr)) stop("For boot_splines test, please specify alpha instead of threshold.")
  } else {
    if (is.null(threshold)) stop("Please specify threshold for this test.")
  }

  # Filter Data:
  if (is.null(aoi)) {
    if (length(unique(data$AOI)) == 1) {
      aoi <- unique(data$AOI)
    }
  } else {
    data <- filter(data, AOI == aoi)
    if (nrow(data)==0) stop("AOI not found in data.")
  }

  # Compute Time Bins:
  the_args <- list(data = data, test = 't.test', formula = 'Prop ~ 1',
                   aoi = aoi, treatment_level = treatment_level, quiet = TRUE)
  the_args[["threshold"]] <- threshold
  for (this_arg in names(dots)) {
    the_args[[this_arg]] <- dots[[this_arg]]$expr
  }
  time_bin_summary <- do.call(analyze_time_against_chance, the_args)
  time_bin_summary <- rename(time_bin_summary, ClusterPos = PositiveRuns, ClusterNeg = NegativeRuns)
  formula <- attr(time_bin_summary, "eyetrackingR")$formula

  # Compute Sum Statistic for each Positive Cluster
  sum_stat_pos <- c()
  for (clust in na.omit(unique(time_bin_summary$ClusterPos))) {
    sum_stat_pos[clust] <- sum(time_bin_summary$Statistic[which(time_bin_summary$ClusterPos==clust)])
  }
  sum_stat_neg <- c()
  for (clust in na.omit(unique(time_bin_summary$ClusterNeg))) {
    sum_stat_neg[clust] <- sum(time_bin_summary$Statistic[which(time_bin_summary$ClusterNeg==clust)])
  }

  # Merge cluster info into original data
  df_timeclust <- left_join(data, time_bin_summary[,c('Time','ClusterNeg','ClusterPos')], by=c('Time'))

  # Collect info about each cluster:
  clusters <-
    data.frame(Cluster = seq_along(c(sum_stat_pos, sum_stat_neg)),
               Direction = unlist(c(rep("Positive", length.out = length(sum_stat_pos)),
                                    rep("Negative", length.out = length(sum_stat_neg)))),
               SumStatistic = c(sum_stat_pos, sum_stat_neg),
               StartTime = unlist(c(sapply(seq_along(sum_stat_pos),
                                           FUN = function(clust) time_bin_summary$Time[first(which(time_bin_summary$ClusterPos==clust))] ),
                                    sapply(seq_along(sum_stat_neg),
                                           FUN = function(clust) time_bin_summary$Time[first(which(time_bin_summary$ClusterNeg==clust))] ) )),
               EndTime = unlist(c(sapply(seq_along(sum_stat_pos),
                                         FUN = function(clust) time_bin_summary$Time[last(which(time_bin_summary$ClusterPos==clust))] ),
                                  sapply(seq_along(sum_stat_neg),
                                         FUN = function(clust) time_bin_summary$Time[last(which(time_bin_summary$ClusterNeg==clust))] ) ))
    )
  clusters$EndTime <-clusters$EndTime + attrs$time_bin_size

  # Output data, add attributes w/ relevant info
  df_timeclust <- as.data.frame(df_timeclust)
  class(df_timeclust) <- unique(c("time_cluster_data", "time_sequence_data", "eyetrackingR_df", class(df_timeclust)))
  attr(df_timeclust, "eyetrackingR") <- c(attrs,
                                          list(clusters = clusters,
                                               predictor_column = predictor_column,
                                               test = test,
                                               threshold = threshold,
                                               alpha = dots$alpha$expr,
                                               formula = formula,
                                               time_bin_summary = time_bin_summary,
                                               the_dots = names(dots),
                                               treatment_level = treatment_level)
  )
  df_timeclust

}

