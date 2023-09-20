#' @title CoxAIPW
#'
#' @description Doubly robust estimation of two-group log hazard ratio under the Cox marginal structural model, works for observational data with informative censoring.
#'
#' @param data A matrix or dataframe object, where the first column is observed event time, second column is event indicator (1 = event, 0 = censored), third
#'     column is a binary group indicator A (numeric 0 or 1), the other columns are baseline covariates Z.
#' @param Tmod working conditional model for failure time T given A,Z, currently supports values ('Cox', 'Spline', 'RSF'), corresponding to Cox PH, hazard
#'     regression, and Random survival forest.
#' @param Cmod working conditional model for censoring time C given A,Z, currently supports values ('Cox', 'Spline', 'RSF'), corresponding to Cox PH, hazard
#'     regression, and Random survival forest.
#' @param PSmod working model for propensity score P(A=1|Z), currently supports values ('logit','RF','twang'), corresponding to logistic regression,
#'     random forest, and the package 'twang' (gradient boosted logistic regression).
#' @param tau the cutoff time of the study, which should be no less than all observed event times. Default to NULL, which sets tau to the largest observed event time.
#' @param k the number of folds used in cross-fitting, default to 5.
#' @param beta0 initial value for beta in the optimization algorithm, default to 0.
#' @param min.S the minimum threshold for estimated survival functions P(C>t|A,Z) and P(T>t|A,Z), any estimated value below min.S is set as min.S to avoid exploding gradient. Default to 0.05
#' @param min.PS the minimum threshold for estimated propensity score, any estimated value below min.PS or above 1-min.PS is set to min.PS or 1-min.PS to avoid exploding gradient.
#' Default to 0.1
#' @param weights a vector of length equal to the samples size that assigns weight to each observation. Also used in Bayesian Bootstrap. Default to NULL (equal weights).
#' @param augmentation a string indicating the type of augmention to consider. 'AIPTCW' performs augmented inverse probability of treatment and censoring weighting, which works for observational data with informative censoring. 'AIPTW' performs augmented inverse probability of treatment weighting, which only works under observational data with random censoring. 'AIPCW' performs augmented inverse probability of censoring weighting, which only works under randomized study with informative censoring. Default to 'AIPTCW'.
#' @param crossFit a logical variable indicating whether to apply cross-fitting. If set to FALSE, all three nuisance function estimators 'Tmod', 'Cmod', and 'PSmod' need to be root-n. Default to TRUE.
#' @return a list object containing 'beta', 'model_se', 'Lambda0', 'survival0', 'survival1' and 'beta_t_plot'. 'beta' is the estimated log hazard ratio, 'model_se' is the model-based asymptotic standard error estimator used for inference, 'Lambda0' is the estimated cumulative baseline hazard function, evaluated on each follow up time, and 'survival0' and 'survival1' are the estimated survival curves, evaluated on each follow up time, for treatment group 0 and 1 respectively. 'beta_t_plot' is a list made up of an x-values vector and a y-values vector, which leads to a plot of time-varying log hazard ratio beta(t) against time t after applying loess or spline smoothing.
#' @examples
#' # run Cox AIPW on the cancer data set from survival package
#' data(cancer, package = 'survival')
#' data = data.frame(
#'   time = cancer$time,
#'   status = cancer$status - 1,
#'   A = cancer$sex - 1,
#'   age = cancer$age)
#' aipw = CoxAIPW(data)
#'
#' # extract beta and model SE.
#' logHR = aipw$beta
#' modelSE = aipw$model_se
#'
#' # extract cumulative baseline hazard function.
#' cumBaseHaz = aipw$Lambda0
#'
#' # extract survival curve for group 0 and survival curve for group 1.
#' surv0 = aipw$survival0
#' surv1 = aipw$survival1
#'
#' # extract x and y values for a plot of time-varying log hazard ratio beta(t).
#' x = aipw$beta_t_plot$x
#' y = aipw$beta_t_plot$y
#' # library(splines)
#' # spline_fit = lm(y ~ ns(x, df=2, intercept =F))
#' # X = seq(0, max(x), length = 40)
#' # plot(X, predict(spline_fit, newdata = data.frame(x = X)), type = 'l')
#' @export
#' @import stats survival randomForestSRC polspline tidyr ranger pracma gbm

#> NULL

CoxAIPW <-
  function(data,
           Tmod = 'Cox',
           Cmod = 'Cox',
           PSmod = 'logit',
           tau = NULL,
           k = 5,
           beta0 = 0,
           min.S = 0.05,
           min.PS = 0.1,
           weights = NULL,
           augmentation = 'AIPTCW',
           crossFit = TRUE) {
    # Estimate PS, S_t and S_c using k-fold cross-fitting
    data = as.data.frame(data)#[sample(nrow(data)),]
    colnames(data)[1:3] = c('X', 'Delta', 'A')
    colnames(data)[-(1:3)] = paste0('Z', 1:(ncol(data) - 3))
    if (is.null(tau))
      tau = max(data$X)
    data$Delta[data$X > tau] = 0
    while (any(duplicated(data$X[data$X < tau])))
      data$X[duplicated(data$X)] = data$X[duplicated(data$X)] + 1e-6
    data$X[data$X > tau] = tau
    X.full = data[, 1]
    X.sort.full = sort(unique(X.full))
    folds = cut(1:nrow(data), breaks = k, labels = F)
    if (is.null(weights))
      weights = rep(1, nrow(data))

    # Assume observational data, perform either AIPTCW or AIPTW
    if (augmentation %in% c('AIPTCW', 'AIPTW')) {
      PS.temp = rep(0, nrow(data))
      S_t.temp1 = S_t.temp0 = S_c.temp1 = S_c.temp0 = matrix(0, nrow = nrow(data), ncol = length(unique(X.full)))
      # sample split functions
      samp.split_S_t = function(dfs,
                                dfn.t,
                                n,
                                X,
                                n_res,
                                Delta,
                                Delta_c,
                                A,
                                Z,
                                p,
                                X.sort,
                                X.full,
                                X.sort.full,
                                Tmod,
                                weights.n) {
        S_t0 = S_t1 = data.frame(matrix(0, nrow = nrow(dfs), ncol = length(unique(X.full))))
        if (Tmod == 'Cox') {
          model = coxph(
            Surv(X, Delta) ~ .,
            timefix = F,
            data = dfn.t,
            weights = weights.n
          )
          Lambda0 = exp(as.matrix(cbind(0, dfs[, -(1:3)])) %*% model$coefficients) %*% t(basehaz(model, centered = F)$hazard)
          S_t0[, X.sort.full %in% X.sort] = exp(-Lambda0)
          Lambda1 = exp(as.matrix(cbind(1, dfs[, -(1:3)])) %*% model$coefficients) %*% t(basehaz(model, centered = F)$hazard)
          S_t1[, X.sort.full %in% X.sort] = exp(-Lambda1)
        }
        if (Tmod == 'Spline') {
          model = hare(X, Delta, cbind(A, Z))
          S_t0 = 1 - sapply(X.sort.full, function(x)
            phare(
              q = x,
              cov = cbind('A' = 0, dfs[, -(1:3)]),
              fit = model
            ))
          S_t0[is.nan(S_t0)] = NA
          S_t1 = 1 - sapply(X.sort.full, function(x)
            phare(
              q = x,
              cov = cbind('A' = 1, dfs[, -(1:3)]),
              fit = model
            ))
          S_t1[is.nan(S_t1)] = NA
        }
        if (Tmod == 'RSF') {
          dfs0 = dfs1 = dfs
          dfs0[, 3] = 0
          dfs1[, 3] = 1
          model = rfsrc(
            Surv(X, Delta) ~ .,
            data = dfn.t,
            ntime = 0,
            splitrule = 'bs.gradient',
            ntree = 100,
            nodesize = 50,
            mtry = 2,
            case.wt = weights.n
          )
          S_t0[, which(X.sort.full %in% model$time.interest)] = predict(model, newdata = dfs0)$survival
          S_t1[, which(X.sort.full %in% model$time.interest)] = predict(model, newdata = dfs1)$survival
        }
        S_t0[, 1] = ifelse(S_t0[, 1] == 0 |
                             is.na(S_t0[, 1]), 1, S_t0[, 1])
        S_t0[S_t0 == 0] = NA
        S_t0 = t(as.matrix(fill(
          data.frame(t(S_t0)), names(data.frame(t(S_t0)))
        )))
        S_t1[, 1] = ifelse(S_t1[, 1] == 0 |
                             is.na(S_t1[, 1]), 1, S_t1[, 1])
        S_t1[S_t1 == 0] = NA
        S_t1 = t(as.matrix(fill(
          data.frame(t(S_t1)), names(data.frame(t(S_t1)))
        )))
        return(list(S_t0, S_t1))
      }

      samp.split_S_c = function(dfs,
                                dfn.c,
                                n,
                                X,
                                n_res,
                                Delta,
                                Delta_c,
                                A,
                                Z,
                                p,
                                X.sort,
                                X.full,
                                X.sort.full,
                                Cmod,
                                weights.n) {
        S_c0 = S_c1 = data.frame(matrix(0, nrow = nrow(dfs), ncol = length(unique(X.full))))
        if (Cmod == 'Cox') {
          model = coxph(
            Surv(X, Delta_c) ~ .,
            timefix = F,
            data = dfn.c,
            weights = weights.n
          )
          Lambda_c0 = exp(as.matrix(cbind(0, dfs[, -(1:3)])) %*% model$coefficients) %*% t(basehaz(model, centered = F)$hazard)
          S_c0[, X.sort.full %in% X.sort] = exp(-Lambda_c0)
          Lambda_c1 = exp(as.matrix(cbind(1, dfs[, -(1:3)])) %*% model$coefficients) %*% t(basehaz(model, centered = F)$hazard)
          S_c1[, X.sort.full %in% X.sort] = exp(-Lambda_c1)
        }
        if (Cmod == 'Spline') {
          model = hare(X, Delta_c, cbind(A, Z))
          S_c0 = 1 - sapply(X.sort.full, function(x)
            phare(
              q = x,
              cov = cbind('A' = 0, dfs[, -(1:3)]),
              fit = model
            ))
          S_c0[is.nan(S_c0)] = NA
          S_c1 = 1 - sapply(X.sort.full, function(x)
            phare(
              q = x,
              cov = cbind('A' = 1, dfs[, -(1:3)]),
              fit = model
            ))
          S_c1[is.nan(S_c1)] = NA
        }
        if (Cmod == 'RSF') {
          dfs0 = dfs1 = dfs
          dfs0[, 3] = 0
          dfs1[, 3] = 1
          model = rfsrc(
            Surv(X, Delta_c) ~ .,
            data = dfn.c,
            ntime = 0,
            splitrule = 'bs.gradient',
            ntree = 100,
            nodesize = 50,
            mtry = 2,
            case.wt = weights.n
          )
          S_c0[, which(X.sort.full %in% model$time.interest)] = predict(model, newdata = dfs0[, -(1:2)])$survival
          S_c1[, which(X.sort.full %in% model$time.interest)] = predict(model, newdata = dfs1[, -(1:2)])$survival
        }
        S_c0[, 1] = ifelse(S_c0[, 1] == 0 |
                             is.na(S_c0[, 1]), 1, S_c0[, 1])
        S_c0[S_c0 == 0] = NA
        S_c0 = t(as.matrix(fill(
          data.frame(t(S_c0)), names(data.frame(t(S_c0)))
        )))
        S_c1[, 1] = ifelse(S_c1[, 1] == 0 |
                             is.na(S_c1[, 1]), 1, S_c1[, 1])
        S_c1[S_c1 == 0] = NA
        S_c1 = t(as.matrix(fill(
          data.frame(t(S_c1)), names(data.frame(t(S_c1)))
        )))
        return(list(S_c0, S_c1))
      }

      samp.split_PS = function(dfs, dfn.ps, n, PSmod, weights.n) {
        if (PSmod == 'logit') {
          model = glm(A ~ .,
                      data = dfn.ps[, -(1:2)],
                      family = quasibinomial,
                      weights = weights.n)
          PS = predict(model, newdata = dfs[, -(1:2)], type = 'response')
        }
        if (PSmod == 'RF') {
          model = ranger(
            A ~ .,
            data = dfn.ps[, -(1:2)],
            probability = TRUE,
            num.trees = 2000
          )
          PS = predict(model, dfs[, -(1:2)])$predictions[, 2]
        }
        if (PSmod == 'twang') {
          N.tr = 200
          model = gbm(
            A ~ .,
            data = dfn.ps[, -(1:2)],
            distribution = 'bernoulli',
            n.trees = N.tr,
            interaction.depth = 1,
            weights = weights.n
          )
          pred.data = data.frame(dfs[, -(1:3)])
          colnames(pred.data) = paste0('Z', 1:(ncol(data) - 3))
          PS = predict(model, pred.data, n.trees = N.tr, type = 'response')
        }
        return(PS)
      }
      if (crossFit) {
        # cross fitting
        for (fold in 1:k) {
          dfs = data[folds == fold, ]
          dfn.t = dfn.c = dfn.ps = data[folds != fold, ]
          n = nrow(dfn.t)
          X = dfn.t[, 1]
          n_res = length(unique(X))
          Delta = dfn.t[, 2]
          Delta_c = (1 - Delta) * (X < tau)
          dfn.c$Delta = Delta_c
          colnames(dfn.c)[2] = 'Delta_c'
          A = dfn.t[, 3]
          Z = as.matrix(dfn.t[, -(1:3)])
          p = ncol(Z)
          X.sort = sort(unique(X))
          S_t.est = samp.split_S_t(
            dfs,
            dfn.t,
            n,
            X,
            n_res,
            Delta,
            Delta_c,
            A,
            Z,
            p,
            X.sort,
            X.full,
            X.sort.full,
            Tmod,
            weights.n = weights[folds != fold]
          )
          S_c.est = samp.split_S_c(
            dfs,
            dfn.c,
            n,
            X,
            n_res,
            Delta,
            Delta_c,
            A,
            Z,
            p,
            X.sort,
            X.full,
            X.sort.full,
            Cmod,
            weights.n = weights[folds != fold]
          )
          PS.est = samp.split_PS(dfs, dfn.ps, n, PSmod, weights.n = weights[folds !=
                                                                              fold])
          S_t.temp0[folds == fold, ] = S_t.est[[1]]
          S_t.temp1[folds == fold, ] = S_t.est[[2]]
          S_c.temp0[folds == fold, ] = S_c.est[[1]]
          S_c.temp1[folds == fold, ] = S_c.est[[2]]
          PS.temp[folds == fold] = PS.est
        }
      }

      # define all parameters and covariates
      n = nrow(data)
      X = data[, 1]
      n_res = length(unique(X))
      Delta = data[, 2]
      Delta_c = (1 - Delta) * (X < tau)
      A = data[, 3]
      Z = as.matrix(data[, -(1:3)])
      p = ncol(Z)
      X.sort = sort(unique(X))
      event_rank = sapply(X, function(x)
        which(x == X.sort))

      if (!crossFit) {
        # obtain PS, S_t and S_c without cross-fitting.
        #S_t
        if (Tmod == 'Cox') {
          S_t.temp0 = S_t.temp1 = data.frame(matrix(0, nrow = nrow(data), ncol = length(unique(X.full))))
          model = coxph(
            Surv(X, Delta) ~ .,
            timefix = F,
            data = data,
            weights = weights
          )
          Lambda0 = exp(as.matrix(cbind(0, data[, -(1:3)])) %*% model$coefficients) %*% t(basehaz(model, centered = F)$hazard)
          S_t.temp0[, X.sort.full %in% X.sort] = exp(-Lambda0)
          Lambda1 = exp(as.matrix(cbind(1, data[, -(1:3)])) %*% model$coefficients) %*% t(basehaz(model, centered = F)$hazard)
          S_t.temp1[, X.sort.full %in% X.sort] = exp(-Lambda1)
        }
        # S_c
        if (Cmod == 'Cox') {
          data$Delta = Delta_c
          colnames(data)[2] = 'Delta_c'
          S_c.temp0 = S_c.temp1 = data.frame(matrix(0, nrow = nrow(data), ncol = length(unique(X.full))))
          model = coxph(
            Surv(X, Delta_c) ~ .,
            timefix = F,
            data = data,
            weights = weights
          )
          Lambda_c0 = exp(as.matrix(cbind(0, data[, -(1:3)])) %*% model$coefficients) %*% t(basehaz(model, centered = F)$hazard)
          S_c.temp0[, X.sort.full %in% X.sort] = exp(-Lambda_c0)
          Lambda_c1 = exp(as.matrix(cbind(1, data[, -(1:3)])) %*% model$coefficients) %*% t(basehaz(model, centered = F)$hazard)
          S_c.temp1[, X.sort.full %in% X.sort] = exp(-Lambda_c1)
        }

        # PS
        if (PSmod == 'logit') {
          model = glm(A ~ .,
                      data = data[, -(1:2)],
                      family = quasibinomial,
                      weights = weights)
          PS.temp = predict(model, newdata = data[, -(1:2)], type = 'response')
        }
      }

      # weight trimming
      S_t0 = as.matrix(S_t.temp0)
      S_t1 = as.matrix(S_t.temp1)
      S_c0 = as.matrix(S_c.temp0)
      S_c1 = as.matrix(S_c.temp1)
      PS = as.vector(PS.temp)

      S_t0[S_t0 < min.S] = min.S
      S_t1[S_t1 < min.S] = min.S
      S_c0[S_c0 < min.S] = min.S
      S_c1[S_c1 < min.S] = min.S
      PS[PS < min.PS] = min.PS
      PS[PS > (1 - min.PS)] = 1 - min.PS

      if (augmentation == 'AIPTW') {
        S_c0[S_c0 >= 0] = 1
        S_c1[S_c1 >= 0] = 1
      }

      # product weight trimming
      S_t0c0 = S_t0 * S_c0
      S_t1c1 = S_t1 * S_c1
      S_t0c0PS0 = (1 - PS) * S_t0c0
      S_t1c1PS1 = PS * S_t1c1
      S_cPS = A * PS * S_c1 + (1 - A) * (1 - PS) * S_c0

      #### beta estimation

      Y = t(sapply(event_rank, function(x)
        c(rep(1, x), rep(0, n_res - x))))   # n x n_res
      dNt = t(sapply(event_rank, function(x)
        tabulate(x, nbins = n_res))) * Delta # n x n_res

      # fixed quantities
      Lambda_c0 = -log(S_c0) # n x n_res  (sample x time)
      lambda_c0 = t(apply(Lambda_c0, 1, function(x)
        c(x[1], diff(x))))  # n x n_res
      Lambda_c1 = -log(S_c1) # n x n_res  (sample x time)
      lambda_c1 = t(apply(Lambda_c1, 1, function(x)
        c(x[1], diff(x))))  # n x n_res
      dNc = t(sapply(event_rank, function(x)
        tabulate(x, nbins = n_res))) * Delta_c # n x n_res
      PS0J0 = t(apply((dNc - Y * lambda_c0) / S_t0c0PS0, 1, cumsum))   # n x n_res  (sample x time)
      PS1J1 = t(apply((dNc - Y * lambda_c1) / S_t1c1PS1, 1, cumsum))   # n x n_res  (sample x time)
      S_t = A * S_t1 + (1 - A) * S_t0
      dS_t = t(apply(S_t, 1, function(x)
        c(x[1] - 1, diff(x))))
      dS_t0 = t(apply(S_t0, 1, function(x)
        c(x[1] - 1, diff(x))))
      dS_t1 = t(apply(S_t1, 1, function(x)
        c(x[1] - 1, diff(x))))
      w = A / PS + (1 - A) / (1 - PS)

      if (augmentation == 'AIPTW') {
        PS0J0[PS0J0 != 0] = 0
        PS1J1[PS1J1 != 0] = 0
      }

      # fixed intermediate quantities
      calc_N_summand_a = function(a, l)
        a ^ l * (1 + a * A * PS1J1 + (1 - a) * (1 - A) * (PS0J0)) * (a * dS_t1 + (1 -
                                                                                    a) * dS_t0)  # n x n_res
      calc_dcalN = function(l)
        A ^ l * (dNt / S_cPS + w * dS_t) - calc_N_summand_a(0, l) - calc_N_summand_a(1, l)  # n x n_res
      dN1 = calc_dcalN(1) * weights
      dN0 = calc_dcalN(0) * weights

      if (crossFit) {
        # A_bar and U
        calc_Gamma_summand_a = function(beta, a, l)
          a ^ l * (1 + a * A * PS1J1 + (1 - a) * (1 - A) * (PS0J0)) * (a * S_t1 + (1 -
                                                                                     a) * S_t0) * exp(beta * a)  # n x n_res
        calc_Gamma = function(beta, l)
          A ^ l * exp(beta * A) * (Y / S_cPS - S_t * w) + calc_Gamma_summand_a(beta, 0, l) + calc_Gamma_summand_a(beta, 1, l) # n x n_res
        calc_A_bar = function(beta) {
          res = matrix(0, n_res, n)
          S1 = calc_Gamma(beta, 1) * weights
          S0 = calc_Gamma(beta, 0) * weights
          for (fold in 1:k) {
            res[, folds == fold] = colMeans(S1[folds == fold, ]) / colMeans(S0[folds == fold, ])
          }
          return(t(res))
        }  # n x n_res
        calc_U = function(beta)
          sum(dN1 - calc_A_bar(beta) * dN0) / n
        calc_dU = function(beta)
          sum((calc_A_bar(beta) ^ 2 - calc_A_bar(beta)) * dN0) / n
        beta = newtonRaphson(calc_U, beta0, dfun = calc_dU)$root

        # model se
        A_bar = calc_A_bar(beta)  # n x n_res
        calc_dLambda0 = function(beta) {
          res = matrix(0, n_res, n)
          S0 = calc_Gamma(beta, 0) * weights
          for (fold in 1:k) {
            res[, folds == fold] = colMeans(dN0[folds == fold,]) / colMeans(S0[folds == fold, ])
          }
          return(t(res))
        }  #  n x n_res
        dLambda0 = calc_dLambda0(beta)
        K = mean(
          rowSums(
            dN1 - weights * calc_Gamma(beta, 1) * dLambda0 - A_bar * dN0 + A_bar * weights *
              calc_Gamma(beta, 0) * dLambda0
          ) ^ 2
        )
        nu = calc_dU(beta)
        model_se = sqrt(K / nu ^ 2 / n)

        # Lambda0 and survival curves of group 0 and 1
        Lambda0 = cumsum(colMeans(dLambda0[!duplicated(dLambda0), ]))
        survival0 = exp(-Lambda0)
        survival1 = exp(-Lambda0 * exp(beta))
        names(Lambda0) = names(survival0) = names(survival1) = X.sort.full

        # Scaled Residuals
        scaled_res = colMeans(
          dN1 - weights * calc_Gamma(beta, 1) * dLambda0 - A_bar * dN0 + A_bar * weights *
            calc_Gamma(beta, 0) * dLambda0
        ) / colMeans((A_bar - A_bar ^ 2) * dN0)
      } else{
        # A_bar and U
        calc_Gamma_summand_a = function(beta, a, l)
          a ^ l * (1 + a * A * PS1J1 + (1 - a) * (1 - A) * (PS0J0)) * (a * S_t1 + (1 -
                                                                                     a) * S_t0) * exp(beta * a)  # n x n_res
        calc_Gamma = function(beta, l)
          A ^ l * exp(beta * A) * (Y / S_cPS - S_t * w) + calc_Gamma_summand_a(beta, 0, l) + calc_Gamma_summand_a(beta, 1, l) # n x n_res
        calc_A_bar = function(beta)
          colMeans(calc_Gamma(beta, 1) * weights) / colMeans(calc_Gamma(beta, 0) *
                                                               weights) # n_res x 1
        calc_U = function(beta)
          sum(dN1 - t(calc_A_bar(beta) * t(dN0))) / n
        calc_dU = function(beta)
          sum(t((
            calc_A_bar(beta) ^ 2 - calc_A_bar(beta)
          ) * t(dN0))) / n
        beta = newtonRaphson(calc_U, beta0, dfun = calc_dU)$root

        # model se
        A_bar = calc_A_bar(beta)  # n_res x 1
        dLambda0 = colMeans(dN0) / colMeans(calc_Gamma(beta, 0) * weights) # n_res x 1
        K = mean(rowSums(dN1 - t(t(
          weights * calc_Gamma(beta, 1)
        ) * dLambda0) - t(A_bar * t(dN0)) + t(
          A_bar * dLambda0 * t(weights * calc_Gamma(beta, 0))
        )) ^ 2)
        nu = calc_dU(beta)
        model_se = sqrt(K / nu ^ 2 / n)

        # Lambda0 and survival curves of group 0 and 1
        Lambda0 = cumsum(dLambda0)
        survival0 = exp(-Lambda0)
        survival1 = exp(-Lambda0 * exp(beta))
        names(Lambda0) = names(survival0) = names(survival1) = X.sort.full

        # Scaled Residuals
        scaled_res = colMeans(dN1 - t(t(weights * calc_Gamma(beta, 1)) * dLambda0) - t(A_bar *
                                                                                         t(dN0)) + t(A_bar * dLambda0 * t(weights * calc_Gamma(beta, 0)))) / colMeans(t((A_bar - A_bar ^
                                                                                                                                                                           2) * t(dN0)))
      }
      # beta(t) plot, first remove time points with no events
      indx = !is.nan(scaled_res)
      y = scaled_res[indx] - mean(scaled_res[indx]) + beta
      x = X.sort.full[indx]
      return(
        list(
          beta = beta,
          model_se = model_se,
          Lambda0 = Lambda0,
          survival0 = survival0,
          survival1 = survival1,
          beta_t_plot = list(y = y, x = x)
        )
      )
    }

    # Assume randomized study, perform AIPCW
    if (augmentation == 'AIPCW') {
      S_t.temp = S_c.temp = matrix(0, nrow = nrow(data), ncol = length(unique(X.full)))
      # sample split functions
      samp.split_S_t = function(dfs,
                                dfn.t,
                                n,
                                X,
                                n_res,
                                Delta,
                                Delta_c,
                                A,
                                Z,
                                p,
                                X.sort,
                                X.full,
                                X.sort.full,
                                Tmod,
                                weights.n) {
        S_t = data.frame(matrix(0, nrow = nrow(dfs), ncol = length(unique(X.full))))
        if (Tmod == 'Cox') {
          model = coxph(
            Surv(X, Delta) ~ .,
            timefix = F,
            data = dfn.t,
            weights = weights.n
          )
          Lambda = exp(as.matrix(dfs[, -(1:2)]) %*% model$coefficients) %*% t(basehaz(model, centered = F)$hazard)
          S_t[, X.sort.full %in% X.sort] = exp(-Lambda)
        }
        if (Tmod == 'RSF') {
          model = rfsrc(
            Surv(X, Delta) ~ .,
            data = dfn.t,
            ntime = 0,
            splitrule = 'bs.gradient',
            ntree = 100,
            nodesize = 50,
            mtry = 2,
            case.wt = weights.n
          )
          S_t[, which(X.sort.full %in% model$time.interest)] = predict(model, newdata = dfs)$survival
        }
        S_t[, 1] = ifelse(S_t[, 1] == 0 | is.na(S_t[, 1]), 1, S_t[, 1])
        S_t[S_t == 0] = NA
        S_t = t(as.matrix(fill(
          data.frame(t(S_t)), names(data.frame(t(S_t)))
        )))
        return(S_t)
      }

      samp.split_S_c = function(dfs,
                                dfn.c,
                                n,
                                X,
                                n_res,
                                Delta,
                                Delta_c,
                                A,
                                Z,
                                p,
                                X.sort,
                                X.full,
                                X.sort.full,
                                Cmod,
                                weights.n) {
        S_c = data.frame(matrix(0, nrow = nrow(dfs), ncol = length(unique(X.full))))
        if (Cmod == 'Cox') {
          model = coxph(
            Surv(X, Delta_c) ~ .,
            timefix = F,
            data = dfn.c[, -3],
            weights = weights.n
          )
          Lambda_c = exp(as.matrix(dfs[, -(1:3)]) %*% model$coefficients) %*% t(basehaz(model, centered = F)$hazard)
          S_c[, X.sort.full %in% X.sort] = exp(-Lambda_c)
        }
        if (Cmod == 'RSF') {
          model = rfsrc(
            Surv(X, Delta_c) ~ .,
            data = dfn.c,
            ntime = 0,
            splitrule = 'bs.gradient',
            ntree = 100,
            nodesize = 50,
            mtry = 2,
            case.wt = weights.n
          )
          S_c[, which(X.sort.full %in% model$time.interest)] = predict(model, newdata = dfs[, -(1:2)])$survival
        }
        S_c[, 1] = ifelse(S_c[, 1] == 0 | is.na(S_c[, 1]), 1, S_c[, 1])
        S_c[S_c == 0] = NA
        S_c = t(as.matrix(fill(
          data.frame(t(S_c)), names(data.frame(t(S_c)))
        )))
        return(S_c)
      }
      if (crossFit) {
        # cross fitting
        for (fold in 1:k) {
          dfs = data[folds == fold, ]
          dfn.t = dfn.c = data[folds != fold, ]
          n = nrow(dfn.t)
          X = dfn.t[, 1]
          n_res = length(unique(X))
          Delta = dfn.t[, 2]
          Delta_c = (1 - Delta) * (X < tau)
          dfn.c$Delta = Delta_c
          colnames(dfn.c)[2] = 'Delta_c'
          A = dfn.t[, 3]
          Z = as.matrix(dfn.t[, -(1:3)])
          p = ncol(Z)
          X.sort = sort(unique(X))
          S_t.temp[folds == fold, ] = samp.split_S_t(
            dfs,
            dfn.t,
            n,
            X,
            n_res,
            Delta,
            Delta_c,
            A,
            Z,
            p,
            X.sort,
            X.full,
            X.sort.full,
            Tmod,
            weights.n = weights[folds != fold]
          )
          S_c.temp[folds == fold, ] = samp.split_S_c(
            dfs,
            dfn.c,
            n,
            X,
            n_res,
            Delta,
            Delta_c,
            A,
            Z,
            p,
            X.sort,
            X.full,
            X.sort.full,
            Cmod,
            weights.n = weights[folds != fold]
          )
        }
        S_t = as.matrix(S_t.temp)
        S_c = as.matrix(S_c.temp)
      }

      # define all parameters and covariates
      n = nrow(data)
      X = data[, 1]
      n_res = length(unique(X))
      Delta = data[, 2]
      Delta_c = (1 - Delta) * (X < tau)
      A = data[, 3]
      Z = as.matrix(data[, -(1:3)])
      p = ncol(Z)
      X.sort = sort(unique(X))
      event_rank = sapply(X, function(x)
        which(x == X.sort))
      Y = t(sapply(event_rank, function(x)
        c(rep(1, x), rep(0, n_res - x))))   # n x n_res  (sample x time)
      dNt = t(sapply(event_rank, function(x)
        tabulate(x, nbins = n_res))) * Delta    # n x n_res  (sample x time)

      # obtain S_t and S_c without cross-fitting.
      if (!crossFit) {
        #S_t
        S_t = data.frame(matrix(0, nrow = nrow(data), ncol = length(unique(X.full))))
        if (Tmod == 'Cox') {
          model = coxph(
            Surv(X, Delta) ~ .,
            timefix = F,
            data = data,
            weights = weights
          )
          Lambda = exp(as.matrix(data[, -(1:2)]) %*% model$coefficients) %*% t(basehaz(model, centered = F)$hazard)
          S_t[, X.sort.full %in% X.sort] = exp(-Lambda)
        }
        S_t[, 1] = ifelse(S_t[, 1] == 0 | is.na(S_t[, 1]), 1, S_t[, 1])
        S_t[S_t == 0] = NA
        S_t = t(as.matrix(fill(
          data.frame(t(S_t)), names(data.frame(t(S_t)))
        )))

        # S_c
        data$Delta = Delta_c
        colnames(data)[2] = 'Delta_c'
        S_c = data.frame(matrix(0, nrow = nrow(data), ncol = length(unique(X.full))))
        if (Cmod == 'Cox') {
          model = coxph(
            Surv(X, Delta_c) ~ .,
            timefix = F,
            data = data,
            weights = weights
          )
          Lambda_c = exp(as.matrix(data[, -(1:2)]) %*% model$coefficients) %*% t(basehaz(model, centered = F)$hazard)
          S_c[, X.sort.full %in% X.sort] = exp(-Lambda_c)
        }
        S_c[, 1] = ifelse(S_c[, 1] == 0 | is.na(S_c[, 1]), 1, S_c[, 1])
        S_c[S_c == 0] = NA
        S_c = t(as.matrix(fill(
          data.frame(t(S_c)), names(data.frame(t(S_c)))
        )))

        S_t = as.matrix(S_t)
        S_c = as.matrix(S_c)
      }
      # Trim weights
      S_t[S_t < min.S] = min.S
      S_c[S_c < min.S] = min.S

      # beta estimation
      Lambda_c = -log(S_c) # n x n_res  (sample x time)
      lambda_c = t(apply(Lambda_c, 1, function(x)
        c(x[1], diff(x))))    # n x n_res  (sample x time)
      # calculate fixed quantities
      dNc = t(sapply(event_rank, function(x)
        tabulate(x, nbins = n_res))) * Delta_c # n x n_res  (sample x time)
      dMc = dNc - Y * lambda_c   # n x n_res  (sample x time)
      J = t(apply(dMc / S_t / S_c, 1, cumsum))   # n x n_res  (sample x time)
      dS_t = t(apply(S_t, 1, function(x)
        c(x[1] - 1, diff(x))))
      dcalN = (dNt / S_c - J * dS_t) * weights  # n x n_res  (sample x time)

      # function for calculating S_l
      if (crossFit) {
        calc_Si = function(beta, l = 0)
          A ^ l * exp(beta * A) * (Y / S_c + J * S_t) * weights # n x n_res
        calc_A_bar = function(beta) {
          res = matrix(0, n_res, n)
          S1 = calc_Si(beta, 1)
          S0 = calc_Si(beta, 0)
          for (fold in 1:k) {
            res[, folds == fold] = colMeans(S1[folds == fold, ]) / colMeans(S0[folds == fold, ])
          }
          return(t(res))
        }  # n x n_res
        calc_U = function(beta)
          sum((A - calc_A_bar(beta)) * dcalN) / n    #  1x1
        calc_dU = function(beta)
          sum((calc_A_bar(beta) ^ 2 - calc_A_bar(beta)) * dcalN) / n   #  1x1
        beta = newtonRaphson(calc_U, beta0, dfun = calc_dU)$root

        # model_se
        A_bar = calc_A_bar(beta)  # n x n_res
        calc_dLambda0 = function(beta) {
          res = matrix(0, n_res, n)
          S0 = calc_Si(beta, 0)
          for (fold in 1:k) {
            res[, folds == fold] = colMeans(dcalN[folds == fold,]) / colMeans(S0[folds == fold, ])
          }
          return(t(res))
        }  #  n x n_res
        dLambda0 = calc_dLambda0(beta)
        K = mean(rowSums((A - calc_A_bar(beta)) * (dcalN - calc_Si(beta, 0) *
                                                     dLambda0)) ^ 2)
        nu = calc_dU(beta)
        model_se = sqrt(K / nu ^ 2 / n)

        # Lambda0 and survival curves of group 0 and 1
        Lambda0 = cumsum(colMeans(dLambda0[!duplicated(dLambda0), ]))
        survival0 = exp(-Lambda0)
        survival1 = exp(-Lambda0 * exp(beta))
        names(Lambda0) = names(survival0) = names(survival1) = X.sort.full

        # scaled residual
        scaled_res = colMeans((A - calc_A_bar(beta)) * (dcalN - calc_Si(beta, 0) *
                                                          dLambda0)) / colMeans((A_bar - A_bar ^ 2) * dcalN)
      } else{
        calc_S = function(beta, l = 0)
          colMeans(A ^ l * exp(beta * A) * (Y / S_c + J * S_t) * weights) # n_res x 1   (time x 1)
        calc_A_bar = function(beta)
          (calc_S(beta, 1) / calc_S(beta, 0)) # n_res x 1   (time x 1)
        calc_U = function(beta)
          sum((dcalN * outer(A, calc_A_bar(beta), '-'))) / n   #  1x1
        calc_dU = function(beta)
          - sum((t(dcalN) * (
            calc_A_bar(beta) - calc_A_bar(beta) ^ 2
          ))) / n   #  1x1
        beta = newtonRaphson(calc_U, beta0, dfun = calc_dU)$root

        # model se
        A_bar = calc_A_bar(beta)  # n_res x 1
        K = mean(rowSums((dcalN - t(
          t(exp(beta * A) * (Y / S_c + J * S_t)) / calc_S(beta) * colMeans(dcalN)
        )) * outer(A, A_bar, '-')) ^ 2)
        nu = calc_dU(beta)
        model_se = sqrt(K / nu ^ 2 / n)

        # Lambda0 and survival curves of group 0 and 1
        Lambda0 = cumsum(colMeans(dcalN) / calc_S(beta, 0))
        survival0 = exp(-Lambda0)
        survival1 = exp(-Lambda0 * exp(beta))
        names(Lambda0) = names(survival0) = names(survival1) = X.sort.full

        # scaled residual
        scaled_res = colMeans((dcalN - t(
          t(exp(beta * A) * (Y / S_c + J * S_t)) / calc_S(beta) * colMeans(dcalN)
        )) * outer(A, A_bar, '-')) / rowMeans(t(dcalN) * (calc_A_bar(beta) - calc_A_bar(beta) ^
                                                            2))
      }
      # beta(t) plot, remove time points with no events
      indx = !is.nan(scaled_res)
      y = scaled_res[indx] - mean(scaled_res[indx]) + beta
      x = X.sort.full[indx]
      return(
        list(
          beta = beta,
          model_se = model_se,
          Lambda0 = Lambda0,
          survival0 = survival0,
          survival1 = survival1,
          beta_t_plot = list(y = y, x = x)
        )
      )
    }
  }

