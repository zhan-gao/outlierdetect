#' Iterative Hard-thresholding (IHT) (Algorithm 1)
#'
#' @param x
#' @param y
#' @param k
#' @param intercept
#' @param b0 Null hypothesis
#' @param b_init
#' @param tol Convergence criteria
#' @param MaxIter
#'
#' @return A list contains estimation results
#'   \item{b} Estimated coefficient
#'   \item{a} Estimated fixed effects
#'   \item{d} = 0 if detected outlier; = 1 if not.
#'
#'
#' @export


iht <-
    function(x,
             y,
             k,
             intercept = TRUE,
             b0 = NULL,
             b_init = NULL,
             tol = 1e-6,
             MaxIter = 2000) {

        # Projected block-coordinate descent

        if (is.null(b_init)) {
            b_init <- lad_reg(x, y, intercept)
        }

        if (intercept) {
            x <- cbind(1, x)
        }

        e_init <- y - x %*% b_init
        alpha_init <- H(e_init, k)
        ind_init <- (alpha_init == 0)

        iter <- 1
        diff <- 100
        while (diff > tol && iter <= MaxIter) {
            b_update <-
                lsfit(x[ind_init,], y[ind_init], intercept = FALSE)$coef
            e_update <- y - x %*% b_update
            alpha_update <- H(e_update, k)
            ind_update <- (alpha_update == 0)

            diff <- max(abs(ind_update - ind_init))
            ind_init <- ind_update
            iter <- iter + 1
        }
        b <- b_update
        a <- alpha_update
        d <- as.numeric(ind_update)

        # t_stat <- post_inference(x, y, b, d, b0)

        return(list(b = b,
                    a = a,
                    d = d))
    }

# ---
#’ Hard-threshold operator
H <- function(d, k){

    # Get the indices of largest k elements of c
    ind <- order(abs(d), decreasing = TRUE)
    # set coefficients not top k to be 0
    d[-ind[1:k]] <- 0

    return(d)
}

# ---

#' Local Combinatorial Search (Algorithm 2)
#'
#' @param x
#' @param y
#' @param k
#' @param l local exactness level
#' @param intercept
#' @param b_init
#' @param tol
#' @param MaxIter
#' @param bd_const
#' @param time_limit
#' @param suppress_msg
#' @param tol_iht
#' @param MaxIter_iht
#' @param int_tol tolerance for integer in gurobi
#' @param MIPGap termination condition for gurobi
#' @param MIPFocus strategies for gurobi
#'
#' @export
#'
#'

local_search <- function(x,
                         y,
                         k,
                         l,
                         intercept = TRUE,
                         b_init = NULL,
                         tol = 1e-6,
                         MaxIter = 100,
                         bd_const = 2,
                         time_limit = 300,
                         suppress_msg = TRUE,
                         tol_iht = 1e-7,
                         MaxIter_iht = 2000,
                         int_tol = 1e-8,
                         MIPGap  = 0.001,
                         MIPFocus = 3) {


    if (is.null(b_init)) {
        b_init <- lad_reg(x, y, intercept)
    }

    if (l == 0) {
        return(iht(x, y, k,
                   intercept = intercept,
                   b_init = b_init,
                   tol = tol_iht,
                   MaxIter = MaxIter_iht))
    } else if ((l < 0) || (l > k)) {
        stop("Incorrect input of local exactness level l.")
    }

    fit_j <- NULL
    for (j in 1:MaxIter) {
        iht_fit <- iht(
            x,
            y,
            k,
            intercept = intercept,
            b_init = b_init,
            tol = tol_iht,
            MaxIter = MaxIter_iht
        )

        start_value <- c(iht_fit$b, iht_fit$a, iht_fit$d)
        mio_fit <-
            local_search_mio(x, y, k, l, intercept, start_value,
                             bd_const, time_limit, suppress_msg,
                             int_tol, MIPGap, MIPFocus)

        dist <- sum(abs(mio_fit$d - iht_fit$d))

        if (dist < tol) {
            mio_fit$b <- iht_fit$b
            mio_fit$a <- iht_fit$a
            fit_j <- mio_fit
            break
        }
        else {
            fit_j <- mio_fit
            b_init <- mio_fit$b
        }
    }

    return(fit_j)
}


#' local combinatorial search MIO problem
#'
#' @param x
#' @param y
#' @param k
#' @param l
#' @param intercept
#' @param start_value (b, a, gamma) solution to be checked
#' @param bd_const parameter for bounds of variables for Gurobi solver
#' @param time_limit
#' @param suppress_msg
#' @param int_tol tolerance for integer in gurobi
#' @param MIPGap termination condition for gurobi
#' @param MIPFocus strategies for gurobi
#'
#'
#' @export

local_search_mio <- function(x,
                             y,
                             k,
                             l,
                             intercept = TRUE,
                             start_value,
                             bd_const = 2,
                             time_limit = 300,
                             suppress_msg = TRUE,
                             int_tol = 1e-8,
                             MIPGap  = 0.001,
                             MIPFocus = 3) {

    if (intercept) {
        x <- cbind(1, x)
    }

    # x: n-by-m matrix
    # y: n-by-1 matrix
    # dimension
    n <- nrow(x)
    m <- ncol(x)

    # construct the optimization problem
    model <- list()
    model$modelsense <- "min"
    model$start <- start_value

    x_tilde <- cbind(x, diag(n))
    # objective
    model$Q <- rbind(cbind(t(x_tilde) %*% x_tilde, Matrix::spMatrix(n + m, n)),
                     Matrix::spMatrix(n, m + 2 * n))
    model$obj <- c(-2 * t(x_tilde) %*% y, rep(0, n))

    # SOS-1 constraints - A list of lists
    sos <- list()
    for (i in 1:n) {
        sos[[i]] <- list(type  = 1,
                         index = c(m + i, m + n + i),
                         weight = c(1, 2))
    }
    model$sos <- sos

    # Linear constraints
    in_ind <- start_value[-(1:(m + n))]
    out_ind <- 1 - in_ind
    model$A <- rbind(c(rep(0, m + n ), rep(1, n)),
                     c(rep(0, m + n ), in_ind),
                     c(rep(0, m + n ), out_ind))
    model$rhs <- c(n - k, n - k - l, l)
    model$sense <- c(">=", ">=", "<=")

    # variable type and bounds
    model$vtype <- c(rep("C", m + n), rep("B", n))

    # bounds for coefficients
    bd <- bd_const * max(abs(y))
    bd_beta <- bd_const * abs(start_value[1:m])
    model$lb <- c(-bd_beta, rep(-bd, n), rep(0, n))
    model$ub <- c(bd_beta, rep(bd, n), rep(1, n))

    # solve the model
    params <- list()
    params$NonConvex <- 2
    params$TimeLimit <- time_limit
    params$IntFeasTol <- int_tol
    params$MIPGap <- MIPGap
    params$MIPFocus <- MIPFocus

    if (suppress_msg) params$OutputFlag <- 0
    gurobi_out <- gurobi::gurobi(model, params = params)


    b <-  gurobi_out$x[1:m]
    a <- gurobi_out$x[(m + 1):(m + n)]
    d <- gurobi_out$x[(m + n + 1):(m + 2 * n)]

    return(list(b = b,
                a = a,
                d = d,
                sol_status = gurobi_out$status,
                objval = gurobi_out$objval,
                mipgap = gurobi_out$mipgap))
}

# ---
#' Robust Linear Regression with L0 regularization via MIO
#'
#' @param x
#' @param y
#' @param k
#' @param b0 Null hypothesis
#' @param bd_const parameter for bounds of variables for Gurobi solver
#' @param time_limit time limit for Gurobi solver
#' @param start_value start value for Gurobi solver; by default using IHT
#' @param suppress_msg
#' @param tol_iht
#' @param MaxIter_iht
#' @param int_tol tolerance for integer in gurobi
#' @param MIPGap termination condition for gurobi
#' @param MIPFocus strategies for gurobi
#'
#'
#' @return A list contains estimation results
#'   \item{b} Estimated coefficient
#'   \item{a} Estimated fixed effects
#'   \item{d} = 0 if detected outlier; = 1 if not.
#'   \item{t}
#'   \item{sol_status}
#'
#' @import gurobi
#'
#' @export

lsfit_robust_mio <- function(x,
                             y,
                             k,
                             intercept = TRUE,
                             b0 = NULL,
                             bd_const = 2,
                             time_limit = 300,
                             start_value = NULL,
                             suppress_msg = TRUE,
                             tol_iht = 1e-6,
                             MaxIter_iht = 2000,
                             int_tol = 1e-8,
                             MIPGap  = 0.001,
                             MIPFocus = 3) {
    # construct the optimization problem
    model <- list()
    model$modelsense <- "min"

    # Warm start
    if (is.null(start_value)) {
        iht_fit <- iht(x,
                       y,
                       k,
                       intercept = intercept,
                       tol = tol_iht,
                       MaxIter = MaxIter_iht)
        start_value <- c(iht_fit$b, iht_fit$a, iht_fit$d)
    }
    model$start <- start_value

    if (intercept) {
        x <- cbind(1, x)
    }

    # x: n-by-m matrix
    # y: n-by-1 matrix
    # dimension
    n <- nrow(x)
    m <- ncol(x)

    x_tilde <- cbind(x, diag(n))

    # objective
    model$Q <- rbind(cbind(t(x_tilde) %*% x_tilde, Matrix::spMatrix(n + m, n)),
                     Matrix::spMatrix(n, m + 2 * n))
    model$obj <- c(-2 * t(x_tilde) %*% y, rep(0, n))
    # SOS-1 constraints - A list of lists
    sos <- list()
    for (i in 1:n) {
        sos[[i]] <- list(type  = 1,
                         index = c(m + i, m + n + i),
                         weight = c(1, 2))
    }
    model$sos <- sos
    # Linear constraints
    model$A <- rbind(c(rep(0, m + n), rep(1, n)))
    model$rhs <- c(n - k)
    model$sense <- c(">=")
    # variable type and bounds
    model$vtype <- c(rep("C", m + n), rep("B", n))
    # bounds for coefficients
    bd_beta <- bd_const * abs(start_value[1:m])
    bd = bd_const * max(abs(y))
    model$lb <- c(-bd_beta, rep(-bd, n), rep(0, n))
    model$ub <- c(bd_beta, rep(bd, n), rep(1, n))
    # solve the model
    params <- list()
    params$NonConvex <- 2
    params$TimeLimit <- time_limit
    if (suppress_msg) params$OutputFlag <- 0
    params$IntFeasTol <- int_tol
    params$MIPGap <- MIPGap
    params$MIPFocus <- MIPFocus

    gurobi_out <- gurobi::gurobi(model, params = params)


    b <-  gurobi_out$x[1:m]
    a <- gurobi_out$x[(m + 1):(m + n)]
    d <- gurobi_out$x[(m + n + 1):(m + 2 * n)]

    # t_stat <- post_inference(x, y, b, d, b0)

    return(list(b = b,
                a = a,
                d = d,
                sol_status = gurobi_out$status,
                objval = gurobi_out$objval,
                mipgap = gurobi_out$mipgap))
}

#' Huber regression
#'
#' @param x
#' @param y
#' @param intercept
#' @param huber_thres
#'
#' @return estimated coefficient
#'
#' @import CVXR
#'
#' @export

huber_reg_cvxr <- function(x, y, intercept = TRUE, huber_thres = 1.345) {
    # huber_thres = 1.345 suggested by Huber (1981)

    if (intercept) {
        x <- cbind(1, x)
    }

    m <- ncol(x)
    beta <- Variable(m)
    obj <- sum(CVXR::huber(y - x %*% beta, huber_thres))
    problem <- Problem(Minimize(obj))
    prob_data <- get_problem_data(problem, solver = "ECOS")
    ECOS_dims <- ECOS.dims_to_solver_dict(prob_data$data[["dims"]])
    solver_output <- ECOSolveR::ECOS_csolve(
        c = prob_data$data[["c"]],
        G = prob_data$data[["G"]],
        h = prob_data$data[["h"]],
        dims = ECOS_dims,
        A = prob_data$data[["A"]],
        b = prob_data$data[["b"]]
    )
    result <- unpack_results(
        problem, solver_output,
        prob_data$chain, prob_data$inverse_data
    )

    b_hat <- c(result$getValue(beta))
    e_hat <- y - x %*% b_hat
    d_hat <- as.numeric(abs(e_hat) < huber_thres)

    a_hat <- rep(NA, length(y))
    a_hat[abs(e_hat) < huber_thres] <- 0
    a_hat[abs(e_hat) >= huber_thres] <- e_hat[abs(e_hat) >= huber_thres]
    # a_hat[e_hat >= huber_thres] <- e_hat[e_hat >= huber_thres] - huber_thres
    # a_hat[e_hat <= -huber_thres] <- e_hat[e_hat <= -huber_thres] + huber_thres

    return(list(b = b_hat,
                a = a_hat,
                d = d_hat))
}

#' Huber regression based on Yi and Huang (2016)
#'
#' Yi, C. and Huang, J. (2016) Semismooth Newton Coordinate Descent Algorithm
#' for Elastic-Net Penalized Huber Loss Regression and Quantile Regression
#'
#'
#' @param x
#' @param y
#' @param intercept
#' @param huber_thres
#'
#' @return estimated coefficient
#'
#' @import hqreg
#'
#' @export
#'
huber_reg <- function(x, y, intercept = TRUE, huber_thres = 1.345) {

    n <- nrow(x)

    if (intercept) {
        huber_fit <- hqreg::hqreg(x, y, lambda = c(1, 0), gamma = huber_thres)
        b_hat <- huber_fit$beta[, 2]
        e_hat <- y - cbind(1, x) %*% b_hat
        d_hat <- as.numeric(abs(e_hat) < huber_thres)

        a_hat <- rep(NA, n)
        a_hat[abs(e_hat) < huber_thres] <- 0
        a_hat[abs(e_hat) >= huber_thres] <- e_hat[abs(e_hat) >= huber_thres]
    } else {
        return(huber_reg_cvxr(x, y, huber_thres = huber_thres))
    }

    return(list(b = b_hat,
                a = a_hat,
                d = d_hat))

}

#' LAD
#'
#' @param x
#' @param y
#' @param intercept
#' @param tau quantile
#'
#' @return estimated coefficient
#'
#' @export
#'
lad_reg <- function(x, y, intercept = TRUE, tau = 0.5) {
    if (intercept) {
        x <- cbind(1, x)
    }
    quantreg::rq.fit(x, y, tau)$coefficients
}



#' Information Criteria
#'
#' @param x
#' @param y
#' @param b
#' @param a
#' @param ic_ind

#' @return The computed information criterion
#'
#' @export
#'

ic <- function(x, y, b, a, ic_ind = 3) {

    if (length(b) > ncol(x)) {
        x <- cbind(1, x)
    }

    n <- nrow(x)
    p <- ncol(x)
    k <- sum(a != 0)

    rss <- sum((y - x %*% b - a)^2)
    ic1 <- n * log(rss / (n - k)) + k * log(n)
    ic2 <- n * log(rss / (n - k)) + (k + p) * (log(n))
    ic3 <- (n - k) * log(rss / (n)) + k * log(n)
    ic4 <- (n - k) * log(rss / (n)) + (k + p) * log(n)
    ic5 <- (n - k) * log(rss / (n - k)) + k * log(n)
    ic6 <- (n - k) * log(rss / (n - k)) + (k + p) * log(n)

    ic_vec <- c(ic1, ic2, ic3, ic4, ic5, ic6)

    return(ic_vec[ic_ind])
}

#' Choose number of outliers k
#'
#' @param x
#' @param y
#' @param l if l = 0, use IHT without local search.
#' @param k_vec candidate k values
#' @param intercept
#' @param ic_ind which IC to be used.
#' @param b_init
#' @param tol_ls
#' @param MaxIter_ls
#' @param bd_const
#' @param time_limit
#' @param suppress_msg
#' @param tol_iht
#' @param MaxIter_iht
#' @param int_tol tolerance for integer in gurobi
#' @param MIPGap termination condition for gurobi
#' @param MIPFocus strategies for gurobi
#'
#'
#'
#' @return selected k (scalar or vector, len = length(ic_ind))
#'
#' @export
#'
choose_k <- function(x, y,  k_vec, l = 0,
                     intercept = TRUE,
                     ic_ind = 3,
                     b_init = NULL,
                     tol_ls = 1e-6,
                     MaxIter_ls = 100,
                     bd_const = 2,
                     time_limit = 600,
                     suppress_msg = TRUE,
                     tol_iht = 1e-6,
                     MaxIter_iht = 2000,
                     int_tol = 1e-8,
                     MIPGap  = 0.001,
                     MIPFocus = 3) {


    if (is.null(b_init)) {
        b_init <- lad_reg(x, y, intercept)
    }

    num_k <- length(k_vec)
    num_ic <- length(ic_ind)

    ic_mat <- matrix(0, num_k, num_ic)
    for (i in 1:num_k) {
        k <- k_vec[i]

        if (l == 0) {
            fit <- iht(x, y, k, intercept, b_init = b_init,
                       tol = tol_iht, MaxIter = MaxIter_iht)
        } else {
            fit <- local_search(x, y, k, l, intercept, b_init,
                                tol = tol_ls,
                                MaxIter = MaxIter_ls,
                                bd_const = bd_const,
                                time_limit = time_limit,
                                suppress_msg = suppress_msg,
                                tol_iht = tol_iht,
                                MaxIter_iht = MaxIter_iht,
                                int_tol = int_tol,
                                MIPGap  = MIPGap,
                                MIPFocus = MIPFocus)
        }

        ic_mat[i, ] <- ic(x, y, fit$b, fit$a, ic_ind)
    }

    k_out <- k_vec[apply(ic_mat, 2, which.min)]

    return(k_out)
}


#' Choose tuning parameters for L1 method
#'
#' @param x
#' @param y
#' @param intercept
#' @param lambda_seq candidate k values
#' @param ic_ind which IC to be used.
#' @param nlambda
#' @param lambda_min_ratio
#'
#' @import doParallel
#' @import foreach
#' @import parallel
#'
#' @export
#'
choose_lambda <- function(x, y, intercept = TRUE,
                          lambda_seq = NULL, ic_ind = 3, parallel = FALSE,
                          nlambda = 100, lambda_min_ratio = 0.0001) {

    if (is.null(lambda_seq)) {
        lambda_seq <- get_lambda_seq(x, y, intercept, nlambda, lambda_min_ratio)
    }

    num_lambda <- length(lambda_seq)
    num_ic <- length(ic_ind)

    if (parallel) {

        if(!foreach::getDoParRegistered()) {
            num_cores <- parallel::detectCores()
            cltr <- parallel::makeCluster(num_cores)
            doParallel::registerDoParallel(cltr, cores = num_cores)
        }

        ic_mat <- foreach(i = 1:num_lambda, .combine = rbind,
                          .packages = c("hqreg"),
                          .export = c("huber_reg",
                                      "huber_reg_cvxr")) %dopar% {
                              lambda <- lambda_seq[i]
                              huber_fit <- huber_reg(x, y, intercept, lambda)
                              ic(x, y, huber_fit$b, huber_fit$a, ic_ind)
                          }
    } else {
        ic_mat <- matrix(0, num_lambda, num_ic)
        for (i in 1:num_lambda) {
            lambda <- lambda_seq[i]
            huber_fit <- huber_reg(x, y, intercept, lambda)
            ic_mat[i, ] <- ic(x, y, huber_fit$b, huber_fit$a, ic_ind)
        }
    }

    lambda_out <- lambda_seq[apply(ic_mat, 2, which.min)]

    return(lambda_out)
}


#' The function generate candidate tuning parameters for L1 method
#'
#'
#'
get_lambda_seq <- function(x, y, intercept,
                           nlambda = 100, lambda_min_ratio = 0.0001) {
    if (intercept) {
        x <- cbind(1, x)
    }

    b_ols <- solve(t(x) %*% x, t(x) %*% y)
    lambda_max <- max(abs(y - x %*% b_ols))
    lambda_min <- lambda_min_ratio * lambda_max
    lambda_seq <- exp(
        seq(log(lambda_max), log(lambda_min), length.out = nlambda)
    )

    return(lambda_seq)

}
