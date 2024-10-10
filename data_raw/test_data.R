dgp <- function(n, b, a, p,
                outlier_fe_type = 1, fe_para = c(0, 3), chisq_df = 2) {

    # number of regressors
    m <- length(b)

    # Innovations
    innov_mat_1 <- matrix(rnorm(n * (chisq_df)), n, chisq_df)
    innov_mat_2 <- matrix(rnorm(n * (m - 1)), n, m - 1)

    # Regressors
    x1 <- (rowSums(innov_mat_1^2) - chisq_df) / sqrt(2 * chisq_df)
    x_mat <- t(apply(cbind(x1, innov_mat_2), 1, cumsum))

    # outlier fixed effects
    if (outlier_fe_type == 1) {
        if(length(fe_para) != 2) stop("Wrong fe_para input.")
        mu <- fe_para[1]
        s <- fe_para[2]
        fe <- rnorm(n, mu, s)
    } else if (outlier_fe_type == 2) {
        if(length(fe_para) != 1) stop("Wrong fe_para input.")
        fe <- fe_para * rowSums(cbind(innov_mat_1, innov_mat_2))
    } else if (outlier_fe_type == 3) {
        if(length(fe_para) != 2) stop("Wrong fe_para input.")
        mu <- fe_para[1]
        s <- fe_para[2]
        delta <- rnorm(n, mu, s)
        fe <- delta * x_mat[, m]
    } else {
        stop("Wrong outlier_fe_type input.")
    }

    # Outlier indexes
    if (p > 1 || p < 0) {
        stop("Fraction of Outliers must lie in the (0,1) interval!")
    }
    k0 <- floor(p * n)
    outlier_ind <- c(rep(1, k0), rep(0, n - k0))

    # generate y
    u <- rnorm(n)
    y <- a + x_mat %*% b + (fe * outlier_ind) + u

    return(
        list(x = x_mat,
             y = y,
             k0 = k0)
    )
}


n <- 200
b <- c(1, 1)
a <- 0.5
p <- 0.1
outlier_fe_type <- 3
fe_para <- c(10, 10)
set.seed(100)
data_dgp <- dgp(n, b, a, p, outlier_fe_type, fe_para)

x <- data_dgp$x
y <- data_dgp$y
k <- data_dgp$k

test_data <- list(x = x, y = y)
usethis::use_data(test_data)


