# outlierdetect

This package is the accompanying package for [Gao, Z. and Moon H.R. (2024) Robust Estimation of Regression Models with Potentially Endogenous Outliers via a Modern Optimization Lens](https://arxiv.org/abs/2302.14380). It implements the local combinatorial search algorithm for robust estimation method with $L_0$-regularization on the case-specific parameters, as detailed in the paper.


## Installation
The package can be installed from the Github repository.
```R
install.packages("devtools")
devtools::install_github("zhan-gao/outlierdetect")
```

## Demo
```R
data("test_data")
# sample size n = 200, number of outliers = 20
x <- test_data$x
y <- test_data$y

k_hat <- choose_k(x, y, k_vec = 1:40, l = 1)
print(k_hat) # Selected k
lcs_out <- local_search(x, y, k_hat, l = 2)
print(lcs_out$b) # estimated beta
print(lcs_out$d) # detected inliers
```

```
> print(k_hat)
[1] 19
> print(lcs_out$b)
                 x1           
0.5613816 1.2078783 0.8584481
> print(lcs_out$d)
  [1] 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1
 [83] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
[165] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
```

