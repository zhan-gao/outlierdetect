library(readxl)
library(zoo)
library(lubridate)
library(tidyverse)

# This script processes the Welch and Goyal (2008) data for return forecasting

m <- 1

# Load the raw data
Monthly.Raw = read_excel("./data_raw/PredictorData2023.xlsx",
                            sheet = "Monthly",
                            col_names = TRUE,
                            na = "NaN")

# Construct the data frame, start from the dates
Monthly.Data = data.frame(date = as.yearmon( ymd(Monthly.Raw$yyyymm, truncated = 2) ))

attach(Monthly.Raw)

# Monthly.Data$ExReturn = log(Index / c(1,Index[-length(Index)])) - log(1 + tbl/12)
# Since tbl is included in the predictors (on the right hand side),
# whether include it in the construction of ExReturn is not important.
# Monthly.Data$ExReturn = log(Index + D12) - log(c(1,Index[-length(Index)]) + c(1,D12[-length(D12)]))
Monthly.Data$ExReturn = log(Index) - log(c(1,Index[-length(Index)]))

# Long-horizon return
T = nrow(Monthly.Data)
LongReturn = rep(NA, T)
for(t in 1:(T-m+1)){
    LongReturn[t] =  sum( Monthly.Data$ExReturn[t:(t+m-1)] )
}

Monthly.Data$LongReturn = LongReturn

Monthly.Data$dp = log(D12) - log(Index)
Monthly.Data$dy = log(D12) - log( c(1, Index[-length(Index)]) )
Monthly.Data$ep = log(E12) - log(Index)
Monthly.Data$tms = lty - tbl
Monthly.Data$dfy = BAA - AAA
Monthly.Data$dfr = corpr - ltr

detach(Monthly.Raw)

Monthly.Data = cbind( Monthly.Data,
                        Monthly.Raw[c("b/m",
                                    "tbl",
                                    "ltr",
                                    # "lty",
                                    "ntis",
                                    "svar",
                                    "infl")] )

Monthly.Data = Monthly.Data[-1, ]

welch_goyal_data <- na.omit(Monthly.Data)
usethis::use_data(welch_goyal_data, overwrite = TRUE)



