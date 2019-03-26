library("devtools")
getwd()
setwd("./fdapoi")

devtools::document() ## makes the documentation using roxygen

## Checking
setwd("..")
# devtools::check("fdapoi")

remove.packages("fdapoi")
devtools::install("fdapoi")


# library("fdapoi")
# help(package = "fdapoi")

# BM.SIM(0,1,10)
# OU.SIM(0,1,10,3)
# 
# ?FUN_PoI_BIC
# ?is.error
