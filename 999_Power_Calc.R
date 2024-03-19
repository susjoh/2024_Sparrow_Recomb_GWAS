bonf <- 0.05/57000

source("r/gwas_power_calc.R")

pow = power_n_hsq(n = 6500, qsq = 0.01, pval = bonf)
pow

power_plot(pow, "x","y")

