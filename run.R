system('Rscript scripts/train_test_historical_forecasts_glmnet.R')

system('Rscript scripts/make_trends_2021.R')

system('Rscript scripts/make_forecasts_2021.R')
