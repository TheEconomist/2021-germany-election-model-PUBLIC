# _The Economist_ Germany federal parliamentary election forecast (Bundestagswahl 2021)

This repository stores R, Python and Stan code to replicate _The Economist_'s polls-based forecasting model for the 2021 Germany federal elections. This model combines structural "prior" information with an aggregate of publicly released opinion polls to predict the percentage breakdown of proportional party votes (sometimes called "second votes") in the Bundestag. 

Details about model methodology can be found [here](https://www.economist.com/graphic-detail/2021/08/10/how-the-economists-german-election-model-works). Daily-updated predictions are located at [this web address](https://www.economist.com/graphic-detail/who-will-succeed-angela-merkel).

# Setup

The model pipeline uses is based almost entirely on Python and R scripts (run in that order), though you will need to do some setup to enable R to run various models with the [Stan](https://mc-stan.org/) computational backend for Bayesian modeling. The model was most recently run with Python v3.9.6 and R v4.1.1. These can be installed from [python.org](https://www.python.org/) and [r-project.org](https://www.r-project.org/), though might work on earlier versions of R.

To install dependencies, navigate to the root of the directory and run `pipenv install`. You can install pipenv in MacOS via running `pip install pipenv` after you have setup python. You can read more about pip [here](https://pip.pypa.io/en/stable/installation/).

For R installation, you will need to install packages in the standard fashion and setup the backend to Stan, one statistical language that is called by the model's R scripts. Run this code in your R GUI of choice (we prefer RStudio):

```r
# package setup:
lapply(c("tidyverse","zoo","lubridate","gridExtra","mgcv","pbapply","brms","mvtnorm","lqmm","imputeTS","glmnet","data.table","caret","doParallel","foreach","gamlss"),install.packages)

# stan
remotes::install_github("stan-dev/cmdstanr")
cmdstanr::install_cmdstan()
```

# Running the model

Once everything is setup, you have your choice of scripts to run. Each is located in the `/scripts/` directory:

1. **`train_test_historical_forecasts_glmnet.R`**: This script trains our set of machine-learning and Bayesian models on historical polling data for German parliamentary elections strecthing from 1953 through 2017. The script takes some time to finish and will require a good amount of RAM — say, 8GB or so to be safe. This script only needs to be run once. It outputs a set of parameters and model objects that are used to make predictions on future data.

2. **`make_trends_2021.R`**: This script downloads the latest polling data for the 2021 election cycle, as aggregated by the excellent academics responsible for `https://www.wahlrecht.de/`. This script should be run every day.

3. **`make_forecasts_2021.R`**: This script takes the current polling data and combines the coefficients produced by `train_test_historical_forecasts_glmnet.R` to produce projections of party vote shares and uncertainty intervals for the 2021 cycle. It will output key summary data to the file `/output-data/site-data/election_day_summaries.csv`.

If you find any bugs


# Licence
This software is published by _[The Economist](https://www.economist.com)_ under the [MIT licence](https://opensource.org/licenses/MIT). The data generated by _The Economist_ are available under the [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

The licences include only the data and the software authored by _The Economist_, and do not cover any _Economist_ content or third-party data or content made available using the software. More information about licensing, syndication and the copyright of _Economist_ content can be found [here](https://www.economist.com/rights/).
© 2021 GitHub, Inc.