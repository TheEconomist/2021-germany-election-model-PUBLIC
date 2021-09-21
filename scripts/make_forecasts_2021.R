# package setup:
# lapply(c("tidyverse","zoo","lubridate","gridExtra","mgcv","pbapply","brms","mvtnorm","lqmm","imputeTS","glmnet","data.table","caret","doParallel","foreach","gamlss"),install.packages)
# remotes::install_github("stan-dev/cmdstanr")
# cmdstanr::install_cmdstan()


library(tidyverse)
library(zoo)
library(lubridate)
library(gridExtra)
library(pbapply)
library(mgcv)
library(brms)
library(cmdstanr)
library(mvtnorm)
library(lqmm)
library(imputeTS)
library(beepr)
library(glmnet)
library(gamlss)

# update the polling average first
# system('Rscript scripts/make_trends_2021.R')
 
# metadata
MASTER_days_until_2021_election = as.numeric(ymd('2021-09-26') - pmin(Sys.Date(),ymd('2021-09-26')))


# read in 2021 polls ------------------------------------------------------
polls_2021 <- read_csv('output-data/polls-with-trends-2021.csv',guess_max = 1e09) 

# now, get forecast prior for 2021 ----------------------------------------
# read prior information
prior <- read_csv('output-data/historical_plus2021_priors.csv',guess_max = 1e09) 

prior <- prior %>%
  mutate(chancellor_party = as.character(chancellor_party),
         gov = as.character(gov),
         prop_poll300_diff_from_lag_vote = abs(polls_300_days_out - lag_voteshare)/lag_voteshare)

plot(prior$prop_poll300_diff_from_lag_vote, abs(prior$lag_voteshare - prior$voteshare))


# model result as function of normal distribution
## mean is predicted using lag vote share, polls at ED-300, and structural party info
## sigma is predicted using polls at ED-300 and lag vote share
student_prior_model <- brm(bf(voteshare | trunc(lb=0,ub=1) ~ lag_voteshare + polls_300_days_out + chancellor_party + chancellor_party:term,
                              sigma | trun(lb=0,ub=0.2) ~ polls_300_days_out + lag_voteshare ), # + chancellor_party + chancellor_party:term),
                           data = prior,
                           family = 'student',
                           prior = c(set_prior('normal(0.01, 0.05)',class='Intercept'),
                                     set_prior('normal(0.23,  0.05)',class='b',coef='lag_voteshare'),
                                     set_prior('normal(0.69,  0.05)',class='b',coef='polls_300_days_out'),
                                     set_prior('normal(0.03, 0.02)',class='b',coef='chancellor_party1'),
                                     set_prior('normal(0, 0.01)',class='b',coef='chancellor_party1:term'),
                                     set_prior('normal(-0.01, 0.03)',class='b',coef='chancellor_party0:term'),
                                     # prior for nu 
                                     set_prior('normal(12,  5)',class='nu'),
                                     # priors for sigma
                                     set_prior('normal(-3.8,  0.15)',class='Intercept',dpar='sigma'),
                                     set_prior('normal(0.6,  0.3)',class='b',coef='polls_300_days_out',dpar='sigma'),
                                     set_prior('normal(0.74,  0.3)',class='b',coef='lag_voteshare',dpar='sigma') #,
                                     #set_prior('normal(0.8, 0.5)',class='b',coef='chancellor_party1',dpar='sigma'),
                                     #set_prior('normal(-0.04, 0.1)',class='b',coef='chancellor_party1:term',dpar='sigma'),
                                     #set_prior('normal(-0.27, 0.15)',class='b',coef='chancellor_party0:term',dpar='sigma')
                           ),
                           iter = 4000,
                           chains = 4, 
                           cores = 4,
                           warmup = 500,
                           backend='cmdstan',
                           control = list(adapt_delta = 0.9,refresh=100)) 


student_prior_model
mcmc_plot(student_prior_model)
plot(student_prior_model,'b_Intercept')

# add 2021 covars, including polls
prior_2021 <-  prior %>% 
  filter(election_year==2021) %>% dplyr::select(-polls_300_days_out) %>%
  left_join(
    polls_2021 %>%
      filter(poll == 'trend') %>%
      group_by(party) %>%
      mutate(election_year=2021,
             days_until_ed = as.numeric(ymd('2021-09-24') - mid),
             most_recent_poll_trend = last(trend)) %>%
      filter(days_until_ed == 300) %>% 
      dplyr::select(election_year,party,polls_300_days_out=trend,most_recent_poll_trend) %>%
      ungroup()
    
  )

preds <- predict(object = student_prior_model,
                 newdata = prior_2021,
                 probs = c(0.025,0.975)
)


prior_2021 <- prior_2021 %>% 
  mutate(est = preds[,1],
         se = preds[,2],
         lower = preds[,3],
         upper = preds[,4]) %>%
  arrange(desc(lag_voteshare))


prior_2021 %>% 
  mutate(party = as.factor(party),
         party = fct_inorder(party)) %>%
  ggplot(.,aes(x=party)) +
  geom_segment(aes(xend=party,y=lower,yend=upper),col='gray30') +
  geom_point(aes(y=lag_voteshare,shape='previous vote'),size=2) +
  geom_point(aes(y=most_recent_poll_trend,shape='latest poll'),size=2) +
  geom_point(aes(y=polls_300_days_out,shape='polls 300d out'),size=2) +
  geom_point(aes(y=est,shape='prior'),size=2) +
  theme_minimal() +
  labs(x='Party',y='Prior forecast, share of second votes',shape='') +
  scale_y_continuous(labels=function(x){x*100})

print('ready to blend with polls')

# now blend the prior with the polls --------------------------------------
historical_trends_and_priors <- read_csv('output-data/historical_trends_and_priors.csv',guess_max = 1e09) 
polls_2021_with_trends_forecast <- read_csv('output-data/raw-for-forecast-polls-with-trends-2021.csv',guess_max = 1e09) 

print('loaded polls')

as.numeric(max(ymd("2021-09-26") - polls_2021_with_trends_forecast$date,na.rm=T))

print('loading loo poll error')

# add das til ed and expected movement in polls
loo_poll_error_function <- read_rds('output-data/loo_poll_error_function.rds')

print('loaded loo poll error')

polls_2021_with_trends_forecast <- polls_2021_with_trends_forecast %>%
  mutate(days_until_ed = as.numeric(ymd('2021-09-26') - mid),
         loess_trend = extrapolated_loess_trend) %>%
  mutate(loess_se = loo_poll_error_function(vote = loess_trend, days = days_until_ed)) %>%
  dplyr::select(mid,days_until_ed,party,loess_trend,loess_se) 


# this runs a separate glmnet model for every day in past campaign
holdout_year <- 2021
historical_coefs_raw <- 
  map_df(0:635,
         function(x){
           # filter data
           train_validate_data <- historical_trends_and_priors %>% 
             filter(election_year != holdout_year, days_until_ed == x)  %>%
             mutate(fold = as.integer(as.factor(election_year)))
           
           # train model on historical data
           todays_mod <- cv.glmnet(x = model.matrix(~ 0 + prior_est + loess_trend + prior_est*loess_trend ,
                                                    train_validate_data),
                                   y = train_validate_data$voteshare,
                                   family = 'gaussian',
                                   type.measure="mse",
                                   intercept = TRUE,
                                   alpha = 0.9,
                                   lambda = seq(0.001,0.1,by = 0.001),
                                   nfolds = length(unique(train_validate_data$election_year)), 
                                   foldid = train_validate_data$fold, 
                                   parallel = FALSE)
           
           # get coefficients
           as.matrix(t(coef(todays_mod,s = 'lambda.min'))) %>% 
             as_tibble() %>%
             mutate(days_until_ed = x) %>%
             rename(intercept = `(Intercept)`,
                    #prior_est_prior_se = `prior_est:prior_se`,
                    #loess_trend_loess_se = `loess_trend:loess_se`,
                    loess_trend_prior_est = `prior_est:loess_trend`) 
           
           
           
         }) %>%
  bind_rows 

historical_coefs_temp <- historical_coefs_raw

names(historical_coefs_temp) <- paste0(names(historical_coefs_temp),'.c')
historical_coefs_temp <- historical_coefs_temp %>% rename(days_until_ed = days_until_ed.c)


# graph coefficients to gauge smoothness
ggplot(historical_coefs_temp,aes(x=days_until_ed)) +
  geom_line(aes(y=intercept.c,col='intercept')) +
  geom_line(aes(y=loess_trend.c,col='poll')) +
  geom_line(aes(y=prior_est.c,col='prior')) +
  #geom_line(aes(y=prior_est_prior_se.c,col='prior*prior_sd')) +
  #geom_line(aes(y=loess_trend_loess_se.c,col='polls*polls_sd')) +
  geom_line(aes(y=loess_trend_prior_est.c,col='polls*prior')) +
  scale_x_reverse(limits=c(300,0),breaks=seq(0,300,30)) +
  labs(y= 'coefficient')

names(historical_coefs_temp)

knots <- 100
basis <- 'ts'

intercept_smooth_coef <- gam(intercept.c ~ s(days_until_ed, k = knots, bs=basis) , data = historical_coefs_temp)

prior_smooth_coef <- gam(prior_est.c ~ s(days_until_ed, k = knots, bs=basis) , data = historical_coefs_temp)
poll_smooth_coef <- gam(loess_trend.c ~ s(days_until_ed, k = knots, bs=basis) , data = historical_coefs_temp)

#poll_se_coef <- gam(loess_se.c ~ s(days_until_ed, k = knots) , data = historical_coefs_temp)
#prior_se_coef <- gam(prior_se.c ~ s(days_until_ed, k = knots) , data = historical_coefs_temp)

prior_poll_smooth_coef <- gam(loess_trend_prior_est.c ~ s(days_until_ed, k = knots, bs=basis) ,data = historical_coefs_temp)

#poll_trend_poll_se.c <- gam(loess_trend_loess_se.c ~ s(days_until_ed, k = 50) , data = historical_coefs_temp)
#prior_trend_prior_se.c <- gam(prior_est_prior_se.c ~ s(days_until_ed, k = 50) , data = historical_coefs_temp)


#plot(predict(intercept_smooth_coef, historical_coefs_temp),historical_coefs_temp$intercept.c)
#plot(predict(poll_trend_poll_se.c, historical_coefs_temp),historical_coefs_temp$loess_trend_loess_se.c)
#plot(predict(prior_poll_smooth_coef, historical_coefs_temp),historical_coefs_temp$loess_trend_prior_est.c)

# use smooth coefficients (commented out = raw coefs)
historical_coefs_temp <- historical_coefs_temp %>%
  mutate(intercept.c = predict(intercept_smooth_coef,.),
         
         prior_est.c = predict(prior_smooth_coef,.),
         loess_trend.c = predict(poll_smooth_coef,.),
         
         #loess_se.c = predict(poll_se_coef, .),
         #prior_se.c = predict(prior_se_coef, .),
         
         #loess_trend_loess_se.c = predict(poll_trend_poll_se.c, .),
         #prior_est_prior_se.c = predict(prior_trend_prior_se.c, .),
         loess_trend_prior_est.c = predict(prior_poll_smooth_coef,.)
         
  )

ggplot(historical_coefs_temp,aes(x=days_until_ed)) +
  geom_line(aes(y=intercept.c,col='intercept')) +
  geom_line(aes(y=loess_trend.c,col='poll')) +
  geom_line(aes(y=prior_est.c,col='prior')) +
  #geom_line(aes(y=prior_est_prior_se.c,col='prior*prior_sd')) +
  #geom_line(aes(y=loess_trend_loess_se.c,col='polls*polls_sd')) +
  geom_line(aes(y=loess_trend_prior_est.c,col='polls*prior')) +
  scale_x_reverse(limits=c(300,0),breaks=seq(0,300,30)) +
  labs(y= 'coefficient')


# generate out-of-sample predictions for the TEST year, based on smooth coefs from validation model
prior_data_posterior <- polls_2021_with_trends_forecast %>%
  # first, predict smoothed coefficients
  mutate(intercept.c = predict(intercept_smooth_coef,.),
         
         prior_est.c = predict(prior_smooth_coef,.),
         loess_trend.c = predict(poll_smooth_coef,.),
         
         #loess_se.c = predict(poll_se_coef, .),
         #prior_se.c = predict(prior_se_coef, .),
         
         #loess_trend_loess_se.c = predict(poll_trend_poll_se.c, .),
         #prior_est_prior_se.c = predict(prior_trend_prior_se.c, .),
         loess_trend_prior_est.c = predict(prior_poll_smooth_coef,.)
         
  ) %>%
  # now combine the predictions
  left_join(prior_2021 %>% dplyr::select(party, prior_est=est, prior_se=se), by = "party") %>%
  mutate(lm_pred = 
           # I
           intercept.c + 
           # 1-way 
           (prior_est*prior_est.c) + 
           (loess_trend*loess_trend.c) + 
           #(prior_se*prior_se.c) + 
           #(loess_se*loess_se.c) + 
           # interactions
           #(loess_trend*loess_se*loess_trend_loess_se.c) +
           #(prior_est*prior_se*prior_est_prior_se.c) +
           (loess_trend*prior_est*loess_trend_prior_est.c)  
  ) %>%
  mutate(election_year = 2021) %>%
  dplyr::select(election_year,party,mid,days_until_ed,loess_trend,loess_se,prior_est,prior_se,lm_pred) %>%
  ungroup()


prior_data_posterior <- prior_data_posterior[!duplicated(prior_data_posterior),]

# filter to election year
prior_data_posterior <- prior_data_posterior %>%
  filter(days_until_ed <= 635)

print ('loading lm error model')

# add expected error for the lm
err_mod <- read_rds('output-data/lm_error_model.rds')
print('loaded error mod, loading training')
err_mod_training_data <- read_rds('output-data/lm_error_model_training_data.rds')
print('loaded training')

prior_data_posterior <- prior_data_posterior %>%
  mutate(party_is_cdu = party=='cdu')  %>%
  left_join(
    historical_trends_and_priors %>%
      group_by(election_year) %>%
      summarise(rmse_poll_error_cycle = sqrt(mean((voteshare - loess_trend)^2))) %>%
      bind_rows(tibble(election_year = 2021,rmse_poll_error_cycle = NA)) %>%
      mutate(rmse_poll_error_last_three_cycles = 
               (lag(rmse_poll_error_cycle) +
                  lag(rmse_poll_error_cycle,2) + 
                  lag(rmse_poll_error_cycle,3) ) / 3) %>% 
      mutate(rmse_poll_error_last_three_cycles = na_locf(rmse_poll_error_last_three_cycles)) %>%
      dplyr::select(c(election_year,rmse_poll_error_last_three_cycles )),
    by = 'election_year'
  ) 
  

test_data_sigmas <- predict(object = err_mod, 
                            data = err_mod_training_data, 
                            newdata = 
                              as.data.frame( prior_data_posterior %>%
                                               dplyr::select(lm_pred,days_until_ed,party_is_cdu,rmse_poll_error_last_three_cycles)), 
                            what = 'sigma',type = 'response')


test_data_nus <- predict(object = err_mod, 
                            data = err_mod_training_data, 
                            newdata = 
                              as.data.frame( prior_data_posterior %>%
                                               dplyr::select(lm_pred,days_until_ed,party_is_cdu,rmse_poll_error_last_three_cycles)), 
                            what = 'nu',type = 'response')

prior_data_posterior <- prior_data_posterior %>%
  ungroup() %>%
  mutate(lm_fitted_error = pmax(0.01,test_data_sigmas),
         lm_fitted_df = pmax(2, test_data_nus)) 


plot(prior_data_posterior[prior_data_posterior$days_until_ed<=100,]$days_until_ed,
     prior_data_posterior[prior_data_posterior$days_until_ed<=100,]$lm_fitted_error)


# remove any duplicated days
prior_data_posterior = prior_data_posterior %>%
  group_by(election_year, days_until_ed, party) %>%
  filter(row_number() == 1) %>%
  ungroup() 

prior_data_posterior %>% group_by(party) %>% mutate(days_elapsed = mid - lag(mid)) %>% 
  pull(days_elapsed) %>% table %>%
  print


# now simulate ------------------------------------------------------------
print('simulating elections')

# first, look at pairs of daily variations
# one round of printing for debugging
prior_data_posterior %>% 
  mutate(election_year_day = paste0(election_year, days_until_ed)) %>%
  dplyr::filter(days_until_ed <= (300 + 365)) %>%  # for all of last 365 days
  dplyr::select(election_year_day,party,loess_trend) %>%
  slice(185:195) %>%
  print


prior_data_posterior %>% 
  mutate(election_year_day = paste0(election_year, days_until_ed)) %>%
  dplyr::filter(days_until_ed <= (300 + 365)) %>%  # for all of last 365 days
  dplyr::select(election_year_day,party,loess_trend) %>%
  distinct() %>%
  spread(party,loess_trend) %>%
  # mutate_at(unique(election_year_data_full$party),
  #           function(x){ pmax(-0.01,pmin(0.01,(x-lag(x)))  ) }) %>%
  # na.omit() %>%
  # mutate_at(unique(election_year_data_full$party),
  #           function(x){ ifelse(abs((x - mean(x))/sd(x) ) > 1.628, NA, x  )   }) %>%
  mutate_at(unique(prior_data_posterior$party),
            function(x){x + rnorm(length(x),0,0.001)}) %>%
  dplyr::select(unique(prior_data_posterior$party)) %>%
  na.omit() %>%
  ggplot(.,aes(x=cdu,y=afd)) +
  geom_point() +
  geom_smooth(method='lm')

# this extracts the cor
master_correlation_between_parties <- prior_data_posterior %>% 
  mutate(election_year_day = paste0(election_year, days_until_ed)) %>%
  filter(days_until_ed <= (300 + 365)) %>%  # for all of last 365 days
  dplyr::select(election_year_day,party,loess_trend) %>%
  distinct() %>%
  spread(party,loess_trend) %>%
  # mutate_at(unique(election_year_data_full$party),
  #           function(x){ pmax(-0.01,pmin(0.01,(x-lag(x)))  ) }) %>%
  # na.omit() %>%
  # mutate_at(unique(election_year_data_full$party),
  #           function(x){ ifelse(abs((x - mean(x))/sd(x) ) > 1.628, NA, x  )   }) %>%
  mutate_at(unique(prior_data_posterior$party),
            function(x){x + rnorm(length(x),0,0.001)}) %>%
  dplyr::select(unique(prior_data_posterior$party)) %>%
  na.omit() %>%
  dplyr::select(afd,cdu,fdp,gru,lin,oth,spd) %>%
  cor

# function for making forecasts for year y at t days before the election
simulate_election_y_from_time_t <- function(t, num_sims = 10000, return_all = FALSE){
  print(sprintf('forecasts for 2021 made at ED-%s',t))
  
  # get correlation between parties UP TO THAT DAY
  election_year_data_full <- prior_data_posterior %>% 
    dplyr::filter(days_until_ed >= t) %>%
    mutate(election_year = 2021)
  
  correlation_between_parties_year <- election_year_data_full %>%
    mutate(election_year_day = paste0(election_year, days_until_ed)) %>%
    filter(days_until_ed <= (t + 365)) %>%  # for all of last 365 days
    arrange(days_until_ed) %>%
    dplyr::select(election_year_day,party,loess_trend) %>%
    distinct() %>%
    spread(party,loess_trend) %>%
    # mutate_at(unique(election_year_data_full$party),
    #           function(x){ pmax(-0.01,pmin(0.01,(x-lag(x)))  ) }) %>%
    # na.omit() %>%
    # mutate_at(unique(election_year_data_full$party),
    #           function(x){ ifelse(abs((x - mean(x))/sd(x) ) > 1.628, NA, x  )   }) %>%
    mutate_at(unique(prior_data_posterior$party),
              function(x){x + rnorm(length(x),0,0.001)}) %>%
    dplyr::select(unique(prior_data_posterior$party)) %>%
    na.omit() %>%
    # ggplot(aes(x=cdu,y=gru)) + geom_point() + geom_smooth(method='lm')
    cor
  
  # correlation_between_parties_year <- correlation_between_parties_year * 0.75
  # correlation_between_parties_year[correlation_between_parties_year > 0.8] <- 0.8
  # correlation_between_parties_year[correlation_between_parties_year < -0.8] <- -0.8
  diag(correlation_between_parties_year) <- 1
  
  correlation_between_parties_year <- lqmm::make.positive.definite(correlation_between_parties_year)
  correlation_between_parties_year
  
  # get data for that day
  election_year_data <- prior_data_posterior  %>%
    mutate(election_year = 2021) %>% 
    filter(days_until_ed == t)
  
  # get mean and sigma
  party_mus_sigmas <- election_year_data %>% 
    dplyr::select(party,lm_pred,lm_fitted_error,lm_fitted_df)
  
  # simulate 10k random walks
  simulated_errors <-  mvtnorm::rmvt(n = (t+1)*num_sims, 
                                     sigma = correlation_between_parties_year[match(rownames(correlation_between_parties_year),party_mus_sigmas$party),],
                                     df = 10) # party_mus_sigmas$lm_fitted_df)
  
  # convert to random walk
  simulated_errors <- simulated_errors %>% 
    as_tibble %>% 
    mutate(trial = floor((row_number() - 1) / (t+1)) +1 ) %>% 
    group_by(trial) %>%
    mutate_at(vars(-group_cols()), cumsum) %>%
    # force to start at 0
    mutate_at(vars(-group_cols()), function(x){x-first(x)}) %>%
    ungroup() %>%
    dplyr::select(-trial) %>%
    as.matrix
  
  # multiply error times sd
  simulated_errors <- simulated_errors * 
    outer(rep.int(1L, nrow(simulated_errors)), 
          party_mus_sigmas$lm_fitted_error / sqrt(t+1))
  
  simulated_errors %>% 
    as_tibble %>% 
    mutate(trial = floor((row_number() - 1) / (t+1)) +1 ) %>% 
    group_by(trial) %>%
    mutate(step = row_number()) %>%
    ungroup() %>%
    filter(trial %in% round(runif(100,0,1000))) %>%
    ggplot(.,aes(x=step,y=V1,col=trial,group=trial)) +
    geom_line()
  
  
  # finally, add the intercept
  simulated_errors <-   simulated_errors + 
    outer(rep.int(1L, nrow(simulated_errors)),
          party_mus_sigmas$lm_pred)
  
  simulated_errors %>% 
    as_tibble %>% 
    mutate(trial = floor((row_number() - 1) / (t+1)) +1 ) %>% 
    group_by(trial) %>%
    mutate(step = row_number()) %>%
    group_by(step) %>%
    summarise(mean = mean(V1),
              upper = quantile(V1,0.975),
              lower = quantile(V1,0.025)) %>%
    ggplot(.,aes(x=step)) +
    geom_line(aes(y=mean)) +
    geom_line(aes(y=upper),linetype=2) +
    geom_line(aes(y=lower),linetype=2) 
  
  
  # clean up and look
  simulated_errors <- apply(simulated_errors,2,function(x){pmin(1,x)})
  simulated_errors <- apply(simulated_errors,2,function(x){pmax(0,x)})
  
  simulated_errors <- simulated_errors %>%
    as_tibble() %>%
    ungroup() 
  
  rownames(simulated_errors) <- NULL
  names(simulated_errors) <- party_mus_sigmas$party
  
  # add trial index
  simulated_errors <- simulated_errors %>%
    mutate(trial = floor((row_number() - 1) / (t+1)) +1 ) %>% 
    group_by(trial) %>%
    mutate(days_until_ed = (t+1) - dplyr::row_number()) %>%
    ungroup() %>%
    mutate(election_year = 2021)
  
  range(simulated_errors$days_until_ed)
  range(simulated_errors$trial)
  
  # forecast generated on? 
  simulated_errors <- simulated_errors %>% 
    mutate(forecast_date = t)
  
  # nothing greater or less than one
  simulated_errors <- simulated_errors %>%
    mutate_at(party_mus_sigmas$party,
              function(x){pmax(0,pmin(1,x))})
  
  # make sure doesn't sum above 1
  row_sums <- rowSums(simulated_errors[1:(ncol(simulated_errors)-4)])
  
  simulated_errors[,1:(ncol(simulated_errors)-4)] <- apply(X = simulated_errors[1:(ncol(simulated_errors)-4)],
                                                           MARGIN = 2,
                                                           FUN = function(x){
                                                             x / row_sums
                                                           }) %>% as_tibble
  
  
  # again, nothign greater or less than one
  simulated_errors <- simulated_errors %>%
    mutate_at(party_mus_sigmas$party,
              function(x){pmax(0,pmin(1,x))})
  
  rowSums(simulated_errors[1:(ncol(simulated_errors)-4)])
  
  # and vote shares can't be above 1 or below 0
  simulated_errors <- simulated_errors %>%
    mutate(trial = floor((row_number() - 1) / (t+1)) +1 ) %>% 
    group_by(trial) %>%
    mutate(days_until_ed = (t+1) - dplyr::row_number()) %>%
    ungroup() %>%
    mutate(election_year = 2021)
  
  range(simulated_errors$days_until_ed)
  range(simulated_errors$trial)
  
  party_mus_sigmas %>%
    left_join(
      simulated_errors %>% 
        filter(days_until_ed == 0) %>% 
        gather(party,vote,1:7) %>% group_by(party) %>% summarise(sd = sd(vote,na.rm=T))
    )
  
  # simulated_errors %>%
  #   dplyr::filter(trial %in% 1:100) %>%
  #   gather(party,simulated_vote,1:(length(.)-4)) %>%
  #   ggplot(.,aes(x=days_until_ed, y=simulated_vote, col=party, group=paste0(party,trial) ),fill=NA) +
  #   geom_line() +
  #   scale_x_reverse()
  
  if(return_all){return(simulated_errors)}else{return(simulated_errors[simulated_errors$days_until_ed==0,])}
  
}

election_day_simulation_data <- map_df(MASTER_days_until_2021_election + c(0,7,14,21,30,60,90,120,150,200),
                                       function(x){
                                         simulate_election_y_from_time_t(t = x,num_sims = 10000)
                                       })

names(election_day_simulation_data)[1:(length(election_day_simulation_data)-4)]



# look at the relationsip between party vote shares in the simulations
party_pairs <- as_tibble(t(combn(names(election_day_simulation_data)[1:(length(election_day_simulation_data)-4)],2))) %>% setNames(.,c('A','B')) %>%
  filter(A != 'oth', B !=' oth')

party_pairs_plots <- vector('list',nrow(party_pairs))

for(row in 1:length(party_pairs_plots)){
  party_pairs_plots[[row]] <- 
    ( 
      election_day_simulation_data %>% 
        mutate(A = election_day_simulation_data[[party_pairs[row,]$A]],
               B = election_day_simulation_data[[party_pairs[row,]$B]]) %>%
        dplyr::filter(days_until_ed == 0) %>% sample_n(5000) %>%
        ggplot(., aes(x = A, y = B)) +
        geom_point(alpha=0.1) +
        theme_minimal() +
        labs(x= party_pairs[row,]$A, 
             y = party_pairs[row,]$B) +
        geom_smooth(method='lm',linetype=2)  +
        coord_cartesian(xlim=c(0,0.6),ylim=c(0,0.6))
    )
}

# do.call("grid.arrange", c(party_pairs_plots, ncol=5))


prior_data_posterior %>%
  left_join(
    election_day_simulation_data %>% 
      gather(party,prediction,1:(length(.)-4)) %>%
      group_by(election_year, party, days_until_ed = forecast_date) %>%
      summarise(yhat_mean = mean(prediction),
                yhat_upper_95 = quantile(prediction,0.975),
                yhat_lower_95 = quantile(prediction,0.025),
                yhat_upper_80 = quantile(prediction,0.9),
                yhat_lower_80 = quantile(prediction,0.1)) %>%
      ungroup(),
    by = c("election_year","party", "days_until_ed")
  ) %>%
  mutate(party = as.factor(party),
         party = fct_reorder(party, -loess_trend, mean)) %>%
  na.omit() %>%
  ggplot(., aes(x=days_until_ed,col=party,fill=party,group=party)) +
  geom_line(aes(y=loess_trend,col='polling average'),linetype=2) +
  geom_line(aes(y=yhat_mean,col='prediction')) +
  geom_ribbon(aes(ymin=yhat_lower_95,ymax=yhat_upper_95),alpha=0.6,col=NA) +
  # geom_point(aes(x=0,y=voteshare)) +
  facet_wrap(~party) +
  scale_x_reverse() +
  theme_minimal() +
  labs(x='forecast generated at x days before election',
       y='predicted vote share',
       subtitle = 'forecasts for 2021')


election_day_simulation_data <- simulate_election_y_from_time_t(t = MASTER_days_until_2021_election, 
                                                                num_sims = 50000, return_all = TRUE)


# output ------------------------------------------------------------------
# output the raw simulations
print('writing sims')
write_rds(election_day_simulation_data,'output-data/site-data/2021_forecast_simulations.rds',compress = 'gz')
print('wrote sims')

print('re-reading sims')
# output the e-day summaries
election_day_simulation_data <- read_rds('output-data/site-data/2021_forecast_simulations.rds')
print('read sims')

## party vote share
party_summaries <- election_day_simulation_data %>% 
  filter(days_until_ed == 0) %>%
  gather(party,prediction,1:(length(.)-4)) %>%
  group_by(election_year, party) %>%
  summarise(yhat_mean = mean(prediction),
            yhat_upper_95 = quantile(prediction,0.975),
            yhat_lower_95 = quantile(prediction,0.025),
            yhat_upper_80 = quantile(prediction,0.9),
            yhat_lower_80 = quantile(prediction,0.1)) %>%
  ungroup() %>%
  left_join(
    prior_data_posterior %>% 
      filter(days_until_ed == min(days_until_ed,na.rm=T)) %>%
      dplyr::select(party,loess_trend,loess_se,prior_est,prior_se)
  ) %>%
  mutate(type = 'party')


## party SEAT share 
party_summaries <- party_summaries %>%
  left_join(
    election_day_simulation_data %>% 
      # transform to seats
      filter(days_until_ed == 0) %>%
      dplyr::select(-oth) %>%
      mutate_at(c('afd','fdp','gru'), # reasonable these would be kicked out if they didn't meet national threshold
                function(x){ifelse(x >= 0.05, x, 0)}) %>%
      mutate(seat_denom = afd + cdu + fdp + gru + lin + spd) %>%
      mutate(afd = afd / seat_denom,
             cdu = cdu / seat_denom,
             fdp = fdp / seat_denom,
             gru = gru / seat_denom,
             lin = lin / seat_denom,
             spd = spd/ seat_denom) %>%
      # summarise as above
      gather(party,seats_won,1:(length(.)-4))  %>%
      group_by(election_year, party) %>%
      summarise(yhat_seat_mean = mean(seats_won),
                yhat_seat_upper_95 = quantile(seats_won,0.975),
                yhat_seat_lower_95 = quantile(seats_won,0.025),
                yhat_seat_upper_80 = quantile(seats_won,0.9),
                yhat_seat_lower_80 = quantile(seats_won,0.1),
                pr_majority = mean(seats_won > 0.5))
  )


# party largest seat probabilities
prob_largest_party = election_day_simulation_data %>% 
  # transform to seats
  filter(days_until_ed == 0) %>%
  dplyr::select(-oth) %>%
  mutate_at(c('afd','fdp','gru'), # reasonable these would be kicked out if they didn't meet national threshold
            function(x){ifelse(x >= 0.05, x, 0)}) %>%
  mutate(seat_denom = afd + cdu + fdp + gru + lin + spd) %>%
  mutate(afd = afd / seat_denom,
         cdu = cdu / seat_denom,
         fdp = fdp / seat_denom,
         gru = gru / seat_denom,
         lin = lin / seat_denom,
         spd = spd/ seat_denom) %>% 
  group_by(days_until_election=forecast_date) %>%
  summarise(spd_largest_party = mean(spd > afd &
                                       spd > cdu &
                                       spd > fdp & 
                                       spd > gru &
                                       spd > lin),
            cdu_largest_party = mean(cdu > afd &
                                       cdu > spd &
                                       cdu > fdp &
                                       cdu > gru &
                                       cdu > lin),
            gru_largest_party = mean(gru > afd &
                                       gru > spd &
                                       gru > fdp &
                                       gru > cdu &
                                       gru > lin),
            fdp_largest_party = mean(fdp > afd &
                                       fdp > spd & 
                                       fdp > gru & 
                                       fdp > cdu &
                                       fdp > lin),
            lin_largest_party = mean(lin > afd &
                                       lin > spd & 
                                       lin > gru & 
                                       lin > cdu &
                                       lin > fdp)) %>% 
  dplyr::select(c(2:ncol(.))) %>%
  set_names(.,c('spd','cdu','gru','fdp','lin')) %>%
  gather(party, pr_largest_party,1:5)

party_summaries = party_summaries %>%
  left_join(prob_largest_party)


# coalitions
## vote and seat share in one function
get_coalition_summary <- function(parties){
  # votes
  votes <- election_day_simulation_data %>% 
    filter(days_until_ed == 0) %>%
    dplyr::select(all_of(c('trial',parties))) %>%
    gather(party,prediction,2:length(.)) %>%
    group_by(trial) %>%
    summarise(prediction = sum(prediction)) %>%
    ungroup() %>%
    summarise(yhat_mean = mean(prediction),
              yhat_upper_95 = quantile(prediction,0.975),
              yhat_lower_95 = quantile(prediction,0.025),
              yhat_upper_80 = quantile(prediction,0.9),
              yhat_lower_80 = quantile(prediction,0.1)) %>%
    ungroup() 
  
  
  seats <- election_day_simulation_data %>% 
    # transform to seats
    filter(days_until_ed == 0) %>%
    dplyr::select(-oth) %>%
    mutate_at(c('afd','fdp','gru'),
              function(x){ifelse(x >= 0.05, x, 0)}) %>%
    mutate(seat_denom = afd + cdu + fdp + gru + lin + spd) %>%
    mutate(afd = afd / seat_denom,
           cdu = cdu / seat_denom,
           fdp = fdp / seat_denom,
           gru = gru / seat_denom,
           lin = lin / seat_denom,
           spd = spd/ seat_denom) %>%
    # choose just select parties
    dplyr::select(all_of(c('trial',parties))) %>%
    # summarise as above
    gather(party,seats_won,2:length(.))  %>%
    group_by(trial) %>%
    summarise(seats_won = sum(seats_won)) %>%
    ungroup() %>%
    summarise(yhat_seat_mean = mean(seats_won),
              yhat_seat_upper_95 = quantile(seats_won,0.975),
              yhat_seat_lower_95 = quantile(seats_won,0.025),
              yhat_seat_upper_80 = quantile(seats_won,0.9),
              yhat_seat_lower_80 = quantile(seats_won,0.1),
              pr_majority = mean(seats_won > 0.5))
  

    votes %>% cbind(seats) %>% return
}

bind_rows(
  # bg = cdu + gru
  get_coalition_summary(c('cdu','gru')) %>% mutate(party = 'bg'),
  # jamaica = cdu + gru + fdp
  get_coalition_summary(c('cdu','gru','fdp')) %>% mutate(party = 'jamaica'),
  # traffic light = spd + gru + fdp
  get_coalition_summary(c('spd','gru','fdp')) %>% mutate(party = 'traffic'),
  # rrg = spd + lin + gru
  get_coalition_summary(c('spd','lin','gru')) %>% mutate(party = 'rrg'),
  # grand = cdu + spd
  get_coalition_summary(c('cdu','spd')) %>% mutate(party = 'grand'),
  # by = cdu + fdp
  get_coalition_summary(c('cdu','fdp')) %>% mutate(party = 'by'),
  # germany = cdu + spd + fdp
  get_coalition_summary(c('cdu','spd','fdp')) %>% mutate(party = 'germany'),
  # kenya = cdu + spd + gru
  get_coalition_summary(c('cdu','spd','gru')) %>% mutate(party = 'kenya'),
  # rg = spd + gru
  get_coalition_summary(c('spd','gru')) %>% mutate(party='rg')
) %>%
  mutate(type='coalition',election_year = 2021) -> coalition_summaries



party_summaries %>%
  bind_rows(coalition_summaries) %>%
  write_csv(.,'output-data/site-data/election_day_summaries.csv')


nsims = nrow(election_day_simulation_data %>% filter(days_until_ed == 0))
seat_sims = election_day_simulation_data %>%
  filter(days_until_ed == 0) %>%
  mutate(
    seat_denom = 
      ifelse(afd >= 0.05, afd, 0) +
      ifelse(cdu >= 0.00, cdu, 0) +
      ifelse(fdp >= 0.05, fdp, 0) +
      ifelse(gru >= 0.05, gru, 0) +
      ifelse(lin >= 0.00, lin, 0) +
      ifelse(spd >= 0.00, spd, 0)
    # seat_denom = sum(sapply(c(afd, cdu, fdp, lin, gru, spd), function(x) { ifelse(x >= 0.05, x, 0)}))
  ) %>%
  pivot_longer(cols=c('afd','cdu','oth','spd','lin','gru','fdp'), names_to='party', values_to='pct') %>%
  mutate(seat_pct = ifelse((party %in% c('cdu','spd','lin') | pct >= 0.05 ) & party != 'oth', pct / seat_denom, 0)) %>%
  mutate(days_until_ed=NULL, election_year=NULL, forecast_date=NULL)



aggregate_coalition_sims = function(data, parties, key) {
  data %>%
    filter(party %in% parties) %>%
    group_by(trial, seat_denom) %>%
    summarise(
      pct = sum(pct),
      seat_pct = sum(seat_pct)
    ) %>%
    mutate(party = key) %>%
    ungroup()
}

seat_sims_all = bind_rows(
  aggregate_coalition_sims(seat_sims, c('cdu','gru'), 'bg'),
  aggregate_coalition_sims(seat_sims, c('cdu','gru','fdp'), 'jamaica'),
  aggregate_coalition_sims(seat_sims, c('spd','gru','fdp'), 'traffic'),
  aggregate_coalition_sims(seat_sims, c('spd','gru','lin'), 'rrg'),
  aggregate_coalition_sims(seat_sims, c('cdu','spd'), 'grand'),
  aggregate_coalition_sims(seat_sims, c('cdu','fdp'), 'by'),
  aggregate_coalition_sims(seat_sims, c('cdu','gru','spd'), 'kenya'),
  aggregate_coalition_sims(seat_sims, c('gru','spd'), 'rg'),
  aggregate_coalition_sims(seat_sims, c('cdu','fdp','spd'), 'germany')
) %>%
  mutate(type='coalition') %>%
  bind_rows(
    seat_sims %>% mutate(type='party')
  )
  

# conditional probabilities
# - probability that no two-way coalition has a majority. (Pretty high, I assume.)
pr_no_two_way_majority = seat_sims_all %>%
  filter(party %in% c('bg','by','grand')) %>%
  group_by(trial) %>%
  summarise(no_two_way_majority = ifelse(all(seat_pct < 0.5), 1, 0)) %>%
  summarise(pr_no_two_way_majority = mean(no_two_way_majority))

# - chances that Jamaica and traffic-light have a majority, and that RRG does not
pr_jamaica_traffic_no_rrg = seat_sims_all %>%
  filter(party %in% c('traffic','jamaica','rrg')) %>%
  dplyr::select(trial,party,seat_pct) %>%
  pivot_wider(names_from = c(party), values_from = c(seat_pct)) %>%
  group_by(trial) %>%
  summarise(jamaica_traffic_no_rrg = ifelse(jamaica > 0.5 & traffic > 0.5 & rrg < 0.5, 1, 0)) %>%
  summarise(jamaica_traffic_no_rrg = mean(jamaica_traffic_no_rrg))

# - chances of SPD, Union, and Greens coming first
prob_largest_party



# some graphs -------------------------------------------------------------
# show prior and forecast point predictiona and ci
prior_2021 %>% 
  mutate(type='a') %>%
  bind_rows(party_summaries %>%
              dplyr::select(party,mean=yhat_mean,upper=yhat_upper_95,lower=yhat_lower_95) %>%
              mutate(type='b')) %>%
  mutate(party = as.factor(party),
         party = fct_inorder(party)) %>%
  # graph
  ggplot(.,aes(x=type)) +
  geom_segment(aes(xend=type,y=lower,yend=upper,group=type),position=position_dodge(),col='gray30') +
  geom_point(aes(y=lag_voteshare,shape='previous vote'),size=2) +
  geom_point(aes(y=most_recent_poll_trend,shape='latest poll'),size=2) +
  geom_point(aes(y=polls_300_days_out,shape='polls 300d out'),size=2) +
  geom_point(aes(y=est,shape='prior'),size=2) +
  geom_point(aes(y=mean,shape='posterior'),size=2) +
  labs(x=' ',y='Sshare of second votes',shape='') +
  scale_y_continuous(labels=function(x){x*100}) +
  facet_grid(cols=vars(party)) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())


# show probabilities for coalitions over time
read_csv('output-data/historical_coalition_shares_probs.csv',guess_max = 1e09) %>%
  bind_rows(
    coalition_summaries %>% mutate(election_year = 2021)  %>% rename(coalition=party)
  ) %>%
  ggplot(., aes(x=election_year,col=coalition)) + 
  geom_line(aes(y=yhat_seat_mean))
