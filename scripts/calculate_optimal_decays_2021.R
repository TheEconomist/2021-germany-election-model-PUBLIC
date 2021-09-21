rm(list=ls())
library(tidyverse)
library(zoo)
library(lubridate)
library(gridExtra)
library(parallel)
library(pbapply)
library(beepr)
library(imputeTS)
library(mgcv)
library(data.table)
library(dtplyr)

# scrape new polls and read them in
system("python3 scripts/fetch-data.py")

polls_2020 <- read_csv('output-data/all-polls.csv') %>% 
  select(-X1) %>%
  mutate(start = as_date(start), end = as_date(end),
         start = as_date(ifelse(start > end, start - 365, start)),
         mid = as_date(start + round(as.numeric(end - start)/2)))

# optimized loess trend ---------------------------------------------------
polls_2020 <- polls_2020 %>%
  group_by(poll) %>%
  mutate(n = ifelse(is.na(n),mean(n,na.rm=T),n)) %>%
  ungroup() %>%
  mutate(n = ifelse(is.na(n),mean(n,na.rm=T),n)) %>%
  ungroup()

mean(polls_2020$n)

# first, overall trends
polls_2020_demo <- polls_2020 
loess_spans <- lapply(c(seq(0.04,1,0.01),seq(1,4,0.1)),
                      function(x){
                        #print(x)
                        rmse_in <- try((polls_2020_demo %>%
                                          dplyr::filter(party == 'FDP') %>%
                                          mutate(mid = as.numeric(mid)) %>%
                                          arrange(mid) %>% 
                                          mutate(loess_trend = loess(share ~ mid, span = x,
                                                                     weights = sqrt(n/1825))$fitted) %>%
                                          mutate(error = loess_trend - lead(share)) %>%
                                          pull(error) %>% na.omit %>% .^2 %>% mean %>% sqrt)
                        )
                        
                        if(!is.numeric(rmse_in)){
                          return(NULL)
                        }else{
                          #print(rmse_in)
                          return(  tibble(span = x, rmse = rmse_in)  )
                        }
                        
                        
                      }) %>% 
  bind_rows %>% na.omit

loess_spans
optimal_span <- loess_spans$span[which(loess_spans$rmse == min(loess_spans$rmse))][1]
optimal_span

# expand DF to include one row for each day * demog var combination
polls_2020_demo <- 
  expand_grid(mid = seq.Date(min(polls_2020_demo$mid),
                             max(polls_2020_demo$mid),1),
              party = unique(polls_2020_demo$party)) %>%
  left_join(polls_2020_demo)

# nested smoothing for each category 
polls_2020_with_trends <-  lapply(unique(polls_2020_demo$party),
                                  function(z){
                                    polls_2020_demo %>% 
                                      dplyr::filter(party == z) %>%
                                      mutate(mid = as.numeric(mid)) %>%
                                      group_by(party) %>%
                                      arrange(mid) %>% 
                                      mutate(loess_span = optimal_span,
                                             loess_trend = predict(loess(share ~ mid, span = optimal_span,
                                                                         weights = sqrt(n/1818),
                                                                         control = loess.control(surface = "direct")),
                                                                   newdata=.),
                                             
                                             loess_trend = case_when(loess_trend > 1 ~ 1,
                                                                     loess_trend < 0 ~ 0,
                                                                     TRUE ~ loess_trend),
                                             loess_se = predict(loess(share ~ mid, span = optimal_span,
                                                                      weights = sqrt(n/1818),
                                                                      control = loess.control(surface = "direct")),
                                                                newdata=.,se=TRUE)$se.fit,
                                             loess_moe = sqrt(mean((loess_trend - share)^2,na.rm=T)) * 2
                                      ) %>%
                                      mutate(mid = as_date(mid))
                                  }) %>%
  bind_rows


# for remainder of time, inflate error by constant and regress towards mean -----------

test_date = ymd('2018-03-15')
t = 100
dat = polls_2020_with_trends
mean_weight_decay = -0.005 # very slow decay from last point + trend to mean
trend_weight_decay = -0.15


predict_from_nested_gam <- function(y, x, knots, weight){
  model = gam(y ~ s(x,k=knots))
  return(predict(model, tibble(x=x)))
}


# first, make a function to calculate the error in projections at a random time in the future
# GIVEN certain decays AND a test date
calculate_projection_error_at_time_t <- function(
  test_date = ymd('2019-06-10'),
  t = 100, 
  dat = polls_2020_with_trends,
  mean_weight_decay = -0.005, # very slow decay from last point + trend to mean
  trend_weight_decay = -0.15){ # very fast decay in share of trend for (last point + trend)
  
  # we're too close to the end of the series, return NULL
  if(isTRUE(ymd(test_date) + t >= ymd(max(dat$mid)))){
    return(NULL)
  }else{
    
    test_dat <- dat %>% 
      filter(ymd(mid) <= ymd(test_date),
             ymd(mid) >= (ymd(test_date) - 730))
    
    # first, expand ts to today
    test_dat <- 
      expand_grid(mid = seq.Date(min(test_dat$mid),
                                 max(test_dat$mid) + t,
                                 1),
                  party = unique(test_dat$party)) %>%
      left_join(test_dat, by = c("mid", "party")) %>%
      mutate(days_since_last_survey = as.numeric(mid - max(test_dat$mid)),
             to_impute = ifelse(days_since_last_survey > 0, TRUE, FALSE),
             mid = as_date(mid))
    
    test_dat <- lazy_dt(test_dat)
    
    # then calculate guess at how trend line would continue
    test_dat <- test_dat %>% 
      ungroup() %>%
      #dplyr::filter(party == z) %>%
      mutate(mid = as.numeric(as_date(mid))) %>%
      group_by(party) %>%
      arrange(party, mid) %>% 
      # as time goes on without polls, put less weight on trend and more on average
      mutate(
        #continued_loess_trend = na_kalman(loess_trend,smooth=T,model='auto.arima'),
        
        continued_loess_trend = predict_from_nested_gam(y=loess_trend, x=mid, knots = 50),
        
        rolling_mean_trend = rollapply(data=loess_trend,width=30*3,align='right',FUN=mean,partial=T),
        
        continued_loess_trend = case_when(continued_loess_trend < 0 ~ 0, 
                                          continued_loess_trend > 1 ~ 1,
                                          TRUE ~ continued_loess_trend),
        rolling_mean_trend = case_when(rolling_mean_trend < 0 ~ 0, 
                                       rolling_mean_trend > 1 ~ 1,
                                       TRUE ~ rolling_mean_trend),
        
        last_loess_trend = last(loess_trend[!is.na(loess_trend)]),
        last_rolling_mean_trend = last(rolling_mean_trend[!is.na(rolling_mean_trend)]),
        last_loess_moe = last(loess_moe[!is.na(loess_moe)]),
        
        weight_on_mean = pmax(0,1 - exp(days_since_last_survey * mean_weight_decay)),
        weight_on_trend = exp(days_since_last_survey * trend_weight_decay),
        weight_on_latest_poll = 1 - (weight_on_mean + weight_on_trend),
        
        weight_on_mean = to_impute * weight_on_mean,
        weight_on_trend = to_impute * weight_on_trend,
        weight_on_latest_poll = to_impute * weight_on_latest_poll,
        
        
        extrapolated_loess_trend = ifelse(to_impute,
                                          (last_rolling_mean_trend * weight_on_mean) + 
                                            (last_loess_trend * weight_on_latest_poll) +
                                            (continued_loess_trend * weight_on_trend),
                                          loess_trend),
        mid = as_date(mid)) %>%
      as.data.table()
    
    # compare the kalman trend and linear extrapolation with actual loess trends in the future
    if(F){
      gg <- test_dat %>% 
        as_tibble() %>%
        mutate(mid = as_date(mid)) %>%
        group_by(party) %>%
        filter(row_number() >= max(row_number()) - 500) %>%
        left_join(dat %>%       
                    filter(mid == max(test_dat$mid)) %>%
                    select(mid,party,actual_loess_trend = loess_trend), by = c("mid", "party")) %>%
        ggplot(., aes(x=mid,col=party)) + 
        geom_point(aes(y=share),alpha=0.2) + 
        geom_point(aes(y=actual_loess_trend),size=3) +
        geom_line(aes(y=loess_trend,linetype='1 actual trend')) +
        geom_line(aes(y=extrapolated_loess_trend,linetype='2 extrapolation')) +
        geom_line(aes(y=continued_loess_trend,linetype='3 continued trend')) +
        geom_line(aes(y=rolling_mean_trend,linetype='4 90-day rolling average'))
      
      print(gg)
    }
    
    # get error for every day in the projection series
    output <- test_dat %>% 
      filter(days_since_last_survey > 0) %>%
      select(-c(loess_trend)) %>%
      left_join(dat %>% select(mid,party,loess_trend), by = c("mid", "party")) %>% 
      ungroup() %>%
      mutate(error_in_extrapolation = extrapolated_loess_trend - loess_trend) %>%
      as.data.table()
    
    return(output %>% as_tibble)
  }
  
}

system.time(calculate_projection_error_at_time_t(test_date = "2021-01-01",t = 100) )


# this is a new function that will calculate summary statistics of errors for a RANDOM set of test dates from 2018 to 2021
calculate_projection_error_at_many_times_t <- function(
  num_tests = 8,
  lead_interval = 200,
  mean_weight_decay_optim = -0.005, # very slow decay from last point + trend to mean
  trend_weight_decay_optim = -0.15){
  
  errors <- pblapply(unique(polls_2020_with_trends$mid)[runif(num_tests,400,1400)],
                     cl = 14,
                     function(test_date_i){
                       calculate_projection_error_at_time_t(test_date = test_date_i,
                                                          t = lead_interval,
                                                          mean_weight_decay = mean_weight_decay_optim,
                                                          trend_weight_decay = trend_weight_decay_optim) 
                     }) %>%
    bind_rows 
  
  # hist(abs(errors$error_in_extrapolation))
  # print(sqrt(mean(errors$error_in_extrapolation^2)))
  return(errors)
  
  
}

system.time(
  simulated_errors <- calculate_projection_error_at_many_times_t(num_tests = 8,
                                                                 mean_weight_decay_optim = -0.005, # very slow decay from last point + trend to mean
                                                                 trend_weight_decay_optim = -0.15)
)

hist(abs(simulated_errors$error_in_extrapolation))
hist(abs(simulated_errors$error_in_extrapolation) / simulated_errors$loess_trend)

ggplot(simulated_errors[simulated_errors$days_since_last_survey != 1,], aes(x=abs(error_in_extrapolation))) +
  geom_density(aes(col=days_since_last_survey,group=days_since_last_survey ))

ggplot(sample_n(simulated_errors,1000), aes(x=loess_trend, y=abs(error_in_extrapolation))) +
  geom_point() + 
  geom_smooth(method='loess')

lm(abs(error_in_extrapolation) ~ days_since_last_survey + sqrt(days_since_last_survey) + loess_trend,
   data = simulated_errors) %>% summary


# now actually optimize weight parameters ---------------------------------

return_mae_given_weights <- function(parlist){
  
  simulated_errors <- calculate_projection_error_at_many_times_t(num_tests = 200,
                                                                 mean_weight_decay_optim = parlist[1],
                                                                 trend_weight_decay_optim = parlist[2])
  
  return(median(abs(simulated_errors$error_in_extrapolation)))
}

return_mae_given_weights(parlist = c(-0.5, -0.9))

# this will take more than a few minutes
if(FALSE){
  system.time(
    optimized_weights_on_trend_and_mean <- optim(par = c(-0.005, -0.15), 
                                                 fn = return_mae_given_weights,
                                                 control = list(maxit = 2000, abstol = 0.0001, # will stop after 500 iterations or gets RMSEs that aren't 0.01 percentage points different between runs
                                                                parscale = c(0.001, 0.01))
    )
  )
  beepr::beep(2)
  
  optimized_weights_on_trend_and_mean
}

optimized_mean_weight_decay = -0.003216095
optimized_trend_weight_decay = -0.150951232


# generate one set of historical predictions with the right weights -------
# so that we can get the margin of error of it
set_of_historical_sims <- simulated_errors <- calculate_projection_error_at_many_times_t(
  num_tests = 256,
  lead_interval = 200,
  mean_weight_decay_optim = optimized_mean_weight_decay, 
  trend_weight_decay_optim = optimized_trend_weight_decay)

set_of_historical_sims %>%
  group_by(days_since_last_survey) %>%
  summarise(rmse = sqrt(mean(error_in_extrapolation^2)),
            upper = quantile(abs(error_in_extrapolation),0.975),
            lower = quantile(abs(error_in_extrapolation),0.025)) %>%
  ggplot(., aes(x=days_since_last_survey)) +
  geom_line(aes(y=rmse)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),fill=NA,linetype=2,col='black')

set_of_historical_sims %>%
  group_by(days_since_last_survey) %>%
  summarise(rmse = sqrt(mean(error_in_extrapolation^2)),
            upper = quantile(abs(error_in_extrapolation),0.975),
            lower = quantile(abs(error_in_extrapolation),0.025)) %>%
  ggplot(., aes(x=days_since_last_survey)) +
  geom_line(aes(y=upper / rmse))

set_of_historical_sims %>%
  filter(days_since_last_survey > 150) %>%
  ggplot(., aes(x=abs(error_in_extrapolation),group=days_since_last_survey)) +
  geom_density()


lm(abs(error_in_extrapolation) ~ sqrt(days_since_last_survey)* loess_trend, data = set_of_historical_sims) %>% summary

get_predicted_sd <- function(days_since_last_survey, est_share){
  pmax(0,
       0.00316150 + 
         (sqrt(days_since_last_survey) * -0.00048204) + 
         (est_share * -0.03417641) +
         (sqrt(days_since_last_survey) *est_share * 0.01442846) 
       
  )
  
}

get_predicted_sd(set_of_historical_sims$days_since_last_survey, set_of_historical_sims$loess_trend)

# actually calculate with real MOE ----------------------------------------

# nested smoothing for each category 
# first, expand ts to today
polls_2020_with_trends_forecast <- 
  expand_grid(mid = seq.Date(min(polls_2020_with_trends$mid),
                             max(polls_2020_with_trends$mid) + as.numeric(ymd('2021-09-26') - Sys.Date()),
                             1),
              party = unique(polls_2020_with_trends$party)) %>%
  left_join(polls_2020_with_trends, by = c("mid", "party")) %>%
  mutate(days_since_last_survey = as.numeric(mid - max(polls_2020_with_trends$mid)),
         to_impute = ifelse(days_since_last_survey > 0, TRUE, FALSE),
         mid = as_date(mid))


# then calculate guess at how trend line would continue
polls_2020_with_trends_forecast <- polls_2020_with_trends_forecast %>% 
  ungroup() %>%
  #dplyr::filter(party == z) %>%
  mutate(mid = as.numeric(as_date(mid))) %>%
  group_by(party) %>%
  arrange(party, mid) %>% 
  # as time goes on without polls, put less weight on trend and more on average
  mutate(
    #continued_loess_trend = na_kalman(loess_trend,smooth=T,model='auto.arima'),
    
    continued_loess_trend = predict_from_nested_gam(y=loess_trend, x=mid, knots = 300),
    
    rolling_mean_trend = rollapply(data=loess_trend,width=30*3,align='right',FUN=mean,partial=T),
    
    continued_loess_trend = case_when(continued_loess_trend < 0 ~ 0, 
                                      continued_loess_trend > 1 ~ 1,
                                      TRUE ~ continued_loess_trend),
    rolling_mean_trend = case_when(rolling_mean_trend < 0 ~ 0, 
                                   rolling_mean_trend > 1 ~ 1,
                                   TRUE ~ rolling_mean_trend),
    
    last_loess_trend = last(loess_trend[!is.na(loess_trend)]),
    last_rolling_mean_trend = last(rolling_mean_trend[!is.na(rolling_mean_trend)]),
    last_loess_moe = last(loess_moe[!is.na(loess_moe)]),
    
    weight_on_mean = pmax(0,1 - exp(days_since_last_survey * optimized_mean_weight_decay)),
    weight_on_trend = exp(days_since_last_survey * optimized_trend_weight_decay),
    weight_on_latest_poll = 1 - (weight_on_mean + weight_on_trend),
    
    weight_on_mean = to_impute * weight_on_mean,
    weight_on_trend = to_impute * weight_on_trend,
    weight_on_latest_poll = to_impute * weight_on_latest_poll,
    
    
    extrapolated_loess_trend = ifelse(to_impute,
                                      (last_rolling_mean_trend * weight_on_mean) + 
                                        (last_loess_trend * weight_on_latest_poll) +
                                        (continued_loess_trend * weight_on_trend),
                                      loess_trend),
    extrapolated_loess_moe = ifelse(to_impute,
                                    last_loess_moe + get_predicted_sd(days_since_last_survey,extrapolated_loess_trend)*3,
                                    loess_moe),
    mid = as_date(mid)) 


# plot
plot(polls_2020_with_trends_forecast$share, polls_2020_with_trends_forecast$extrapolated_loess_trend)

gg1 <- polls_2020_with_trends_forecast %>%
  filter(party != 'Sonstige') %>%
  #filter(mid >= ymd("2020-01-01")) %>%
  ggplot(., aes(x = mid, col=party)) +
  geom_point(aes(y=share),alpha=0.1) +
  geom_line(aes(y=loess_trend))  +
  geom_line(aes(y=rolling_mean_trend),linetype=2) +
  labs(x='',y='',col='Party') +
  theme_minimal() +
  theme(legend.position = 'top')

gg2 <- polls_2020_with_trends_forecast %>%
  filter(party %in% c("CDU/CSU","GRÜNE","AfD","SPD","FDP")) %>%
  filter(mid >= ymd("2021-01-01")) %>%
  ggplot(., aes(x = mid, col=party)) +
  geom_point(aes(y=share),alpha=0.1) +
  geom_line(aes(y=rolling_mean_trend),linetype=2) +
  # geom_line(aes(y=continued_loess_trend),linetype=3) +
  geom_line(aes(y=extrapolated_loess_trend)) +
  geom_ribbon(aes(ymin=extrapolated_loess_trend - extrapolated_loess_moe, 
                  ymax=extrapolated_loess_trend + extrapolated_loess_moe,fill=party),
              col=NA,alpha=0.2) +
  labs(x='2021',y='') +
  theme_minimal() +
  theme(legend.position = 'none') +
  scale_x_date(date_breaks = '1 month',date_labels ='%b')

grid.arrange(gg1,gg2)
