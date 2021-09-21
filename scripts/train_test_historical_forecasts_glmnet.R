# package setup:
# lapply(c("tidyverse","zoo","lubridate","gridExtra","mgcv","pbapply","brms","mvtnorm","lqmm","imputeTS","glmnet","data.table","caret","doParallel","foreach","gamlss"),install.packages)
# remotes::install_github("stan-dev/cmdstanr")
# cmdstanr::install_cmdstan()

options(scipen = 999)
library(tidyverse)
library(zoo)
library(lubridate)
library(gridExtra)
library(mgcv)
library(pbapply)
library(brms)
library(cmdstanr)
library(mvtnorm)
library(lqmm)
library(imputeTS)
library(glmnet)
library(data.table)
library(caret)
library(doParallel)
library(foreach)
library(gamlss)
library(stringr)


redo_models_flag <- F

# 1. Generate polling trends for all years --------------------------------
# read the data
historical_polling <- read_csv('data/historical_polls.csv') %>%
  mutate(
    start = as_date(start)-1, end = as_date(end),
    start = as_date(ifelse(start > end, start - 365, start)),
    mid = as_date(start + round(as.numeric(end - start) / 2)),
    election_year = year(election_date),
    share = share / 100
  ) 


historical_polling <- historical_polling %>%
  group_by(election_year,poll) %>%
  mutate(n = ifelse(is.na(n), mean(n, na.rm = T), n)) %>%
  group_by(election_year) %>%
  mutate(n = ifelse(is.na(n), mean(n, na.rm = T), n)) %>%
  ungroup() %>%
  group_by(poll) %>%
  mutate(n = ifelse(is.na(n), mean(n, na.rm = T), n)) %>%
  ungroup() %>%
  mutate(n = ifelse(is.na(n), mean(n, na.rm = T), n)) 

mean(historical_polling$n)

# remove polls for parties at NA vote share
historical_polling <- historical_polling %>% filter(!is.na(share))

# expand DF to include one row for each day * party var combination
historical_polling_demo <-
  lapply(unique(historical_polling$election_year),
         function(elec_year){
           print(elec_year)
           expand_grid(
             mid = seq.Date(
               min(historical_polling[historical_polling$election_year == elec_year,]$mid),
               unique(historical_polling[historical_polling$election_year == elec_year,]$election_date), 1
             ),
             party = unique(historical_polling$party),
             election_year = elec_year,
             election_date = unique(historical_polling[historical_polling$election_year == elec_year,]$election_date), 
           ) %>%
             left_join(historical_polling)
         }) %>% bind_rows

# delete rows before each party's formation date
historical_polling_demo <- historical_polling_demo %>%
  filter(!(party == 'gru' & ymd(mid) < ymd('1980-01-01')),
         !(party == 'lin' & ymd(mid) < ymd('1990-01-01')),
         !(party == 'afd' & ymd(mid) < ymd('2013-01-01')),
  )



# OLD CODE for averaging:
if(F){
  historical_polling_with_trends <-
    pblapply(unique(historical_polling_demo$election_year),
             cl=16,
             function(y){
               
               # for each party...
               lapply(unique(historical_polling_demo[historical_polling_demo$election_year == y,]$party),
                      function(z) {
                        print(sprintf('%s - %s',y,z))
                        
                        # check if enough observations
                        historical_polling_year_demo <- historical_polling_demo %>%
                          dplyr::filter(election_year == y, party == z) %>%
                          mutate(mid = as.numeric(mid)) %>%
                          group_by(party) %>%
                          arrange(mid)
                        
                        # if empty....
                        if(all(is.na(historical_polling_year_demo$share))){
                          return(NULL)
                        }else if(length(na.omit(historical_polling_year_demo$share)) < 10){ # if not enough rows....
                          historical_polling_year_demo %>%
                            mutate(loess_trend = rep(mean(historical_polling_year_demo$share,na.rm=T),
                                                     length(historical_polling_year_demo$share)),
                                   loess_trend = case_when(loess_trend > 1 ~ 1,
                                                           loess_trend < 0 ~ 0,
                                                           TRUE ~ loess_trend),
                                   loess_moe = sqrt(mean((loess_trend - share)^2, na.rm = T)) * 2,
                                   mid = as_date(mid),
                                   k = NA) %>%
                            return
                          
                          
                        }else{# else!
                          # get optimal knots
                          gam_knots <- lapply(
                            seq(4,pmin(150, nrow(na.omit(historical_polling_year_demo))),2),
                            function(x) {
                              print(x)
                              rmse_in <- try((historical_polling_year_demo %>%
                                                dplyr::filter(election_year == y, party == z) %>%
                                                mutate(mid = as.numeric(mid)) %>%
                                                arrange(mid) %>%
                                                ungroup() %>%
                                                mutate(loess_trend = predict(gam(share ~ s(mid,k=x),
                                                                                 weights = sqrt(n / 1917),
                                                ),newdata = .)) %>%
                                                mutate(loess_trend = case_when(loess_trend > 1 ~ 1,
                                                                               loess_trend < 0 ~ 0,
                                                                               TRUE ~ loess_trend),
                                                       error = loess_trend - lead(share)) %>%
                                                pull(error) %>% na.omit() %>% .^2 %>% mean() %>% sqrt() ))
                              
                              if (!is.numeric(rmse_in)) {
                                return(NULL)
                              } else {
                                # print(rmse_in)
                                return(tibble(knots = x, rmse = rmse_in))
                              }
                            }
                          ) %>%
                            bind_rows() %>%
                            na.omit()
                          
                          gam_knots
                          optimal_knots <- gam_knots$knots[which(gam_knots$rmse == min(gam_knots$rmse))][1]
                          optimal_knots
                          
                          # and then draw trends
                          historical_polling_year_demo %>%
                            mutate(
                              loess_trend = predict(gam(share ~ s(mid,k=optimal_knots),
                                                        weights = sqrt(n / 1917)),
                                                    newdata = .),
                              loess_trend = case_when(loess_trend > 1 ~ 1,
                                                      loess_trend < 0 ~ 0,
                                                      TRUE ~ loess_trend),
                              loess_moe = sqrt(mean((loess_trend - share)^2, na.rm = T)) * 2) %>%
                            mutate(mid = as_date(mid),
                                   k = optimal_knots) %>%
                            return()
                        }
                        
                        
                      }) %>% bind_rows
             }) %>%
    bind_rows  %>%
    mutate(mid = as_date(mid)) %>%
    ungroup()
  
  # save
  historical_polling_with_trends <- historical_polling_with_trends %>%
    mutate(days_until_ed = as.numeric(election_date - mid)) %>%
    ungroup()
  
  # extract just the trends now
  historical_trends <- historical_polling_with_trends %>%
    dplyr::select(election_year, election_date, mid, days_until_ed, party, loess_trend, loess_moe,k)
  
  sum(any(duplicated(historical_trends)))
  
  historical_trends <- historical_trends[!duplicated(historical_trends),]
  
  historical_trends %>%
    write_csv('data/historical_trends.csv')
}else{
  historical_trends <- read_csv('data/historical_trends.csv')
}


# NEW CODE -- now parallelize by day
# but only do this if we haven't already
# get indexing vars
historical_polling_demo$days_until_ed = as.integer(historical_polling_demo$election_date - historical_polling_demo$mid)
historical_polling_demo$weeks_until_ed = floor(historical_polling_demo$days_until_ed / 7)

if(!any(grepl('historical_trends.csv',list.files('data')))){
  iter_dataset <- lapply(unique(historical_polling_demo$election_year),
                         function(y){
                           lapply(unique(historical_polling_demo[historical_polling_demo$election_year == y,]$party),
                                  function(z){
                                    days <- unique(historical_polling_demo[historical_polling_demo$election_year == y & historical_polling_demo$party == z,]$days_until_ed) 
                                    
                                    # only save data for two-ish years before election
                                    days = days[days <= 635]
                                    
                                    tibble(days_until_ed = days) %>%
                                      mutate(party = z, election_year = y) %>%
                                      return
                                  }) %>% bind_rows %>% return
                         }) %>% bind_rows %>% return
  
  nrow(iter_dataset)
  
  iter_dataset <- iter_dataset %>% mutate(rowid = paste0(election_year,party,days_until_ed))
  
  # only parallelize parties we haven't done yet
  if(length(list.files('data/historical_samples/')) == 0){
    iter_dataset
  }else{
    ids_done <- map_df(list.files('data/historical_samples/',full.names = T), read_csv) %>% 
      pull(rowid)
    nrow(iter_dataset)
    iter_dataset <- iter_dataset %>% filter(!rowid %in% ids_done)
    nrow(iter_dataset)
  }
  
  iter_rows <- sample(1:nrow(iter_dataset), replace = F) # , size = 5000) 
  
  
  # ...and on each day (in parallel):
  workers <- makeCluster(detectCores() - 2, outfile = '')
  registerDoParallel(workers)
  set.seed(1843)
  
  # chunk it
  num_chunks <- ceiling(nrow(iter_dataset) / 100)
  length_of_chunk <- round(nrow(iter_dataset)/num_chunks)
  for( chunk in 1:num_chunks-1){
    
    # here the data gets split up
    iter_chunk <- iter_rows[(chunk*length_of_chunk):(chunk*length_of_chunk+length_of_chunk)]
    
    #foreach with those 5 datarows
    # will take another few hour and a half
    chunk_results <- foreach(idx = 1:length(iter_chunk), 
                             .export = c('num_chunks','length_of_chunk','chunk','iter_chunk','iter_rows','iter_dataset','historical_polling_demo'),
                             .packages = c('mgcv','dplyr','lubridate'),
                             .combine = 'rbind'
    ) %dopar%
      {
        print("--------------")
        print(sprintf('chunk %s of %s (%s%%), idx %s/%s = %s%%',
                      chunk+1, num_chunks,round((chunk+1)/num_chunks*100),
                      idx,length(iter_chunk),round(idx/length(iter_chunk)*100)))
        
        row <- iter_chunk[idx]
        y <- iter_dataset[row,]$election_year
        z <- iter_dataset[row,]$party
        d <- iter_dataset[row,]$days_until_ed
        master_rowid <- iter_dataset[row,]$rowid
        
        print(sprintf('%s - %s - ED-%s days',y,z,d))
        
        
        # check if enough observations
        historical_polling_year_demo <- historical_polling_demo %>%
          dplyr::filter(election_year == y, party == z, days_until_ed >= d) %>%
          mutate(mid = as.numeric(mid)) %>%
          group_by(party) %>%
          arrange(mid)
        
        # if empty....
        if(all(is.na(historical_polling_year_demo$share))){
          output <- NULL
        }else if(length(na.omit(historical_polling_year_demo$share)) < 10){ # if not enough rows....
          output <- historical_polling_year_demo %>% 
            mutate(loess_trend = rep(mean(historical_polling_year_demo$share,na.rm=T),
                                     length(historical_polling_year_demo$share)),
                   loess_trend = case_when(loess_trend > 1 ~ 1, 
                                           loess_trend < 0 ~ 0,
                                           TRUE ~ loess_trend),
                   loess_moe = sqrt(mean((loess_trend - share)^2, na.rm = T)) * 2,
                   mid = as_date(mid),
                   k = NA_real_) %>%
            filter(days_until_ed == d) %>%
            dplyr::select(mid, party, election_year, election_date, poll, start, end, share, 
                          days_until_ed, weeks_until_ed, loess_trend, loess_moe, k) %>%
            mutate(rowid = master_rowid)
          
          
        }else{# else!
          # get optimal knots
          ks_to_iter_over <- round(seq(pmax(5, pmin(70,nrow(na.omit(historical_polling_year_demo)))-10),
                                       pmin(120, nrow(na.omit(historical_polling_year_demo))),
                                       
                                       (pmin(120, nrow(na.omit(historical_polling_year_demo))) - 
                                          pmax(5, pmin(70,nrow(na.omit(historical_polling_year_demo)))-10)) / 5))
          
          gam_knots <- lapply(ks_to_iter_over,
                              function(x) {
                                # print(x)
                                rmse_in <- try(silent = T,
                                               expr = (historical_polling_year_demo %>%
                                                         mutate(mid = as.numeric(mid)) %>%
                                                         arrange(mid) %>%
                                                         ungroup() %>% 
                                                         mutate(loess_trend = predict(gam(share ~ s(mid,k=x),
                                                                                          weights = sqrt(n / 1917),
                                                         ),newdata = .)) %>%
                                                         mutate(loess_trend = case_when(loess_trend > 1 ~ 1, 
                                                                                        loess_trend < 0 ~ 0,
                                                                                        TRUE ~ loess_trend),
                                                                error = loess_trend - lead(share)) %>%
                                                         pull(error) %>% na.omit() %>% .^2 %>% mean() %>% sqrt() ))
                                if (!is.numeric(rmse_in)) {
                                  return(NULL)
                                } else {
                                  # print(rmse_in)
                                  return(tibble(knots = x, rmse = rmse_in))
                                }
                              }
          ) %>%
            bind_rows() %>%
            na.omit()
          
          if(nrow(gam_knots) == 0){
            optimal_knots = pmin(100, nrow(na.omit(historical_polling_year_demo))/2)
          }else{
            gam_knots
            optimal_knots <- mean(gam_knots$knots[which(gam_knots$rmse == min(gam_knots$rmse))],na.rm=T)
            optimal_knots
          }
          
          optimal_knots
          
          # and then draw trends
          historical_polling_year_demo_trends <- historical_polling_year_demo %>% 
            mutate(
              loess_trend = predict(gam(share ~ s(mid,k=optimal_knots),
                                        weights = sqrt(n / 1917)),
                                    newdata = .),
              loess_trend = case_when(loess_trend > 1 ~ 1, 
                                      loess_trend < 0 ~ 0,
                                      TRUE ~ loess_trend),
              loess_moe = sqrt(mean((loess_trend - share)^2, na.rm = T)) * 2) %>%
            mutate(mid = as_date(mid),
                   k = optimal_knots) 
          
          
          # return that day's value
          output <- historical_polling_year_demo_trends %>% 
            filter(days_until_ed == d) %>%
            dplyr::select(mid, party, election_year, election_date, poll, start, end, share, 
                          days_until_ed, weeks_until_ed, loess_trend, loess_moe, k) %>%
            mutate(rowid = master_rowid)
        }
        
        output
        
        return(output)
        
      }
    
    
    # save your foreach results and then begin again
    if (length(list.files('data/historical_samples/')) == 0){
      write_csv(chunk_results, file= "data/historical_samples/chunked_historical_data.csv")
    }else{
      write_csv(chunk_results, file="data/historical_samples/chunked_historical_data.csv", append=TRUE, col_names = FALSE)
    }
    
  }
  
  stopCluster(workers)
  beepr::beep(2)
  
  # read in the full data
  daily_averages_knots <- read_csv('data/historical_samples/chunked_historical_data.csv')
  
  historical_polling_with_trends <- daily_averages_knots %>%
    arrange(election_year, party, mid)
  
  # extract just the trends now
  historical_trends <- historical_polling_with_trends %>%
    dplyr::select(election_year, election_date, mid, days_until_ed, party, loess_trend, loess_moe, k)
  
  sum(any(duplicated(historical_trends)))
  
  historical_trends <- historical_trends[!duplicated(historical_trends),]
  
  historical_trends %>%
    write_csv('data/historical_trends.csv')
}else{
  historical_trends <- read_csv('data/historical_trends.csv')
}


# look at the averages for all years
pdf('polls_all_years.pdf',width=14,height=10)
ggplot(historical_trends, aes(x=days_until_ed, col=as.factor(party))) +
  geom_line(aes(y = loess_trend)) + 
  geom_point(data = historical_polling %>% mutate(days_until_ed = election_date - date), aes(y=share),
             alpha=0.5,size=0.2) +
  facet_wrap(~election_year,scales = 'free_x')  +
  theme_minimal() +
  geom_vline(xintercept=300,col='gray40') +
  geom_vline(xintercept=200,col='gray60') +
  scale_x_reverse(limits=c(365,0)) 
dev.off()

# 2. Loop generate LOOCV structural prior for all years -------------------------
# FORMULA IS vote in year t ~ vote share [t-1] + polls in january + PMs party dummy

# attach actual vote shares
prior <- read_csv('data/prior.csv') %>%
  rename(election_year = year) %>%
  mutate(voteshare = voteshare / 100)

# get lagged vote shares
prior <- prior %>% 
  group_by(party) %>%
  mutate(lag_voteshare = lag(voteshare)) %>%
  ungroup() 

# add polls at 200 days before the election
prior <- prior %>%
  left_join(historical_trends[historical_trends$days_until_ed==300,] %>%
              dplyr::select(election_year, party, polls_300_days_out = loess_trend)) %>%
  left_join(historical_trends[historical_trends$days_until_ed==200,] %>%
              dplyr::select(election_year, party, polls_200_days_out = loess_trend))

sqrt(mean((prior$voteshare - prior$polls_300_days_out)^2,na.rm=T))
sqrt(mean((prior$voteshare - prior$polls_200_days_out)^2,na.rm=T))

# handle new parties
prior <- prior %>%
  # handle new parties
  mutate(lag_voteshare = case_when(party == 'afd' & election_year == 2013 ~ 0.001,
                                   party == 'gru' & election_year == 1980 ~ 0.001,
                                   party == 'lin' & election_year == 1990 ~ 0.001,
                                   TRUE ~ lag_voteshare),
         polls_300_days_out = case_when(party == 'afd' & election_year == 2013 ~ 0.001,
                                        party == 'gru' & election_year == 1980 ~ 0.001,
                                        party == 'lin' & election_year == 1990 ~ 0.001,
                                        TRUE ~ polls_300_days_out),
         polls_200_days_out = case_when(party == 'afd' & election_year == 2013 ~ 0.001,
                                        party == 'gru' & election_year == 1980 ~ 0.001,
                                        party == 'lin' & election_year == 1990 ~ 0.001,
                                        TRUE ~ polls_200_days_out))



# EDA for model covariates
gg1 <- ggplot(prior, aes(x=lag_voteshare*100,y=voteshare*100)) +
  geom_point() +
  geom_abline() +
  geom_smooth(method='lm',linetype=2,se=F) +
  theme_minimal()  +
  labs(subtitle='Germany, vote shares by party, 1953-2017',
       x='Share of second votes in previous election, %',
       y='Share of second votes, %')


gg2 <- ggplot(prior, aes(x=polls_300_days_out*100,y=voteshare*100)) +
  geom_point() +
  geom_abline() +
  geom_smooth(method='lm',linetype=2,se=F) +
  theme_minimal() +
  labs(subtitle=' ',
       x='Share of second votes in polls 300 days before election, %',
       y='Share of second votes, %')


residualized_swing_model <- lm(I(voteshare - lag_voteshare) ~ lag_voteshare + polls_300_days_out,
                               data = prior)

prior %>%
  mutate(baseline_swing = predict(residualized_swing_model,.),
         residual_error = (voteshare - lag_voteshare) - baseline_swing) %>%
  lm(residual_error ~ chancellor_party*term, data =.) %>% summary


gg3 <- prior %>%
  mutate(baseline_swing = predict(residualized_swing_model,.)) %>%
  ggplot(., aes(x=chancellor_party,y=((voteshare - lag_voteshare) - baseline_swing)*100)) +
  geom_point() +
  geom_smooth(method='lm',linetype=2,se=F) +
  scale_x_continuous(expand = c(0.3,0.3),breaks=c(0,1)) +
  theme_minimal() +
  labs(subtitle=' ',
       x="Indicator for whether party\nis the chancellor's party, %",
       y='Swing in share of second votes\nResidualized, percentage points')

gg4 <-  prior %>%
  mutate(baseline_swing = predict(residualized_swing_model,.)) %>%
  filter(gov != 0) %>%
  ggplot(., aes(x=gov*term,y=((voteshare - lag_voteshare) - baseline_swing)*100)) +
  geom_point() +
  geom_smooth(method='lm',linetype=2,se=F) +
  theme_minimal() +
  labs(subtitle=' ',
       x='Number of years party has been\npart of governing coalition, %',
       y='Swing in share of second votes\nResidualized, percentage points')


grid.arrange(gg1,gg2,gg3,gg4,ncol=2)


# looping to train stan models for heldout predictions
validation_years <- c(1953, 1957, 1961, 1965, 1969, 1972, 1976, 1980,
                      1983, 1987, 1990, 1994, 1998, 2002, 2005, 2009,
                      2013, 2017)

prior <- prior %>%
  mutate(chancellor_party = as.character(chancellor_party),
         gov = as.character(gov),
         prop_poll300_diff_from_lag_vote = abs(polls_300_days_out - lag_voteshare)/lag_voteshare)

plot(prior$prop_poll300_diff_from_lag_vote, abs(prior$lag_voteshare - prior$voteshare))

# model result as function of normal distribution
## mean is predicted using lag vote share, polls at ED-300, and structural party info
## sigma is predicted using polls at ED-300 and lag vote share
normal_prior_model <- brm(bf(voteshare | trunc(lb=0,ub=1) ~ lag_voteshare + polls_300_days_out + chancellor_party + chancellor_party:term,
                             sigma | trun(lb=0,ub=0.2) ~ polls_300_days_out + lag_voteshare ), # + chancellor_party + chancellor_party:term),
                          data = prior,
                          family = 'normal',
                          prior = c(set_prior('normal(0.01, 0.05)',class='Intercept'),
                                    set_prior('normal(0.23,  0.05)',class='b',coef='lag_voteshare'),
                                    set_prior('normal(0.69,  0.05)',class='b',coef='polls_300_days_out'),
                                    set_prior('normal(0.03, 0.02)',class='b',coef='chancellor_party1'),
                                    set_prior('normal(0, 0.01)',class='b',coef='chancellor_party1:term'),
                                    set_prior('normal(-0.01, 0.03)',class='b',coef='chancellor_party0:term'),
                                    # prior for nu 
                                    # set_prior('normal(15,  2)',class='nu'),
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


loo(normal_prior_model,student_prior_model)

mcmc_plot(student_prior_model)
plot(student_prior_model,'b_Intercept')
pp_check(student_prior_model,nsamples = 100)


# look at in-sample fit
prior %>% 
  na.omit() %>%
  mutate(normal = predict(student_prior_model,newdaata=.,probs = c(0.025,0.975))[,1],
         lower = predict(student_prior_model,newdaata=.,probs = c(0.025,0.975))[,3],
         upper = predict(student_prior_model,newdaata=.,probs = c(0.025,0.975))[,4]) %>%
  summarise(rmse = sqrt(mean((normal - voteshare)^2)),
            pct_below_lower = mean(voteshare < lower),
            pct_above_upper = mean(voteshare > upper))

# this loop trains a bunch of models
if(redo_models_flag){
  lapply(validation_years,
         function(holdout_year){
           print(sprintf('running holdout predictions for %s (max year %s)',holdout_year,max(validation_years)))
           prior_model_holdout <- brm(bf(voteshare | trunc(lb=0,ub=1) ~ lag_voteshare + polls_300_days_out + chancellor_party + chancellor_party:term,
                                         sigma | trun(lb=0,ub=0.2) ~ polls_300_days_out + lag_voteshare ), #+ chancellor_party + chancellor_party:term),
                                      data = prior[prior$election_year != holdout_year,],
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
           
           
           write_rds(prior_model_holdout,sprintf('data/training_models/prior_%s.rds',holdout_year),compress = 'gz')
         })
}


# this loop uses the models trained above to make predictions
prior_with_predictions <- lapply(validation_years,
                                 function(holdout_year){
                                   print(holdout_year)
                                   preds <- predict(object = read_rds(sprintf('data/training_models/prior_%s.rds',holdout_year)), 
                                                    newdata = prior %>% filter(election_year == holdout_year),
                                                    probs = c(0.025,0.975)
                                   )
                                   
                                   prior %>% 
                                     filter(election_year == holdout_year) %>%
                                     mutate(est = preds[,1],
                                            se = preds[,2],
                                            lower = preds[,3],
                                            upper = preds[,4]) %>%
                                     return
                                   
                                 }) %>% bind_rows

# errors normally distributed or no?
ggplot(prior_with_predictions, aes(x=voteshare-est)) +
  geom_histogram(binwidth=0.001) +
  theme_minimal()

prior_with_predictions %>% # this code leads me to believe, yes, they are fat-failed
  mutate(overall_loo_prior_rmse = sqrt(mean((est - voteshare)^2)),
         outside_95ci = 
           ( voteshare > (est + overall_loo_prior_rmse*1.96) ) |
           ( voteshare < (est - overall_loo_prior_rmse*1.96) ),
         outside_99ci = 
           ( voteshare > (est + overall_loo_prior_rmse*2.576) ) |
           ( voteshare < (est - overall_loo_prior_rmse*2.576) ),
         outside_99.5ci = 
           ( voteshare > (est + overall_loo_prior_rmse*2.807) ) |
           ( voteshare < (est - overall_loo_prior_rmse*2.807) ),
         outside_99.9ci = 
           ( voteshare > (est + overall_loo_prior_rmse*3.291) ) |
           ( voteshare < (est - overall_loo_prior_rmse*3.291) )
  ) %>%
  summarise(outside_95ci = mean(outside_95ci),
            outside_99ci = mean(outside_99ci),
            outside_99.5ci = mean(outside_99.5ci),
            outside_99.9ci = mean(outside_99.5ci)
  )


prior_with_predictions %>% # actual 95% ci probs
  mutate(outside_95ci = 
           ( voteshare > upper ) |
           ( voteshare < lower )
  ) %>%
  summarise(outside_95ci = mean(outside_95ci)
  )


plot(prior_with_predictions$est, (prior_with_predictions$upper - prior_with_predictions$est) / prior_with_predictions$est)


# let's look at model accuracy
ggplot(prior_with_predictions, aes(x=est,y=voteshare,col=party,group=1)) +
  geom_point() +
  geom_segment(aes(yend=voteshare,x=lower,xend=upper)) +
  geom_abline() +
  geom_smooth(method='lm',linetype=2) +
  theme_minimal()

prior_with_predictions %>%
  group_by(party) %>%
  summarise(avg_voteshare = mean(voteshare),
            rmse = sqrt(mean((est - voteshare)^2))) %>%
  arrange(desc(rmse))


priors_to_attach <- prior_with_predictions %>% 
  dplyr::select(election_year,party,prior_est=est,prior_se=se,voteshare)


# 3. Loop train and test daily models excluding validation year ---------------------
# Result ~ prior + polls

# look
historical_trends %>% 
  filter(days_until_ed <= 250) %>% 
  ggplot(., aes(x=days_until_ed, col=as.factor(party))) +
  geom_line(aes(y = loess_trend)) + 
  facet_wrap(~election_year,scales = 'free_x')  +
  theme_minimal() +
  geom_vline(xintercept=300,col='gray40') +
  geom_vline(xintercept=200,col='gray60') +
  scale_x_reverse()

# get historical polling error at day T
# with adjustments for party share of the vote
poll_errors_data <- historical_trends %>% 
  filter(days_until_ed <= 250) %>%
  left_join(priors_to_attach, by = c("election_year", "party")) %>% 
  mutate(abs_error = abs(loess_trend - voteshare)) %>%
  group_by(vote_bucket = cut(loess_trend,c(0,0.1,0.2,0.3,0.4,1)), days_until_ed) %>%
  summarise(loess_trend = mean(loess_trend),
            rmse = sqrt(mean(abs_error^2))) %>%
  group_by(vote_bucket) %>%
  mutate(rmse = pmax(rmse,first(rmse)))


loo_poll_error_function <- function(vote, days){
  0.016 + sqrt(days)*0.0005 + pmin((vote/0.4),1)*0.005 + sqrt(days)*0.0015*pmin((vote/0.4),1)
}

ggplot(poll_errors_data, aes(x=loess_trend,y=rmse,col=days_until_ed)) + geom_point() + geom_smooth(method='lm',aes(group=1),formula=y~sqrt(x),fullrange=T)
ggplot(poll_errors_data, aes(x=days_until_ed,y = rmse,col=vote_bucket)) + geom_line() + geom_smooth(method='lm',formula=y~sqrt(x),fullrange=T) +
  geom_line(dat = tibble(days_until_ed = 0:250,
                         rmse = loo_poll_error_function(0.4,0:250)),
            aes(col='test'),col='black')

# first, combine trends and priors (for election-year dates only) 
##  be sure to get the sd from the historical polling average
historical_trends_and_priors <- historical_trends %>% 
  dplyr::select(-loess_moe) %>%
  left_join(priors_to_attach, by = c("election_year", "party"))%>% 
  filter(days_until_ed <= 1095) %>%
  mutate(loess_se = loo_poll_error_function(vote = loess_trend, days = days_until_ed)) %>%
  dplyr::select(election_year, election_date, mid, days_until_ed, party, voteshare,
                loess_trend, loess_se, prior_est, prior_se)

ggplot(historical_trends_and_priors , aes(x=days_until_ed,loess_se,col=party,group=paste(party,election_year))) +
  geom_line() +
  scale_y_continuous(breaks=seq(0,1,0.005)) +
  labs(subtitle='estimate standard deviation of polling average, year-party pairs')


# add standard deviation
historical_trends_and_priors <- historical_trends_and_priors %>%
  group_by(party,election_year) %>%
  arrange(election_year,party,desc(days_until_ed)) %>%
  mutate(sd_polls_last_year = rollapply(data=loess_trend,width=365,FUN=sd,partial=T,fill=NA,align='right'),
         sd_polls_last_year = pmax(0.005,sd_polls_last_year)) %>%
  # fill in gaps
  arrange(election_year,party,desc(days_until_ed)) %>%
  mutate(sd_polls_last_year = imputeTS::na_locf(sd_polls_last_year)) %>%
  arrange(election_year,party,(days_until_ed)) %>%
  mutate(sd_polls_last_year = imputeTS::na_locf(sd_polls_last_year)) %>%
  ungroup() 


# test predictions for holdout years, with a cross-validated ridge regression (cv.glnet)
if(any(grepl('historical_prior_data_posterior_l',list.files('data/')))){
  historical_prior_data_posterior_l <- read_rds('data/historical_prior_data_posterior_l.rds')
} else{
  # rm(historical_prior_data_posterior_l)
  workers <- makeCluster(detectCores() /2, outfile = '')
  registerDoParallel(workers)
  set.seed(1843)
  
  historical_prior_data_posterior_l <- foreach(holdout_year = unique(historical_trends_and_priors$election_year),
                                               .packages = c('glmnet','dplyr','purrr','mgcv','ggplot2'),
                                               .combine = 'rbind') %dopar%
    {
      # holdout_year <- 2017; x <- 13;
      
      print(sprintf('blending polls and prior for %s...',holdout_year))
      
      # this runs a separate glmnet model for every day in past campaign
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
      
      Sys.sleep(2)
      
      
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
      
      # intercept_smooth_coef <- gam(intercept.c ~ s(days_until_ed, k = 50), data = historical_coefs_temp)
      # 
      # prior_smooth_coef <- gam(prior_est.c ~ intercept.c + s(days_until_ed, k = 50), data = historical_coefs_temp)
      # poll_smooth_coef <- gam(loess_trend.c ~ intercept.c + prior_est.c, data = historical_coefs_temp) 
      # 
      # poll_se_coef <- gam(loess_se.c ~ intercept.c + prior_est.c + loess_trend.c + s(days_until_ed, k = 50), data = historical_coefs_temp) 
      # prior_se_coef <- gam(prior_se.c ~ intercept.c + prior_est.c + loess_trend.c + loess_se.c, data = historical_coefs_temp) 
      # 
      # poll_trend_poll_se.c <- gam(loess_trend_loess_se.c ~ intercept.c + prior_est.c + prior_se.c + loess_trend.c + loess_se.c, data = historical_coefs_temp) 
      # prior_trend_prior_se.c <- gam(prior_est_prior_se.c ~ intercept.c + prior_est.c + prior_se.c + loess_trend.c + loess_se.c + loess_trend_loess_se.c, data = historical_coefs_temp) 
      # prior_poll_smooth_coef <- gam(loess_trend_prior_est.c ~ intercept.c + prior_est.c + loess_trend.c,data = historical_coefs_temp) 
      
      
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
      test_data <- historical_trends_and_priors %>%
        filter(election_year == holdout_year) %>%
        group_by(election_year, party, days_until_ed) %>%
        summarise(prior_se = unique(prior_se),
                  loess_se = unique(loess_se),
                  prior_est = unique(prior_est),
                  loess_trend = unique(loess_trend),
                  voteshare = unique(voteshare),
                  sd_polls_last_year = unique(sd_polls_last_year)) %>%
        ungroup() 
      
      
      test_data <- test_data %>%
        left_join(historical_coefs_temp, by = c('days_until_ed')) %>%
        mutate(bayes_pred = ((prior_se^2/(prior_se^2+loess_se^2))*prior_est) + ((loess_se^2/(prior_se^2+loess_se^2))*loess_trend),
               bayes_pred_se = sqrt(((1/(prior_se^2)) + (1/(loess_se^2)))^-1) ,
               lm_pred = 
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
        ungroup()
      
      test_data %>% 
        group_by(days_until_ed) %>%
        summarise(rmse = sqrt(mean(c(voteshare - lm_pred)^2)))
      
      return(test_data)
      
    }
  
  beepr::beep(2)
  stopCluster(workers)
  
  write_rds(historical_prior_data_posterior_l,'data/historical_prior_data_posterior_l.rds',compress = 'gz')
}



# bind the data from individual years
historical_prior_data_posterior <- historical_prior_data_posterior_l

as.data.frame(head(historical_prior_data_posterior)) 


# add some variables on context in polling, then trim to 730 days til ED (2 years out)
historical_prior_data_posterior <- historical_prior_data_posterior %>%
  filter(days_until_ed <= 635) %>%
  group_by(party,election_year) %>%
  mutate(error_lm_pred = abs(voteshare - lm_pred)) %>%
  arrange(election_year,party,desc(days_until_ed)) %>%
  ungroup()

# check which model is better. (the linear model is better)
historical_prior_data_posterior %>% 
  filter(days_until_ed %in% c(0,14,30,60,90,200,300)) %>%
  group_by(days_until_ed) %>%
  summarise(bayes_rmse = sqrt(mean((bayes_pred - voteshare)^2)),
            lm_rmse = sqrt(mean((lm_pred - voteshare)^2)),
            prior_rmse = sqrt(mean((prior_est - voteshare)^2)),
            loess_rmse = sqrt(mean((loess_trend - voteshare)^2))) 


# graph rmse by day for all models
historical_prior_data_posterior %>% 
  group_by(days_until_ed) %>%
  summarise(bayes_rmse = sqrt(mean((bayes_pred - voteshare)^2)),
            lm_rmse = sqrt(mean((lm_pred - voteshare)^2)),
            prior_rmse = sqrt(mean((prior_est - voteshare)^2)),
            loess_rmse = sqrt(mean((loess_trend - voteshare)^2))) %>%
  ggplot(.,aes(x=days_until_ed)) +
  geom_line(aes(y=bayes_rmse,col='bayes'),linetype=1) +
  geom_line(aes(y=lm_rmse,col='lm'),linetype=1) +
  geom_line(aes(y=loess_rmse,col='poll'),linetype=2) +
  geom_line(aes(y=prior_rmse,col='prior'),linetype=2) +
  scale_x_reverse() +
  scale_y_continuous(limits=c(0.01,0.05),breaks=seq(0,1,0.01)) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())  +
  labs(subtitle='errors in predictions')


# look at different projections and results in eg 2017
historical_prior_data_posterior %>% 
  filter(election_year == 2017) %>%
  mutate(party = as.factor(party),
         party = fct_reorder(party, -voteshare, mean)) %>%
  ggplot(.,aes(x=days_until_ed)) +
  geom_line(aes(y=bayes_pred,col='bayes_pred'),linetype=1) +
  geom_line(aes(y=lm_pred,col='lm_pred'),linetype=1) +
  geom_line(aes(y=loess_trend,col='poll'),linetype=2) +
  geom_line(aes(y=prior_est,col='prior'),linetype=2) +
  scale_x_reverse() +
  geom_point(aes(x=0,y=voteshare),size=2) +
  scale_y_continuous(expand=c(0.02,0.02),breaks=seq(0,1,0.02)) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  facet_wrap(~party,scales='free_y')


# graph rmse by party

# graph rmse by day for all models
historical_prior_data_posterior %>% 
  group_by(days_until_ed,party) %>%
  summarise(lm_rmse = sqrt(mean((lm_pred - voteshare)^2))) %>%
  ggplot(.,aes(x=days_until_ed)) +
  geom_line(aes(y=lm_rmse,col=party),linetype=1) +
  scale_x_reverse() +
  scale_y_continuous(limits=c(0.01,0.06),breaks=seq(0,1,0.01)) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())  +
  labs(subtitle='errors in predictions')



# 4. Fit sigma of error distribution, with mu = glmnet prediction ------------
# train one master model for saving later

historical_prior_data_posterior <- historical_prior_data_posterior %>%
  left_join(
    historical_prior_data_posterior %>%
      group_by(election_year) %>%
      summarise(rmse_poll_error_cycle = sqrt(mean((voteshare - loess_trend)^2))) %>%
      mutate(rmse_poll_error_last_three_cycles = 
               (lag(rmse_poll_error_cycle) +
                  lag(rmse_poll_error_cycle,2) + 
                  lag(rmse_poll_error_cycle,3) ) / 3) %>% 
      mutate(rmse_poll_error_last_three_cycles = na_locf(rmse_poll_error_last_three_cycles)) 
  ) %>% 
  mutate(party_is_cdu = party == 'cdu')


# train one model for 2021 later
err_mod_training_data <- historical_prior_data_posterior 
err_mod <- gamlss(formula = voteshare ~ lm_pred, 
                  sigma.formula =  ~lm_pred*poly(days_until_ed, 2) + party_is_cdu*poly(days_until_ed, 2) + rmse_poll_error_last_three_cycles,
                  nu.formula =  ~lm_pred*poly(days_until_ed,2) + party_is_cdu*poly(days_until_ed,2),
                  family = "TF", # TK - test ST1
                  data = err_mod_training_data,
                  control = gamlss.control(n.cyc = 500))


historical_prior_data_posterior %>% 
  mutate(lm_fitted_error = fitted(err_mod,'sigma')) %>%
  group_by(days_until_ed, party) %>%
  summarise(lm_rmse = sqrt(mean((lm_pred - voteshare)^2)),
            lm_fitted_error = mean(lm_fitted_error)) %>%
  ungroup() %>%
  # graph
  ggplot(.,aes(x=days_until_ed,col=party)) +
  geom_line(aes(y=lm_rmse,linetype='historical rmse')) +
  geom_line(aes(y=lm_fitted_error,linetype='predicted error')) +
  scale_x_reverse() +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())  +
  labs(subtitle='errors in predictions')


historical_prior_data_posterior %>% 
  mutate(lm_fitted_df = fitted(err_mod,'nu')) %>%
  group_by(days_until_ed) %>% 
  summarise(lm_fitted_df = mean(lm_fitted_df)) %>% 
  ggplot(.,aes(x=days_until_ed, y=lm_fitted_df)) + geom_line() + scale_y_log10(limits=c(1,2000))


# now we want to parallelize over each year, witholding it for calculating test-set prediction
if(any(grepl('historical_prior_data_posterior_m',list.files('data/')))){
  historical_prior_data_posterior_m <- read_rds('data/historical_prior_data_posterior_m.rds')
} else{
  
  workers <- makeCluster(detectCores() /2, outfile = '')
  registerDoParallel(workers)
  set.seed(1843)
  
  historical_prior_data_posterior_m <- foreach(test_year = unique(historical_prior_data_posterior$election_year),
                                               .packages = c('gamlss','dplyr','purrr','mgcv','ggplot2'),
                                               .combine = 'rbind') %dopar% 
    {
      # test_year <- 2013
      print(sprintf("Fitting LOO sigma models to predict on %s...",test_year))
      
      train_validate_data <- historical_prior_data_posterior %>%
        filter(election_year != test_year) %>%
        mutate(party_is_cdu = party=='cdu')
      
      # TK- optimize the gamlss with a loop for AIC, we want it to be as parsimonious as possible
      temp_error_mod <- gamlss(formula = voteshare ~ lm_pred, 
                               sigma.formula =  ~lm_pred*poly(days_until_ed, 2) + party_is_cdu*poly(days_until_ed, 2) + rmse_poll_error_last_three_cycles,
                               nu.formula =  ~lm_pred*poly(days_until_ed,2) + party_is_cdu*poly(days_until_ed,2),
                               family = "TF", # TK - test ST1
                               data = train_validate_data,
                               control = gamlss.control(n.cyc = 500))
      
      
      # predict on holdout year:
      test_data <- historical_prior_data_posterior %>%
        filter(election_year == test_year) %>%
        ungroup() %>%
        mutate(lm_rmse = sqrt(mean((lm_pred - voteshare)^2)),
               party_is_cdu = party=='cdu')
      
      ## sigmas
      test_data_sigmas <- predict(object = temp_error_mod, 
                                  data = train_validate_data, 
                                  newdata = as.data.frame(test_data %>% dplyr::select(lm_pred,days_until_ed,party_is_cdu,rmse_poll_error_last_three_cycles)), 
                                  what = 'sigma',type = 'response')
      
      ## nus
      test_data_nus <- predict(object = temp_error_mod, 
                               data = train_validate_data, 
                               newdata = as.data.frame(test_data %>% dplyr::select(lm_pred,days_until_ed,party_is_cdu,rmse_poll_error_last_three_cycles)), 
                               what = 'nu',type = 'response')
      
      # add
      test_data %>%
        ungroup() %>%
        mutate(lm_fitted_error = pmax(0.01,test_data_sigmas),
               lm_fitted_df = pmax(2, test_data_nus)) %>%
        return
      
      
    }
  
  stopCluster(workers)
  beepr::beep(2)
  
  write_rds(historical_prior_data_posterior_m,'data/historical_prior_data_posterior_m.rds',compress = 'gz')
  
}

# bind together
if(exists('historical_prior_data_posterior_param')){rm(historical_prior_data_posterior_param)}
historical_prior_data_posterior_param <- historical_prior_data_posterior_m
historical_prior_data_posterior_param 


# check accuracy of sigma prediction
historical_prior_data_posterior_param %>% 
  group_by(days_until_ed, party) %>%
  summarise(lm_rmse = sqrt(mean((lm_pred - voteshare)^2)),
            lm_fitted_error = mean(lm_fitted_error)) %>%
  ungroup() %>%
  # graph
  ggplot(.,aes(x=days_until_ed,col=party)) +
  geom_line(aes(y=lm_rmse,linetype='historical rmse')) +
  geom_line(aes(y=lm_fitted_error,linetype='predicted error')) +
  scale_x_reverse() +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())  +
  labs(subtitle='errors in predictions')

# check nu
historical_prior_data_posterior_m %>% 
  group_by(days_until_ed) %>% 
  summarise(lm_fitted_df = mean(lm_fitted_df))


# OK, we are done. now check how often the result is inside the prediction interval
historical_prior_data_posterior_param %>%
  summarise(
    inside_99_ci = mean(voteshare < (lm_pred + lm_fitted_error*2.576) & 
                          voteshare > (lm_pred - lm_fitted_error*2.576),na.rm=T),
    inside_95_ci = mean(voteshare < (lm_pred + lm_fitted_error*1.96) & 
                          voteshare > (lm_pred - lm_fitted_error*1.96),na.rm=T),
    inside_80_ci = mean(voteshare < (lm_pred + lm_fitted_error*1.282) & 
                          voteshare > (lm_pred - lm_fitted_error*1.282),na.rm=T)
  )



# 5. Extract correlations and generate predictions ------------------------

# function for making forecasts for year y at t days before the election
simulate_election_y_from_time_t <- function(y, t, num_sims = 10000, return_all = FALSE){
  # y <- 2017; t <- 100; num_sims <- 10000; return_all <- FALSE
  print(sprintf('forecasts for %s made at ED-%s',y,t))
  
  # get correlation between parties UP TO THAT DAY
  election_year_data_full <- historical_prior_data_posterior_param %>% 
    dplyr::filter(election_year == y, days_until_ed >= t)
  
  election_year_data_full %>%
    filter(party %in% c('afd','cdu','spd','gru')) %>%
    ggplot(., aes(x=days_until_ed, y=loess_trend,col=party)) + geom_line() + coord_cartesian(xlim=c(0,635)) + facet_wrap(~party,scales = 'free_y') 
  
  
  election_year_data_full %>%
    mutate(election_year_day = paste0(election_year, days_until_ed)) %>%
    arrange(days_until_ed) %>%
    dplyr::select(election_year_day,party,loess_trend) %>%
    spread(party,loess_trend) %>%
    # mutate_at(unique(election_year_data_full$party),
    #           function(x){ pmax(-0.01,pmin(0.01,(x-lag(x)))  ) }) %>%
    # na.omit() %>%
    # mutate_at(unique(election_year_data_full$party),
    #           function(x){ 
    #             case_when(sd(x) == 0 | is.na(sd(x)) ~ x,
    #                       abs((x - mean(x))/sd(x) ) > 1.628 ~ NA_real_,
    #                       TRUE ~ x  )
    #             }) %>%
    mutate_at(unique(election_year_data_full$party),
              function(x){x + rnorm(length(x),0,0.001)}) %>%
    na.omit() %>%
    filter(row_number() <= 365) %>% 
    ggplot(.,aes(x=fdp,y=spd)) + geom_point() + geom_smooth(method='lm')
  
  
  correlation_between_parties_year <- election_year_data_full %>%
    mutate(election_year_day = paste0(election_year, days_until_ed)) %>%
    arrange(days_until_ed) %>%
    dplyr::select(election_year_day,party,loess_trend) %>%
    spread(party,loess_trend) %>%
    # mutate_at(unique(election_year_data_full$party),
    #           function(x){ pmax(-0.01,pmin(0.01,(x-lag(x)))  ) }) %>%
    # na.omit() %>%
    # mutate_at(unique(election_year_data_full$party),
    #           function(x){ ifelse(abs((x - mean(x))/sd(x) ) > 1.628, NA, x  )   }) %>%
    mutate_at(unique(election_year_data_full$party),
              function(x){x + rnorm(length(x),0,0.001)}) %>%
    dplyr::select(unique(election_year_data_full$party)) %>%
    na.omit() %>%
    filter(row_number() <= 365) %>% 
    cor
  
  # correlation_between_parties_year <- correlation_between_parties_year * 0.75
  # correlation_between_parties_year[correlation_between_parties_year > 0.8] <- 0.8
  # correlation_between_parties_year[correlation_between_parties_year < -0.8] <- -0.8
  # diag(correlation_between_parties_year) <- 1
  
  correlation_between_parties_year <- lqmm::make.positive.definite(correlation_between_parties_year)
  correlation_between_parties_year
  
  # get data for that day
  election_year_data <- historical_prior_data_posterior_param %>% filter(election_year == y, days_until_ed == t)
  
  # get mean and sigma
  party_mus_sigmas <- election_year_data %>% 
    dplyr::select(party,lm_pred,lm_fitted_error,lm_fitted_df)
  
  # simulate 10k random walks, allowing forecast to change over time
  simulated_errors <-  mvtnorm::rmvt(n = (t+1)*num_sims,
                                     sigma = correlation_between_parties_year[match(rownames(correlation_between_parties_year),party_mus_sigmas$party),],
                                     df = 10) # party_mus_sigmas$lm_fitted_df)
  
  # convert to random walk
  simulated_errors <- simulated_errors %>% 
    as_tibble %>% 
    mutate(trial = floor((row_number() - 1) / (t+1)) +1 ) %>% 
    group_by(trial) %>%
    mutate_at(vars(-group_cols()), cumsum) %>%
    ungroup() %>%
    dplyr::select(-trial) %>%
    as.matrix
  
  # multiply error times sd, add mean
  simulated_errors <- simulated_errors * 
    outer(rep.int(1L, nrow(simulated_errors)), 
          party_mus_sigmas$lm_fitted_error / sqrt(t+1)) # one added step, bc ed=0 == 1 step
  
  simulated_errors <- simulated_errors + 
    outer(rep.int(1L, nrow(simulated_errors)),
          party_mus_sigmas$lm_pred)
  
  # glance
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
    as_tibble() %>%
    ungroup() %>%
    mutate(trial = floor((row_number() - 1) / (t+1)) +1 ) %>% 
    group_by(trial) %>%
    mutate(days_until_ed = (t+1) - dplyr::row_number()) %>%
    ungroup() %>%
    mutate(election_year = y)
  
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
    mutate(election_year = y)
  
  range(simulated_errors$days_until_ed)
  range(simulated_errors$trial)
  
  # check if simulated sd matches target
  apply(simulated_errors[simulated_errors$days_until_ed==0,][,1:(ncol(simulated_errors)-4)], 2, sd)
  party_mus_sigmas$lm_fitted_error
  
  apply(simulated_errors[simulated_errors$days_until_ed==0,][,1:(ncol(simulated_errors)-4)], 2, mean)
  party_mus_sigmas$lm_pred
  
  # glance again
  # simulated_errors %>% sample_n(pmin(1000,max(simulated_errors$trial)/5)) %>% ggplot(.,aes(cdu,afd)) + geom_point() + geom_smooth(method='lm')
  simulated_errors %>%
    dplyr::filter(trial %in% 1:100) %>%
    gather(party,simulated_vote,1:(length(.)-4)) %>%
    ggplot(.,aes(x=days_until_ed, y=simulated_vote, col=party, group=paste0(party,trial) ),fill=NA) +
    geom_line() +
    scale_x_reverse()
  
  if(return_all){return(simulated_errors)}else{return(simulated_errors[simulated_errors$days_until_ed==0,])}
  
}

simulate_election_y_from_time_t(y=1953,t=30)
simulate_election_y_from_time_t(y=1980,t=30)


# test confidence interval coverage for various time horizons for all years
test_forecasts <- map_df(validation_years[validation_years >= 1953],
                         function(xi){
                           map_df(c(0,14,30,60,90,120),
                                  function(xj){
                                    simulate_election_y_from_time_t(xi,xj,num_sims = 20000)
                                  })
                         })

beepr::beep(2)

test_forecasts <- test_forecasts %>% relocate(cdu,fdp,oth,spd,gru,lin,afd)

# replicate the figure above, for all years
pdf('test_forecasts.pdf',width=10,height=20)
test_forecasts %>% 
  dplyr::filter(days_until_ed == 0,election_year >= 2002) %>% 
  dplyr::select(-oth) %>%
  gather(party,prediction,1:(length(.)-4)) %>%
  group_by(election_year, party,forecast_date) %>%
  summarise(yhat_mean = mean(prediction),
            yhat_upper_95 = quantile(prediction,0.975,na.rm=T),
            yhat_lower_95 = quantile(prediction,0.025,na.rm=T),
            yhat_upper_80 = quantile(prediction,0.9,na.rm=T),
            yhat_lower_80 = quantile(prediction,0.1,na.rm=T)) %>%
  left_join(priors_to_attach %>% dplyr::select(election_year, party, prior_est, voteshare)) %>%
  left_join(prior %>% dplyr::select(election_year, party, lag_voteshare)) %>%
  arrange(desc(forecast_date)) %>%
  ungroup() %>%
  mutate(party = as.factor(party),
         party = fct_reorder(party, -voteshare, mean),
         forecast_date = as.factor(forecast_date),
         forecast_date = fct_inorder(forecast_date)) %>%
  ggplot(.,aes(x=forecast_date)) +
  geom_point(aes(y=yhat_mean,col='prediction')) +
  geom_line(aes(x=as.numeric(forecast_date),y=voteshare,col='vote',linetype='vote')) +
  geom_line(aes(x=as.numeric(forecast_date),y=prior_est,col='prior',linetype='prior')) +
  geom_line(aes(x=as.numeric(forecast_date),y=lag_voteshare,col='lag vote',linetype='lag vote')) +
  geom_segment(aes(y=yhat_lower_95,yend=yhat_upper_95,xend=(forecast_date)),alpha=0.4) +
  geom_segment(aes(y=yhat_lower_80,yend=yhat_upper_80,xend=(forecast_date)),alpha=0.2) +
  theme_minimal() +
  labs(x='forecast generated at x weeks before election',
       y='predicted vote share',
       subtitle = 'forecasts for select parties (columns) in select years (rows)')  +
  facet_grid(rows = vars(election_year), cols = vars(party)) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
  coord_cartesian(ylim=c(0,0.51)) +
  scale_color_manual(values=c('vote' ='gray40','prior'='gray40','lag vote'='gray60','prediction'='black')) +
  scale_linetype_manual(values=c('vote' =1,'prior'=2,'lag vote'=3)) 
dev.off()


# look at confidence interval coverage
test_forecasts %>% 
  filter(days_until_ed == 0) %>% 
  gather(party,prediction,1:(length(.)-4)) %>%
  na.omit() %>%
  group_by(election_year, party, forecast_date) %>%
  summarise(yhat_mean = mean(prediction),
            yhat_upper_99 = quantile(prediction,0.995),
            yhat_lower_99 = quantile(prediction,0.005),
            yhat_upper_95 = quantile(prediction,0.975),
            yhat_lower_95 = quantile(prediction,0.025),
            yhat_upper_80 = quantile(prediction,0.9),
            yhat_lower_80 = quantile(prediction,0.1)) %>%
  left_join(priors_to_attach %>% dplyr::select(election_year, party, voteshare)) %>%
  ungroup() %>%
  mutate(party = as.factor(party),
         party = fct_reorder(party, -voteshare, mean)) %>%
  summarise(inside_99_ci = mean(voteshare < yhat_upper_99 & voteshare > yhat_lower_99),
            inside_95_ci = mean(voteshare < yhat_upper_95 & voteshare > yhat_lower_95),
            inside_80_ci = mean(voteshare < yhat_upper_80 & voteshare > yhat_lower_80)) 


forecast_coverage <- test_forecasts %>% 
  filter(days_until_ed == 0) %>% 
  gather(party,prediction,1:(length(.)-4)) %>%
  na.omit() %>%
  group_by(election_year, party, forecast_date) %>%
  summarise(yhat_mean = mean(prediction),
            yhat_upper_99 = quantile(prediction,0.995),
            yhat_lower_99 = quantile(prediction,0.005),
            yhat_upper_95 = quantile(prediction,0.975),
            yhat_lower_95 = quantile(prediction,0.025),
            yhat_upper_80 = quantile(prediction,0.9),
            yhat_lower_80 = quantile(prediction,0.1)) %>%
  left_join(priors_to_attach %>% dplyr::select(election_year, party, voteshare)) %>%
  ungroup() %>%
  mutate(party = as.factor(party),
         party = fct_reorder(party, -voteshare, mean)) %>%
  group_by(forecast_date) %>%
  summarise(inside_99_ci = mean(voteshare < yhat_upper_99 & voteshare > yhat_lower_99),
            inside_95_ci = mean(voteshare < yhat_upper_95 & voteshare > yhat_lower_95),
            inside_80_ci = mean(voteshare < yhat_upper_80 & voteshare > yhat_lower_80)) 


ggplot(forecast_coverage, aes(x=forecast_date)) +
  geom_line(aes(y=inside_99_ci,linetype='99% CI coverage')) + 
  geom_line(aes(y=inside_95_ci,linetype='95% CI coverage')) + 
  geom_line(aes(y=inside_80_ci,linetype='80% CI coverage')) + 
  scale_x_reverse() +
  scale_y_continuous(limits=c(0.7,1),breaks=c(0.7, 0.8, 0.95, 0.99)) +
  labs(subtitle = 'how often are results inside the CI?\nabove y benchmark indicates underfitting -- eg simulating too much error')


# which party-years are the biggest source of error?
test_forecasts %>% 
  filter(days_until_ed == 0) %>% 
  gather(party,prediction,1:(length(.)-4)) %>%
  na.omit() %>%
  group_by(election_year, party, forecast_date) %>%
  summarise(yhat_mean = mean(prediction),
            yhat_upper_99 = quantile(prediction,0.995),
            yhat_lower_99 = quantile(prediction,0.005),
            yhat_upper_95 = quantile(prediction,0.975),
            yhat_lower_95 = quantile(prediction,0.025),
            yhat_upper_80 = quantile(prediction,0.9),
            yhat_lower_80 = quantile(prediction,0.1))  %>%
  left_join(priors_to_attach %>% dplyr::select(election_year, party, voteshare)) %>%
  mutate(outside_99_ci = (voteshare > yhat_upper_99) | (voteshare < yhat_lower_99)) %>%
  filter(outside_99_ci) %>% 
  dplyr::select(election_year, party, days_until_ed = forecast_date,
                yhat_mean, yhat_upper_99, yhat_lower_99, voteshare,
                outside_99_ci) %>% as.data.frame %>%
  print(.,digits=3)

test_forecasts %>% 
  filter(days_until_ed == 0) %>% 
  gather(party,prediction,1:(length(.)-4)) %>%
  na.omit() %>%
  group_by(election_year, party, forecast_date) %>%
  summarise(yhat_mean = mean(prediction),
            yhat_upper_99 = quantile(prediction,0.995),
            yhat_lower_99 = quantile(prediction,0.005),
            yhat_upper_95 = quantile(prediction,0.975),
            yhat_lower_95 = quantile(prediction,0.025),
            yhat_upper_80 = quantile(prediction,0.9),
            yhat_lower_80 = quantile(prediction,0.1))  %>%
  left_join(priors_to_attach %>% dplyr::select(election_year, party, voteshare)) %>%
  mutate(outside_95_ci = (voteshare > yhat_upper_95) | (voteshare < yhat_lower_95)) %>%
  filter(outside_95_ci) %>% 
  dplyr::select(election_year, party, days_until_ed = forecast_date,
                yhat_mean, yhat_upper_95, yhat_lower_95, voteshare,
                outside_95_ci) %>% as.data.frame %>%
  print(.,digits=3)


# tell me each coalition's chances on election day for each election
get_coalition_summary <- function(dat, parties){
   dat %>% 
    # transform to seats
    filter(days_until_ed == 0) %>%
    dplyr::select(-oth) %>%
    mutate_at(c('afd','cdu','fdp','gru','lin','spd'),
              function(x){case_when(is.na(x) ~ 0,
                                    x >= 0.05 ~ x, 
                                    TRUE ~ 0)}) %>%
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
    summarise(seats_won = sum(seats_won,na.rm=T)) %>%
    ungroup() %>%
    summarise(yhat_seat_mean = mean(seats_won,na.rm=T),
              yhat_seat_upper_95 = quantile(seats_won,0.975),
              yhat_seat_lower_95 = quantile(seats_won,0.025),
              yhat_seat_upper_80 = quantile(seats_won,0.9),
              yhat_seat_lower_80 = quantile(seats_won,0.1),
              pr_majority = mean(seats_won > 0.5)) %>% 
    return
}


coalition_probs <- map_df(validation_years[validation_years>=1980],
    function(x){
      print(x)
      get_coalition_summary(
        test_forecasts %>% 
          dplyr::filter(election_year == x, days_until_ed == 0,forecast_date==0) ,
        c('cdu','gru')
      ) %>%
        mutate(coalition = 'bg',election_year = x) %>%
        bind_rows(
          get_coalition_summary(
            test_forecasts %>% 
              dplyr::filter(election_year == x, days_until_ed == 0,forecast_date==0) ,
            c('cdu','gru','fdp')
          ) %>%
            mutate(coalition = 'jamaica',election_year = x)
        ) %>%
        bind_rows(
          get_coalition_summary(
            test_forecasts %>% 
              dplyr::filter(election_year == x, days_until_ed == 0,forecast_date==0) ,
            c('spd','gru','fdp')
          ) %>%
            mutate(coalition = 'traffic',election_year = x)
        ) %>%
        bind_rows(
          get_coalition_summary(
            test_forecasts %>% 
              dplyr::filter(election_year == x, days_until_ed == 0,forecast_date==0) ,
            c('spd','lin','gru')
          ) %>%
            mutate(coalition = 'rrg',election_year = x)
        ) %>%
        bind_rows(
          get_coalition_summary(
            test_forecasts %>% 
              dplyr::filter(election_year == x, days_until_ed == 0,forecast_date==0) ,
            c('cdu','spd')
          ) %>%
            mutate(coalition = 'grand',election_year = x)
        ) %>%
        bind_rows(
          get_coalition_summary(
            test_forecasts %>% 
              dplyr::filter(election_year == x, days_until_ed == 0,forecast_date==0) ,
            c('cdu','fdp')
          ) %>%
            mutate(coalition = 'by',election_year = x)
        )
    })

ggplot(coalition_probs,aes(x=election_year,col=coalition)) + 
  geom_line(aes(y=yhat_seat_mean))


# 6. Save outputs ---------------------------------------------------------
# Historical polls, priors, predictions and margins of error
# Functions to be used to turn 2021 polls and priors into predictions


# `prior` - the dataset of covariates used in historical priors, and subsequent predictions
write_csv(prior,'output-data/historical_plus2021_priors.csv')

# `historical_trends_and_priors` - the dataset used to train a lm to blend historical priors and trends
write_csv(historical_trends_and_priors,'output-data/historical_trends_and_priors.csv')

# `historical_prior_data_posterior_param` - the dataset of predictions and fitted lm, rmse values
write_csv(historical_prior_data_posterior_param, 'output-data/historical_prior_data_posterior_param.csv')

# the several functions used to figure out how much error to simulate
write_rds(loo_poll_error_function,'output-data/loo_poll_error_function.rds')
write_rds(err_mod,'output-data/lm_error_model.rds')
write_rds(err_mod_training_data,'output-data/lm_error_model_training_data.rds')

# `coalition_probs` historical coalition probabilities
write_csv(coalition_probs,'output-data/historical_coalition_shares_probs.csv')
