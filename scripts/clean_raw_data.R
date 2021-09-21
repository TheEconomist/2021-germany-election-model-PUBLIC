library(tidyverse)
library(readstata13)

# clean historical election polls -----------------------------------------
dat <- readstata13::read.dta13('raw-data/polls_btw_wide.dta')

dat <- dat %>% filter(!is.na(election_date))


# elecdate, date, poll, start, end, n, party, share,

dat %>% 
  arrange(election_date, date) %>%
  gather(party,share,c('cdu','spd','fdp','gru','lin','afd','oth')) %>%
  mutate(start = date, end = date) %>%
  select(election_date, date,poll=institute,start,end,n=sample_size,party,share) %>%
  write_csv('data/historical_polls.csv')



# clean structural prior data ---------------------------------------------



readstata13::read.dta13('raw-data/ger_model_df.dta') %>% 
  select(year, party, voteshare, chancellor,chancellor_party,
         gov, term, unemp, cdsu_gov, spd_gov, fdp_gov, gru_gov, lin_gov) %>% 
  write_csv('data/prior.csv')