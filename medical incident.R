knitr::opts_chunk$set(echo = TRUE)

pacman::p_load(tidyverse,dplyr,ggplot2,lubridate,brms, data.table, cmdstanr, dlnm, splines, survival)

#if you want to use cmdstanr to run code without having to use C++, is supposed to help run MCMC chains faster 
# install.packages("remotes")
# remotes::install_github("stan-dev/cmdstanr")


#loading in data sets 
#medical incidents 
medical_incidents <- read.csv("Fire_Department_and_Emergency_Medical_Services_Dispatched_Calls_for_Service_20260309.csv")

daymet_current <- read.csv("daymet_current.csv")

daymet_p95_byZIP <- read.csv("daymet_p95_byZIP.csv")


per_day <- medical_incidents |>
  mutate(Call.Date = mdy(Call.Date)) |>
  count(Zipcode.of.Incident, Call.Date)

per_day <- per_day |>
  mutate(year = year(Call.Date),
         month = month(Call.Date)) |>
  relocate(year, .after = "Call.Date") |> 
  filter(month >= 6 & month < 11) 

#seasonal bias to ONLY track assumed heatwave related days, need to find paper to support this 


# first make a date column like in previous dataframe and add zero infront of site
daymet_current$site <- paste0("0", daymet_current$site)
daymet_p95_byZIP$site <- as.character(daymet_p95_byZIP$site)


# calculate AT from temp and vp 
daymet_current <- daymet_current %>%
  mutate(date = make_date(year = year) + days(yday - 1),
         month = month(date),
         apparent_tmin = -1.3 + 0.92*tmin + 2.2*(vp/1000),
         apparent_tmax = -1.3 + 0.92*tmax + 2.2*(vp/1000),
         apparent_tmean = (apparent_tmin+apparent_tmax)/2
         )


# merge data sets
daymet_current_join <- left_join(daymet_current, daymet_p95, by = "zip_code")

# remove duplicate observations
# Deduplicate by zip
daymet_heatwaves <- daymet_current_join %>%
  distinct(zip_code, date, .keep_all = TRUE) 

#tally heatwave days (above 95th percentile)
daymet_heatwaves <- daymet_heatwaves %>%
  filter(month >= 6 & month < 11) %>% #restrict to SF summer dates to avoid seasonal bias
  mutate(hw = apparent_tmean > p95_tmean) #trying mean now
daymet_heatwaves$zip_code <- as.character(daymet_heatwaves$zip_code)

daymet_heatwaves <- daymet_heatwaves %>%
  group_by(site.x, zip_code) %>%
  mutate(event_id = rleid(hw))

#summarise TRUE runs and 2+ day heatwaves
hw_events_2day <- daymet_heatwaves %>%
  filter(hw) %>% #keeps only heatwave days
  group_by(site.x, event_id, zip_code, year) %>%
  summarise(
    duration = n(),
    start_date = min(date),
    end_date = max(date),
    mean_tmean = mean(apparent_tmean),
    .groups = "drop"
  ) %>%
  filter(duration >= 2) # enforce 2+ days

#summarise TRUE runs and 3+ day heatwaves
hw_events_3day <- daymet_heatwaves %>%
  filter(hw) %>% #keeps only heatwave days
  group_by(site.x, event_id, zip_code, year) %>%
  summarise(
    duration = n(), #adjust to just "duration" if 
    start_date = min(date),
    end_date = max(date),
    mean_tmean = mean(apparent_tmean),
    .groups = "drop"
  ) %>%
  filter(duration >= 3) # enforce 3+ days



per_day <- per_day %>% 
  mutate(Zipcode.of.Incident = as.character(Zipcode.of.Incident))

daymet_heatwaves <- daymet_heatwaves %>% 
  mutate(zip_code = as.character(zip_code))
  
overall_current <- left_join(
  per_day,
  daymet_heatwaves,
  by = c("Call.Date" = "date", "Zipcode.of.Incident" = "zip_code", "year"= "year", "month" = "month"),
  na_matches = "never"
)

overall_current <- overall_current %>% 
  mutate(day_abbr = format(as.Date(Call.Date), "%a"),
         date = Call.Date,
         )

library(cmdstanr)
#install_cmdstan()
set_cmdstan_path(path = NULL)

fit <- brm(
  n ~ hw + (1 | Zipcode.of.Incident), #the outcome ~ the exposure, grouping by zipcode 
  data = overall_current, 
  family = negbinomial(), #negbinomial is for true/false counts   
   prior = c( #priors were set after running test runs of the model 
    prior(normal(0,1), class = "b"), #percentage affect of hw 
    prior(exponential(1), class = "sd"), #standard deviation
    prior(exponential(1), class = "shape"), #spatial distribution 
    prior(normal(2,1), class = "Intercept")), #what is your baseline calls on a false day     
  chains = 4, #how many times running MC 
  cores = 4, #how fast ur machine can run the model CPU 
  threads = threading(2), #forces iterations to be run at the same time much faster, MIGHT mess with coefficients  
  iter = 3000, #how many samples are you taking 
  warmup = 1000, #"throwaway" sample size 
  backend = "cmdstanr", 
  control = list(adapt_delta = 0.95),
  seed = 2149, #random but not really
  sample_prior = TRUE , #sampling prior while also fitting the model 
  ) 

#saving the model after it's been run for future analysis
saveRDS(fit, "050626_bayesian.rds") 

#bayesian <- readRDS("post_bayesian.rds") 


# bayesian <- readRDS("911_bayesian.rds") 
summary(fit)
plot(fit) #i am not entirely sure how to read these graphs 



knitr::purl(input = "medical incidents.Rmd", output = "medical incident.R",documentation = 0)

