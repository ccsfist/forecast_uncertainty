
library(dplyr)
library(RJSONIO)
library(curl)
library(abind)
library(sp)
library(tidyr)
library(httr)
library(ggplot2)


month_names <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")

### SETUP - CHANGE STUFF HERE ###

## set what data you'd like to download

# country / project

country <- c("niger")

# freqs (0-100)

freqs <- seq(5,95,by=5)

# months (numbered 0-11)

months <- c(0:5) 

months_df <- data.frame("month" = c(0:11),"name" = month_names)

prev_months_df <- data.frame("month" = c(1:11,0),"name" = month_names)

# obs dataset (see config for list of possible options)

obs_dataset <- c("chirps-precip-jas")

# server (e.g. production or test)

server <- "https://iridl.ldeo.columbia.edu/fbfmaproom2/"

## get key guide

# admin 0

admin_0_key <- bind_rows(fromJSON(paste0(server,"regions?country=",country,"&level=0")))

## function to set query paths 

fun_download_path <- function(country,level,region,month,freq,predictor,predictand){
  
  query_path <- paste0(server,country,"/export?",
                       "&mode=",level,
                       "&region=",region,
                       "&season=season1&issue_month0=",month,
                       "&freq=",freq,
                       "&predictor=",predictor,
                       "&predictand=",predictand) 
  return(query_path)
  
}

### DOWNLOADS ###

## admin 0

admin_0_downloads <- list()

for (freq in freqs) {
  
  admin_0_downloads[[as.character(freq)]] <- list()

  for (month in months) {
    
    api_download <- fromJSON(fun_download_path(country,0,admin_0_key[1,1],month,freq,"pnep",obs_dataset[1]),nullValue = NA)
    admin_0_downloads[[as.character(freq)]][[as.character(month)]] <- api_download
    

  }
  
}

admin_0_bad_years <- fromJSON(fun_download_path(country,0,admin_0_key[1,1],0,freq,"pnep",'bad-years-v2'),nullValue = NA)


bootstrap <- function(obj,n_omit) {
  
  years <- unique(obj$year)
  
  omit_index <- round(runif(n_omit,min=1,max=length(years)))
  
  years_select <- years[-omit_index]
  
  obj_select <- obj %>% filter(year %in% years_select) %>%
    group_by(month,severity) %>%
    summarise(new_threshold = quantile(forecast,1-(max(freq,na.rm=T)/100),na.rm=T),old_threshold = max(threshold,na.rm=T),
              freq = max(freq,na.rm=T)) %>%
    ungroup()
  
  return(obj_select)
  
}

bootstrap_results <- lapply(c(1:100),function(x) bootstrap(admin_0_triggers,5))
bootstrap_results <- bind_rows(bootstrap_results)

ggplot(bootstrap_results,aes(y=new_threshold,x=month)) + 
  facet_grid(cols=vars(severity)) +
  geom_boxplot()  +
  geom_point(data=admin_0_triggers %>% filter(year == max(year)),aes(y=forecast,x=month),color='blue',alpha=0.7) +
  geom_text(data=admin_0_triggers %>% filter(year == max(year)),aes(y=forecast,x=month,label=round(forecast,2)),color='blue',alpha=0.7,vjust=-0.5) +
  geom_point(aes(y=old_threshold,x=month),color='red') + 
  geom_text(data = bootstrap_results %>% group_by(month,freq) %>% summarise_all(max), aes(y=old_threshold,x=month,label=round(old_threshold,2)),color='red',alpha=0.7,vjust=-0.5) +
  geom_hline(aes(yintercept = freq),linetype='dashed') +
  theme_bw() +
  xlab("Issue month") +
  ylab("Prob. non-exceedence") +
  ggtitle("Trigger uncertainty range",subtitle = "AA threshold in red; current forecast in blue; clim. odds dotted line")

