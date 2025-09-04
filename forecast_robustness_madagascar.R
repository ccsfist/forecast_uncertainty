

## This file generates estimates of forecast uncertainty for the BAMS paper 
# "Quantifying Uncertainty in Forecast-Based Anticipatory Action Triggers Using Bootstrapping Methods"
## code written by Max Mauerman (max.mauerman@columbia.edu)


source("~/iri-fist/aa_dashboard_common_code/aa_functions.R")
setwd("~/forecast_uncertainty")

options(timeout=1000000000000)

month_names <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")

### SETUP  ###

## set what data you'd like to download

# country / project

countries <- c("madagascar")

# freqs (0-100)

freqs <- c(15)

# season (1 or 2)

season <- 'season1'

# months (numbered 0-11)

months_list <- list(c(8, 9, 10)) 
months_df <- data.frame("month" = c(0:11),"name" = month_names)
prev_months_df <- data.frame("month" = c(2:11,0,1),"name" = month_names)

# define dataset names

forecast_list <- c("pnep-v3")
bad_years_list <- c("enacts-precip-all") 

# server (e.g. production or test)

server <- "https://iridl.ldeo.columbia.edu/fbfmaproom2/"

i <- 1
for (country in countries) {

  ## get key file
  
  admin_0_key <- get_key(server,country,0)
  
  ### DOWNLOADS ###
  
  admin_0_rain <- download_data_national(server=server,country=country,key=admin_0_key[1,1],season=season,months=months_list[[i]],
                                         freqs=freqs,predictor=forecast_list[i],predictand=bad_years_list[i])
  
  admin_0_triggers <- compile_data_national(forecast=admin_0_rain,badyears=admin_0_rain,ndvi=NA,prevseas=NA)
  
  # select_data <- admin_0_triggers %>% filter(month == max(month), freq == max(freq)) 
  select_data <- admin_0_triggers
  bootstrap <- trigger_uncertainty(select_data,5) 
  
  colors <- c("Cumulative Forecasts, 1983-present"="black",
              "Empirical Trigger Threshold"="red",
              "Bootstrapped Probability of Passing Threshold" = "blue",
              "Driest Years"="darkgreen")
  
  years_tohighlight <- select_data %>% 
    group_by(month) %>%
    mutate(pctile=ecdf(forecast)(forecast)) %>%
    filter(badyear==1)
  
  freq_shown <- max(admin_0_triggers$freq)
                                           
  plot <- ggplot(select_data,aes(x=forecast)) + stat_ecdf(aes(color="Cumulative Forecasts, 1983-present")) + 
    facet_wrap(~month) +
    stat_ecdf(data=bootstrap,aes(x=new_threshold,color='Bootstrapped Probability of Passing Threshold')) +
    geom_vline(data=bootstrap,aes(xintercept = old_threshold,color='Empirical Trigger Threshold'),linetype='dashed') +
    geom_rect(xmin=quantile(bootstrap$new_threshold,0.10),xmax=quantile(bootstrap$new_threshold,0.90),ymin=Inf,ymax=-Inf,fill='blue',alpha=0.01) +
    geom_point(data=years_tohighlight,aes(x=forecast,y=pctile,color='Driest Years')) +
    geom_text(data=years_tohighlight,aes(x=forecast,y=pctile,label=year,color='Driest Years'),size=4,hjust=1,vjust=-1) +
    scale_color_manual(values=colors) +
    labs(x = "Forecast Strength", y = "Cumulative Probability", color = "Legend") +
    scale_y_continuous(labels = scales::percent) +
    ggtitle(paste0("Forecast Trigger Uncertainty Range For ",country),subtitle = paste0("for ",freq_shown,"th percentile drought event")) +
    xlim(c(0,100)) +
    theme(legend.position = 'bottom')
  
  print(plot)

  i <- i + 1
}




