

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

countries <- c("madagascar","ethiopia","niger","guatemala-jja","djibouti-ond")

# freqs (0-100)

freqs <- c(20,30)

# season (1 or 2)

season <- 'season1'

# months (numbered 0-11)

months_list <- list(c(8, 9, 10),c(11, 0, 1),c(3, 4, 5),c(2, 3, 4),c(6, 7, 8)) 
months_df <- data.frame("month" = c(0:11),"name" = month_names)
prev_months_df <- data.frame("month" = c(2:11,0,1),"name" = month_names)

# define dataset names

forecast_list <- c("pnep-v3","pnep-v5","prcp-v5","pnep-v6","pnep-v2")
bad_years_list <- c("bad-years-v3","bad-years","bad-years-v3","bad-years","bad-years-rank-v2") 

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
  
  select_data <- admin_0_triggers %>% filter(month == max(month), freq == max(freq)) ## at some point we need to update this for all freqs, months
  bootstrap <- trigger_uncertainty_new(select_data,5) ## at some point we need to make this the default uncertainty function
  
  colors <- c("Cumulative Forecasts, 1983-present"="black",
              "Empirical Trigger Threshold"="red",
              "Bootstrapped Probability of Passing Threshold" = "blue")
  
  plot <- ggplot(select_data,aes(x=forecast)) + stat_ecdf(aes(color="Cumulative Forecasts, 1983-present")) + 
    stat_ecdf(data=bootstrap,aes(x=new_threshold,color='Bootstrapped Probability of Passing Threshold')) +
    geom_vline(data=bootstrap,aes(xintercept = old_threshold,color='Empirical Trigger Threshold'),linetype='dashed') +
    scale_color_manual(values=colors) +
    labs(x = "Forecast Strength", y = "Probability", color = "Legend") +
    scale_y_continuous(labels = scales::percent) +
    ggtitle(paste0("Forecast Trigger Uncertainty Range For ",country))
  
  png(file=paste0(getwd(),"/figures/",country,".png"),width=650,height=350)
  print(plot)
  dev.off()
  
  i <- i + 1
}




