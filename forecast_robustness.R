

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
  
  bootstrap_results <- trigger_uncertainty(admin_0_triggers,5)
  
  plot <- ggplot(bootstrap_results,aes(y=new_threshold,x=month)) + 
    facet_grid(cols=vars(freq)) +
    geom_boxplot()  +
    geom_point(aes(y=old_threshold,x=month),color='red') + 
    geom_text(data = bootstrap_results %>% group_by(month,freq) %>% summarise_all(max), aes(y=old_threshold,x=month,label=round(old_threshold,2)),color='red',alpha=0.7,vjust=-0.5) +
    geom_hline(aes(yintercept = as.numeric(freq)),linetype='dashed') +
    theme_bw() +
    ylim(c(0,50)) +
    xlab("Issue month") +
    ylab("Prob. non-exceedence") +
    ggtitle(paste0(country," trigger uncertainty range"),subtitle = "AA threshold in red; clim. odds dotted line")
  
  png(file=paste0(getwd(),"/figures/",country,".png"),width=650,height=400)
  print(plot)
  dev.off()
  
  i <- i + 1
}