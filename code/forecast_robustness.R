

## This file generates estimates of forecast uncertainty for the BAMS paper 
# "Quantifying Uncertainty in Forecast-Based Anticipatory Action Triggers Using Bootstrapping Methods"
## code written by Max Mauerman (max.mauerman@columbia.edu)


source("~/iri-fist/aa_dashboard_common_code/aa_functions.R")
setwd("~/forecast_uncertainty")
library(ggpubr)
library(rnaturalearth)
library(rnaturalearthdata)
library(cowplot)

options(timeout=1000000000000)

month_names <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")

### SETUP  ###

## set what data you'd like to download

# country / project
# replace guatemala with lesotho, add another west african country 

countries <- c("madagascar","ethiopia","niger","guatemala-jja","djibouti-ond") # make display names nicer

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

plots <- list()
ranges <- list()
i <- 1
for (country in countries) {

  print(paste0("Running for ",country))
  
  ## get key file
  
  admin_0_key <- get_key(server,country,0)
  
  ### DOWNLOADS ###
  
  admin_0_rain <- download_data_national(server=server,country=country,key=admin_0_key[1,1],season=season,months=months_list[[i]],
                                         freqs=freqs,predictor=forecast_list[i],predictand=bad_years_list[i])
  
  admin_0_triggers <- compile_data_national(forecast=admin_0_rain,badyears=admin_0_rain,ndvi=NA,prevseas=NA)
  
  select_data <- admin_0_triggers %>% filter(month == max(month), freq == max(freq)) 
  
  write.csv(select_data,paste0(getwd(),"/1. data/",country,".csv"),row.names = FALSE)
  
  bootstrap <- trigger_uncertainty(select_data,5) 
  
  colors <- c("Forecasts, 1983-present"="black",
              "Empirical Trigger Threshold"="red",
              "Bootstrapped Probability of Passing Threshold" = "blue",
              "Current Forecast"="green")
  
  current_forecast <- select_data$forecast[select_data$year == max(select_data$year)]
  forecast_cdf <- ecdf(select_data$forecast)
  current_pctile = forecast_cdf(current_forecast)
  freq_shown <- max(admin_0_triggers$freq)
                                           
  plots[[i]] <- ggplot(select_data,aes(x=forecast)) + stat_ecdf(aes(color="Forecasts, 1983-present")) + 
    stat_ecdf(data=bootstrap,aes(x=new_threshold,color='Bootstrapped Probability of Passing Threshold')) +
    geom_vline(data=bootstrap,aes(xintercept = old_threshold,color='Empirical Trigger Threshold'),linetype='dashed') +
    geom_rect(xmin=quantile(bootstrap$new_threshold,0.10),xmax=quantile(bootstrap$new_threshold,0.90),ymin=Inf,ymax=-Inf,fill='blue',alpha=0.01) +
    geom_point(x=current_forecast,y=current_pctile,size=3,aes(color='Current Forecast')) +
    scale_color_manual(values=colors,drop=FALSE) +
    labs(x = "Forecast Strength", y = "Cumulative Probability", color = "Legend") +
    scale_y_continuous(labels = scales::percent) +
    ggtitle(paste0("Forecast Trigger Uncertainty Range For ",country),subtitle = paste0("for ",freq_shown,"th percentile drought event")) +
    theme_bw() +
    scale_x_continuous(breaks=seq(0,100,by=5), limits = c(0,100)) +
    theme(legend.position = 'bottom') +
    guides(color = guide_legend(nrow=2))
  
  png(file=paste0(getwd(),"/2. figures/",country,".png"),width=650,height=350)
  print(plots[[i]])
  dev.off()
  
  write.csv(bootstrap,paste0(getwd(),"/3. results/",country,".csv"),row.names = FALSE)
  
  ranges[[i]] <- bootstrap %>% summarise(p10=quantile(new_threshold,0.1),
                                         p90=quantile(new_threshold,0.9),country=country)
  
  i <- i + 1
}

## fig 1 for paper - example of threshold calculation with obs data

plain <- theme(
  axis.text = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.title = element_blank(),
  panel.background = element_rect(fill = "white"),
  plot.title = element_text(hjust = 0.5)
)

admin_0_key <- get_key(server,countries[1],0)

admin_0_rain <- download_data_national(server=server,country=countries[1],key=admin_0_key[1,1],season=season,months=months_list[[1]],
                                       freqs=freqs,predictor=forecast_list[1],predictand='enacts-precip-mon')

admin_0_triggers <- compile_data_national(forecast=admin_0_rain,badyears=admin_0_rain,ndvi=NA,prevseas=NA)

fig1_data <- admin_0_triggers %>% filter(month == max(month), freq == max(freq)) %>% mutate(forecast = rain)

bootstrap <- trigger_uncertainty(fig1_data,5) 

current_forecast <- fig1_data$forecast[fig1_data$year == 2023]
forecast_cdf <- ecdf(fig1_data$forecast)
current_pctile = forecast_cdf(current_forecast)
freq_shown <- max(admin_0_triggers$freq)

# example of intermediate bootstrap calculations

bootstrap_fun <- function(obj,n_block) {
  
  years <- unique(obj$year)
  years_touse <- years[c(1:(length(years)-(n_block-1)))]
  
  years_block <- lapply(sample(c(1:length(years_touse)),round(length(years)/n_block)),function(x){
    data.frame(years[seq(x,x+(n_block-1) )])
  })
  
  return(years_block)
  
}

bootstrap_example <- bootstrap_fun(fig1_data,5)

bootstrap_example <- bind_rows( bootstrap_example, .id = "id" ) %>%
  rename_at(c(2),~c("year"))

fig1a <- ggplot(bootstrap_example,aes(y=id,x=factor(year))) + geom_text(aes(label=year),size=3,color='blue') + plain +
  ggtitle("Example of single bootstrap sample")

# histogram

fig1_data_select <- bootstrap_example %>% left_join(fig1_data,by="year")

new_threshold <- fig1_data_select %>%
  mutate(freq_inverted = 1-(as.numeric(freq)/100)) %>%
  group_by(freq,month) %>%
  mutate(threshold_new = quantile(forecast,freq_inverted,na.rm=T)) %>%
  summarise(new_threshold = max(threshold_new)) %>%
  ungroup()

fig1b <- ggplot(fig1_data,aes(x=forecast)) + geom_histogram(fill='lightgrey',binwidth = 25,aes(color="Forecasts, 1983-present", y = after_stat(count / sum(count)))) + theme_bw() +
  geom_histogram(data=fig1_data_select, fill=NA,binwidth = 25,color='blue',alpha=0.3,aes(y = after_stat(count / sum(count)))) +
  scale_x_continuous(breaks=seq(0,300,by=50), limits = c(0,300)) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  geom_vline(data=bootstrap,aes(xintercept = old_threshold,color='Empirical Trigger Threshold'),linetype='dashed') +
  geom_vline(data=new_threshold,aes(xintercept = new_threshold),color='blue',linetype='dashed') +
  geom_point(x=current_forecast,y=current_pctile,size=3,aes(color='Current Forecast')) +
  labs(x = "Cumulative Rainfall (mm)", y = "Proportion of Years", color = "Legend") +
  ggtitle(paste0("Historical Distribution of Rainfall For ",countries[1],", 1992-2024"),
          subtitle = "Result of 1 bootstrap in blue") +
  theme(legend.position = 'none') +
  scale_color_manual(values=colors,drop=FALSE) +
  guides(color = guide_legend(nrow=2))
  
# uncertainty plot

fig1c <- ggplot(fig1_data,aes(x=forecast)) + stat_ecdf(aes(color="Forecasts, 1983-present")) + 
  stat_ecdf(data=bootstrap,aes(x=new_threshold,color='Bootstrapped Probability of Passing Threshold')) +
  geom_vline(data=bootstrap,aes(xintercept = old_threshold,color='Empirical Trigger Threshold'),linetype='dashed') +
  geom_rect(xmin=quantile(bootstrap$new_threshold,0.10),xmax=quantile(bootstrap$new_threshold,0.90),ymin=Inf,ymax=-Inf,fill='blue',alpha=0.01) +
  geom_point(x=current_forecast,y=current_pctile,size=3,aes(color='Current Forecast')) +
  scale_color_manual(values=colors) +
  labs(x = "Cumulative Rainfall (mm)", y = "Cumulative Probability", color = "Legend") +
  scale_y_continuous(labels = scales::percent) +
  ggtitle(paste0("Trigger Uncertainty Range For ",countries[1]),subtitle = "After 1,000 bootstrap replicates") +
  theme_bw() +
  scale_x_continuous(breaks=seq(0,300,by=50), limits = c(0,300)) +
  theme(legend.position = 'bottom') +
  scale_color_manual(values=colors,drop=FALSE) +
  guides(color = guide_legend(nrow=2))

# put it all together

png(file=paste0(getwd(),"/2. figures/fig1.png"),width=600,height=600)
ggarrange(fig1a,fig1b,fig1c,nrow=3,heights = c(1,1,2))
dev.off()


## fig 2 for paper - map of WFP countries colored by degree of forecast uncertainty, w/ inset plots

uncertainty_summary <- bind_rows(ranges) %>% mutate(iqr = p90-p10)

world <- ne_countries(scale="medium",returnclass = "sf")

world_select <- world %>%
  mutate(selected = ifelse(admin %in% c("Madagascar","Ethiopia","Niger","Guatemala","Djibouti"),1,0)) %>%
  filter(selected == 1) %>%
  mutate(iqr = NA)

world_select$iqr[world_select$admin=="Madagascar"] <- uncertainty_summary$iqr[1]
world_select$iqr[world_select$admin=="Ethiopia"] <- uncertainty_summary$iqr[2]
world_select$iqr[world_select$admin=="Niger"] <- uncertainty_summary$iqr[3]
world_select$iqr[world_select$admin=="Guatemala"] <- uncertainty_summary$iqr[4]
world_select$iqr[world_select$admin=="Djibouti"] <- uncertainty_summary$iqr[5]

ethiopia_bbox <- st_bbox(world %>% filter(admin=="Ethiopia"))
niger_bbox <- st_bbox(world %>% filter(admin=="Niger"))
africa_bbox <- st_bbox(world %>% filter(continent=="Africa"))

fig2a <- ggplot(data = world) + 
  geom_sf(fill=NA,color='lightgrey') +
  geom_sf(data=world_select,aes(fill=iqr),color='black') +
  geom_sf_text(data=world_select,aes(label=admin),nudge_x = -5, nudge_y = 5,size=4) +
  plain +
  theme(legend.position = 'top') +
  labs(fill="IQR of Trigger Uncertainty") +
  coord_sf(xlim = c(africa_bbox[1],africa_bbox[3]), ylim=c(africa_bbox[2],africa_bbox[4])) +
  geom_rect(xmin = ethiopia_bbox[1],
            ymin = ethiopia_bbox[2],
            xmax = ethiopia_bbox[3],
            ymax = ethiopia_bbox[4],
            fill=NA,color='black',size=0.4) +
  geom_rect(xmin = niger_bbox[1],
            ymin = niger_bbox[2],
            xmax = niger_bbox[3],
            ymax = niger_bbox[4],
            fill=NA,color='black',size=0.4)

fig2b <- ggarrange(plots[[2]],plots[[3]],nrow=2,common.legend = TRUE) # show sharpness ratio?

png(file=paste0(getwd(),"/2. figures/fig2.png"),width=1000,height=500)
ggarrange(fig2a,fig2b,ncol=2)
dev.off()
