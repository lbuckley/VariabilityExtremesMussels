#load libraries
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)

library(ggplot2)

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/ICBClimateBiology/analysis/")
source("TempcyclesAnalysis.R")

#ROBOMUSSEL ANALYSIS

#SITES
# WaOr Tatoosh Island, WA 1660.1 48.39 124.74
# WaOr Boiler Bay, OR 1260.7 44.83 124.05
# WaOr Strawberry Hill, OR 1196 44.25 124.12
# CenCal Hopkins, CA 327.1 36.62 121.90
# CenCal Piedras Blancas, CA 208.11 35.66 121.28
# CenCal Cambria, CA 185.66 35.54 121.10
# SoCal Lompoc, CA 84.175 34.72 120.61
# SoCal Jalama, CA 57.722 34.50 120.50
# SoCal Alegria, CA 37.284 34.47 120.28
# SoCal Coal Oil Point (COP), CA 0 34.41 119.88

#-----------------
#Robomussel data
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/ICBClimateBiology/data/")
site.dat= read.csv("musselREADME.csv")

#Load robomussel data
te.max <- readRDS("tedat.rds")
#te.max= read.csv("tedat.csv")

#ungroup due to formatting error 
te.max= ungroup(te.max)

#Fix duplicate CP in WA
te.max$lat= as.character(te.max$lat)
te.max$site= as.character(te.max$site)
te.max[which(te.max$lat==48.45135),"site"]<-"CPWA"
te.max$site= as.factor(te.max$site)

#drop 1 of two close sites
te.max<- te.max[-which(te.max$site=="LB"),]

#subset sites to Boiler Bay
te.max= subset(te.max, te.max$site %in% c("BB") ) 
#te.max2= subset(te.max, te.max$site %in% c("SD","BB","PD") ) 
#te.max2= subset(te.max, te.max$lat %in% c(48.39137,44.83064,35.66582,34.46717) )

#----------------------
#add weather station data 
# Boiler Bay, 44.8336° N, 124.0624° W

#NEWPORT, OR US (USW00024285.csv)
#Or use Coos, OR for full data

library(rnoaa)
library(TrenchR)

#map: https://ghcn-leaflet.herokuapp.com/
#NEWPORT
#tmax=try( ghcnd_search("USW00024285", var = "TMAX") )
#COOS BAY
#tmax=try( ghcnd_search("USW00004141", var = "TMAX") )
#OTIS 2NE USC00356366
tmax=try( ghcnd_search("USC00356366", var = "TMAX") )

dat= tmax$tmax

#split date
date= as.Date(dat$date, "%Y-%m-$d")
dat$year=as.numeric(format(date, "%Y"))
dat$month=as.numeric(format(date, "%m"))
dat$day=as.numeric(format(date, "%d"))
dat$j= unlist(lapply(date, FUN="day_of_year"))

dat1=subset(dat, dat$year %in% 2000:2015)
dat1$doy= dat1$j

#divide temp by 10 to account for format
dat1$tmax = dat1$tmax/10

#----------------------
#TIME SERIES PLOTS

#time series
te.max1= subset(te.max2, te.max2$year==2002)

#restrict to summer
#May 1 through September: 121:273 
te.max1= subset(te.max1, te.max1$doy>120 & te.max1$doy<274)

#FIG 1A, time series
fig2a.bb<- ggplot(data=te.max1[which(te.max1$site=="BB" & te.max1$subsite %in% c(2,13,10) ),], aes(x=doy, y = MaxTemp_C))+
  geom_line(aes(color=subsite), alpha=0.8) +theme_classic()+
  guides(color=FALSE)+labs(x = "Day of year",y="Maximum daily temperature (°C)")
#missing site

#add weather station
fig.timeseries= fig2a.bb+ geom_line(data=dat1[which(dat1$year==2002 & dat1$j %in% 120:274),], aes(x=j, y =tmax), alpha=0.8) +theme_classic()

#------------------
#FREQUENCY ANALYSIS
# https://github.com/georgebiogeekwang/tempcycles/

#power spectrum
#x: frequency (1/days)
#y: log amplitude

fseq= exp(seq(log(0.001), log(1), length.out = 200))

sites= c("BB") #levels(te.max2$site)
subsites=  levels(te.max2$subsite)

pow.out= array(NA, dim=c(length(sites),length(subsites),length(fseq) ) )

for(site.k in 1:length(sites))
{
  te.dat= te.max2[which(te.max2$site==sites[site.k]),]
  subsites1= levels(te.dat$subsite)
  
  for(subsite.k in  1:length(subsites)) {
    te.dat1= te.dat[which(te.dat$subsite==subsites1[subsite.k]),]
    
    pow.out[site.k, subsite.k,] <- spec_lomb_phase(te.dat1$MaxTemp_C, te.dat1$j, freq=fseq)$cyc_range
  }
}

dimnames(pow.out)[[1]]<- sites
dimnames(pow.out)[[2]]<- 1:75

#----------------
#Te frequency
#fill data
library(zoo)
dat1$tmax= na.fill(dat1$tmax, fill="extend")

pow.out.te <- spec_lomb_phase(dat1$tmax, dat1$j, freq=fseq)$cyc_range

#-----------
#to long format
for(site.k in 1:length(sites)){
  pow1= pow.out[site.k,,]
  pow1= na.omit(pow1)
  pow1m= melt(pow1)
  pow1m$site= sites[site.k]
  
  if(site.k==1)pow=pow1m
  if(site.k>1)pow=rbind(pow,pow1m)
}
colnames(pow)[1:3]=c("subsite","freq","cyc_range")

#correct freq values
pow$freq= fseq[pow$freq]

#subset to intertidal heights
#pow1= subset(pow, pow$subsite %in% c("2","13","10") ) #"10","3","4"
#13: High, #10 and 2 low
#pow1= subset(pow, pow$site=="SD" & pow$subsite %in% c("9","12","10") ) 

#pow Te
pow.te= melt(pow.out.te)
pow.te$site="Te"
pow.te$subsite="Te"
pow.te$freq= fseq
colnames(pow.te)[1]="cyc_range"
pow.te= pow.te[,c("subsite","freq","cyc_range","site")]

#bind
pow=rbind(pow,pow.te)

#sort by frequency
pow= pow[order(pow$site, pow$subsite, pow$freq),]
pow$subsite= factor(pow$subsite)

#Add heights
subs= c(2, 13, 10, "Te")
height= c("low","high","low")
pow1$Height= height[match(pow1$subsite, subs)]

table(pow1$subsite)

#subsites
pow1= pow[pow$subsite %in% subs,]

fig2b<- ggplot(data=pow1, aes(x=log(freq), y = log(cyc_range/2) ))+geom_line(alpha=0.8, aes(color=subsite)) +theme_classic()+ 
  #guides(color=FALSE, size=FALSE)+
  geom_vline(xintercept=-2.639, color="gray")+geom_vline(xintercept=-1.946, color="gray")+geom_vline(xintercept=-3.40, color="gray")+geom_vline(xintercept=-5.9, color="gray")+
  labs(x = "log (frequency) (1/days)",y="log (amplitude)")+
  annotate(geom="text", x=-2.1, y=-5.5, label="1 week", size=3, color="black",angle=90)+ 
  annotate(geom="text", x=-2.8, y=-5.5, label="2 weeks", size=3, color="black",angle=90)+ 
  annotate(geom="text", x=-3.6, y=-5.5, label="1 month", size=3, color="black",angle=90)+ 
  annotate(geom="text", x=-6.2, y=-5.5, label="1 year", size=3, color="black",angle=90)+
  ylim(range(-5.8,2.05))

#--------------------
#average over different time periods and show performance implications
#Jensen's inequality

#pick subsites

#=================================
#combine
library(patchwork)

#Figures
#TPC through seasons? TPCs with variability?

#boiler bay
fig2a.bb+ fig2b + plot_layout(ncol = 1)

#count by site, subsite, year
te.count = te.max %>% group_by(year,site, subsite, zone, tidal.height..m.) %>% summarise( count=length(MaxTemp_C)  )
te.count= as.data.frame(te.count)
#find sub sites to use




