#load libraries
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(zoo)

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
te.max2= subset(te.max, te.max$lat %in% c(48.39137,44.83064,35.66582,34.46717) )

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

vir.cols=c("#5ec962","#21918c","#3b528b","#440154")

#time series
te.max1= subset(te.max2, te.max2$year==2007)

#restrict to summer
#May 1 through September: 121:273 
te.max1= subset(te.max1, te.max1$doy>120 & te.max1$doy<274)

#subset
te.sub= te.max1[which(te.max1$site=="BB" & te.max1$subsite %in% c(5,2,13) ),]
#Add heights
subs= c(5,2, 13, "Te")
height= c("lower-mid","mid","high", "Tair")
te.sub$Height= height[match(te.sub$subsite, subs)]
#order factor
te.sub$Height= ordered(te.sub$Height, levels=c("lower-mid","mid","high", "Tair") )

#FIG 1A, time series
fig2a.bb<- ggplot(data=te.sub, aes(x=doy, y = MaxTemp_C))+
  geom_line(aes(color=Height), alpha=0.8) +
  labs(x = "day of year",y="maximum daily temperature (°C)")+
  scale_color_manual(values=vir.cols)+
  guides(color=FALSE)

#add weather station
fig.timeseries= fig2a.bb+ geom_line(data=dat1[which(dat1$year==2007 & dat1$j %in% 120:274),], aes(x=j, y =tmax, col=vir.cols[4]), alpha=0.8) +theme_classic(base_size = 14)

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
#Tair weather station frequency
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
subs= c(5,2, 13, "Te")
height= c("lower-mid","mid","high", "Tair")
pow$Height= height[match(pow$subsite, subs)]
#order factor
pow$Height= ordered(pow$Height, levels=c("lower-mid","mid","high", "Tair") )

table(pow$subsite)

#subsites
pow1= pow[pow$subsite %in% subs,]

fig2b<- ggplot(data=pow1, aes(x=log(freq), y = log(cyc_range/2) ))+geom_line(alpha=0.8, aes(color=Height)) +theme_classic(base_size = 14)+ 
  #guides(color=FALSE, size=FALSE)+
  geom_vline(xintercept=-2.639, color="gray")+geom_vline(xintercept=-1.946, color="gray")+geom_vline(xintercept=-3.40, color="gray")+geom_vline(xintercept=-5.9, color="gray")+
  labs(x = "log (frequency) (1/days)",y="log (amplitude)")+
  scale_color_manual("height", values=vir.cols)+
  annotate(geom="text", x=-2.1, y=-5, label="1 week", size=3, color="black",angle=90)+ 
  annotate(geom="text", x=-2.8, y=-5, label="2 weeks", size=3, color="black",angle=90)+ 
  annotate(geom="text", x=-3.6, y=-5, label="1 month", size=3, color="black",angle=90)+ 
  annotate(geom="text", x=-6.2, y=-5, label="1 year", size=3, color="black",angle=90)+
  ylim(range(-5.8,2.05))+
  theme(legend.position="bottom")

#count by site, subsite, year
te.count = te.max %>% group_by(year,site, subsite, zone, tidal.height..m.) %>% summarise( count=length(MaxTemp_C)  )
te.count= as.data.frame(te.count)
#find sub sites to use

#---------------
#estimate performance
#show effects of different time averaging
#From MusselTPCs_fits
#M. trossulus summer

tperf= function(x) {
  perf= beta_2012(x, tpc.beta[1], tpc.beta[2], tpc.beta[3], tpc.beta[4], tpc.beta[5])
  #set performance to zero outside TPC
  perf[which(!is.finite(perf))]=0
  return(perf)
}

#temporal averaging
te.sub= te.max2[which(te.max2$site=="BB" & te.max2$subsite %in% c(5,2,13) & te.max2$year==2007),]

#add Tair data
dat1.sub= dat1[which(dat1$year==2007),]
te.tair= matrix(NA, nrow=nrow(dat1.sub),ncol=ncol(te.sub))
te.tair= as.data.frame(te.tair)
colnames(te.tair)= colnames(te.sub)
te.tair$doy= dat1.sub$doy
te.tair$month= dat1.sub$month
te.tair$MaxTemp_C= dat1.sub$tmax
te.tair$site= "Tair"
#combine
te.sub= rbind(te.sub, te.tair)

#rolling averages
doys= 1:365

ts5= as.data.frame(cbind(doys, site=rep(5, length(doys)), temps=rep(NA, length(doys)) ))
match1= match(doys, te.sub[which(te.sub$subsite==5),]$doy)
inds= which(!is.na(match1))
ts5$temps[inds]= te.sub[which(te.sub$subsite==5),]$MaxTemp_C[match1[inds]]
ts5$t.wk= rollmean(ts5$temps, k=7, fill=NA)
ts5$t.2wk= rollmean(ts5$temps, k=14, fill=NA)
ts5$t.mo= rollmean(ts5$temps, k=30, fill=NA)
ts5$t.3mo= rollmean(ts5$temps, k=90, fill=NA)

ts2= as.data.frame(cbind(doys, site=rep(2, length(doys)), temps=rep(NA, length(doys)) ))
match1= match(doys, te.sub[which(te.sub$subsite==2),]$doy)
inds= which(!is.na(match1))
ts2$temps[inds]= te.sub[which(te.sub$subsite==2),]$MaxTemp_C[match1[inds]]
ts2$t.wk= rollmean(ts2$temps, k=7, fill=NA)
ts2$t.2wk= rollmean(ts2$temps, k=14, fill=NA)
ts2$t.mo= rollmean(ts2$temps, k=30, fill=NA)
ts2$t.3mo= rollmean(ts2$temps, k=90, fill=NA)

ts13= as.data.frame(cbind(doys, site=rep(13, length(doys)), temps=rep(NA, length(doys)) ))
match1= match(doys, te.sub[which(te.sub$subsite==13),]$doy)
inds= which(!is.na(match1))
ts13$temps[inds]= te.sub[which(te.sub$subsite==13),]$MaxTemp_C[match1[inds]]
ts13$t.wk= rollmean(ts13$temps, k=7, fill=NA)
ts13$t.2wk= rollmean(ts13$temps, k=14, fill=NA)
ts13$t.mo= rollmean(ts13$temps, k=30, fill=NA)
ts13$t.3mo= rollmean(ts13$temps, k=90, fill=NA)

#combine
ts= rbind(ts5, ts2, ts13)

# #add month
# te.mo= aggregate(te.sub[,c("MaxTemp_C","doy")], by=list(te.sub$subsite, te.sub$year, te.sub$month), FUN=mean, na.rm=TRUE)
# names(te.mo)[1:3]=c("subsite","year","month")
# 
# ts$ds=paste(ts$doys, ts$site, sep="_")
# ds=paste(round(te.mo$doy), te.mo$subsite, sep="_")
# ts$t.mo=NA
# match1= match(ts$ds, ds)
# inds= which(!is.na(match1))
# ts$t.mo[inds]= te.mo$MaxTemp_C[match1[inds]]

#estimate for robomussels at different heights and Tair
#restrict to time underwater?

#estimate performance
ts$p.d= tperf(ts$temps)
ts$p.wk= tperf(ts$t.wk)
ts$p.2wk= tperf(ts$t.2wk)
ts$p.mo= tperf(ts$t.mo)
ts$p.3mo= tperf(ts$t.3mo)

summary(ts)

#melt
ts.l= melt(ts, id.vars=c("doys","site"))
#label temp vs perf
ts.l$type="temperature"
ts.l$type[which(ts.l$variable %in% c("p.d","p.wk","p.2wk","p.mo","p.3mo"))]="performance"
#label timescale
ts.l$timescale= "day"
ts.l$timescale[which(ts.l$variable %in% c("p.wk","t.wk"))]= "week"
ts.l$timescale[which(ts.l$variable %in% c("p.2wk","t.2wk"))]= "2 week"
ts.l$timescale[which(ts.l$variable %in% c("p.mo","t.mo"))]= "month"
ts.l$timescale[which(ts.l$variable %in% c("p.3mo","t.3mo"))]= "3 month"
ts.l$timescale= ordered(ts.l$timescale, levels=c("day","week","2 week","month","3 month"))

#time series plot
fig.perf= ggplot(data=ts.l[ts.l$type %in% "performance" & ts.l$site==5,], aes(x=doys, y =value, color=timescale))+
  #facet_grid(site~., scales="free_y")+
  geom_line() +theme_classic(base_size = 14)+
  scale_color_viridis_d()+theme(legend.position="bottom")+
  ylab("clearance rate at lower-mid site (ml/min)")+xlab("day of year")+
  xlim(120,274)

#density plot
ggplot(data=ts.l[ts.l$type %in% "performance",], aes(x =value, color=timescale))+
  facet_grid(type~site, scales="free")+
  geom_density()

#=================================
#combine #boiler bay figures
library(patchwork)

#Figures
#TPC through seasons? TPCs with variability?

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PLoSextremes/figures/")
pdf("Fig2_musselTempPerf.pdf",height = 10, width = 10)
fig.timeseries +fig2b +fig.perf +plot_layout(ncol = 1)+ 
  plot_annotation(tag_levels = 'A')
dev.off()





