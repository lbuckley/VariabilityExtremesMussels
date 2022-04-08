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
#Site data
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

#----------------------
#PLOTS

#Time series
#clim2 = clim2 %>% group_by(Year,Site) %>% summarise(Min= mean(Min, na.rm=TRUE),Max= mean(Max, na.rm=TRUE),Mean= mean(Mean, na.rm=TRUE) )

#count by site, subsite, year
te.count = te.max %>% group_by(year,site, subsite) %>% summarise( count=length(MaxTemp_C)  )
te.count= as.data.frame(te.count)

#subset sites
te.max2= subset(te.max, te.max$site %in% c("SD","BB","PD") ) # "HS",
#te.max2= subset(te.max, te.max$lat %in% c(48.39137,44.83064,35.66582,34.46717) )

# USWACC	48.5494	-123.0059667	Colins Cove
# USWACP	48.45135	-122.9617833	Cattle Point
#* USWASD	48.39136667	-124.7383667	Strawberry Point

#* USORBB	44.83064	-124.06005	Boiler Bay

#* USCAPD	35.66581667	-121.2867167	Piedras
# USCAAG	34.46716667	-120.2770333	Alegria


#time series
te.max1= subset(te.max2, te.max2$year==2002)

#restrict to summer
#May 1 through September: 121:273 
te.max1= subset(te.max1, te.max1$doy>120 & te.max1$doy<274)

#ggplot(data=te.max1, aes(x=doy, y = MaxTemp_C, color=subsite ))+geom_line() +theme_bw()+facet_wrap(~site)
#by tidal height
#ggplot(data=te.max1, aes(x=doy, y = MaxTemp_C, color=height ))+geom_line() +theme_bw()+facet_wrap(~site)

#update labels
te.max1$labs= as.character(te.max1$site)
te.max1$labs[which(te.max1$labs=="BB")]<- "44.8° Boiler Bay, OR"
te.max1$labs[which(te.max1$labs=="PD")]<- "35.7° Piedras Blancas, CA"
te.max1$labs[which(te.max1$labs=="SD")]<- "48.4° Strawberry Point, WA"

#FIG 1A, time series
#by lat
fig2a<- ggplot(data=te.max1, aes(x=doy, y = MaxTemp_C, color=subsite ))+geom_line(alpha=0.8) +theme_classic()+
  facet_wrap(~labs, nrow=1)+ guides(color=FALSE)+labs(x = "Day of year",y="Maximum daily temperature (°C)")

#------------------
#FREQUENCY
# https://github.com/georgebiogeekwang/tempcycles/

#power spectrum
#x: frequency (1/days)
#y: log amplitude

fseq= exp(seq(log(0.001), log(1), length.out = 200))

sites= c("SD","BB","PD") #levels(te.max2$site)
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

#sort by frequency
pow= pow[order(pow$site, pow$subsite, pow$freq),]
pow$subsite= factor(pow$subsite)

#freq, amp plot

#add latitude
site.dat1=  te.max %>% group_by(site) %>% summarise( lat=lat[1],zone=zone[1],tidal.height..m.=tidal.height..m.[1],substrate=substrate[1] )
match1= match(pow$site, site.dat1$site)
pow$lat= site.dat1$lat[match1]

#update labels
pow$labs= as.character(pow$site)
pow$labs[which(pow$labs=="BB")]<- "44.8° Boiler Bay, OR"
pow$labs[which(pow$labs=="PD")]<- "35.7° Piedras Blancas, CA"
pow$labs[which(pow$labs=="SD")]<- "48.4° Strawberry Point, WA"

fig2b<- ggplot(data=pow, aes(x=log(freq), y = log(cyc_range/2) ))+geom_line(alpha=0.8, aes(color=subsite)) +theme_classic()+facet_wrap(~labs, nrow=1)+ guides(color=FALSE, size=FALSE)+
  geom_vline(xintercept=-2.639, color="gray")+geom_vline(xintercept=-1.946, color="gray")+geom_vline(xintercept=-3.40, color="gray")+geom_vline(xintercept=-5.9, color="gray")+
  labs(x = "log (frequency) (1/days)",y="log (amplitude)")+
  annotate(geom="text", x=-2.1, y=-5.5, label="1 week", size=3, color="black",angle=90)+ 
  annotate(geom="text", x=-2.8, y=-5.5, label="2 weeks", size=3, color="black",angle=90)+ 
  annotate(geom="text", x=-3.6, y=-5.5, label="1 month", size=3, color="black",angle=90)+ 
  annotate(geom="text", x=-6.2, y=-5.5, label="1 year", size=3, color="black",angle=90)+ylim(range(-5.8,1.2))
#add lines for 1 week, 2 week, month, year

# #save plots
# pdf("Fig2a.pdf",height = 4, width = 8)
# fig2a
# dev.off()
# 
# pdf("Fig2b.pdf",height = 4, width = 8)
# fig2b
# dev.off()

#----
# setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ICBClimateBiology\\figures\\") 
# pdf("Fig2.pdf",height = 8, width = 10)
# grid.newpage()
# pushViewport(viewport(layout=grid.layout(2,1)))
# vplayout<-function(x,y)
#   viewport(layout.pos.row=x,layout.pos.col=y)
# print(fig2a,vp=vplayout(1,1))
# print(fig2b,vp=vplayout(2,1))
# dev.off()

#========================================
#PLOT HEIGHTS
#subset to intertidal heights

fig2a<- ggplot(data=te.max1, aes(x=doy, y = MaxTemp_C, color=subsite ))+geom_line(alpha=0.8) +theme_classic()+
  facet_wrap(~labs, nrow=1)+ guides(color=FALSE)+labs(x = "Day of year",y="Maximum daily temperature (°C)")

fig2b<- ggplot(data=pow, aes(x=log(freq), y = log(cyc_range/2) ))+geom_line(alpha=0.8, aes(color=subsite)) +theme_classic()+facet_wrap(~labs, nrow=1)+ guides(color=FALSE, size=FALSE)+
  geom_vline(xintercept=-2.639, color="gray")+geom_vline(xintercept=-1.946, color="gray")+geom_vline(xintercept=-3.40, color="gray")+geom_vline(xintercept=-5.9, color="gray")+
  labs(x = "log (frequency) (1/days)",y="log (amplitude)")+
  annotate(geom="text", x=-2.1, y=-5.5, label="1 week", size=3, color="black",angle=90)+ 
  annotate(geom="text", x=-2.8, y=-5.5, label="2 weeks", size=3, color="black",angle=90)+ 
  annotate(geom="text", x=-3.6, y=-5.5, label="1 month", size=3, color="black",angle=90)+ 
  annotate(geom="text", x=-6.2, y=-5.5, label="1 year", size=3, color="black",angle=90)+ylim(range(-5.8,1.2))
#add lines for 1 week, 2 week, month, year

#subset
#BB, PD, SD
#subsites 10,12,15,21,3
te.max1$zone2= factor(te.max1$zone, levels=c("Low","Mid","High") )
tes= subset(te.max1, te.max1$site=="BB" & te.max1$subsite %in% c(2,13,10) )
summary(tes)

#tes= subset(te.max1, te.max1$site=="SD" & te.max1$subsite %in% c(9,12,10) )

table(tes$subsite, tes$zone)

fig2a<- ggplot(data=tes, aes(x=doy, y = MaxTemp_C, color=zone ))+geom_line(alpha=0.8) +theme_classic()+
  guides(color=FALSE)+labs(x = "Day of year",y="Maximum daily temperature (°C)")

#---
pow1= subset(pow, pow$site=="BB" & pow$subsite %in% c("2","13","10") ) #"10","3","4"
#13: High, #10 and 2 low

#pow1= subset(pow, pow$site=="SD" & pow$subsite %in% c("9","12","10") ) 

#Add heights
subs= c(2, 13, 10)
height= c("low","high","low")
pow1$Height= height[match(pow1$subsite, subs)]

table(pow1$subsite)

fig2b<- ggplot(data=pow1, aes(x=log(freq), y = log(cyc_range/2) ))+geom_line(alpha=0.8, aes(color=subsite)) +theme_classic()+ 
  guides(color=FALSE, size=FALSE)+
  geom_vline(xintercept=-2.639, color="gray")+geom_vline(xintercept=-1.946, color="gray")+geom_vline(xintercept=-3.40, color="gray")+geom_vline(xintercept=-5.9, color="gray")+
  labs(x = "log (frequency) (1/days)",y="log (amplitude)")+
  annotate(geom="text", x=-2.1, y=-5.5, label="1 week", size=3, color="black",angle=90)+ 
  annotate(geom="text", x=-2.8, y=-5.5, label="2 weeks", size=3, color="black",angle=90)+ 
  annotate(geom="text", x=-3.6, y=-5.5, label="1 month", size=3, color="black",angle=90)+ 
  annotate(geom="text", x=-6.2, y=-5.5, label="1 year", size=3, color="black",angle=90)+
  ylim(range(-5.8,2.05))

#===================================================
#Grasshopper

#BuckleyUnsynch> FinishedProjects>FinishedUNC>Te

#Load NSRDB station list
setwd("/Volumes/GoogleDrive/My Drive/BuckleyUnsynch/FinishedProjects/FinishedUNCProjects/Te/data/NSRDB/")
stats= read.csv("rockymntsNSRDB_select.csv")
stats=subset(stats, stats$Select==1) 

#CO sites: 
#Limon 1695m
#Eagle 1980m
#Alamosa 2296m
#Leadville 3026m

ids= stats[c(2,4,5,8),"USAF"]
stat.names= c("Limon","Eagle","Alamosa","Leadville")

#read daily data
for(st in 1:length(ids)){
  id=ids[st]
  dir=paste("/Volumes/GoogleDrive/My Drive/BuckleyUnsynch/FinishedProjects/FinishedUNCProjects/Te/data/NSRDB/", id,"/", sep="")
  
  setwd(dir)
  filename=paste(id,"daily_Grasshopper.csv", sep="")
  dat.day= read.csv(filename)
  dat.day$stat= id
  dat.day$station= stat.names[st]
  
  if(st==1) dat.all= dat.day
  if(st>1) dat.all= rbind(dat.all, dat.day)
}

#------------------
#FREQUENCY
# https://github.com/georgebiogeekwang/tempcycles/

#power spectrum
#x: frequency (1/days)
#y: log amplitude

fseq= exp(seq(log(0.001), log(1), length.out = 200))

sites= stat.names

pow.out= matrix(NA, length(sites),length(fseq) )

for(site.k in 1:length(sites))
{
  te.dat1= dat.all[which(dat.all$station==sites[site.k]),]
  pow.out[site.k,] <- spec_lomb_phase(te.dat1$Te.max, te.dat1$Julian, freq=fseq)$cyc_range
  
}

rownames(pow.out)<- sites
colnames(pow.out)<- 1:200

#to long format
pow= melt(pow.out)
colnames(pow)[1:3]=c("site","freq","cyc_range")

#correct freq values
pow$freq= fseq[pow$freq]

#sort by frequency
pow= pow[order(pow$site, pow$freq),]

#add elevation
elevs= c("1695m","1980m","2296m","3026m")
pow$Elevation= elevs[match(pow$site, sites)]

#CO sites: 
#Limon 1695m
#Eagle 1980m
#Alamosa 2296m
#Leadville 3026m

#freq, amp plot
fig2b_CO<- ggplot(data=pow, aes(x=log(freq), y = log(cyc_range/2) ))+geom_line(alpha=0.8, aes(color=Elevation)) +theme_classic()+
  guides(color=FALSE, size=FALSE)+
  geom_vline(xintercept=-2.639, color="gray")+geom_vline(xintercept=-1.946, color="gray")+geom_vline(xintercept=-3.40, color="gray")+geom_vline(xintercept=-5.9, color="gray")+
  labs(x = "log (frequency) (1/days)",y="log (amplitude)")+
  annotate(geom="text", x=-2.1, y=-5.5, label="1 week", size=3, color="black",angle=90)+ 
  annotate(geom="text", x=-2.8, y=-5.5, label="2 weeks", size=3, color="black",angle=90)+ 
  annotate(geom="text", x=-3.6, y=-5.5, label="1 month", size=3, color="black",angle=90)+ 
  annotate(geom="text", x=-6.2, y=-5.5, label="1 year", size=3, color="black",angle=90)+
  ylim(range(-5.8,2.05))
#add lines for 1 week, 2 week, month, year

#---------------------------------
#Load NRCS station list
setwd("/Volumes/GoogleDrive/My Drive/BuckleyUnsynch/FinishedProjects/FinishedUNCProjects/Te/data/NRCS/")
stats= read.csv("NRCS_Stations.csv")
stats=stats[stats$Select==1,]
stats=stats[order(stats$State, stats$Elev_m),]
lab=paste(stats$SiteNum," ", stats$State,", ", round(stats$Elev_m),"m", sep="")

ids= stats$SiteNum
stat.names= stats$Elev_m

#read daily data
for(st in 1:length(ids)){
  id=ids[st]
  
  setwd("/Volumes/GoogleDrive/My Drive/BuckleyUnsynch/FinishedProjects/FinishedUNCProjects/Te/data/NRCS/")
  filename=paste(id,"daily_Grasshopper.csv", sep="")
  dat.day= read.csv(filename)
  dat.day$stat= id
  dat.day$station= stat.names[st]
  
  if(st==1) dat.all= dat.day
  if(st>1) dat.all= rbind(dat.all, dat.day)
}

#------------------
#FREQUENCY
# https://github.com/georgebiogeekwang/tempcycles/

#power spectrum
#x: frequency (1/days)
#y: log amplitude

fseq= exp(seq(log(0.001), log(1), length.out = 200))

sites= stat.names

pow.out= matrix(NA, length(sites),length(fseq) )

for(site.k in 1:length(sites))
{
  te.dat1= dat.all[which(dat.all$station==sites[site.k]),]
  pow.out[site.k,] <- spec_lomb_phase(te.dat1$Te.max, te.dat1$Julian, freq=fseq)$cyc_range
  
}

rownames(pow.out)<- sites
colnames(pow.out)<- 1:200

#to long format
pow= melt(pow.out)
colnames(pow)[1:3]=c("site","freq","cyc_range")

#correct freq values
pow$freq= fseq[pow$freq]

#sort by frequency
pow= pow[order(pow$site, pow$freq),]
pow= na.omit(pow)

#add site info
pow$state= stats$State[match(pow$site, stats$Elev_m)]
pow$site.name= stats$SCANsite[match(pow$site, stats$Elev_m)]

#freq, amp plot
fig2b_trop<- ggplot(data=pow[pow$state=="HI",], aes(x=log(freq), y = log(cyc_range/2) ))+geom_line(alpha=0.8, aes(color=factor(site))) +theme_classic()+
  guides(color=FALSE, size=FALSE)+
  geom_vline(xintercept=-2.639, color="gray")+geom_vline(xintercept=-1.946, color="gray")+geom_vline(xintercept=-3.40, color="gray")+geom_vline(xintercept=-5.9, color="gray")+
  labs(x = "log (frequency) (1/days)",y="log (amplitude)")+
  annotate(geom="text", x=-2.1, y=-5.5, label="1 week", size=3, color="black",angle=90)+ 
  annotate(geom="text", x=-2.8, y=-5.5, label="2 weeks", size=3, color="black",angle=90)+ 
  annotate(geom="text", x=-3.6, y=-5.5, label="1 month", size=3, color="black",angle=90)+ 
  annotate(geom="text", x=-6.2, y=-5.5, label="1 year", size=3, color="black",angle=90)+
  ylim(range(-5.8,2.05))
#add lines for 1 week, 2 week, month, year

#=================================
#combine
library(patchwork)

#add titles
fig2b= fig2b +ggtitle('a. OR Intertidal Heights')
fig2b_CO= fig2b_CO +ggtitle('b. CO Elevations (1695-3026m)')
fig2b_trop= fig2b_trop +ggtitle('c. HI Elevations (353-926m)')

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/Proposals/NSF_BoCP/figures/")
pdf("Fig_freq.pdf",height = 6, width = 10)
fig2b + fig2b_CO + fig2b_trop + plot_layout(ncol = 3)
dev.off()

#boiler bay

#--------------------------
#Figures
#TPC through seasons? TPCs with variability?




