
############################################################################################################################################################################
# ONE STATION METABOLISM CODE for Methods in Stream Ecology, 3rd edition (Hall & Hotchkiss)

#This code provides all of the functions necessary to estimate single station open-channel stream metaboism from diel O2.  
#It gives you light data if you need it.  
#It also allows calculating gas exchange from nighttime regression.  
#The methods to estimate metabolism are via non-linear minimization of the -log likelihood of the modelled O2 vs the O2 data.  
#There are two appproaches: one with unknown gas exchange (K; function rivermetabK), and one with known K (function rivermetab). Most will use the unknown K function, with the caveat that solving for 3 parameters (GPP, ER, K) may give an over parameterized model or multiple solutions for ER and K, and thus uncertain estimates of ER (GPP is fairly robust to uncertainty in K).
#We also include 2 data files (S4 & S5)

#Bob Hall 7 Oct 2015.
#updated annotation, Erin Hotckhiss 2016-10-31
############################################################################################################################################################################


###########################################
# [1] DATA IMPORT AND MANAGEMENT
###########################################	

## Call and name the O2 logger data. Use your own path specific to your computer.
spring<-read.csv("/Users/bobhall/Rivers/spring14.csv")  		##data from Spring Creek WY, across the street from Bob's house
french<-read.csv("/Users/mpeipoch/Dropbox/b.mscripts/biofilmMARC/Datasets/Metabolism Data/meta model/Online Supplement Chapter 34_Example Data File_French Creek WY.csv")	##data from French Creek WY, Hotchkiss & Hall (2015) Ecology

marc<-read.table("clipboard",header=TRUE, sep="\t")

## for this code to work you need your data looking like the example datasets below:
	## temp = water temperature in Celcius
	## oxy = dissolved oxygen in mg/L
head(spring)  	##time is in seconds since 1970 and in UTC
head(french)	##time is date and time in separate columns. Follow the date and time format exactly (i.e. seconds on the time as hh:mm:ss; mm/dd/yy for date)
head(marc)

## TO USE DATES & TIMES IN R: First make make date and time a chron object. Load package chron from the R website and look online or via "help(chron)" for how it works.  
## How you make a chron object depends on your data. PME MiniDOTS give data in seconds since 1970 which make for a really easy translation into a chron object.  Hydrolabs are a bit tougher.  Below are two ways of doing it.
## If not already installed, first enter command "install.packages("chron")"

## Load chron package to work with dates/times in R
install.packages("chron")
library(chron)

## PME Minidot: make time correction (versus UTC) and chron object
## Minidots report time in UTC and give option to adjust time with miniDOT Plot software <or> can do it here in R with original CAT file (we prefer to do it in R).
spring$dtime<-chron(spring$time/86400)-(6/24)	##the subtraction of 6/24 (0.25 d) is because local time in WY is 6 h behind UTC in summer.  

## The below makes a chron object for dates and times where each is in its own column
french$dtime<-chron(dates=as.character(french $date), times=as.character(french $time))
marc$dtime<-chron(dates=as.character(marc $date), times=as.character(marc $time))

## ONLY IF NEEDED: if date and time are in the same column; split the string (under the header 'datetime'). Make sure your datetime is to the second!  
## Format is mm/dd/yy hh:mm:ss  The function strsplit will split at the space bewteen date and time. 
## Change "YOURDATA" to whatever you named your data file (as we did for "spring" or "french")
dtparts<-t(as.data.frame(strsplit(as.character(YOURDATA$datetime),' ')))
row.names(dtparts)=NULL
## make chron object
YOURDATA$dtime<-chron(dates=dtparts[,1], times=dtparts[,2])

## Check new chron object dtime and O2 data
## here is where you can target and remove outliers if you want to check before running models 
plot(spring$dtime, spring$oxy)
plot(french$dtime, french$oxy)
plot(marc$dtime, marc$oxy)

####### END data important and management ####### 

####### 


###########################################
# [2] LOAD O2 SATURATION FUNCTION
###########################################	

## LOAD function to calculate oxygen saturation given water temperature and barometric pressure.

## oxygen saturation. From García and Gordonn (1992).
	## calculates mL gas / dm^3 water; the below equation converts to mg/L.    
	## BP in mm Hg. 
	## 1.42905 is 'correct' and accounts for the fact that O2 is a non-ideal gas.  
	## u is the vapor pressure of water

osat<- function(temp, bp) {
	
	tstd<-log((298.15-temp) / (273.15 + temp))
	
	a0<-2.00907
	a1<-3.22014
	a2<-4.0501
	a3<-4.94457
	a4<- -0.256847
	a5<- 3.88767
	
	u<-10^(8.10765-(1750.286/(235+temp)))
	
	sato<-(exp(a0 + a1*tstd + a2*tstd^2 + a3*tstd^3 + a4*tstd^4+a5*tstd^5))*((bp-u)/(760-u))*1.42905
	sato
}

#end of function

####### END loading O2 SAT calc ####### 


###########################################
# [3] LOAD BAROMETRIC PRESSURE FUNCTION
###########################################	

## Function to correct barometric pressure for altitude. From Colt (2012).
## This function gives bp in mmHg for altitude given nearby measurement of standardized barometric pressure. 
	## temp is degC 
	## alt is m
	## bpst is in inches of Hg and the sort of data you will find from U.S. weather stations.  Delete the *25.4 if you in the rest of the world where pressure is given in metric units

####function returns mm of Hg
bpcalc<- function(bpst, alt) {
bpst*25.4*exp((-9.80665*0.0289644*alt)/(8.31447*(273.15+15)))
}
#end of function


####### END loading BP calc function ####### 


###########################################
# [4] LOAD GAS EXCHANGE (K) FUNCTIONS
###########################################	

# NOTE: The functions you use will depend on the methods used to estimate K (O2, propane, SF6). The temperature correction (Kcor) is embedded in the models below, but the Kcor function must be loaded into R before running the model.

# UNITS are day^(-1)
	
## This code does the opposite of Kcor below; it estimates K600 for KO2 at a given temperature. From Wanninkhof (1992).
K600fromO2<-function (temp, KO2) {
    ((600/(1800.6 - (120.1 * temp) + (3.7818 * temp^2) - (0.047608 * temp^3)))^-0.5) * KO2
}
#end of function

## This calculates K600 from K measured in a propane addition. From Raymond et al. (2012).
K600frompropane<-function (temp, Kpropane) {
    ((600/(2864 - (154.14 * temp) + (3.791 * temp^2) - (0.0379 * temp^3)))^-0.5) * Kpropane
}
#end of function

## This calculates K600 from K measured in a SF6 addition. From Raymond et al. (2012).
K600fromSF6<-function(temp, KSF6) {
	((600/(3255.3-(217.13*temp)+(6.837*temp^2)-(0.08607*temp^3)))^-0.5)*KSF6
}
#end of function

## This function calculates KO2 at any given tempertaure from K600. via schmidt number scaling.  The scaling equation if From Jähne et al. (1987), Schmidt number conversions from Wanninkhof et al. 1992.
Kcor<-function (temp,K600) {
	K600/(600/(1800.6-(temp*120.1)+(3.7818*temp^2)-(0.047608*temp^3)))^-0.5
}
#end of function

####### END loading K functions ####### 


###########################################
# [5] NIGHTTIME REGRESSION to estimate K -- OPTIONAL
###########################################	
	
### nighttime regression code to estimate K. Approach as in Hornberger and Kelly (1975).
	## o2 file is your oxygen data (defined in subsetting)
	## bp is barometric pressure in mm Hg for your site, 
	## ts is the time step in MINUTES (not days as in metabolism code below)
	
nightreg<-function(o2file, bp, ts){

	temp<-o2file$temp
	oxy<-o2file$oxy

# moving average on oxy data
oxyf1<-filter(o2file$oxy,rep(1/3,3), sides=2)

# trim the ends of the oxy data
oxyf2<- oxyf1[c(-1,-length(oxyf1))]

# calculate delO/delt; convert to units of days by ts in min-1*1440
deltaO2<-((oxyf2[-1]-oxyf2[-length(oxyf2)])/ts)*1440

# Trim the first two and last one from the temp data to match the filter oxy data
temptrim<-temp[c(-2:-1,-length(temp))]

# calc the dodef
satdef<-osat(temptrim,bp)-oxyf2[-1]

# fit linear model and plot using linear model fitting (lm) and abline functions in R stats package
nreg<-lm(deltaO2~satdef)
plot(satdef,deltaO2)
abline(nreg)

# use coef function in R stats package to get lm coefficients
coeff<-coef(nreg)

# output gives lm coeff and K600 (converted from nighttime regression estimate of KO2)
out<-list(coeff, K600fromO2(mean(temp), coeff[2]))
out

}
#end of function

# NOTE: this approach works better for some sites/dates than others; always check that your model fit is good and that your K600 estimate makes sense!

#Call as:  The first argument in the function defines when to pull data.  In this case on 10/27/204 (for spring creek) between 18:05 and 23:00
	
	nightreg(spring[spring$dtime>=as.numeric(chron(dates="10/27/14", times="18:05:00")) & spring$dtime<=as.numeric(chron(dates="10/27/14", times="23:00:00")), ], bp=595, ts=10)
	
	nightreg(french[french$dtime>=as.numeric(chron(dates="09/21/12", times="19:40:00")) & french$dtime<=as.numeric(chron(dates="09/21/12", times="23:00:00")), ], bp=523, ts=5)
	
	nightreg(french[french$dtime>=as.numeric(chron(dates="09/22/12", times="19:40:00")) & french$dtime<=as.numeric(chron(dates="09/22/12", times="23:00:00")), ], bp=523, ts=5)

####### END nighttime regression calculation of K ####### 


###########################################
# [6] LIGHT MODEL -- OPTIONAL
###########################################

## now make up light data if you don't have it

## From Yard et al. (1995) Ecological Modelling.  Remember your trig?  
## calculate light as umol photon m-2 s-1.
## Arguments are:  
	## time = a date and time input (i.e. a chron object) 
	## lat = latitude of field site
	## longobs = longitude of field site
	## longstd = standard longitude of the field site (NOTE: watch daylight savings time!!!). For PST, longstd is be 120 degrees. But during PDT it is 105 degrees. MST is 105 deg. MDT is 90. 
	## year = the year for which you collected data and is entered as "2013-01-01"

# convert degrees to radians
radi<-function(degrees){(degrees*pi/180)}

# function to estimate light
lightest<- function (time, lat, longobs, longstd, year ) {
	
	jday<-as.numeric(trunc(time)-as.numeric(as.Date(year)))
	E<- 9.87*sin(radi((720*(jday-81))/365)) - 7.53*cos(radi((360*(jday-81))/365)) - 1.5*sin(radi((360*(jday-81))/365))
	LST<-as.numeric (time-trunc(time))
	ST<-LST+(3.989/1440)*(longstd-longobs)+E/1440
	solardel<- 23.439*sin(radi(360*((283+jday)/365)))
	hourangle<-(0.5-ST)*360
	theta<- acos( sin(radi(solardel)) * sin(radi(lat)) + cos(radi(solardel)) * cos(radi(lat)) * cos(radi(hourangle)) )
	suncos<-ifelse(cos(theta)<0, 0, cos(theta))
GI<- suncos*2326
GI	
	
}
#end of function

## onespecific time just to see if function works
lightest(time=chron(15500.5), lat=41.33,  longobs=105.7, longstd= 105, year="2013-01-01")

## estimate light for spring creek and french creek using lightest function
spring$light<- lightest(time=spring$dtime, lat=41.33,  longobs=105.6, longstd= 90, year="2014-01-01")
french$light<- lightest(time=french$dtime, lat=41.33,  longobs=106.3, longstd= 90, year="2012-01-01")
marc$light<- lightest(time=marc$dtime, lat=45.5,  longobs=113.4, longstd= 90, year="2014-01-01")

# check that modeled light makes sense
plot(french$dtime, french$light)
plot(marc$dtime, marc$light)

####### END light model ####### 


###########################################
# [7] MODEL v1 - RIVERMETABK; solves for GPP, ER, AND K
##
# [7a] LOAD RIVERMETABK functions
###########################################

## VERSION 1 of possible metabolism model

##### now estimate metabolism

## parameter names (and units) for BOTH rivermetab functions (solving for or fixing K600):
	## oxy.mod = modeled O2 (mg/L)
	## oxy = O2 data (mg/L)
	## MET[1] = GPP = estimated daily gross primary production (g O2 m-2 d-1); + flux of O2 production
	## MET[2] = ER = estimated daily ecosystem respiration (g O2 m-2 d-1); - flux of O2 consumption 
	## MET[3] = K = estimated K600 (d-1)
	## z = estimated site depth (m)
	## Kcor = KO2 corrected for temperature, calculated from K600 (d-1)
	## bp = barometric pressure (mmHg)
	## temp = water temperature (C)
	## light = PAR data (or modlight from above light function)
	## ts = time step of O2 data from logger (10 min intervals --> units are day, so 10/1440=0.06944)

# Model to calculate GPP, ER and K simultaneously
# This model is advantageous in the sense that one can estimate K from the data, but beware that it may be an overparameterized model.  
# Note that you have the opportunity below to use your light data or modeled light (here as data$light estimated for spring creek and french creek with the light model function above). You decide and modify the function and data names accordingly.

# rivermetabK is the master function that calls the MLE function (onestationmleK) and plotting function (onestationplot) below
rivermetabK<-function(o2file, z, bp, ts){
	
	##pull data out of loaded o2file (subset in CALL function) and give it a name. 
	temp<-o2file$temp
	oxy<-o2file$oxy
	light<-o2file$light

	##This calls onestationmleK (below) to calculate metabolism by non-linear minimization. We use nlm() function in R stats package; for more information see "help(nlm)"  The first argument is the function to be minimized, the second defines the starting values.  The function that is minimized (onestationmleK, see below) always has as its first argument, a vector of parameters (MET), and second argument a vector of starting values for those parameters, p=c(3,-5, 10).
	
	river.mle<-nlm(onestationmleK, p=c(3,-5, 10), oxy=oxy, z=z,temp=temp,light=light, bp=bp, ts=ts)
	
	##plot modeled and measaured O2 given MLE estimates of GPP, ER, and K600.  It calls a function below onestationplot()
	
	onestationplot(GPP=river.mle$estimate[1], ER=river.mle$estimate[2],oxy=oxy,z=z,temp=temp,light=light, K=river.mle$estimate[3], bp=bp, ts=ts)
	
	##return GPP, ER, K600, and MLE values
	b <- list(GPP=river.mle$estimate[1], ER=river.mle$estimate[2], K600=river.mle$estimate[3], neglogL=river.mle$minimum[1])
	b
}
# end of function


# This function returns the negative log likelihood value given O2 data and estimates of GPP, ER, and K (which is vector MET); is included in master rivermetabK function above
onestationmleK<-function(MET,temp, oxy, light, z, bp, ts) {

	# create new vector for modeled O2
	oxy.mod<-numeric(length(data))
	# give starting value from oxygen data; this is the only time O2 data is used to model GPP and ER
	oxy.mod[1]<-oxy[1]

	# this is the metabolism equation as in Van de Bogert et al 2007 L&OMethods
	for (i in 2:length(oxy)) {oxy.mod[i]<-oxy.mod[i-1]+((MET[1]/z)*(light[i]/sum(light)))+ MET[2]*ts/z+(Kcor(temp[i],MET[3]))*ts*(osat(temp[i],bp)-oxy.mod[i-1]) }
	
	##below is MLE calculation; output is -log likelihood
	# diff in measured and modeled O2
	sqdiff<-(oxy-oxy.mod)^2 
	# likelihood function
	L <- length(oxy) * (log(((sum(sqdiff)/length(oxy))^0.5))+0.5*log(6.28)) + ((2*sum(sqdiff)/length(oxy))^-1) * sum(sqdiff)
	L
}
# end of function

# this function plots modeled O2 and O2 data from estimates of daily GPP, ER, and K; is included in master rivermetabK function above
# Calls same metabolism equation as in mle function, but plots modeled O2 as a function of GPP, ER, and K estimates from mle
# use this to visually assess your model estimates (should be good agreement between modeled and measured O2)
onestationplot<-function(GPP, ER, oxy, z, temp, K, light, bp, ts) {
	
	oxy.mod<-numeric(length(oxy))
	oxy.mod[1]<-oxy[1]
	
	# this is the metabolism equation as in Van de Bogert et al (2007) L&OMethods
	for (i in 2:length(oxy)) { oxy.mod[i]<-oxy.mod[i-1]+((GPP/z)*(light[i]/sum(light)))+ ER*ts/z+(Kcor(temp[i],K))*ts*(osat(temp[i],bp)-oxy.mod[i-1]) }

	plot(seq(1:length(oxy)),oxy.mod, type="l",xlab="Time", ylab="Dissolved oxygen  (mg/L)", cex.lab=1.5, cex.axis=1.5, lwd=2 )
	points(seq(1:length(oxy)),oxy)
}
# end of function

####### END LOADING rivermetabK function ####### 


###########################################
# [7b] CALL RIVERMETABK function
###########################################

## Call as:  z is depth in m, bp is im mmHg, ts is time steps in days.

# for spring creek data; 10/28/14
rivermetabK(o2file=spring[ spring $dtime>=as.numeric(chron(dates="10/27/14", times="22:00:00")) & spring $dtime<=as.numeric(chron(dates="10/29/14", times="06:00:00")), ], z=0.18, bp=595, ts=0.006944)

## for french creek data; 09/21/12
rivermetabK(o2file=french[french$dtime>=as.numeric(chron(dates="09/20/12", times="22:00:00")) & french $dtime<=as.numeric(chron(dates="09/22/12", times="06:00:00")), ], z=0.18, bp=526, ts=0.003422)

## OPTIONAL: update to run an additional day from the french creek or your own dataset.... 

####### END example CALLs of rivermetabK function ####### 


###########################################
# [8] MODEL v2 - RIVERMETAB; solves for GPP, ER
##
# [8a] LOAD RIVERMETAB function w/ FIXED K
###########################################

## VERSION 2 of possible metabolism model 
## Here is same function but for fixed and known K600 from propane, SF6, nighttime regression, or other method
rivermetab<-function(o2file, z, bp, ts, K){

	##pull data out of loaded o2file (subset in CALL function) and give it a name. 
	temp<-o2file$temp
	oxy<-o2file$oxy
	light<-o2file$light
	
	##calculate metabolism by non linear minimization of MLE function (below)
	river.mle<-nlm(onestationmle, p=c(3,-5), oxy=oxy, z=z,temp=temp,light=light, bp=bp, ts=ts, K=K)
	
	##plot data; uses same plot function as given for rivermetabK above (relisted below to keep each model as own unit)
	onestationplot(GPP=river.mle$estimate[1], ER=river.mle$estimate[2],oxy=oxy,z=z,temp=temp,light=light, K=K, bp=bp, ts=ts)
	
	##return GPP, ER, and MLE value
	b<-list(GPP= river.mle$estimate[1], ER= river.mle$estimate[2],  neglogL= river.mle$minimum[1])
	b
}
# end of function

# function returns the likelihood value given O2 data and estimates of GPP, ER (which is vector MET); is included in master rivermetab function above; K600 is fixed
onestationmle<-function(MET,temp, oxy, light, z, bp, ts, K) {

	oxy.mod<-numeric(length(data))
	oxy.mod[1]<-oxy[1]

	# this is the metabolism equation as in Van de Bogert et al 2007 L&OMethods
	for (i in 2:length(oxy)) {oxy.mod[i]<-oxy.mod[i-1]+((MET[1]/z)*(light[i]/sum(light)))+ MET[2]*ts/z+(Kcor(temp[i],K))*ts*(osat(temp[i],bp)-oxy.mod[i-1]) }

	## below is MLE calculation; output is likelihood
	# diff in measured and modeled O2
	sqdiff<-(oxy-oxy.mod)^2 
	# likelihood function
	L <- length(oxy)*(log(((sum(sqdiff)/length(oxy))^0.5)) +0.5*log(6.28)) + ((2*sum(sqdiff)/length(oxy))^-1)*sum(sqdiff)
	L
}
#end of function

# (As in rivermetabK)
# this function plots modeled O2 and O2 data from estimates of daily GPP and ER; is included in master rivermetab function above
# Calls same metabolism equation as in mle function, but plots modeled O2 as a function of GPP, ER estimates from mle
# use this to visually assess your model estimates (should be good agreement between modeled and measured O2)
onestationplot<-function(GPP, ER, oxy, z, temp, K, light, bp, ts) {
	
	oxy.mod<-numeric(length(oxy))
	oxy.mod[1]<-oxy[1]
	
	# this is the metabolism equation as in Van de Bogert et al (2007) L&OMethods
	for (i in 2:length(oxy)) { oxy.mod[i]<-oxy.mod[i-1]+((GPP/z)*(light[i]/sum(light)))+ ER*ts/z+(Kcor(temp[i],K))*ts*(osat(temp[i],bp)-oxy.mod[i-1]) }

	plot(seq(1:length(oxy)),oxy.mod, type="l",xlab="Time", ylab="Dissolved oxygen  (mg/L)", cex.lab=1.5, cex.axis=1.5, lwd=2 )
	points(seq(1:length(oxy)),oxy)
}
# end of function


####### END loading rivermetab function ####### 


###########################################
# [8b] CALL RIVERMETAB function
###########################################

#Call as:

## for spring creek data; 10/28/14
rivermetab(o2file=spring[ spring $dtime>=as.numeric(chron(dates="10/27/14", times="22:00:00")) & spring $dtime<=as.numeric(chron(dates="10/29/14", times="06:00:00")), ], z=0.18, bp=595, K=29, ts=0.006944)

## for french creek data; 09/21/12
rivermetab(o2file=french[french$dtime>=as.numeric(chron(dates="09/20/12", times="22:00:00")) & french $dtime<=as.numeric(chron(dates="09/22/12", times="06:00:00")), ], z=0.18, bp=526, ts=0.003422, K=35)


##K and pb are not updated
rivermetab(o2file=marc[marc$dtime>=as.numeric(chron(dates="06/05/14", times="10:00:00")) & marc $dtime<=as.numeric(chron(dates="06/06/14", times="10:00:00")), ], z=0.54, bp=526, ts=0.003422, K=5.25
          )

## OPTIONAL:update to run an additional day from the french creek or your own dataset.... 

## NOTE: your choice of K=... will depend on previous estimates from propane, SF6, nighttime R, or other method

####### END example CALLs of rivermetab function ####### 


############################################################################################################################################################################


# So...what's next? 
# Please refer to Hotchkiss & Hall (Methods in Stream Ecology 3rd edition book chapter) for some guidelines on assessing model output!

#Our code that we supply here is useful for one-off metabolism estimates and for classroom instruction.  Once you undertsand what is 'under the hood' of this model and you want to conduct research on stream metabolism, we recommend that you apply streamMetabolizer.  This package is a flexible way to fit long or short time series of oxygen data.  Computational routines  in this package will evolve. Because it is an R pacakge, it will be possible to create replicable work with this tool.
#See
#Alison P. Appling, Robert O. Hall, Maite Arroita and Charles B. Yackulic (2016). streamMetabolizer:Models for Estimating Aquatic Photosynthesis and Respiration. R package version 0.9.19.  https://github.com/USGS-R/streamMetabolizer


############################################################################################################################################################################


###########################################
# [9] REFERENCES CITED 
###########################################

# Colt, J. 2012 Dissolved Gas Concentration in Water: Computation as Functions of Temperature, Salinity and Pressure. Elsevier.

# García, H.E., and L. I. Gordon. 1992. Oxygen solubility in seawater: Better fitting equations. Limnology & Oceanography 37: 1307-1312.

# Hornberger, G. M., and M. G. Kelly. 1975. Atmospheric reaeration in a river using productivity analysis. Journal of the Environmental Engineering Division 101: 729-739.

# Hotchkiss, E.R., and R. O. Hall. 2015. Whole-stream 13C tracer addition reveals disctint fates of newly fixed carbon. Ecology 96: 403-416.

# Jähne, B., et al. 1987. On the parameters influencing air-water gas exchange. Journal of Geophysical Research 92: 1937-1949.

# Raymond, P. A., et al. 2012. Scaling the gas transfer velocity and hydraulic geometry in streams and small rivers. Limnology & Oceanography: Fluids & Environments 2: 41–53.

# Van de Bogert, M. C., S. R. Carpenter, J. J. Cole, and M. L. Pace. 2007. Assessing pelagic and benthic metabolism using free water measurements. Limnology & Oceanography: Methods 5: 145-155.

# Wannikhof, R. 1992. Relationship between wind speed and gas exchange over the ocean. Journal of Geophysical Research 92: 7373-7387.

# Yard, M., et al. 2005. Influence of topographic complexity on solar insolation estimates for the Colorado River, Grand Canyon, AZ. Ecological Modelling 183: 157-172.

####### END references cited #######


############################################################################################################################################################################
	
	
	