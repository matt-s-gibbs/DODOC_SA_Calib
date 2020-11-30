library(zoo)
library(hydroGOF)
library(tidyverse)
library(optimx)

################
#
# MAKE SURE PATHS IN SERVER.BAT AND RUN.BAT ARE CORRECT, AND THEN RUN SERVER.BAT. OUTPUT FILE PATH SHOULD BE THE SAME AS THE R WORKING DIRECTORY (TYPE getwd() )
#
##############

WriteInputSet<-function(x)
{
  DOCConsumptionrate<-0.2
  Reaerationrate<-x[1]
  RedGumRate<-x[2]
  
  #must match order in files in VegetationAreas folder:
  #blackbox, high shrub,lignum,low shrub, redgum
  #kg/ha, not g/m2 (*10)
  accumulationrates<-c(0.8,0.1,0.4,0.1,1)*RedGumRate #these fractions come from field data, average load under different veg types as a proporation of red gum.
  
  ##########
  #
  # Apply accumulation rates and write Source Input Set
  #
  #########
  
  storagenodes<-list.files("ChowVegetationAreas/",".csv")
  
  #read in data
  vegdata<-lapply(storagenodes,function(x) read.csv(paste0("ChowVegetationAreas/",x)))
  
  #remove .csv from storagenodes
  storagenodes<-sapply(storagenodes,function(x) substr(x,1,nchar(x)-4))
  
  con<-file("../Chowilla/InputSet.txt","w")
  
  writeLines(paste0("Functions.Functions.$Blackwater.f_DOCConsumption.Expression=",DOCConsumptionrate),con)
  writeLines(paste0("Functions.Functions.$Blackwater.f_Reaeration.Expression=",Reaerationrate),con)
  
  for(i in 1:length(vegdata))
  {
    areas<-vegdata[[i]][,-1]
    levels<-vegdata[[i]][,1]
    
    accumulation<-apply(areas,1,function(x) sum(x*accumulationrates)/sum(x))
    sep<-c(rep("][",length(accumulation)-1),"]]")
    writeLines(paste0(storagenodes[i],".Constituent Data.DOC.Processing Model.Leaf A=[[",paste0(levels," ",round(accumulation,2),sep,collapse="")),con)
    
  }
  close(con)
}

PlotResults<-function(file,index)
{
  
  
  #outputs to compare. These need to match the observed files in the ObservedData folder, and the node names in the model
  DOCsites<-data.frame(Obs=c("River Murray Lock 5","River Murray Renmark Sample Pump","River Murray 11km Downstream Lock 6","River Murray 8km downstream Lock 6"),
                       Sim=c("Lock 5  426512","Renmark  426628","River Murray DS Chowilla Creek","Confluence 29"),stringsAsFactors=FALSE)
  
  DOsites<-data.frame(Obs=c("A4261160","A4261107","A4261224","A4261168","A4260704","A4260703"),
                      Sim=c("Gum Flat 281","Chowilla Regulator 254","Chowilla Regulator 254","River Murray DS Chowilla Creek","River Murray DS Chowilla Creek","Renmark  426628"),stringsAsFactors=FALSE)
  
  
  #########
  #
  # Read in observed data
  #
  ########
  
  DOdat<-read.csv("ChowObservedData/DO_Obs.csv",as.is=TRUE)
  DOdat<-zoo(DOdat[,-1],as.Date(DOdat[,1],"%d/%m/%Y"))
  
  #limit by saturation
  DOsat<-read.csv("ChowObservedData/DOSat.csv",as.is=TRUE)
  DOsat<-zoo(DOsat[,-1],as.Date(DOsat[,1],"%d/%m/%Y"))
  
  DOsat<-window(DOsat,start=start(DOdat),end="2017-10-01")
  DOdat<-window(DOdat,end=end(DOsat))
  
  X<-apply(DOdat,2,function(x) pmin(x,DOsat))
  DOdat<-zoo(X,as.Date(rownames(X)))
  
  DOCdat<-read.csv("ChowObservedData/DOCLongData.csv",as.is=TRUE)
  DOCdat<-data.frame(Index=DOCdat$SampledDate,Series=DOCdat$DetailedSamplingPointDescription,Value=DOCdat$ResultValue)
  DOCdat$Index<-as.Date(DOCdat$Index,"%d/%m/%Y")
  
  
  outCols<-as.numeric(read_csv(file,skip=9,n_max=1,col_names=FALSE))
  Cols<-read_csv(file,skip=10,col_names = FALSE,n_max=outCols)
  
  Sim<-read_csv(file,skip=12+outCols,col_names=FALSE) #start row untested, should search for EOH
  Sim<-Sim %>% setNames(c("Date",paste0(Cols$X7,Cols$X11)))
  
  #Source auto output outputs in kg/m3, need to multiply by 1000
  DOSim<-Sim %>% select(contains("Constituents@DO@Downstream")) %>% mutate_all(funs(.*1000))
  DOCSim<-Sim %>% select(contains("Constituents@DOC@Downstream")) %>% mutate_all(funs(.*1000))
  RefractSim<-Sim %>% select(contains("Constituents@Refract@Downstream")) %>% mutate_all(funs(.*1000))
  
  #########
  #
  # Plots
  #
  ########
  XAll<-NULL
  
  error<-0
  
  for(i in 1:nrow(DOsites))
  {
    if(DOsites[i,]$Obs=="A4260704"||DOsites[i,]$Obs=="A4261107") next
    
    weight<-switch(DOsites[i,]$Obs,"A4260703"=0.3,"A4261160"=0.2,"A4261168"=0.5,"A4261224"=1)
    print(weight)
    
    
      obs<-na.trim(DOdat[,as.character(DOsites[i,]$Obs)])
      
      simcol<-grep(gsub(" ",".",DOsites[i,]$Sim),colnames(DOSim))
      sim<-zoo(DOSim[,simcol],Sim$Date)
      sim<-window(sim,start=start(obs),end=end(obs))
      dat<-cbind(obs=obs,sim=sim)
      dat<-na.trim(dat)
      colnames(dat)<-c("Observed","Simulated")
      
      error<-error+(1-NSE(dat$Simulated,dat$Observed))*weight
      
      X<-fortify(dat,melt=TRUE)
      X$Site<-DOsites[i,]$Obs
      XAll<-rbind(XAll,X)

  }
  
  if(index)
  {
  
  p<-ggplot(XAll)+geom_line(aes(Index,Value,colour=Series))+facet_grid(Site~.)+theme(legend.position = "top")+
    ylab("Dissolved Oxygen (mg/L)")+xlab("Date")
  ggsave(paste0("ChowPlots/DOResults",index,".png"),p,width=15,height=22,units="cm")
  
  simall<-NULL
  obsall<-NULL
  
  for(i in 1:nrow(DOCsites))
  {
    simcol<-grep(gsub(" ",".",DOCsites[i,]$Sim),colnames(DOCSim))
    sim<-fortify(zoo(DOCSim[,simcol],Sim$Date),melt=TRUE)
    
    simcol<-grep(gsub(" ",".",DOCsites[i,]$Sim),colnames(RefractSim))
    refract<-fortify(zoo(RefractSim[,simcol],Sim$Date),melt=TRUE)
    
    sim$Value<-sim$Value+refract$Value
    sim$Series<-DOCsites[i,]$Obs
    
    obs<-DOCdat[DOCdat$Series==DOCsites[i,]$Obs,]
    obs<-obs[obs$Index>=min(sim$Index),]
    sim<-sim[sim$Index>=min(obs$Index),]
    
    obsall<-rbind(obsall,obs)
    simall<-rbind(simall,sim)
    
  }

      obsall$Series<-gsub("River Murray","",obsall$Series)
  simall$Series<-gsub("River Murray","",simall$Series)
  
  p<-ggplot()+geom_line(data=simall,aes(Index,Value),col="blue")+geom_point(data=obsall,aes(Index,Value),size=0.5)+
    facet_grid(Series~.)+
    ylab ("Dissolved Organic Carbon (mg/L)")+xlab("Date")
  ggsave(paste0("ChowPlots/DOCResults_",index,".png"),width=15,height=22,units="cm")
  }
  
  return(error)
}


RunModel<-function(x)
{
  WriteInputSet(x)
  system("RunModel.bat")
  error<-PlotResults("Chowilla_4.7_NoWeirs - Scenario 1.res.csv",FALSE)
  return(error)
}

#values are DOC consumption rate, reaeration rate, redgum accumulation rate (kg/ha)
  init<-c(0.03,10)
  lower<-c(0.01,1)
  upper<-c(0.08,20)

  methods<-c("bobyqa","Rcgmin","spg") #try a few different optimisation algorithms
  outs<-list()
  for(method in methods)
  {
    X<-optimx(init,RunModel,lower=lower,upper=upper,method=method)
    save(X,file=paste0(method,".Rdata"))
    outs[[method]]<-X
  }
