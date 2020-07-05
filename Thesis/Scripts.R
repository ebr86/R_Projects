
######################
##R scripts
######################



###Solving Horton equation 
# Solving Horton formula for Fc
horton <- function(fc,f0,k,t){
  f = fc + (f0 -fc)*exp(-k*t)
  return(f)
}

#Fc equation
fcSolver <- function(f,k,tt){
  fc = f / (1 + 1.5*(1-exp(-k*tt))/(k*tt))
  return(fc)
}

#based on given reference Fo=2.5*Fc
ff0 <- function(fc){
  return(2.5*fc)
}

tmax = 7200  # sec, simulation time
t = seq(0,tmax,1) 
i0 = 0.001 # mm/sec infiltration 

k = 1.0/3600 # converting 1/T to 1/s
fc = fcSolver(i0,k,tmax)
f0 = ff0(fc)

png(filename="Horton Equation.png",width = 270, height = 160, units = "mm" ,res = 800)
plot(t,horton(fc,f0,k,t)/i0,t="l",ylim=c(0.5,2.5), ylab = "Horton Equation / i0", xlab = "Time (s)")
lines(t,rep(i0,length(t))/i0,col="green") ## to show Horton infiltration is below constant infiltration of 1
legend( "topright", legend=c( "rain intensity","k = 2","k = 1","k = 0.5", "i0 = 0.001")  ,col=c("blue","green","black","red","purple") , lty=c(2,1,1,1,1) , bty="n",  ncol = 1,lwd = 2 ) 
lines(t,rep(i0,length(t))/i0,col="purple") ## constant rate

## Rain amount
rainDur = 1800
rt = seq(1:rainDur)
rain = 7.5/3600
lines(rt,rep(rain,length(rt))/i0,col="blue", lty= "dashed")

k = 1/3600  # example for k=1 
fc = fcSolver(i0,k,tmax)
f0 = ff0(fc)
lines(t,horton(fc,f0,k,t)/i0,col="black")

dev.off() 

#### Reading Background cases

allslopes=c(0.0000, 0.0010, 0.0020, 0.0050, 0.0080, 0.0100, 0.0120, 0.0140,
            0.0160, 0.0180, 0.0190, 0.0200, 0.0300, 0.0400, 0.0500, 0.0600, 0.0700, 0.0800, 0.0900, 0.1000)    # 20  slopes
allwavelenghts=c(0.1500, 0.4143, 0.6786, 0.9429, 1.2071,
                 1.4714, 1.7357, 2.0000)      # 8   wave-lenghts
allamplitudes=c(0.0100, 0.0200, 0.0300, 0.0400,
                0.0500, 0.0600, 0.0700, 0.0800, 0.0900) # 9   amplitudes


allslopesString = format(allslopes,nsmall=4) # text
allwavelenghtsString = format(allwavelenghts,nsmall=4)
allamplitudesString = format(allamplitudes,nsmall=4)


Ns=length(allslopes) # number of slopes
Nwl=length(allwavelenghts) #8
Na=length(allamplitudes) #9


simTime = c(0,7200) # simulation time
myColors = c("orange","red","blue","green","yellow","violet","pink","grey","brown")
BCmat=array(data.frame(),dim = c(Ns))



for(i in 1:Ns){ 
  slopeName = paste("Slope",allslopesString[i],sep="") 
  allfiles= file.path("D:","BackgroundCases", slopeName, "MT_out.mass")
  BCmat[[i]]= read.table (allfiles, header=T)
}

## Ranges should be adjusted
## plot Discharge
for (i in 1:Ns){
  Fname=paste("QQQ_S",allslopesString[i],sep="_") 
  plotpath = file.path("D:","allplots","BCplots", paste(Fname,".JPG",sep = ""))
  png(file=plotpath, width=4, height=4, units="in", res=600)
  df = BCmat[[i]]
  plot(df$time_.s.,abs(df$bc1_.m3.s.)*1000,path=plotpath, type = "l", col= "blue", xlim = simTime, ylim = qrange ,xlab= "Time [s]", ylab= "Runoff [L/s]", 
       main=paste("Slope",allslopesString[i])) 
}
## Accum Infiltration
for (i in 1:Ns){
  Mname=paste("AccInf_S-Tot",allslopesString[i],sep="_")
  plotpath = file.path("D:","allplots","BCplots", paste(Mname,".JPG",sep = ""))
  png(file=plotpath, width=4, height=4, units="in", res=600)
  Frange = c(0,max(BCmat[[1]]$accumInf_.m3.))
  df = BCmat[[i]]
  plot(df$time_.s.,df$accumInf_.m3.+df$massDomain_.m3.,path=plotpath, type = "l", col= "blue", xlim = simTime, ylim = Frange ,xlab= "Time [s]", ylab= "Total Accumulation Infiltration [m3/s]", 
       main=paste("Slope",allslopesString[i])) 
}

## Infiltration rate

for (i in 1:Ns){
  Uname=paste("InfilRate_S",allslopesString[i],sep="_")
  plotpath = file.path("D:","allplots","BCplots", paste(Uname,".JPG",sep = ""))
  png(file=plotpath, width=4, height=4, units="in", res=600)
  Rrange = c(0,max(BCmat[[1]]$inf_.m3.s.))
  df = BCmat[[i]]
  plot(df$time_.s.,df$inf_.m3.s.,path=plotpath, type = "l", col= "blue", xlim = simTime, ylim = Rrange ,xlab= "Time [s]", ylab= "Infiltration Rate [m3/s]",
       main=paste("Slope",allslopesString[i]))
}

###Reading microtopography cases
## change the filepath for reading constant infiltration setup

mat=array(data.frame(),dim = c(Ns,Nwl,Na))
for(i in 1:Ns){
  for (j in 1:Nwl){
    for (k in 1:Na){
      slopeName = paste("Slope",allslopesString[i],sep="")  
      wlName = paste("WL",allwavelenghtsString[j],sep="")
      aName = paste("A",allamplitudesString[k],sep="")
      allfiles=file.path("D:","MassFiles","R7.5new",slopeName,wlName,aName,"MT_out.mass")
      mat[[i,j,k]]= read.table (allfiles, header=T)
    }
  }
}


#### Plotting Phase
#ranges should be adjusted

simTime = c(0,7200) 
myColors = c("orange","red","blue","green","yellow","violet","pink","grey","brown")

## Discharge or Runoff 
qrange = c(0,0.05) 

for (i in 1:Ns){
  for (j in 1:Nwl){
    
    Fname=paste("QQQ_S",allslopesString[i],"W",allwavelenghtsString[j],sep="_") 
    plotpath = file.path("D:","allplots","Runoff-Discharge", paste(Fname,".png",sep = ""))
    png(file=plotpath, width=4, height=4, units="in", res=600)
    plot(NA, xlim = simTime, ylim = qrange ,xlab= "Time [s]", ylab= "Runoff [L/s]", 
         main=paste("Slope",allslopesString[i],"-","WL",allwavelenghtsString[j])) 
    for (k in 1:Na){k
      df = mat[[i,j,k]]
      lines(df$time_.s.,abs(as.numeric(df$bc1_.m3.s.))*1000,col=myColors[k])
    }
    dev.off()
  }
}

for (i in 1:Ns){
  for (k in 1:Na){
    Gname=paste("QQQ_S",allslopesString[i],"A",allamplitudesString[k],sep="_") 
    plotpath = file.path("D:","allplots","Runoff-Discharge", paste(Gname,".png",sep = ""))
    png(file=plotpath, width=4, height=4, units="in", res=600)
    
    plot(NA, xlim = simTime, ylim = qrange , xlab="Time [s]", ylab="Runoff [L/s]", 
         main= paste("Slope",allslopesString[i],"-","AMP",allamplitudesString[k]))   
    for (j in 1:Nwl){
      df = mat[[i,j,k]]
      lines(df$time_.s.,abs(as.numeric(df$bc1_.m3.s.))*1000,col=myColors[j]) 
    }
    dev.off()
  }
}

for (j in 1:Nwl){
  for(k in 1:Na){
    Hname=paste("QQQ_W",allwavelenghtsString[j],"A",allamplitudesString[k],sep="_") #
    plotpath = file.path("D:","allplots","Runoff-Discharge",paste(Hname,".png",sep = ""))
    png(file=plotpath, width=4, height=4, units="in", res=600)
    plot(NA, xlim = simTime, ylim = qrange , xlab="Time [s]", ylab="Runoff [L/s]",
         main= paste("WL",allwavelenghtsString[j],"-","AMP",allamplitudesString[k]))    
    for (i in 1:Ns){
      df = mat[[i,j,k]]
      lines(df$time_.s.,abs(as.numeric(df$bc1_.m3.s.))*1000,col=myColors[i])      
    }
    dev.off()
  }
}

## Accumulation Infiltration
for (i in 1:Ns){
  for (j in 1:Nwl){
    
    Mname=paste("AccInf_S",allslopesString[i],"W",allwavelenghtsString[j],sep="_")
    plotpath = file.path("D:","allplots","Accum-Infiltration", paste(Mname,".png",sep = ""))
    png(file=plotpath, width=4, height=4, units="in", res=600)
    
    Frange = c(0,0.1) 
    plot(NA, xlim = simTime, ylim = Frange , xlab="Time [s]", ylab="Accumulation Infiltration [m3/s]",
         main=paste("Slope",allslopesString[i],"-","WL",allwavelenghtsString[j])) 
    for (k in 1:Na){k
      df = mat[[i,j,k]]
      lines(df$time_.s.,df$accumInf_.m3.+df$massDomain_.m3. , col=myColors[k]) 
    }
    dev.off()
  }
}
for (i in 1:Ns){
  for (k in 1:Na){
    
    Nname=paste("AccInf_S",allslopesString[i],"A",allamplitudesString[k],sep="_") 
    plotpath = file.path("D:","allplots","Accum-Infiltration", paste(Nname,".png",sep = ""))
    png(file=plotpath, width=4, height=4, units="in", res=600)
    plot(NA, xlim = simTime, ylim = Frange , xlab="Time [s]", ylab="Accumulation Infiltration [m3/s]", 
         main= paste("Slope",allslopesString[i],"-","AMP",allamplitudesString[k]))    
    for (j in 1:Nwl){
      df = mat[[i,j,k]]
      lines(df$time_.s.,df$accumInf_.m3.+df$massDomain_.m3., col=myColors[j])      
    }
    dev.off()
  }
}
for (j in 1:Nwl){
  for(k in 1:Na){
    Oname=paste("AccInf_W",allwavelenghtsString[j],"A",allamplitudesString[k],sep="_")     
    plotpath = file.path("D:","allplots","Accum-Infiltration", paste(Oname,".png",sep = ""))
    png(file=plotpath, width=4, height=4, units="in", res=600)
    
    plot(NA, xlim = simTime, ylim = Frange , xlab="Time [s]", ylab="Acc Infiltration [m3/s]", 
         main =paste("WL",allwavelenghtsString[j],"-","AMP",allamplitudesString[k]))    
    for (i in 1:Ns){
      df = mat[[i,j,k]]
      lines(df$time_.s.,df$accumInf_.m3.+ df$massDomain_.m3., col=myColors[i])      
    }
    dev.off()
  }
}
## Infiltration Rate
for (i in 1:Ns){
  for (j in 1:Nwl){
    Uname=paste("InfilRate_S",allslopesString[i],"W",allwavelenghtsString[j],sep="_")
    plotpath = file.path("D:", "allplots","Infil-Rate", paste(Uname,".png",sep = ""))
    png(file=plotpath, width=4, height=4, units="in", res=600)
    
    Rrange = c(0,0.01)   #max(mat[[1,1,1]]$inf_.m3.s.))
    plot(NA, xlim = simTime, ylim = Rrange , xlab="Time [s]", ylab="Infiltration Rate [L/s]", 
         main=paste("Slope",allslopesString[i],"-","WL",allwavelenghtsString[j])) 
    
    for (k in 1:Na){k
      df = mat[[i,j,k]]
      lines(df$time_.s.,(df$inf_.m3.s.)*1000, col=myColors[k]) 
    }
    dev.off()
  }
}


for (i in 1:Ns){
  for (k in 1:Na){
    
    Vname=paste("InfilRate_S",allslopesString[i],"A",allamplitudesString[k],sep="_") 
    plotpath = file.path("D:","allplots","Infil-Rate", paste(Vname,".png",sep = ""))
    png(file=plotpath, width=4, height=4, units="in", res=600)
    
    plot(NA, xlim = simTime, ylim = Rrange , xlab="Time [s]", ylab="Infiltration Rate [L/s]", 
         main=paste("Slope",allslopesString[i],"-","AMP",allamplitudesString[k])) 
    for (j in 1:Nwl){
      df = mat[[i,j,k]]
      lines(df$time_.s.,(df$inf_.m3.s.)*1000,col=myColors[j]) 
    }
    dev.off()
  }
}
for (j in 1:Nwl){
  for(k in 1:Na){
    Wname=paste("InfilRate_W",allwavelenghtsString[j],"A",allamplitudesString[k],sep="_") 
    plotpath = file.path("D:","allplots","Infil-Rate", paste(Wname,".png",sep = ""))
    png(file=plotpath, width=4, height=4, units="in", res=600)
    plot(NA, xlim = simTime, ylim = Rrange , xlab="Time", ylab="Infiltration Rate [L/s]", 
         main=paste("WL",allwavelenghtsString[j],"-","AMP",allamplitudesString[k])) 
    for (i in 1:Ns){
      df = mat[[i,j,k]]
      lines(df$time_.s.,(df$inf_.m3.s.)*1000, col=myColors[i])  
    }
    dev.off()
  }
}


### Runoff Comparison
## sample for specific S, WL and Amp: Slope0.0000 + WL1.4714 + A0.0100, A0.0300, A0.090

bc0= file.path("D:","BackgroundCases","Slope0.0000", "MT_out.mass")
dbc0 = read.table(bc0,header=T)
md00= file.path("D:", "MassFiles","R7.5new","Slope0.0000","WL1.4714","A0.0100","MT_out.mass")
df00 = read.table(md00,header=T)
md01= file.path("D:","MassFiles","R7.5new","Slope0.0000","WL1.4714","A0.0300","MT_out.mass")
df01= read.table(md01,header=T)
md02= file.path("D:","MassFiles","R7.5new","Slope0.0000","WL1.4714","A0.0900","MT_out.mass")
df02= read.table(md02,header=T)

# plotting runoff
ylab11 = parse(text = "Runoff ~~ group('[', m^3/s, ']')")  # Runoff [m3/s]
png(filename="runoff_s0-0_wl1-4714.png",width = 270, height = 160, units = "mm" ,res = 800)
plot(dbc0$time_.s.,  abs(dbc0$bc1_.m3.s.), t= "l",col= "black",lty="dashed", lwd= "2" ,xlim = simTime,xlab= "Time [s]", ylab= ylab11)
lines(df00$time_.s., abs(df00$bc1_.m3.s.), t= "l",col= "green",  xlim = simTime, xlab= "Time [s]", ylab= ylab11)
lines(df01$time_.s., abs(df01$bc1_.m3.s.), t= "l",col= "blue",   xlim = simTime, xlab= "Time [s]", ylab= ylab11)
lines(df02$time_.s., abs(df02$bc1_.m3.s.), t= "l",col= "purple", xlim = simTime, xlab= "Time [s]", ylab= ylab11)
legend(title= "Amplitude", "topright", legend=c( "BC","0.0100","0.0300","0.0900")  ,col=c("black","green","blue","purple") , lty=c(2,1,1,1) , bty="n",  ncol = 1,lwd = 2 ) 
dev.off()


# IR-plot
# reading Background & Microtopography cases as before

TinfB=array(999,dim = c(Ns)) # empty array

for(i in 1:Ns){ # total infiltration =  accumulation water + ponded water
  BCmat[[i]]$totinf.bc = BCmat[[i]]$accumInf_.m3. + BCmat[[i]]$massDomain_.m3.
  TinfB[[i]] = BCmat[[i]]$totinf.bc[nrow(BCmat[[i]])]  # last value
}

TinfM=array(999,dim = c(Ns,Nwl,Na)) 
for(i in 1:Ns){   
  for (j in 1:Nwl){
    for (k in 1:Na){
      Mmat[[i,j,k]]$totinf.mc = Mmat[[i,j,k]]$accumInf_.m3. +  Mmat[[i,j,k]]$massDomain_.m3.
      TinfM[[i,j,k]] = Mmat[[i,j,k]]$totinf.mc[nrow(Mmat[[i,j,k]])] 
    }
  }
}

mDb=array(-1,dim = c(Ns,Nwl,Na))
for(i in 1:Ns){   
  for (j in 1:Nwl){
    for (k in 1:Na){
      mDb[[i,j,k]]= TinfM[[i,j,k]] / TinfB[[i]]
    }
  }
}



library(reshape2)
mDb.melt= melt(mDb, varnames = c("S","WL","Amp") , id.vars = c("S","WL","Amp"), measure.vars = "value") # ,  )
View(mDb.melt)
mDb.melt.2 = mDb.melt  

for (i in 1:Ns) {  # changing numbers to cases values
  for (j in 1:Nwl) {
    for (k in 1:Na) {
      mDb.melt.2$S[which(mDb.melt.2$S == i)]= allslopes[i]
      mDb.melt.2$WL[which(mDb.melt.2$WL == j)]= allwavelengths[j]
      mDb.melt.2$Amp[which(mDb.melt.2$Amp == k)]= allamplitudes[k]
    }
  }
}


mDb.melt.final = mDb.melt.2
colnames(mDb.melt.final)[4] = "IR"
View(mDb.melt.final)


# sample plot
library(ggplot2)
library(RColorBrewer)
library("viridis")
my_colors= c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")

## IR for all cases
sensIR= ggplot(mDb.melt.final, aes(x= as.factor(Amp), y= as.factor(WL),  fill=as.numeric(IR)))+ 
  geom_tile(color='black')+
  theme_bw() +
  facet_wrap(. ~ as.factor(S))+
  scale_fill_gradientn(colours= my_colors) + 
  scale_y_discrete (name="Wavelength (m)", labels=c("0.15", "0.41", "0.67", "0.94", "1.20", "1.47", "1.73", "2.00")) + 
  scale_x_discrete(name="Amplitude (m)",  labels=c("0.01" , "0.02",  "0.03", "0.04", "0.05", "0.06", "0.07", "0.08", "0.09")) +  
  labs(label=11, title = "Infiltration Ratio", subtitle= "for all slopes" , fill="IR") 
sensIR

####MR vs IR scatterplot

#reading same as 2.5
#change for specific slope

mDb.melt.final$MRTs = round(mDb.melt.final$WL/mDb.melt.final$Amp, digits=2) # calculating & adding MR
# sample plot
scatT= ggplot(mDb.melt.final , aes(x= MRTs, y= as.numeric(IR) , col= IR)) +
  geom_point(size=2) +
  theme_bw() +
  scale_x_log10() +  # with logarithmic scale
  scale_color_gradientn(colours = my_colors) + 
  facet_wrap(. ~ as.factor(S)) +
  labs(title=" Scatterplot", subtitle = "All Slopes" , x="Microtopography Ratio", y="Infiltration Ratio" )
scatT

#####I1 , I 

# calculating IR for constant infiltration rate same as section 2.5
#  I ^  = Horton IR / constant IR
# I1 = total infiltration of Horton microtopography cases / total infiltration of constant microtopography cases
# plotting similar to 2.5


