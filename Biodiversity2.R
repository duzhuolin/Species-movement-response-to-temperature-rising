setwd("D:/RStudio/ABM code/Biodiversity/")
directory=getwd()
out.dir=paste(directory,"/output",sep='')
source(paste(directory, "/source/Comparetemp.R", sep=''))

#parameters
landscape=200 #number of patches=landscape^2
elev.range=c(0,2000) #range of elevation within the landscape
init.lat=20 #latitude of the southernmost patches within the landscape
incre.lat=0.05 #increment of latitude between two adjacent patches on latitudinal direction
init.temp=30 #temperature at sea level on southernmost latitude of the landscape
mean.temp=0.5 #mean of increment of temperature per year
var.temp=0.1 #standard deviation of temperature increment
num.spe=100 #number of species simulated
num.years=100 #simulation lasting time of this model
replica=10 #number of replicates to run the model

for(rep in 1:replica){
 print(rep)
 pdf(paste(out.dir, "/SpeciesMove_", rep,".pdf", sep=""))

 #initialize the landscape with latitude
 land.lat=matrix(nrow=landscape,ncol=landscape)
 land.lat[1,]=init.lat
 for(r in 2:landscape){
   land.lat[r,]=incre.lat+land.lat[r-1,]
 }

 #initialize the landscape with elevation
 land.elev=matrix(nrow=landscape,ncol=landscape)
 xpeak = sample(10:(landscape-10), 1)
 ypeak = sample(10:(landscape-10), 1)
 land.elev[xpeak, ypeak] = elev.range[2]
 land.elev[xpeak, 1:(ypeak-1)] = round(seq(elev.range[1], elev.range[2], (elev.range[2]-elev.range[1])/(ypeak-2)) + rnorm((ypeak-1),0,2), 0)
 land.elev[xpeak, (ypeak+1):landscape] = round(rev(seq(elev.range[1], elev.range[2], (elev.range[2]-elev.range[1])/(landscape-ypeak-1)) + rnorm((landscape-ypeak),5,2)), 0)
 for(e in (xpeak-1):1){
   land.elev[e,] = land.elev[(e+1),] - round(rnorm(landscape, 5, 2), 0)
 }
 for(e in (xpeak+1):landscape){
   land.elev[e,] = land.elev[(e-1),] - round(rnorm(landscape, 5, 2), 0)
 }

 #determine the initial temperature of patches
 #Temperature falls by 0.8¡æ for 1¡ã of latitudinal increase
 #Temperature falls by 0.65¡æ for 100 meters of altitudinal increase
 land.temp=matrix(nrow=landscape,ncol=landscape)
 land.temp=round(init.temp-(land.lat-init.lat)*0.8-land.elev*0.65/100,2)
 image(land.temp)

 #determine the initial location of species
 x=sample(1:landscape,num.spe)
 y=sample(1:landscape,num.spe)
 loc.spe=t(rbind(x,y))
 points(loc.spe[,1]/200, loc.spe[,2]/200, pch=19, cex=0.1)

 #Assign different temperature tolerance variations to species
 temp.var=sample(abs(round(rnorm(num.spe,1,0.25),2)),num.spe)
 upper.spe=matrix(nrow=1,ncol=num.spe)
 for(i in 1:num.spe) {
   #Initial temperature of patch is assigned to species 
   #located in this patch as its median of temperature tolerance range
   #Upper limit of species' temperature tolerance range is calculated
   #by addition of median and variation
   upper.spe[1,i]=land.temp[loc.spe[i,1],loc.spe[i,2]]+temp.var[i]
 }
 upper.spe=round(upper.spe,2)

 #Species movement driven by temperature rising
 pathways=NULL
 incre.temp=round(rnorm(num.years,mean.temp,var.temp),2)
 for(s in 1:num.spe){
   mvmt=loc.spe[s,,drop=F]
   loc.x=mvmt[1]
   loc.y=mvmt[2]
   move.temp=land.temp
   for(n in 1:num.years){
     move.temp=move.temp+incre.temp[n]
     if(move.temp[loc.x,loc.y]<=upper.spe[s]){
       nxtpt=c(loc.x,loc.y)
     }
     else{
       if(loc.x==landscape | loc.y==landscape | loc.x==1 | loc.y==1){
         mvmt=c(mvmt, rep(c(loc.x,loc.y), (num.years-((length(mvmt)/2)-1))))
         break
       }
       else{
         lwpt=which(move.temp[(loc.x-1):(loc.x+1),(loc.y-1):(loc.y+1)]==min(move.temp[(loc.x-1):(loc.x+1),(loc.y-1):(loc.y+1)]), arr.ind = TRUE)
         if(nrow(lwpt)>1){ 
           p=sample(nrow(lwpt),1)
           a=lwpt[p,1]
           b=lwpt[p,2]
           nxtpt=Comparetemp(a,b,loc.x,loc.y)
           loc.x=nxtpt[1]
           loc.y=nxtpt[2]
         }
         else{
           a=lwpt[1]
           b=lwpt[2]
           nxtpt=Comparetemp(a,b,loc.x,loc.y)
           loc.x=nxtpt[1]
           loc.y=nxtpt[2]
         }
       }
     }
    mvmt=c(mvmt,nxtpt)
   }
   
   lines(mvmt[seq(1,length(mvmt), 2)]/200, mvmt[seq(2,length(mvmt), 2)]/200, lwd=0.1)
   pathways=rbind(pathways,mvmt)
 }
 rownames(pathways) = seq(1,num.spe,1)
 write.table(pathways, paste(out.dir, "/SpeciesMove_", rep, ".csv", sep=""), row.names=F, col.names=F, sep=",")
 dev.off()
}