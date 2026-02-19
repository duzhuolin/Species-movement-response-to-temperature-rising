setwd("D:/RStudio/ABM code/Biodiversity/")
directory=getwd()
out.dir=paste(directory,"/output",sep='')
source(paste(directory, "/source/Comparetemp.R", sep=''))
library(raster)

#parameters
landscape=500 #number of patches=landscape^2,virtually representing a 5°x5° area in real world
coor.nwp=c(100,30) #coordinate(longitude,latitude) of the most northwestern point in the landscape
init.temp=25 #temperature at sea level on the latitude of landscape boundary near the equator
mean.temp=0.05 #mean of increment of temperature per year
var.temp=0.05 #standard deviation of temperature increment
num.spe=500 #number of species simulated
num.years=300 #simulation lasting time of this model
replica=5 #number of replicates to run the model
move.lkl=c(0.85) #likelihood of species moving to coolest neighbor patch
var.stp=c(1:3) #options for step sizes that one species can make each time

for(m in 1:length(move.lkl)){
  for(rep in 1:replica){
   print(rep)
   pdf(paste(out.dir, "/SpeciesMoveInMountain_", rep,".pdf", sep=""))

   #initialize the landscape with latitude
   land.lat=matrix(nrow=landscape,ncol=landscape)
   land.lat[landscape,]=coor.nwp[2]
   #latitudinal span from northmost to southmost of the landscape is 5°
   incre.lat=round(5/(landscape-1),2) #increment of latitude between two adjacent patches on latitudinal direction
   for(r in (landscape-1):1){
     land.lat[r,]=land.lat[r+1,]-incre.lat
   }

   #initialize the landscape with elevation
   dem=getData('SRTM', lon=coor.nwp[1], lat=coor.nwp[2])
   df=as.matrix(dem) #df is a 6000x6000 matrix
   rowarray=rev(round(seq(1,nrow(df),nrow(df)/landscape),0)) #upscale from the data matrix provided by SRTM
   land.elev=matrix(nrow=landscape,ncol=landscape)
   for(row in 1:landscape){
      land.elev[row,]=df[rowarray[row],round(seq(1,ncol(df),ncol(df)/landscape),0)] #upscale from the data matrix provided by SRTM
   }

   #determine the initial temperature of patches
   #Temperature falls by 0.8° for 1° of latitudinal increase
   #Temperature falls by 0.65° for 100 meters of altitudinal increase
   land.temp=matrix(nrow=landscape,ncol=landscape)
   land.temp=t(round(init.temp-(land.lat-land.lat[1,1])*0.8-land.elev*0.65/100,2))
   image(land.temp)

   #determine the initial location of species
   x=sample(1:landscape,num.spe,replace=T)
   y=sample(1:landscape,num.spe,replace=T)
   loc.spe=t(rbind(x,y))
   points(loc.spe[,1]/landscape, loc.spe[,2]/landscape, pch=19, cex=0.3)

   #Assign different temperature tolerance variations to species
   temp.var=sample(abs(round(rnorm(num.spe,2,0.5),2)),num.spe)
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
             stpsz=sample(var.stp,1) #choose maximal step size for current step
             if(loc.x>(landscape-stpsz) | loc.y>(landscape-stpsz) | loc.x<=stpsz | loc.y<=stpsz){
                nxtpt=c(loc.x,loc.y)
                mvmt=c(mvmt,nxtpt)
                next
             }
             else{
                movecool=sample(c(0,1), 1, prob=c((1-move.lkl[m]), move.lkl[m]))
                if(movecool==0){
                  loc.x = sample(c(-stpsz:stpsz), 1) + loc.x
                  loc.y = sample(c(-stpsz:stpsz), 1) + loc.y
                  nxtpt=c(loc.x,loc.y)
                }
                if(movecool==1){
                  lwpt=which(move.temp[(loc.x-stpsz):(loc.x+stpsz),(loc.y-stpsz):(loc.y+stpsz)]==min(move.temp[(loc.x-stpsz):(loc.x+stpsz),(loc.y-stpsz):(loc.y+stpsz)]), arr.ind = TRUE)
                  if(nrow(lwpt)>1){ 
                    p=sample(nrow(lwpt),1)
                    a=lwpt[p,1]
                    b=lwpt[p,2]
                    nxtpt=Comparetemp(a,b,loc.x,loc.y,stpsz)
                    loc.x=nxtpt[1]
                    loc.y=nxtpt[2]
                  }
                  else{
                    a=lwpt[1]
                    b=lwpt[2]
                    nxtpt=Comparetemp(a,b,loc.x,loc.y,stpsz)
                    loc.x=nxtpt[1]
                    loc.y=nxtpt[2]
                  }
                }
             }
          }
       }
       mvmt=c(mvmt,nxtpt)
     }
     lines(mvmt[seq(1,length(mvmt), 2)]/landscape, mvmt[seq(2,length(mvmt), 2)]/landscape, lwd=0.6)
     pathways=rbind(pathways,mvmt)

     #jdjd
   }
   rownames(pathways) = seq(1,num.spe,1)
   write.table(pathways, paste(out.dir, "/SpeciesMoveInMountain_", rep, ".csv", sep=""), row.names=F, col.names=F, sep=",")
   dev.off()
  }
}

