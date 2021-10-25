dir.create("D:/RStudio/ABM code/Biodiversity/")
setwd("D:/RStudio/ABM code/Biodiversity/")
directory=getwd()
dir.create("D:/RStudio/ABM code/Biodiversity/output")
out.dir=paste(directory,"/output",sep='')
dir.create("D:/RStudio/ABM code/Biodiversity/source")
source(paste(directory, "/source/Comparetemp.R", sep=''))

#parameters
landscape=200 #number of patches=landscape^2
elev.range=c(-50,2000) #range of elevtion within the landscape
init.lat=20 #latitude of the southernmost patches within the landscape
incre.lat=0.05 #increment of latitude between two adjecent patches on latitudinal direction
init.temp=30 #temperature at sea level on southernmost latitude of the landscape
incre.temp=0.1 #increment of temperature per year
num.spe=100 #number of species simulated
num.years=100 #simulation lasting time of this model
replica=10 #number of reolicates to run the model

for(rep in 1:replica){
pdf(paste(out.dir, "/SpeciesMove_", rep,".pdf", sep=""))

#initialize the landscape with latitude
land.lat=matrix(nrow=landscape,ncol=landscape)
land.lat[1,]=init.lat
for(r in 2:landscape){
  land.lat[r,]=incre.lat+land.lat[r-1,]
}

#initialize the landscape with elevation
land.elev=matrix(nrow=landscape,ncol=landscape)

#determine the initial temperature of patches
#Temperature falls by 0.8¡æ for 1¡ã of latitudinal increase
#Temperature falls by 0.65¡æ for 100 meters of altitudinal increase
land.temp=matrix(nrow=landscape,ncol=landscape)
land.temp=init.temp-(land.lat-init.lat)*0.8-land.elev*0.65/100
image(land.temp)

#determine the initial location of species
x=sample(1:landscape,num.spe)
y=sample(1:landscape,num.spe)
loc.spe=t(rbind(x,y))

#Assign different temperature tolerance variations to species
temp.var=abs(round(rnorm(num.spe,1,0.25),2))
upper.spe=matrix(nrow=1,ncol=num.spe)
for(i in 1:num.spe) {
  #Initial temperature of patch is assigned to species 
  #located in this patch as its median of temperature tolerance range
  #Upper limit of species' temperature tolerance range is calculated
  #by addition of median and variation
  upper.spe[1,i]=land.temp[loc.spe[i,1],loc.spe[i,2]]+sample(temp.var,1)
}

#Species movement driven by temperature rising
pathways=NULL
for(s in 1:num.spe){
  mvmt=loc.spe[s,,drop=F]
  loc.x=mvmt[1]
  loc.y=mvmt[2]
  move.temp=land.temp
  for(n in 1:num.years){
    move.temp=move.temp+incre.temp
    lwpt=which(move.temp[(loc.x-1):(loc.x+1),(loc.y-1):(loc.y+1)]==min(move.temp[(loc.x-1):(loc.x+1),(loc.y-1):(loc.y+1)]), arr.ind = TRUE)
    if(nrow(lwpt)>1){
      p=sample(nrow(lwpt),1)
      a=lwpt[p,1]
      b=lwpt[p,2]
      nxtpt=Comparetemp(a,b,loc.x,loc.y)
    }
    else{
      a=lwpt[1]
      b=lwpt[2]
      nxtpt=Comparetemp(a,b,loc.x,loc.y)
    }
    mvmt=c(mvmt,nxtpt)
  }
  lines(mvmt[seq(1,length(mvmt), 2)], mvmt[seq(2,length(mvmt), 2)], lwd=2)
  pathways=rbind(pathways,mvmt)
}
rownames(pathways) = seq(1,num.spe,1)
write.table(pathways, paste(out.dir, "/SpeciesMove", ".csv", sep=""), row.names=F, col.names=F, sep=",")
dev.off()
}
