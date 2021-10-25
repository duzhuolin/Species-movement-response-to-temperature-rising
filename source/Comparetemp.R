Comparetemp=function(a,b,loc.x,loc.y){
  if(a==1){loc.x = loc.x - 1 }
  if(a==2){loc.x = loc.x }
  if(a==3){loc.x = loc.x + 1 }
  if(b==1){loc.y = loc.y - 1 }
  if(b==2){loc.y = loc.y }
  if(b==3){loc.y = loc.y + 1 }
  point=c(loc.x,loc.y)
  return(point)
}