prepare.dates<- function(DT,Obs,dyn.type,unit="day"){
  switch (unit,
          "day" = unit<- 1,
          "week" = unit<- 7,
          "month" = unit<- 30.5,
          stop("Incorrect unit. Accepted units are: day, week, and month")
  )
  switch(dyn.type,
         "sta"={
           Obs.Season<-Obs[,.(dinf=min(date),dsup=max(date)),by=.(Season)];
           DT$dinf<- Obs.Season[match(DT$Season,Obs.Season$Season),]$dinf;DT$dsup<- Obs.Season[match(DT$Season,Obs.Season$Season),]$dsup;
           return(DT)
         },
         "dyn"={
           duration<- 4;extension<-0;
         },
         "RTWS"={
           duration<- 4;extension<-2;}
         ,
         stop("Incorrect static/dynamic type. Accepted types are: static, dynamic, and RTWS"))
  add.prepatent(DT,mean.lag = 0,duration = duration,extension = extension,unit=unit)
}
