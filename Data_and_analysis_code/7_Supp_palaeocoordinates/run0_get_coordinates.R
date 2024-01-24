# 2023 04 02 I.Zliobaite
# updated 2023 08 10
# get palaeocoordinates for NOW localities

data_all <- read.csv('data_raw/now_pub_20230818.csv', header = TRUE, sep = "\t", na.string = '\\N')

#palaeocord <- c('PALEOMAP','SETON2012','MULLER2016','MULLER2022','MATTHEWS2016_mantle_ref')
#palaeocord <- c('MULLER2022','MATTHEWS2016_mantle_ref')
#palaeocord <- c('MULLER2019','MATTHEWS2016_pmag_ref','MERDITH2021')
palaeocord <- c('MULLER2016')

#"SETON2012" (Seton et al., 2012) for coastlines and topological plate polygons (0-200 Ma).
#"MULLER2016" (Muller et al., 2016) for coastlines and topological plate polygons (0-230Ma).
#"GOLONKA" (Wright et al. 2013) for coastlines only (0-550 Ma).
#"PALEOMAP" (Scotese, 2016) for coastlines only (0-1100 Ma).
#"MATTHEWS2016_mantle_ref" (Matthews et al., 2016) for coastlines and topological plate polygons (0-410 Ma).
# "MATTHEWS2016_pmag_ref" (Matthews et al., 2016) for coastlines and topological plate polygons (0-410 Ma).
# "MULLER2019" (Müller et al., 2019) for coastlines and static plate polygons. (0-250 Ma).
# "MERDITH2021" (Merdith et al., 2021, default) for coastlines and static plate polygons (0-1000 Ma).
# "MULLER2022" (Müller et al., 2022) for coastlines and static plate polygons (0-1000 Ma).

# assigning localities to discreat ages of millions of years (as per rgplates requirement)
midage <- (data_all[,'MAX_AGE'] + data_all[,'MIN_AGE'])/2
roundmidage <- round(midage)
data_all <- cbind(data_all,midage,roundmidage)

un_age <- unique(roundmidage)
un_age <- un_age[order(un_age)]

library(rgplates)
#xy <-cbind(long=c(95,142), lat=c(54, -33))
#reconstruct(xy, 140, model=NULL)


for (sk3 in 1:length(palaeocord)){
  
  data_cords <- c()
  
  for (sk in 1:length(un_age)){
    age_now <- un_age[sk]
    ind <- which(data_all[,'roundmidage']==age_now)
    data_now <- data_all[ind,]
    un_loc <- unique(data_now[,'LIDNUM'])
    xy <- c()
    for (sk2 in 1:length(un_loc)){
      loc_now <- un_loc[sk2]
      ind2 <- which(data_now[,'LIDNUM']==loc_now)
      xy <- rbind(xy,c(data_now[ind2[1],'LONG'],data_now[ind2[1],'LAT']))
    }
    colnames(xy) <- c('long','lat')
    #rr <- reconstruct(xy, age_now,model = 'SETON2012')
    #rr <- reconstruct(xy, age_now,model = 'PALEOMAP')
    rr <- reconstruct(xy, age_now,model = palaeocord[sk3])
    agevec <- rep(age_now,length(un_loc),1)
    data_cords <- rbind(data_cords,cbind(un_loc,xy,rr,agevec))
  }
  
  colnames(data_cords) <- c('LIDNUM','long','lat','paleolong','paleolat','roundage')
  
  #write.table(data_cords, file = "data_cords.csv",col.names = TRUE,row.names = FALSE, sep = '\t')   
  
  data_cords_all <- c()
  
  for (sk in 1:dim(data_all)[1]){
    loc_now <- data_all[sk,'LIDNUM']
    ind <- which(data_cords[,'LIDNUM'] == loc_now)
    data_cords_all <- rbind(data_cords_all,data_cords[ind,])
  }
  
  #diff <- sqrt((long-paleolong)^2 + (lat-palaeolat)^2)
  
  file_name <- paste('data_working/data_cords_',palaeocord[sk3],'.csv',sep='')
  
  write.table(data_cords_all, file = file_name,col.names = TRUE,row.names = FALSE, sep = '\t')  
  
  
  
  
}




