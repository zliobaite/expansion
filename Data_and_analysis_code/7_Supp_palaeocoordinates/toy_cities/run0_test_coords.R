# 2023 0810 I.Zliobaite

# simulation where is Helsinki

data_cities <- read.csv('coords_cities.csv', header = TRUE, sep = '\t')

meths_all <- c('SETON2012','MERDITH2021','MULLER2016','MULLER2019','MULLER2022','MATTHEWS2016_mantle_ref','MATTHEWS2016_pmag_ref','PALEOMAP')
meths_names <- c('SET12','MER21','MUL16','MUL19','MUL22','MAT16m','MAT16p','PAL')

all_years <- c(0:50)

city_now <- 'Berlin'

ind_city <- which(data_cities[,1]==city_now)

all_longs <- c()
all_lats <-c()

library(rgplates)

for (sk2 in 1:length(meths_all)){
  
  meth_now <- meths_all[sk2]
  print(meth_now)
  
  xy <- c(data_cities[ind_city,'Long'],data_cities[ind_city,'Lat'])
  
  rr <- xy
  
  for (sk in all_years){
    
    rr <- rbind(rr, reconstruct(xy, sk, model = meth_now))
    
  }
  
  all_longs <- cbind(all_longs,rr[,1])
  all_lats <- cbind(all_lats,rr[,2])
}

all_years <- c('0',all_years)

all_longs <- cbind(all_longs,all_years)
all_lats <- cbind(all_lats,all_years)

colnames(all_longs) <- c(meths_all,'age')
colnames(all_lats) <- c(meths_all,'age')


file_long <- paste('all_longs_',city_now,'.csv',sep='')
file_lat <- paste('all_lats_',city_now,'.csv',sep='')

write.table(all_longs, file = file_long, col.names = TRUE,row.names = FALSE, sep = '\t')   
write.table(all_lats, file = file_lat, col.names = TRUE,row.names = FALSE, sep = '\t')   


for (sk in 1:length(all_years)){
  year_now <- all_years[sk]
  
  data_now <- c()
  
  for (sk2 in 1:dim(all_longs)[2]){
    data_now <- rbind(data_now,c(all_longs[sk,sk2],all_lats[sk,sk2]))  
  }
  
  plot_name <- paste('plots/',city_now,'/plot_',city_now,'_',year_now,'.pdf',sep = '')
  pdf(plot_name,height = 4, width = 3.5)
  plot(NA,NA,xlim = c(data_cities[ind_city,'Long']-8,data_cities[ind_city,'Long']+8),ylim = c(data_cities[ind_city,'Lat']-8,data_cities[ind_city,'Lat']+8), xlab = 'Long now',ylab = 'Lat now', main = paste(city_now,year_now,'Ma'))
  points(data_now,pch = 16)
  points(data_cities[ind_city,'Long'],data_cities[ind_city,'Lat'],pch = 16, col='red')
  text(data_now,meths_names,pos = 4, cex = 0.4)
  text(data_cities[ind_city,'Long'],data_cities[ind_city,'Lat'],city_now,pos = 4, cex = 0.4, col='red')
  dev.off()
}