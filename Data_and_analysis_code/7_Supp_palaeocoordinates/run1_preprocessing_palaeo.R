# 2023 03 16 I.Zliobaite

data_all <- read.csv('data_raw/now_pub_20230818.csv', header = TRUE, sep = "\t", na.string = '\\N')

palaeocord <- 'MULLER2016'

file_name <- paste('data_working/data_cords_',palaeocord[1],'.csv',sep='')

data_palaeocord <- read.csv(file_name, header = TRUE, sep = "\t", na.string = c('\\N','NA'))

data_all <- cbind(data_all,data_palaeocord)

ind <- which(!is.na(data_all[, "paleolat"]))
data_all <- data_all[ind,]

age_threshold <- 3 #localities with larger than this time intervals are excluded

by_species <- TRUE
# FALSE = by genera

ind <- which(data_all[,'ORDER']!='Sirenia')
data_all <- data_all[ind,]
ind <- which(data_all[,'ORDER']!='Cetacea')
data_all <- data_all[ind,]
ind <- which(data_all[,'ORDER']!='incertae sedis')
data_all <- data_all[ind,]


#data_all[,'LAT'] <- round(data_all[,'LAT'],digits = 1)
#data_all[,'LONG'] <- round(data_all[,'LONG'],digits = 1)
#data_all[,'paleolat'] <- round(data_all[,'paleolat'],digits = 1)
#data_all[,'paleolong'] <- round(data_all[,'paleolong'],digits = 1)

if (by_species){
  genusspecies <- paste( data_all[,'GENUS'], data_all[,'SPECIES']) 
  print('by species')
}else{
  genusspecies <- data_all[,'GENUS']  
  print('by genus')
}

midage <- (data_all[,'MAX_AGE'] + data_all[,'MIN_AGE'])/2
midage <- round(midage,digits = 2)
diffage <- (data_all[,'MAX_AGE'] - data_all[,'MIN_AGE'])
data_all <- cbind(data_all,genusspecies,midage,diffage)

print(dim(data_all))
ind <- which(data_all[,'diffage'] <= age_threshold)
data_all <- data_all[ind,]

print(dim(data_all))
indet_set <- c('Gen.', 'gen.','indet.','sp.','Indet','indet','Gen','gen','Sp','sp','incertae sedis')
ind <- which(!(data_all[,'GENUS'] %in% indet_set))
data_all <- data_all[ind,]
if (by_species){
  ind <- which(!(data_all[,'SPECIES'] %in% indet_set))
  data_all <- data_all[ind,]
  print(dim(data_all))
}

# getting rid of genera or species that have alternative identificiations
ind <- which(!grepl( '/', data_all[,'GENUS'], fixed = TRUE))
data_all <- data_all[ind,]
if (by_species){
  ind <- which(!grepl( '/', data_all[,'SPECIES'], fixed = TRUE))
  data_all <- data_all[ind,]
}
print(dim(data_all))

un_sp <- unique(data_all[,'genusspecies'])
print('No. unique taxa')
print(length(un_sp))

data_sum <- c()

for (sk in 1:length(un_sp)){
  sp_now <- un_sp[sk]
  
  ind_sp_loc <- which(data_all[,'genusspecies']==sp_now)
  
  no_loc <- length(ind_sp_loc)
  no_countries <- length(unique(data_all[ind_sp_loc,'COUNTRY']))
  min_mid_age <- min(data_all[ind_sp_loc,'midage'])
  max_mid_age <- max(data_all[ind_sp_loc,'midage'])
  sp_duration <- round(max_mid_age - min_mid_age,digits = 2)
  
  min_lat <- min(data_all[ind_sp_loc,'LAT'])
  ind_geo <- which(data_all[ind_sp_loc,'LAT']==min_lat)
  max_lat <- max(data_all[ind_sp_loc,'LAT'])
  ind_geo <- c(ind_geo,which(data_all[ind_sp_loc,'LAT']==max_lat))
  min_long <- min(data_all[ind_sp_loc,'LONG'])
  ind_geo <- c(ind_geo,which(data_all[ind_sp_loc,'LONG']==min_long))
  max_long <- max(data_all[ind_sp_loc,'LONG'])
  ind_geo <- c(ind_geo,which(data_all[ind_sp_loc,'LONG']==max_long))
  age_peak <- mean(data_all[ind_sp_loc[ind_geo],'midage'])
  duration_peak <- age_peak - min_mid_age
  
  diff_lat <- round(max_lat - min_lat,digits = 2)
  diff_long <- round(max_long - min_long,digits = 2)
  sp_area_el <- 3.14*(diff_lat/2)*(diff_long/2)*111*111/1000000
  sp_area_square <- diff_lat*diff_long*111*111/1000000
  sp_range_width <- sqrt(sp_area_square)
  
  min_pal_lat <- min(data_all[ind_sp_loc,'paleolat'])
  max_pal_lat <- max(data_all[ind_sp_loc,'paleolat'])
  min_pal_long <- min(data_all[ind_sp_loc,'paleolong'])
  max_pal_long <- max(data_all[ind_sp_loc,'paleolong'])
  
  diff_pal_lat <- round(max_pal_lat - min_pal_lat,digits = 2)
  diff_pal_long <- round(max_pal_long - min_pal_long,digits = 2)
  
  sp_pal_area_sq <- diff_pal_lat*diff_pal_long*111*111/1000000

  data_sp <- data_all[ind_sp_loc[1],c('SIDNUM','ORDER','FAMILY','SUBFAMILY','GENUS','SPECIES','genusspecies','BODYMASS','TCRWNHT','CROWNTYP')]
  data_sum <- rbind(data_sum, c(data_sp,sp_now,no_loc,no_countries,sp_duration,max_mid_age,min_mid_age,sp_area_square,sp_area_el,sp_range_width,max_lat,min_lat,max_long,min_long,duration_peak,sp_pal_area_sq))
}

colnames(data_sum) <- c('SIDNUM','ORDER','FAMILY','SUBFAMILY','GENUS','SPECIES','genusspecies','BODYMASS','TCRWNHT','CROWNTYP','sp_now','no_loc','no_countries','sp_duration','max_mid_age','min_mid_age','sp_area_square','sp_area_el','sp_range_width','max_lat','min_lat','max_long','min_long','duration_peak','sp_pal_area_sq')

write.table(data_sum, file = "data_working/data_sum.csv",col.names = TRUE,row.names = FALSE, sep = '\t')   

plot_name <- paste('plots/plot_dist_',palaeocord,'.pdf',sep = '')
pdf(plot_name,height = 6,width = 6)
plot(NA,NA,xlim = c(0,50), ylim = c(0,50),xlab = 'Area in present day coords, km2', ylab = 'Are in palaeocoords, km2',main = palaeocord)
points(data_sum[,'sp_area_square'],data_sum[,'sp_pal_area_sq'],pch = 16)
dev.off()

sp_test = 27190

ind <- which(data_all[,'SIDNUM']==sp_test)
print(data_all[ind,])
data_test <- data_all[ind,]
ind25 <- which(data_test[,'roundage']==25)
ind26 <- which(data_test[,'roundage']==26)
ind28 <- which(data_test[,'roundage']==28)
ind30 <- which(data_test[,'roundage']==30)

pdf('plots/plot_test_species_NOW.pdf',height = 4.1, width = 3.6)
plot(data_test[,'LONG'],data_test[,'LAT'],pch = 16,main = paste(data_all[ind[1],'GENUS'],data_all[ind[1],'SPECIES'],'\n NOW'),xlim = c(0,10),ylim = c(40,50),xlab = 'Long now',ylab = 'Lat now')
points(data_test[ind25,'LONG'],data_test[ind25,'LAT'],pch = 16,col = "#E69F00",cex  = 0.7)
points(data_test[ind26,'LONG'],data_test[ind26,'LAT'],pch = 16,col = "#009E73",cex  = 0.7)
points(data_test[ind28,'LONG'],data_test[ind28,'LAT'],pch = 16,col = "#CC79A7",cex  = 0.7)
points(data_test[ind30,'LONG'],data_test[ind30,'LAT'],pch = 16,col = "#0072B2",cex  = 0.7)
dev.off()

plot_name_pal <- paste('plots/plot_test_species_',palaeocord,'.pdf',sep = '')
pdf(plot_name_pal,height = 4.1, width = 3.6)
plot(data_test[,'paleolong'],data_test[,'paleolat'],pch = 16,main = paste(data_all[ind[1],'GENUS'],data_all[ind[1],'SPECIES'],'\n',palaeocord),xlim = c(0,10),ylim = c(40,50),xlab = 'Long now',ylab = 'Lat now')
points(data_test[ind25,'paleolong'],data_test[ind25,'paleolat'],pch = 16,col = "#E69F00",cex  = 0.7)
points(data_test[ind26,'paleolong'],data_test[ind26,'paleolat'],pch = 16,col = "#009E73",cex  = 0.7)
points(data_test[ind28,'paleolong'],data_test[ind28,'paleolat'],pch = 16,col = "#CC79A7",cex  = 0.7)
points(data_test[ind30,'paleolong'],data_test[ind30,'paleolat'],pch = 16,col = "#0072B2",cex  = 0.7)
dev.off()

write.table(data_test, file = "data_working/data_test.csv",col.names = TRUE,row.names = FALSE, sep = '\t')   

sp_test = 20243

ind <- which(data_all[,'SIDNUM']==sp_test)
data_test <- data_all[ind,]

pdf('plots/plot_test_species_NOW.pdf',height = 4.1, width = 3.6)
plot(data_test[,'LONG'],data_test[,'LAT'],pch = 16,main = paste(data_all[ind[1],'GENUS'],data_all[ind[1],'SPECIES'],'\n NOW'),xlim = c(-20,60),ylim = c(-0,80),xlab = 'Long now',ylab = 'Lat now')
dev.off()

plot_name_pal <- paste('plots/plot_test_species_',palaeocord,'.pdf',sep = '')
pdf(plot_name_pal,height = 4.1, width = 3.6)
plot(data_test[,'paleolong'],data_test[,'paleolat'],pch = 16,main = paste(data_all[ind[1],'GENUS'],data_all[ind[1],'SPECIES'],'\n',palaeocord),xlim = c(-20,60),ylim = c(0,80),xlab = 'Long now',ylab = 'Lat now')
dev.off()

write.table(data_test, file = "data_working/data_test.csv",col.names = TRUE,row.names = FALSE, sep = '\t')   
