# 2023 03 16 I.Zliobaite
# 2024 01 20 adding hominin cleaning

data_all <- read.csv('data_raw/now_pub_20230818.csv', header = TRUE, sep = "\t", na.string = '\\N')

data_hom <- read.csv('data_working/data_sum_hominini_man.csv', header = TRUE, sep = "\t", na.string = c('\\N','NA',''))

do_lumping <- FALSE
do_bootstrap <- FALSE

if (do_bootstrap){
  
  set.seed(1982)
  
  boo <- sample(1:dim(data_all)[1],replace = TRUE)
  boo <- unique(boo)    
  
  data_all <- data_all[boo,]
  
  print('bootstrapping yes')
}



ind_hominini <- which(data_all[,'SUBFAMILY']=='Hominini')

ind <- intersect(ind_hominini,which(data_all[,'SPECIES']=='rhodesiensis'))
data_all[ind,'SPECIES'] <- 'heidelbergensis'



if (do_lumping){
  for (sk in 1:dim(data_hom)[1]){
    sp_now <- data_hom[sk,'SPECIES']
    ind <- intersect(ind_hominini,which(data_all[,'SPECIES']==sp_now))
    if (is.na(data_hom[sk,'Synsp'])){
      data_all[ind,'SUBFAMILY'] <- 'Hominini0'
      print(sp_now)
    }else{
      data_all[ind,'SPECIES'] <- data_hom[sk,'Synsp']
    }
  }
}

age_threshold <- 3 #localities with larger than this time intervals are excluded

by_species <- TRUE
# FALSE = by genera

ind <- which(data_all[,'ORDER']!='Sirenia')
data_all <- data_all[ind,]
ind <- which(data_all[,'ORDER']!='Cetacea')
data_all <- data_all[ind,]
ind <- which(data_all[,'ORDER']!='incertae sedis')
data_all <- data_all[ind,]


data_all[,'LAT'] <- round(data_all[,'LAT'],digits = 1)
data_all[,'LONG'] <- round(data_all[,'LONG'],digits = 1)

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
  #sp_range_width<- sqrt(sp_area_el)
  sp_range_width <- sqrt(sp_area_square)

  data_sp <- data_all[ind_sp_loc[1],c('SIDNUM','ORDER','FAMILY','SUBFAMILY','GENUS','SPECIES','genusspecies','BODYMASS','TCRWNHT','CROWNTYP')]
  data_sum <- rbind(data_sum, c(data_sp,sp_now,no_loc,no_countries,sp_duration,max_mid_age,min_mid_age,sp_area_square,sp_area_el,sp_range_width,max_lat,min_lat,max_long,min_long,duration_peak))
}

colnames(data_sum) <- c('SIDNUM','ORDER','FAMILY','SUBFAMILY','GENUS','SPECIES','genusspecies','BODYMASS','TCRWNHT','CROWNTYP','sp_now','no_loc','no_countries','sp_duration','max_mid_age','min_mid_age','sp_area_square','sp_area_el','sp_range_width','max_lat','min_lat','max_long','min_long','duration_peak')

if (do_lump){
  if (do_bootstrap){
    write.table(data_sum, file = "data_working/data_sum_lump_boo.csv",col.names = TRUE,row.names = FALSE, sep = '\t')       
  }else{
    write.table(data_sum, file = "data_working/data_sum_lump.csv",col.names = TRUE,row.names = FALSE, sep = '\t')       
  }
  }else{
    if (do_bootstrap){
      write.table(data_sum, file = "data_working/data_sum_split_boo.csv",col.names = TRUE,row.names = FALSE, sep = '\t')       
    }else{
      write.table(data_sum, file = "data_working/data_sum_split.csv",col.names = TRUE,row.names = FALSE, sep = '\t')       
    }
}


