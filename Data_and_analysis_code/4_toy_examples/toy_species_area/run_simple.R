# 2023 08 19 I.Zliobaite


n_species <- 100
n_area_per_species <- 10
n_areas <- 1000

which_survivorship <- 'powerlaw'
#  c('simple','exponential_area','exponential_width')

data_species_area <- matrix(0, n_species, n_areas)

library("poweRlaw")

for (sk in 1:n_species){
  if (which_survivorship=='simple'){
    file_name1 <- 'fig_species_area_simple.pdf'
    file_name2 <- 'fig_species_area_simple_log.pdf'
  }
  if (which_survivorship=='exponential_area'){
    n_area_per_species <- round(rexp(1, rate = 0.1))
    file_name1 <- 'fig_species_area_exp_area.pdf'
    file_name2 <- 'fig_species_area_exp_area_log.pdf'
  }
  if (which_survivorship=='exponential_width'){
    n_area_per_species <- round(rexp(1, rate = 1/sqrt(10))^2)
    file_name1 <- 'fig_species_area_exp_width.pdf'
    file_name2 <- 'fig_species_area_exp_width_log.pdf'
  }
  if (which_survivorship=='powerlaw'){
    n_area_per_species <- round(rexp(1, rate = 1/sqrt(10))^2)
    file_name1 <- 'fig_species_area_pow.pdf'
    file_name2 <- 'fig_species_area_pow_log.pdf'
  }
  ind_where <- sample(1:n_areas, n_area_per_species, replace=FALSE) 
  data_species_area[sk,ind_where] <- 1
}

data_power <- c(1,mean(apply(data_species_area,2,sum)))

for (sk in 2:1000){
  print(sk)
  sample <- c()
  for (sk2 in 1:20){
    ind_where <- sample(1:n_areas, sk, replace=FALSE) 
    su <- apply(data_species_area[,ind_where],1,sum)
    ind <- which(su>0)
    su[ind] <- 1
    sample <- c(sample,sum(su))
  }
  data_power <- rbind(data_power,c(sk,mean(sample)))
}

pdf(file_name1,width = 4.5,height = 5)
plot(data_power[,1],data_power[,2])
dev.off()

pdf(file_name2,width = 4.5,height = 5)
plot(log10(data_power[,1]),log10(data_power[,2]))
dev.off()
