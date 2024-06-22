# 2021 01 21 I.Zliobaite

# random seed to reproduce the plots exactly
set.seed(1981)

# parameters
n_years <- 12000 #also number of species thousands of years, one species per 2x thousand of years
step_years <- 2000
tsteps <- 100
fossilization_rate <- 1/10^11
lifespan <- 14 #years
pop_density <- 10 #individuals per km2
lambdat <- 0.3
lambdapeak <- 1.44

spcies_statistics <- c()
species_alive <- rep(0,n_years)
species_fossils <- rep(0,n_years)
species_fossils_known <- rep(0,n_years)

species_alive_tsteps <- rep(0,n_years/tsteps)
species_fossils_tsteps <- rep(0,n_years/tsteps)
species_fossils_known_tsteps <- rep(0,n_years/tsteps)

# 75 kg
# 14 years
# 33 individuals per km2

return_big_time <- function(sk_now,vector_now,t_step,t_total){
  species_now <- rep(0,t_total)
  t_vec <- length(vector_now)
  species_now[sk_now:(sk_now+t_vec-1)] <- vector_now
  n_steps <- t_total/t_step # must be integer
  species_big <- rep(0,n_steps)
  for (sk_cik in 1:n_steps){
    ind_start_now <- (sk_cik-1)*t_step+1
    ind_end_now <- sk_cik*t_step
    sum_within <- sum(species_now[ind_start_now:ind_end_now])
    if (sum_within>0){
      species_big[sk_cik] <- 1
    }
    
  }
  return(species_big)
}

for (sk in 1:n_years){
  if (sk %% 1000 == 0){
    print(sk)
  }
  #ex <- rexp(10000, rate = 1.85)
  
  dur0 <- rexp(1, rate = lambdat)*1000000
  dur <- ceiling(dur0/step_years)
  peak <- rexp(1, rate = lambdapeak)*1000 #km
  peak <- peak*peak
  abu <- round(peak*pop_density*step_years/lifespan)
  
  if (dur == 1){
    vector_alive <- abu
  }else{
    if (dur %% 2 == 0){
      half_life <- dur/2
      slope <- abu/half_life
      first_half_years <- 1:half_life
      last_half_years  <- (half_life+1):dur
      first_abundances <- round(first_half_years*slope)
      last_abundances <- rev(first_abundances)
      vector_alive <- c(first_abundances,last_abundances)
    }else{
      half_life <- ceiling(dur/2)
      slope <- abu/half_life
      first_half_years <- 1:half_life
      last_half_years  <- (half_life):dur
      first_abundances <- round(first_half_years*slope)
      last_abundances <- rev(first_abundances)
      last_abundances <- last_abundances[2:length(last_abundances)]
      vector_alive <- c(first_abundances,last_abundances)
    }
  }
  
  # 1 - dbinom(0,1000,0.001) <- probability to fossilize per species
  
  vector_fossils <- c()
  for (sk2 in 1:length(vector_alive)){
    vector_fossils <- c(vector_fossils,rbinom(1,vector_alive[sk2],fossilization_rate))
  }
  ind <- which(vector_fossils>0)
  vector_fossils[ind] <- 1
  if (length(ind)>0){
    ind_start <- min(ind)
    ind_end <- max(ind)  
    vector_fossils_known <- vector_fossils
    vector_fossils_known[ind_start:ind_end] <- 1
  }else{
    vector_fossils_known <- vector_fossils
  }
  ind <- which(vector_alive>0)
  vector_alive[ind] <- 1
  
  
  
  n_end <- sk-1+length(vector_alive)

  if (n_end>n_years){
    diff <-  n_end - n_years
    vector_alive0 <- vector_alive[1:(length(vector_alive)-diff)]
    vector_fossils0 <- vector_fossils[1:(length(vector_fossils)-diff)]
    vector_fossils_known0 <- vector_fossils_known[1:(length(vector_fossils_known)-diff)]
  }else{
    vector_alive0 <- vector_alive
    vector_fossils0 <- vector_fossils
    vector_fossils_known0 <- vector_fossils_known
  }
  
  t_alive <- length(vector_alive0)
  species_alive[sk:(sk+t_alive-1)] <- species_alive[sk:(sk+t_alive-1)] + vector_alive0
  species_fossils[sk:(sk+t_alive-1)] <- species_fossils[sk:(sk+t_alive-1)] + vector_fossils0
  species_fossils_known[sk:(sk+t_alive-1)] <- species_fossils_known[sk:(sk+t_alive-1)] + vector_fossils_known0
  
  species_alive_big_now <- return_big_time(sk,vector_alive0,tsteps,n_years)
  species_alive_tsteps <- species_alive_tsteps + species_alive_big_now
  
  species_fossils_big_now <- return_big_time(sk,vector_fossils0,tsteps,n_years)
  species_fossils_tsteps <- species_fossils_tsteps + species_fossils_big_now
  
  species_fossils_known_big_now <- return_big_time(sk,vector_fossils_known0,tsteps,n_years)
  species_fossils_known_tsteps <- species_fossils_known_tsteps + species_fossils_known_big_now
  
  spcies_statistics <- rbind(spcies_statistics,c(sk,dur0,dur,peak,abu,sum(vector_alive),sum(vector_fossils),sum(vector_fossils_known),sum(vector_alive0),sum(vector_fossils0),sum(vector_fossils_known0),sum(species_alive_big_now),sum(species_fossils_big_now)))
}

colnames(spcies_statistics) <- c('t','duration_y','duration_steps','peak_km2','peak_abu','duration_live','time_fos','duration_fos','duration_live_short','time_fos_short','duration_fos_short','duration_live_big','duration_fos_big')
  
pdf('fig_species_alive.pdf',height = 6, width = 10)
plot(species_alive, type = 'l')
points(species_fossils, type = 'l',col='darkgreen')
points(species_fossils_known, type = 'l',col='darkblue')
dev.off()

pdf('fig_species_alive_big.pdf',height = 6, width = 10)
plot(species_alive_tsteps, type = 'l')
points(species_fossils_tsteps, type = 'l',col='darkgreen')
points(species_fossils_known_tsteps, type = 'l',col='darkblue')
dev.off()


spcies_statistics <- spcies_statistics[2001:dim(spcies_statistics)[1],]

pdf('fig_durations_alive_vs_fossils.pdf',height = 7, width = 7)
plot(spcies_statistics[,'duration_live_short'],spcies_statistics[,'duration_fos_short'],pch = 16)
dev.off()

ind <- which(spcies_statistics[,'duration_fos_short']!=0)

pdf('fig_durations_vs_abu_fossils_log.pdf',height = 7, width = 7)
plot(log10(spcies_statistics[,'duration_steps']),log10(spcies_statistics[,'peak_abu']),pch = 16)
points(log10(spcies_statistics[ind,'duration_steps']),log10(spcies_statistics[ind,'peak_abu']),pch = 16,col='darkred')
dev.off()

pdf('fig_durations_vs_abu_fossils.pdf',height = 7, width = 7)
plot(spcies_statistics[,'duration_steps'],spcies_statistics[,'peak_abu'],pch = 16)
points(spcies_statistics[ind,'duration_steps'],spcies_statistics[ind,'peak_abu'],pch = 16,col='darkred')
dev.off()


library(MASS)

fit1 <- fitdistr(spcies_statistics[ind,'duration_fos_short'], "exponential")
#tt1 <- ks.test(spcies_statistics[ind,'duration_fos_short'], "pexp", fit1$estimate)

pdf('fig_hist.pdf')
hist(spcies_statistics[ind,'duration_fos_short'], freq = FALSE, breaks = 50)
curve(dexp(x, rate = fit1$estimate), from = 0, col = "black", add = TRUE)
dev.off()


fit2 <- fitdistr(spcies_statistics[ind,'duration_live_short'], "exponential")
#tt2 <- ks.test(spcies_statistics[ind,'duration_live_short'], "pexp", fit2$estimate)

pdf('fig_hist_live.pdf')
hist(spcies_statistics[ind,'duration_live_short'], freq = FALSE, breaks = 50)
curve(dexp(x, rate = fit2$estimate), from = 0, col = "black", add = TRUE)
dev.off()


law_law_live <- c()
law_law_livefos <- c()
law_law_fos <- c()

dur_live_all <- spcies_statistics[,'duration_live_short']*2000/1000000
dur_live_all <- dur_live_all[order(dur_live_all)]
un_dur_live_all <- unique(dur_live_all)
n <- length(dur_live_all)
for (sk in 1:length(un_dur_live_all)){
  ind1<- which(dur_live_all>un_dur_live_all[sk])
  law_law_live <- rbind(law_law_live,c(un_dur_live_all[sk],length(ind1)/n))
}

dur_livefos <- spcies_statistics[ind,'duration_live_short']*2000/1000000
dur_livefos <- dur_livefos[order(dur_livefos)]
un_dur_livefos <- unique(dur_livefos)
n_livefos <- length(dur_livefos)
for (sk in 1:length(un_dur_livefos)){
  ind2<- which(dur_livefos>un_dur_livefos[sk])
  law_law_livefos <- rbind(law_law_livefos,c(un_dur_livefos[sk],length(ind2)/n_livefos))
}


dur_fos <- spcies_statistics[ind,'duration_fos_short']*2000/1000000
dur_fos <- dur_fos[order(dur_fos)]
un_dur_fos <- unique(dur_fos)
n_fos <- length(dur_fos)
for (sk in 1:length(un_dur_fos)){
  ind2<- which(dur_fos>un_dur_fos[sk])
  law_law_fos <- rbind(law_law_fos,c(un_dur_fos[sk],length(ind2)/n_fos))
}

law_law_live <- cbind(law_law_live,log10(law_law_live[,2]))
colnames(law_law_live) <- c('duration','psurvive','log10p')
nnn <- dim(law_law_live)[1]
law_law_live <- law_law_live[1:(nnn-1),]

law_law_livefos <- cbind(law_law_livefos,log10(law_law_livefos[,2]))
colnames(law_law_livefos) <- c('duration','psurvive','log10p')
nnn <- dim(law_law_livefos)[1]
law_law_livefos <- law_law_livefos[1:(nnn-1),]

law_law_fos <- cbind(law_law_fos,log10(law_law_fos[,2]))
colnames(law_law_fos) <- c('duration','psurvive','log10p')
nnn <- dim(law_law_fos)[1]
law_law_fos <- law_law_fos[1:(nnn-1),]

write.table(law_law_live, file = "law_law_live.csv",col.names = TRUE,row.names = FALSE, sep = '\t')   
law_law_live <- read.csv('law_law_live.csv', header = TRUE, sep = "\t")
fit_law_live <- lm(log10p ~ duration,data = law_law_live)
cf <- coef(fit_law_live)
intercept_live <- paste('intercept',round(cf[1],digits = 3))
slope_live <- paste('slope',round(cf[2],digits = 3))

write.table(law_law_livefos, file = "law_law_livefos.csv",col.names = TRUE,row.names = FALSE, sep = '\t')   
law_law_livefos <- read.csv('law_law_livefos.csv', header = TRUE, sep = "\t")
fit_law_livefos <- lm(log10p ~ duration,data = law_law_livefos)
cf <- coef(fit_law_livefos)
intercept_livefos <- paste('intercept',round(cf[1],digits = 3))
slope_livefos <- paste('slope',round(cf[2],digits = 3))

write.table(law_law_fos, file = "law_law_fos.csv",col.names = TRUE,row.names = FALSE, sep = '\t')   
law_law_fos <- read.csv('law_law_fos.csv', header = TRUE, sep = "\t")
fit_law_fos <- lm(log10p ~ duration,data = law_law_fos)
cf <- coef(fit_law_fos)
intercept_fos <- paste('intercept',round(cf[1],digits = 3))
slope_fos <- paste('slope',round(cf[2],digits = 3))


pdf('fig_law_simulation1.pdf', width = 4, height = 4.5)
plot(law_law_live[,1],log10(law_law_live[,2]),pch=16,cex = 0.7,xlab = 'time, Myr',ylab = 'log10 (Proportion of species surviving)',main = 'Simulation1',xlim = c(0,18),ylim = c(-4,0),col='darkgrey')
#points(law_law_livefos[,1],log10(law_law_livefos[,2]),pch=16,cex = 0.7,col='grey')
points(law_law_fos[,1],log10(law_law_fos[,2]),pch=16,cex = 0.7,col='black')
abline(fit_law_live,col='darkgrey')
#abline(fit_law_livefos,col='grey')
abline(fit_law_fos,col='black')
text(14,-0.2,'Simulated living',col = 'darkgrey',cex = 0.7)
text(15.2,-0.5,intercept_live,col = 'darkgrey',cex = 0.7)
text(16,-0.8,slope_live,col = 'darkgrey',cex = 0.7)
text(14,-1.1,'Simulated fossilized',col = 'black',cex=0.7)
text(15.2,-1.4,intercept_fos,col = 'black',cex = 0.7)
text(16,-1.7,slope_fos,col = 'black',cex = 0.7)
dev.off()

