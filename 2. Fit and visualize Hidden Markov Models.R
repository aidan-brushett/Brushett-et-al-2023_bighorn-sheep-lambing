##### SET-UP #####
library(ggplot2)
library(depmixS4)
library(tidyr)
library(cowplot)
library(zoo)
library(dplyr)
library(lubridate)
library(Rsolnp)

select <- dplyr::select
filter <- dplyr::filter
graph.theme <- theme_classic() + 
  theme(axis.line.x.bottom=element_blank(),
        axis.line.y.left =element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))
xaxis_plotlabel <- ggdraw() + 
  draw_label("Date",
             hjust = 0.5,
             size=12)
setwd("C:/_sheep")


info <- read.csv("data/validated-parturient-sheep.csv") %>% mutate(
  windowstart=as.Date(windowstart),
  windowend=as.Date(windowend),
  id=paste(name, year, sep="-")) %>%
  rowwise %>%
  mutate(
    monthstart=mean.Date(c(windowstart, windowend))-days(15),
    monthend=mean.Date(c(windowstart, windowend))+days(15))
info2 <- read.csv("data/validated-non-parturient-sheep.csv") %>% mutate(id=paste(name, year, sep="-"))

##### READ IN THE DATA #####
#data8 <- data7 %>% 
#  mutate(year=year(date), month=month(date), day=day(date))
# OR
data8 <- read.csv("data/BHS-female_averaged-movement-terrain_clean.csv") %>% 
  mutate(date=as.POSIXct(date), year=year(date), month=month(date), day=day(date), id=gsub("_", "-", id))





set.seed(69420)
data9 <- data8 %>% filter(id%in%info$id, month==5&day>=15 | month==6 | month==7&day<=15) %>%
  mutate(speed=ifelse(speed==0, 0.0001, speed)) %>% # a lazy fix to avoid needing to use zero-inflation
  select(-X, -X.1) %>%
  drop_na() %>%
  mutate(hr.48=hr.48/100)



data9 <- data9 %>% mutate(hr.48=hr.48/hrmax,
                          rt100=rt100/rtmax,
                          speed=speed/speedmax)


##### FIT THE HIDDEN MARKOV MODEL #####






unique(data9$id)
# These are all inputs/parameters for the depmixS4 functions
K <- 3
ntimes <- numeric(length=length(unique(data9$id)))
for (i in 1:length(unique(data9$id))) {
  temp <- data9 %>% filter(id==unique(data9$id)[i], month %in% 5:7)
  ntimes[i] <- nrow(temp)
  rm(temp)
  }

# Fit the model for the specified animals and time frame:
mod <- depmix(list(speed~1, rt100~1, hr.48~1), data=data9, nstates=K, 
              family=list(Gamma(link="log"), Gamma(link="log"), Gamma(link="log")),
              ntimes=ntimes)
mod <- fit(mod, emcontrol=em.control(maxit=1000, tol=0.01, crit="relative"))
#mod <- fit(mod, verbose=T, method="rsolnp")

data10 <- cbind(data9, posterior(mod, type="viterbi"))

print(mod)


summary(mod)







fill=c("3"="red", "2"="grey90", "1"="grey60")
# Plot the viterbi behavioural states against the raw data ()
{pdf("TESTWITHTOL/HMM-parturient_viterbistates_May15-Jul15.pdf", width=6, height=5.5)
  
  # This is just for assigning the correct colours to each movement state

  for(i in 1:length(unique(data10$id))){
    temp <- data10%>%filter(id==unique(data10$id)[i])
    plotA <- ggplot(data=temp, aes(date, rt100.48)) + 
      geom_rect(data=temp, color=NA, aes(ymin=-Inf, ymax=Inf, xmin=date, 
                                         xmax=lead(date), fill=factor(state, levels=c("1", "2", "3"))), alpha=0.4) +
      scale_fill_manual(values=fill) +
      geom_line() +
      graph.theme +
      scale_x_datetime(expand=c(0, 0)) +
      theme(legend.position="none") +
      theme(axis.title.x=element_blank()) +
      labs(y="RT100 (h)", x="Date") +
      theme(strip.text=element_text(face="bold"), strip.background=element_rect(colour=NA))
    
    plotB <- ggplot(data=temp, aes(date, dist.48)) + 
      geom_rect(data=temp, color=NA, aes(ymin=-Inf, ymax=Inf, xmin=date, 
                                          xmax=lead(date), fill=factor(state)), alpha=0.4) +
      scale_fill_manual(values=fill) +
      geom_line() +
      graph.theme +
      scale_x_datetime(expand=c(0, 0)) +
      theme(legend.position="none") +
      theme(axis.title.x=element_blank()) +
      labs(y="DIST (m)", x="Date") +
      theme(strip.text=element_text(face="bold"), strip.background=element_rect(colour=NA))
    
    plotC <- ggplot(data=temp, aes(date, hr.48)) + 
      geom_rect(data=temp, color=NA, aes(ymin=-Inf, ymax=Inf, xmin=date, 
                                          xmax=lead(date), fill=factor(state)), alpha=0.4) +
      scale_fill_manual(values=fill) +
      geom_line() +
      graph.theme +
      scale_x_datetime(expand=c(0, 0)) +
      theme(legend.position="none") +
      theme(axis.title.x=element_blank()) +
      labs(y="HR (Ha)", x="Date") +
      theme(strip.text=element_text(face="bold"), strip.background=element_rect(colour=NA))
    
    title <- ggdraw() + 
      draw_label(
        paste(unique(data10$id)[i]),
        fontface = 'bold',
        x = 0,
        hjust =-0.5)
    
    plotABC <- plot_grid(title, plotA, NULL, plotB, NULL, plotC, xaxis_plotlabel, 
                         ncol=1, align="v", rel_heights=c(0.1, 1, -0.12, 1, -0.12, 1, 0.08))
    
    print(plotABC)
    rm(temp, title)
  }
dev.off()
rm(middlelevel, upperlevel, lowerlevel)
}

# Plot the viterbi behavioural states against the Pr[State3]
{pdf("TESTWITHTOL/HMM-parturient_state-probability_May15-Jul15.pdf", width=7, height=3)

middlelevel <- data10[which.max(data10$rt100.48),]$state
lowerlevel <- ifelse(middlelevel==3, 1, 3)
upperlevel <- ifelse(middlelevel==2, 1, 2)

for(i in 1:length(unique(data10$id))){
  temp <- data10%>%filter(id==unique(data10$id)[i])
  plotD <- ggplot(data=temp, aes(date, S3)) + 
    geom_rect(data=temp, color=NA, aes(ymin=-Inf, ymax=Inf, xmin=date, 
                                       xmax=lead(date), fill=factor(state)), alpha=0.4) +
    scale_fill_manual(values=fill) +
    geom_line(size=0.2) +
    graph.theme +
    scale_x_datetime(expand=c(0, 0)) +
    theme(legend.position="none") +
    labs(y="Probability of low-movement state", x="Date") +
    theme(strip.text=element_text(face="bold"), strip.background=element_rect(colour=NA))
  
  title <- ggdraw() + 
    draw_label(
      paste(unique(data10$id)[i]),
      fontface = 'bold',
      x = 0,
      hjust =-0.5)
  
  plotD <- plot_grid(title, plotD, ncol=1, align="v", rel_heights=c(0.1, 1))
  
  print(plotD)
  rm(temp, title)
  }
dev.off()
rm(middlelevel, upperlevel, lowerlevel)
}


##### Summarize the data and determine inferred parturition dates -----
data11 <- data10 %>% mutate(parturition=ifelse(state==3, 1, 0)) %>%
  group_by(id) %>%
  mutate(sumpart.7=rollsum(parturition, 84, align="left", na.pad=T)/84,
         sumpart.1=rollsum(parturition, 12, align="left", na.pad=T)/12,
         sumpart.2=rollsum(parturition, 24, align="left", na.pad=T)/24,
         sumpart.3=rollsum(parturition, 36, align="left", na.pad=T)/36,
         sumpart.4=rollsum(parturition, 48, align="left", na.pad=T)/48) %>%
  ungroup()
write.csv(data11, "data/BHS-female_averaged-movement-terrain-hmm-states-viterbi-rollingavg_clean.csv")

parturitionsummary <- data11 %>% group_by(id) %>%
  summarise(maxsumpart.2=max(sumpart.2, na.rm=T),
            meansumpart.2=mean(sumpart.2, na.rm=T),
            sdsumpart.2=sd(sumpart.2, na.rm=T),
            part.date=date[which.max(sumpart.2)],
            part.x=x[which.max(sumpart.2)],
            part.y=y[which.max(sumpart.2)])
parturitionsummary <- left_join(parturitionsummary, info, by="id")
write.csv(parturitionsummary, "data/BHS-female_inferred-parturition-sites-dates.csv")

# Plot the proportion of time in a given state in a rolling time window, AND the state. 
{pdf("TESTWITHTOL/HMM-parturient_2day-state-proportion_May15-Jul15.pdf", width=7, height=3)
  

  i=2
  for(i in 1:length(unique(data11$id))){
    temp <- data11%>%filter(id==parturitionsummary$id[i])
    plotE <- ggplot(data=temp, aes(date, sumpart.2)) + 
      geom_rect(data=temp, color=NA, aes(ymin=-Inf, ymax=Inf, xmin=date, 
                                         xmax=lead(date), fill=factor(state)), alpha=0.4) +
      scale_fill_manual(values=fill) +
      geom_line(size=0.04) +
      geom_vline(xintercept=parturitionsummary$part.date[i], color="blue") +
      graph.theme +
      scale_x_datetime(expand=c(0, 0)) +
      theme(legend.position="none") +
      labs(y="Proportion of 48-hr period \nspent in low-movement state", x="Date") +
      theme(strip.text=element_text(face="bold"), strip.background=element_rect(colour=NA))

    title <- ggdraw() + 
      draw_label(
        paste(unique(data10$id)[i]),
        fontface = 'bold',
        x = 0,
        hjust =-0.5)
    
    plotE <- plot_grid(title, plotE, ncol=1, align="v", rel_heights=c(0.1, 1))
    
    print(plotE)
    rm(temp, title)
  }
  dev.off()
  rm(middlelevel, upperlevel, lowerlevel)
  }







# Plot the viterbi low-movement state against the raw data and add the inferred parturition date
{pdf("TESTWITHTOL/RESIZED-HMM-parturient_viterbi-lambingdate-FINAL_May15-Jul15.pdf", width=6.69, height=5.7)
  fill=c("3"="lightpink1", "2"="grey85", "1"="white")
  
  # This is just for assigning the correct colours to each movement state
  # middlelevel <- data10[which.max(data10$rt100.48),]$state
  # parturitionstate <- middlelevel
  # lowerlevel <- ifelse(middlelevel==3, 1, 3)
  # upperlevel <- ifelse(middlelevel==2, 1, 2)
  
  for(i in 1:length(unique(data10$id))){
    
    temp <- data10%>%filter(id==parturitionsummary$id[i]) %>%
      drop_na() %>%
      mutate(state=as.factor(state),
             group=cumsum(state!=lag(state, default=first(state))),
             hr.48=hr.48*100, 
             hr.48=ifelse(hr.48>=100, 100, hr.48),
             dist=ifelse(dist>=500, 500, dist))
    
    
    plotA <- ggplot(data=temp, aes(date, rt100)) + 
      geom_rect(data=temp, color=NA, aes(ymin=-Inf, ymax=Inf, xmin=date, 
                                         xmax=lead(date), fill=factor(state)), alpha=1, group=temp$state) +
      scale_fill_manual(values=fill) +
      geom_vline(xintercept=parturitionsummary$part.date[i], color="red", size=1) +
      geom_line() +
      graph.theme +
      scale_x_datetime(expand=c(0, 0)) +
      theme(legend.position="none") +
      theme(axis.title.x=element_blank()) +
      labs(y="Residence Time (h)", x="Date") +
      theme(strip.text=element_text(face="bold"), strip.background=element_rect(colour=NA))


  plotB <- ggplot(data=temp, aes(date, dist)) + 
    geom_rect(data=temp, color=NA, aes(ymin=-Inf, ymax=Inf, xmin=date, 
                                       xmax=lead(date), fill=factor(state)), alpha=1, group=temp$state) +
    scale_fill_manual(values=fill) +    
    geom_vline(xintercept=parturitionsummary$part.date[i], color="red", size=1) +
    geom_line() +
    graph.theme +
    scale_x_datetime(expand=c(0, 0)) +
    theme(legend.position="none") +
    theme(axis.title.x=element_blank()) +
    labs(y="Step Length (m)", x="Date") +
    theme(strip.text=element_text(face="bold"), strip.background=element_rect(colour=NA))


  plotC <- ggplot(data=temp, aes(date, hr.48)) + 
    geom_rect(data=temp, color=NA, aes(ymin=-Inf, ymax=Inf, xmin=date, 
                                       xmax=lead(date), fill=factor(state)), alpha=1, group=temp$state) +
    scale_fill_manual(values=fill) +
    geom_vline(xintercept=parturitionsummary$part.date[i], color="red", size=1) +
    geom_line() +
    graph.theme +
    scale_x_datetime(expand=c(0, 0)) +
    theme(legend.position="none") +
    theme(axis.title.x=element_blank()) +
    labs(y="Home Range (Ha)", x="Date") +
    theme(strip.text=element_text(face="bold"), strip.background=element_rect(colour=NA))

  
  
    title <- ggdraw() + 
      draw_label(
        paste(unique(data10$id)[i], "(Parturient)", sep=" "),
        fontface = 'bold',
        x = 0,
        hjust =-0.5)
    
    plotABC <- plot_grid(title, plotA, NULL, plotB, NULL, plotC, xaxis_plotlabel, ncol=1, align="v", rel_heights=c(0.1, 1, -0.12, 1, -0.12, 1, 0.08))
    
    print(plotABC)
    rm(temp, title)
  }
  dev.off()
  rm(middlelevel, upperlevel, lowerlevel)
}















##### K-FOLD/LEAVE-ONE-OUT CROSS VALIDATION ON PARTURIENT SHEEP #####

K <- 3

for(i in 1:length(unique(data9$id))){
  id_test <- unique(data9$id)[i]
  data9k <- data9 %>% filter(id!=id_test)
  print(id_test)
  
  ntimes <- numeric(length=length(unique(data9k$id)))
  for (j in 1:length(unique(data9k$id))) {
    temp <- data9k %>% filter(id==unique(data9k$id)[j], month %in% 5:7)
    ntimes[j] <- nrow(temp)
    rm(temp)
    }
  
  modk <- depmix(list(speed~1, rt100~1, hr.48~1), data=data9k, nstates=K, 
                family=list(Gamma(link="log"), Gamma(link="log"), Gamma(link="log")),
                ntimes=ntimes)
  modk <- fit(mod, emcontrol=em.control(maxit=1000, tol=0.01, crit="relative"))
  print(modk)
  summary(modk)
  
  
  data9kb <- data9 %>% filter(id==id_test)
  modkNew <- depmix(list(speed~1, rt100~1, hr.48~1), data=data9kb, nstates=K, 
                family=list(Gamma(link="log"), Gamma(link="log"), Gamma(link="log")),
                ntimes=nrow(data9kb))
  modkNew <- setpars(modkNew,getpars(modk))
  modkNew <- fit(modkNew, emcontrol=em.control(maxit=1000, tol=0.01, crit="relative"))
  print(modkNew)
  summary(modkNew)
  
  data10kb <- cbind(data9kb, posterior(modkNew, type='viterbi'))
  
  partstate <- data10kb[which.max(data10kb$rt100.48),]$state
  faststate <- data10kb[which.max(data10kb$hr.48),]$state# parturitionstate <- middlelevel
  slowstate <- setdiff(c(1, 2, 3), c(partstate, faststate))
  
  data10kb <- data10kb %>%
    mutate(state=case_when(
      state==partstate ~ 1,
      state==faststate ~ 2,
      state==slowstate ~ 3,
      
    ))
  
  ifelse(i==1, data10k<-data10kb, data10k <- rbind(data10k, data10kb))
}  


data11k <- data10k %>% mutate(parturition=ifelse(state==1, 1, 0)) %>%
  group_by(id) %>%
  mutate(sumpart.7=rollsum(parturition, 84, align="left", na.pad=T)/84,
         sumpart.1=rollsum(parturition, 12, align="left", na.pad=T)/12,
         sumpart.2=rollsum(parturition, 24, align="left", na.pad=T)/24,
         sumpart.3=rollsum(parturition, 36, align="left", na.pad=T)/36,
         sumpart.4=rollsum(parturition, 48, align="left", na.pad=T)/48) %>%
  ungroup()

parturitionsummaryk <- data11k %>% group_by(id) %>%
  summarise(maxsumpart.2=max(sumpart.2, na.rm=T),
            meansumpart.2=mean(sumpart.2, na.rm=T),
            sdsumpart.2=sd(sumpart.2, na.rm=T),
            part.date=date[which.max(sumpart.2)],
            part.x=x[which.max(sumpart.2)],
            part.y=y[which.max(sumpart.2)])
parturitionsummaryk <- left_join(parturitionsummaryk, info, by="id")
write.csv(parturitionsummaryk, "kfold/BHS-kfold-female_inferred-parturition-sites-dates.csv")





# Plot the viterbi low-movement state against the raw data and add the inferred parturition date
{pdf("kfold/RESIZED-HMM-kfold-parturient_viterbi-lambingdate-FINAL_May15-Jul15.pdf", width=6.69, height=5.7)
  fill=c("1"="lightpink1", "3"="grey85", "2"="grey90")
  
  # This is just for assigning the correct colours to each movement state
  # middlelevel <- data10[which.max(data10$rt100.48),]$state
  # parturitionstate <- middlelevel
  # lowerlevel <- ifelse(middlelevel==3, 1, 3)
  # upperlevel <- ifelse(middlelevel==2, 1, 2)
  
  for(i in 1:length(unique(data10k$id))){
    
    temp <- data10k%>%filter(id==parturitionsummaryk$id[i]) %>%
      drop_na() %>%
      mutate(state=as.factor(state),
             group=cumsum(state!=lag(state, default=first(state))),
             hr.48=hr.48*100, 
             hr.48=ifelse(hr.48>=100, 100, hr.48),
             dist=ifelse(dist>=500, 500, dist))
    
    
    plotA <- ggplot(data=temp, aes(date, rt100)) + 
      geom_rect(data=temp, color=NA, aes(ymin=-Inf, ymax=Inf, xmin=date, 
                                         xmax=lead(date), fill=factor(state)), alpha=1, group=temp$state) +
      scale_fill_manual(values=fill) +
      geom_vline(xintercept=parturitionsummaryk$part.date[i], color="red", size=1) +
      geom_line() +
      graph.theme +
      scale_x_datetime(expand=c(0, 0)) +
      theme(legend.position="none") +
      theme(axis.title.x=element_blank()) +
      labs(y="Residence Time (h)", x="Date") +
      theme(strip.text=element_text(face="bold"), strip.background=element_rect(colour=NA))
    
    
    plotB <- ggplot(data=temp, aes(date, dist)) + 
      geom_rect(data=temp, color=NA, aes(ymin=-Inf, ymax=Inf, xmin=date, 
                                         xmax=lead(date), fill=factor(state)), alpha=1, group=temp$state) +
      scale_fill_manual(values=fill) +    
      geom_vline(xintercept=parturitionsummaryk$part.date[i], color="red", size=1) +
      geom_line() +
      graph.theme +
      scale_x_datetime(expand=c(0, 0)) +
      theme(legend.position="none") +
      theme(axis.title.x=element_blank()) +
      labs(y="Step Length (m)", x="Date") +
      theme(strip.text=element_text(face="bold"), strip.background=element_rect(colour=NA))
    
    
    plotC <- ggplot(data=temp, aes(date, hr.48)) + 
      geom_rect(data=temp, color=NA, aes(ymin=-Inf, ymax=Inf, xmin=date, 
                                         xmax=lead(date), fill=factor(state)), alpha=1, group=temp$state) +
      scale_fill_manual(values=fill) +
      geom_vline(xintercept=parturitionsummaryk$part.date[i], color="red", size=1) +
      geom_line() +
      graph.theme +
      scale_x_datetime(expand=c(0, 0)) +
      theme(legend.position="none") +
      theme(axis.title.x=element_blank()) +
      labs(y="Home Range (Ha)", x="Date") +
      theme(strip.text=element_text(face="bold"), strip.background=element_rect(colour=NA))
    
    
    
    title <- ggdraw() + 
      draw_label(
        paste(unique(data10$id)[i], "(Parturient)", sep=" "),
        fontface = 'bold',
        x = 0,
        hjust =-0.5)
    
    plotABC <- plot_grid(NULL, plotA, NULL, plotB, NULL, plotC, xaxis_plotlabel, ncol=1, align="v", rel_heights=c(0.1, 1, -0.12, 1, -0.12, 1, 0.08))
    
    print(plotABC)
    rm(temp, title)
  }
  dev.off()
  rm(middlelevel, upperlevel, lowerlevel)
}

print("hello world")



##### FIT MY HMM TO TEST DATA (aka non-parturient sheep) #####

data9b <- data8 %>% filter(id%in%info2$id, month==5&day>=15 | month==6 | month==7&day<=15) %>%
  mutate(speed=ifelse(speed==0, 0.0001, speed)) %>%
  select(-X, -X.1) %>% # a lazy fix to avoid needing to use zero-inflated gamma
  drop_na() %>%
  mutate(hr.48=hr.48/100)

K <- 3
ntimesb <- numeric(length=length(unique(data9b$id)))
for (i in 1:length(unique(data9b$id))) {
  temp <- data9b %>% filter(id==unique(data9b$id)[i], month %in% 5:7)
  ntimesb[i] <- nrow(temp)
  rm(temp)
}
modNew <- depmix(list(speed~1, rt100~1, hr.48~1), data=data9b, nstates=K, 
              family=list(Gamma(link="log"), Gamma(link="log"), Gamma(link="log")),
              ntimes=ntimesb)
modNew <- setpars(modNew,getpars(mod))
modNew <- fit(modNew, emcontrol=em.control(maxit=1000, tol=0.01, crit="relative"))


data10b <- cbind(data9b, posterior(modNew, type='viterbi'))


fill=c("3"="grey85", "2"="grey90", "1"="lightpink1")
parturitionstate=1


data11b <- data10b %>% mutate(parturition=ifelse(state==parturitionstate, 1, 0)) %>%
  group_by(id) %>%
  mutate(sumpart.7=rollsum(parturition, 84, align="left", na.pad=T)/84,
         sumpart.1=rollsum(parturition, 12, align="left", na.pad=T)/12,
         sumpart.2=rollsum(parturition, 24, align="left", na.pad=T)/24,
         sumpart.3=rollsum(parturition, 36, align="left", na.pad=T)/36,
         sumpart.4=rollsum(parturition, 48, align="left", na.pad=T)/48) %>%
  ungroup()
write.csv(data11b, "data/BHS-female-NONPARTURIENT_averaged-movement-terrain-hmm-states-viterbi-rollingavg_clean.csv")

nonparturitionsummary <- data11b %>% 
  filter(id!="B18-2022"|date>"2022-06-10") %>% # remove the data from B18 that's all messed up from pre-collaring
  group_by(id) %>%
  summarise(maxsumpart.2=max(sumpart.2, na.rm=T),
            meansumpart.2=mean(sumpart.2, na.rm=T),
            sdsumpart.2=sd(sumpart.2, na.rm=T),
            part.date=date[which.max(sumpart.2)],
            part.x=x[which.max(sumpart.2)],
            part.y=y[which.max(sumpart.2)],
            part.date=as.POSIXct(ifelse(maxsumpart.2>0.5, part.date, NA), origin="1970-01-01"))
nonparturitionsummary <- left_join(nonparturitionsummary, info2, by="id")
write.csv(nonparturitionsummary, "data/BHS-female-NONPARTURIENT_inferred-parturition-sites-dates.csv")


# Plot the viterbi states and estimated parturition dates, if there are any
{pdf("TESTWITHTOL/RESIZED-HMM-nonparturient_viterbi-lambingdate_May15-Jul15.pdf", width=6.69, height=5.8)
  fill=c("3"="white", "2"="grey85", "1"="lightpink1")
  # This is just for assigning the correct colours to each movement state
  for(i in 1:length(nonparturitionsummary$id)){
    temp <- data10b%>%filter(id==nonparturitionsummary$id[i]) %>%
      filter(id!="B18-2022"|date>"2022-06-10") %>% # remove the data from B18 that's all messed up from pre-collaring
      mutate(state=as.factor(state),
             group=cumsum(state!=lag(state, default=first(state))),
             hr.48=hr.48*100, 
             hr.48=ifelse(hr.48>=100, 100, hr.48),
             dist=ifelse(dist>=500, 500, dist))
      
    plotAb <- ggplot(data=temp, aes(date, rt100)) + 
      geom_rect(data=temp, color=NA, aes(ymin=-Inf, ymax=Inf, xmin=date, 
                                         xmax=lead(date), fill=factor(state)), alpha=1, group=temp$state) +
      scale_fill_manual(values=fill) +
      geom_vline(xintercept=nonparturitionsummary$part.date[i], color="red", size=1) +
      geom_line() +
      graph.theme +
      scale_x_datetime(expand=c(0, 0)) +
      theme(legend.position="none") +
      theme(axis.title.x=element_blank()) +
      labs(y="Residence Time (h)", x="Date") +
      theme(strip.text=element_text(face="bold"), strip.background=element_rect(colour=NA))
    
    plotBb <- ggplot(data=temp, aes(date, dist)) + 
      geom_rect(data=temp, color=NA, aes(ymin=-Inf, ymax=Inf, xmin=date, 
                                         xmax=lead(date), fill=factor(state)), alpha=1, group=temp$state) +
      scale_fill_manual(values=fill) +
      geom_vline(xintercept=nonparturitionsummary$part.date[i], color="red", size=1) +
      geom_line() +
      graph.theme +
      scale_x_datetime(expand=c(0, 0)) +
      theme(legend.position="none") +
      theme(axis.title.x=element_blank()) +
      labs(y="Step Length (m)", x="Date") +
      theme(strip.text=element_text(face="bold"), strip.background=element_rect(colour=NA))
    
    plotCb <- ggplot(data=temp, aes(date, hr.48)) + 
      geom_rect(data=temp, color=NA, aes(ymin=-Inf, ymax=Inf, xmin=date, 
                                         xmax=lead(date), fill=factor(state)), alpha=1, group=temp$state) +
      scale_fill_manual(values=fill) +
      geom_vline(xintercept=nonparturitionsummary$part.date[i], color="red", size=1) +
      geom_line() +
      graph.theme +
      scale_x_datetime(expand=c(0, 0)) +
      theme(legend.position="none") +
      theme(axis.title.x=element_blank()) +
      labs(y="Home Range (Ha)", x="Date") +
      theme(strip.text=element_text(face="bold"), strip.background=element_rect(colour=NA))
    
    title <- ggdraw() + 
      draw_label(
        paste(nonparturitionsummary$id[i], ifelse(nonparturitionsummary$confident[i]==1, "(Nonparturient)", "(Uncertain)"), sep=" "),
        fontface = 'bold',
        x = 0,
        hjust =-0.5)
    
    plotABCb <- plot_grid(title, plotAb, NULL, plotBb, NULL, plotCb, xaxis_plotlabel, 
                         ncol=1, align="v", rel_heights=c(0.1, 1, -0.12, 1, -0.12, 1, 0.08))
    print(plotABCb)
    rm(temp, title)
  }
  dev.off()
}

# Plot of viterbi state vs. Pr[low movement state]
{pdf("TESTWITHTOL/HMM-nonparturient_viterbi-state-probability_May15-Jul15.pdf", width=7, height=3)
  
  for(i in 1:length(unique(data10b$id))){
    temp <- data10b%>%filter(id==unique(data10b$id)[i])
    plotDb <- ggplot(data=temp, aes(date, S1)) + 
      geom_rect(data=temp, color=NA, aes(ymin=-Inf, ymax=Inf, xmin=date, 
                                         xmax=lead(date), fill=factor(state)), alpha=0.4) +
      scale_fill_manual(values=fill) +
      geom_line(size=0.2) +
      graph.theme +
      scale_x_datetime(expand=c(0, 0)) +
      theme(legend.position="none") +
      labs(y="Probability of low-movement state", x="Date") +
      theme(strip.text=element_text(face="bold"), strip.background=element_rect(colour=NA))
    
    title <- ggdraw() + 
      draw_label(
        paste(nonparturitionsummary$id[i]),
        fontface = 'bold',
        x = 0,
        hjust =-0.5)
    
    plotDb <- plot_grid(title, plotDb, ncol=1, align="v", rel_heights=c(0.1, 1))
    
    print(plotDb)
    rm(temp, title)
  }
  dev.off()
  rm(middlelevel, upperlevel, lowerlevel)
}


##### CREATE FIGURE OF PARTURITION AGGREGATE MOVEMENT TRENDS -----

confidence_interval <- function(vector, interval) {
  vec_sd <- sd(vector)
  n <- length(vector)
  vec_mean <- mean(vector)
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  return(error)
}

data12 <- data11 %>% 
  filter(id %in% parturitionsummary$id) %>%
  left_join(parturitionsummary %>% select(part.date, id), by="id") %>%
  mutate(timetopart=round(as.numeric(date-part.date)/(60*60))/24,
         withlamb=ifelse(date>=part.date, 1, 0)) %>%
  filter(timetopart>=-15 & timetopart<=15) %>%
  select(id, year, x, y, date, timetopart, hr.48, rt100.48, dist.48)
unique(data12$id)
write.csv(data12%>%filter(timetopart>=-15&timetopart<=15), "maps/mappingdata-updated-with-outliers-15day.csv")

data13 <- data12 %>%
  group_by(timetopart) %>%
  summarise(meandist = mean(dist.48),
            meanhr = mean(hr.48),
            meanrt100 = mean(rt100.48),
            errordist = confidence_interval(dist.48, 0.95),
            errorhr = confidence_interval(hr.48, 0.95),           
            errorrt100 = confidence_interval(rt100.48, 0.95)) %>%
  mutate(timetopart=timetopart-1) %>%
  filter(timetopart >-10 & timetopart <10)

plotF <- ggplot(data=data13, aes(x=timetopart, y=meandist)) + 
  geom_ribbon(aes(ymin=meandist-errordist, ymax=meandist+errordist), fill="grey80") +  
  geom_vline(xintercept=0, color="red", size=0.8) + 
  geom_line(aes(y=meandist), color="black") +
  graph.theme + 
  labs(x=NULL, y="Step Length (m)") +
  scale_x_continuous(breaks=seq(-10,10,2),expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0))
print(plotF)

plotG <- ggplot(data=data13, aes(x=timetopart, y=meanrt100)) + 
  geom_ribbon(aes(ymin=meanrt100-errorrt100, ymax=meanrt100+errorrt100), fill="grey80") +
  geom_vline(xintercept=0, color="red", size=0.8) + 
  geom_line(aes(y=meanrt100), color="black") +
  graph.theme + 
  labs(x=NULL, y="Residence Time (h)") +
  scale_x_continuous(breaks=seq(-10,10,2),expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0))
print(plotG)

plotH <- ggplot(data=data13, aes(x=timetopart, y=meanhr*100)) + 
  geom_ribbon(aes(ymin=ifelse(meanhr*100-errorhr*100>0, meanhr*100-errorhr*100, 0), ymax=meanhr*100+errorhr*100), fill="grey80") +
  geom_vline(xintercept=0, color="red", size=0.8) + 
  geom_line(aes(y=meanhr*100), color="black") +
  graph.theme + 
  labs(x=NULL, y="Home Range (Ha)") +
  scale_x_continuous(breaks=seq(-10,10,2),expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0))
print(plotH)

xaxis_plotlabel2 <- ggdraw() + 
  draw_label("Days Since Parturition",
             hjust = 0.5,
             size=11)
plotFGH <- plot_grid(plotF, NULL, plotG, NULL, plotH, xaxis_plotlabel2, ncol=1, align = "v", rel_heights=c(1, -0.1, 1, -0.1, 1, 0.08))
print(plotFGH)
{pdf("TESTWITHTOL/RESIZED-HMM-aggregate-postpartum-movement-95CI-leftaligned.pdf", width=5.7, height=5.5)
print(plotFGH)
dev.off()}
