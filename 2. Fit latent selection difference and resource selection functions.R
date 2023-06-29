## SET-UP #####
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(raster)
library(sp)
library(sf)
library(amt)
library(MASS)
library(reshape2)
library(lme4)
library(arm)
library(optimx)
library(survival)
library(mclogit) # I used this for mixed-effect conditional logistic regression. But I could use SURVIVAL instead.
library(repr)
library(adehabitatHR)


setwd("C:/_sheep/")
extract <- raster::extract
select <- dplyr::select
filter <- dplyr::filter
set.seed(69420) # For the same random points every time. 
decay.value <- 10


#Read in the information about sheep, create columns as needed
covariatenames <- read.csv("data/covariate-full-names.csv")
info <- read.csv("data/BHS-female_inferred-parturition-sites-dates.csv") %>%
  mutate(part.date=as.POSIXct(part.date))
info2all <- read.csv("data/validated-non-parturient-sheep.csv") %>%
  mutate(id=paste(name, year, sep="-"))
info2 <- info2all %>% filter(confident==1)

#Read in the data created in script 1) 
data12 <- read.csv("data/BHS-female_basic-movement-terrain_clean.csv") %>%
  mutate(date=as.POSIXct(date),
         year=year(date), 
         month=month(date), 
         day=day(date),
         id=gsub("_", "-", id))
# Read in the aligned rasters and put them in a stack
covariates <- stack(list.files(path="C:/_sheep/rasters", pattern="tif", full.name=T))


## MIXED EFFECT LATENT SELECTION DIFFERENCE -----
###  An  LSD  function  contrasts  the  habitat  use  of  two  classes  of  locations.  
###  The  two  classes  of  concern  for  our  analysis  were  locations  prelambing and postlambing
# Need to scale the variables for this to work. 

data13 <- data12 %>% 
  filter(id %in% info$id) %>%
  left_join(info %>% select(part.date, id), by="id") %>%
  mutate(timetopart=as.numeric(date-part.date)/(60*60*24),
         withlamb=ifelse(date>=part.date, 1, 0)) %>%
  filter(timetopart>=-15 & timetopart<=15) %>%
  select(id, year, x, y, date, withlamb, timetopart, dist, rel.angle)
unique(data13$id)

lsd <- extract_covariates(make_track(data13, x, y, id=id, year=year, withlamb=withlamb, crs=sp::CRS("+init=epsg:26911")), covariates)
lsd_unscaled <- lsd
lsd <- lsd_unscaled #Resets lsd to before I scaled the variables. 
# Scale the variables
lsd <- lsd %>% 
  mutate(dbarren = scale(dbarren), descp = scale(descp), dforest = scale(dforest), dgrassforb = scale(dgrassforb), 
         dgrassforbshrub = scale(dgrassforbshrub), dwater = scale(dwater), 
         elevation = scale(elevation), snowdepth = scale(snowdepth), slope = scale(slope), heatload = scale(heatload), 
         vrm = scale(vrm), percentcrown = scale(percentcrown),
         droad = scale(1 - exp(-decay.value * droad * 0.001)),
         dtrail = scale(1-exp(-decay.value * dtrail * 0.001))
  )

mod2 <- glmer(withlamb ~ heatload + descp + dbarren +
              dgrassforbshrub + droad + dtrail + dwater + elevation + 
              percentcrown + snowdepth + vrm + (1|id),
              data=lsd, family=binomial)
summary(mod2)

binnedplot(fitted(mod2), 
           residuals(mod2, type = "response"), 
           nclass = NULL, 
           xlab = "Expected Values", 
           ylab = "Average residual", 
           main = "Binned residual plot", 
           cex.pts = 0.8, 
           col.pts = 1, 
           col.int = "gray")

mod2_summary <- data.frame(summary(mod2)$coefficients)
mod2_summary$variables <- rownames(mod2_summary)
mod2_CI <- confint(mod2,parm="beta_", method="Wald")
colnames(mod2_CI) <- c("CI.min", "CI.max")
mod2_summary <- cbind(mod2_summary, mod2_CI)
mod2_summary <- mod2_summary %>% filter(variables!="(Intercept)") %>%
  mutate(AIC = summary(mod2)$AICtab[1],
         Deviance = summary(mod2)$AICtab[4],
         df.resid = summary(mod2)$AICtab[5]) %>%
  left_join(covariatenames, by="variables")

LSDplot4 <- ggplot(data=mod2_summary, aes(x=variable_full, y=Estimate)) +
  geom_hline(yintercept=0, color="grey") +
  geom_errorbar(aes(ymin=CI.min, ymax=CI.max), width=0) + 
  geom_point(shape=22, size=3, fill="black") +
  xlab("Landscape Covariate") + ylab("Relative Usage") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, h=1), 
        axis.line = element_line(colour = 'black', size = 0.5),
        strip.background = element_rect(fill="antiquewhite", colour="black", size=0.5),
        strip.text = element_text(size=8, colour="black"),
        axis.title.x = element_text(vjust=7))
print(LSDplot4)

{pdf("peerreview/LSD-randomintercept.pdf", width=5, height=4)
  print(LSDplot4)
  dev.off()}
write.csv(mod2_summary, "data/LSD-randomintercept-summary.csv")


## DON'T USE THIS. TOO SMALL A SAMPLE SIZE TO RELIABLY ESTIMATE RANDOM SLOPE PARAMETERS
#mod3 <- glmer(withlamb ~ heatload + descp + dbarren +
#              dgrassforbshrub + droad + dtrail + dwater + elevation + 
#              percentcrown + snowdepth + vrm +
#              (heatload|id) + (descp|id) + (dbarren|id) +
#              (dgrassforbshrub|id) + (droad|id) + (dtrail|id) + (dwater|id) + (elevation|id) +
#              (percentcrown|id) + (snowdepth|id) + (vrm|id),
#           data=lsd, family=binomial, glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
#summary(mod3)
#
#mod3_summary <- data.frame(summary(mod3)$coefficients)
#mod3_summary$variables <- rownames(mod3_summary)
#mod3_CI <- confint(mod3,parm="beta_", method="Wald")
#colnames(mod3_CI) <- c("CI.min", "CI.max")
#mod3_summary <- cbind(mod3_summary, mod3_CI)
#mod3_summary <- mod3_summary %>% filter(variables!="(Intercept)")
#mod3_summary$AIC <- summary(mod3)$AICtab[1]
#mod3_summary$Deviance <- summary(mod3)$AICtab[4]
#mod3_summary$df.resid <- summary(mod3)$AICtab[5]
#
#binnedplot(fitted(mod3), 
#           residuals(mod3, type = "response"), 
#           nclass = NULL, 
#           xlab = "Expected Values", 
#           ylab = "Average residual", 
#           main = "Binned residual plot", 
#           cex.pts = 0.8, 
#           col.pts = 1, 
#           col.int = "gray")
#
#LSDplot5 <- ggplot(data=mod3_summary, aes(x=variables, y=Estimate)) +
#  geom_hline(yintercept=0, color="grey") +
#  geom_errorbar(aes(ymin=CI.min, ymax=CI.max), width=0) + 
#  geom_point(shape=21, size=3, fill="black") +
#  xlab("Variable") + ylab("Selection Coefficient (95% CI)") +
#  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#        panel.background = element_blank(),
#        axis.text.x = element_text(angle = 90), 
#        axis.line = element_line(colour = 'black', size = 0.5),
#        strip.background = element_rect(fill="antiquewhite", colour="black", size=0.5),
#        strip.text = element_text(size=8, colour="black"))
#print(LSDplot5)
#
#{pdf("C:/sheep_redo/LSDmixedeffect-plots-randomintercept+slope.pdf", width=5, height=5)
#  print(LSDplot5)
#  dev.off()}
#write.csv(model3_summary, "C:/sheep_redo/LSDmixedeffect-summary-scaled-randomintercept+slope.csv")



## RESOURCE SELECTION FUNCTION SEPARATE FOR EACH REPRO STATUS -----
rsf_part <- data14 %>% filter(parturient==1) %>% select(id, x_="x", y_="y", parturient) %>% mutate(case_=TRUE)
id_part <- unique(rsf_part$id)
hr_part <- hr_mcp(make_track(data14, x, y, parturient=parturient, id=id, crs=sp::CRS("+init=epsg:26911")), levels=0.99)
rsf_part <- rbind(rsf_part,
                  random_points(hr_part, n=15*nrow(rsf_part)) %>% mutate(parturient="parturient", id=NA))
rsf_part$id[is.na(rsf_part$id)]=sample(id_part, length(is.na(rsf_part$id)), replace=T)
rsf_part <- extract_covariates(make_track(rsf_part, x_, y_, id=id, parturient=parturient, case_=case_, crs=sp::CRS("+init=epsg:26911")), covariates)
rsf_part_unscaled <- rsf_part
rsf_part <- rsf_part_unscaled #Resets rsf to before I scaled the variables. 
rsf_part <- rsf_part %>% 
  mutate(dbarren = scale(dbarren), descp = scale(descp), dforest = scale(dforest), dgrassforb = scale(dgrassforb), 
         dgrassforbshrub = scale(dgrassforbshrub), dwater = scale(dwater), 
         elevation = scale(elevation), snowdepth = scale(snowdepth), slope = scale(slope), heatload = scale(heatload), 
         vrm = scale(vrm), percentcrown = scale(percentcrown),
         droad = scale(1 - exp(-decay.value * droad * 0.001)),
         dtrail = scale(1-exp(-decay.value * dtrail * 0.001)))

rsf_nonpart <- data14 %>% filter(parturient==0) %>% select(id, x_="x", y_="y", parturient) %>% mutate(case_=TRUE)
id_nonpart <- unique(rsf_nonpart$id)
hr_nonpart <- hr_mcp(make_track(data14, x, y, parturient=parturient, id=id, crs=sp::CRS("+init=epsg:26911")), levels=0.99)
rsf_nonpart <- rbind(rsf_nonpart,
                  random_points(hr_nonpart, n=15*nrow(rsf_nonpart)) %>% mutate(parturient="nonparturient", id=NA))
rsf_nonpart$id[is.na(rsf_nonpart$id)]=sample(id_nonpart, length(is.na(rsf_nonpart$id)), replace=T)
rsf_nonpart <- extract_covariates(make_track(rsf_nonpart, x_, y_, id=id, parturient=parturient, case_=case_, crs=sp::CRS("+init=epsg:26911")), covariates)
rsf_nonpart_unscaled <- rsf_nonpart
rsf_nonpart <- rsf_nonpart_unscaled #Resets rsf to before I scaled the variables. 
rsf_nonpart <- rsf_nonpart %>%
  mutate(dbarren = scale(dbarren), descp = scale(descp), dforest = scale(dforest), dgrassforb = scale(dgrassforb), 
         dgrassforbshrub = scale(dgrassforbshrub), dwater = scale(dwater), 
         elevation = scale(elevation), snowdepth = scale(snowdepth), slope = scale(slope), heatload = scale(heatload), 
         vrm = scale(vrm), percentcrown = scale(percentcrown),
         droad = scale(1 - exp(-decay.value * droad * 0.001)),
         dtrail = scale(1-exp(-decay.value * dtrail * 0.001)))

id_nonpart
id_part

mod9 <- glmer(case_ ~ heatload + descp + dbarren +
                dgrassforbshrub + droad + dtrail + dwater + elevation + 
                percentcrown + snowdepth + vrm + (1|id),
              data=rsf_part, family=binomial)
mod10 <- glmer(case_ ~ heatload + descp + dbarren +
                dgrassforbshrub + droad + dtrail + dwater + elevation + 
                percentcrown + snowdepth + vrm + (1|id),
              data=rsf_nonpart, family=binomial)
summary(mod9)
summary(mod10)

mod9_summary <- data.frame(summary(mod9)$coefficients) %>%
  mutate(variables=row.names(summary(mod9)$coefficients),
         AIC=AIC(mod9),
         n=summary(mod9)$N,
         CI.min=confint(mod9,parm="beta_", method="Wald")[,1],
         CI.max=confint(mod9,parm="beta_", method="Wald")[,2],
         parturient="parturient")
mod10_summary <- data.frame(summary(mod10)$coefficients) %>%
  mutate(variables=row.names(summary(mod10)$coefficients),
         AIC=AIC(mod10),
         n=summary(mod10)$N,
         CI.min=confint(mod10,parm="beta_", method="Wald")[,1],
         CI.max=confint(mod10,parm="beta_", method="Wald")[,2],
         parturient="nonparturient")
mod910_summary <- rbind(mod9_summary, mod10_summary) %>% filter(variables!="(Intercept)") %>%
  left_join(covariatenames, by="variables")


RSFplot1 <- ggplot(data=mod910_summary, aes(x=variable_full, y=Estimate, fill=factor(parturient))) +
  geom_hline(yintercept=0, color="grey") +
  geom_errorbar(aes(x=variable_full, ymin=CI.min, ymax=CI.max), width=0, position=position_dodge(width=0.5)) + 
  geom_point(shape=22, size=3, position=position_dodge(width=0.5)) +
  xlab("Landscape Covariate") + ylab("Selection Coefficient") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45,  
                                   hjust = 1), 
        axis.line = element_line(colour = 'black', size = 0.5),
        strip.background = element_rect(fill="antiquewhite", colour="black", size=0.5),
        strip.text = element_text(size=8, colour="black"),
        axis.title.x = element_text(vjust=8)) +
  scale_fill_manual(values = c("white", "black"), labels=c("Non-parturient  ", "Parturient")) +
  theme(legend.background = element_rect(fill="grey90"),
        legend.title=element_blank(),
        legend.position = c(0.89, 0.9),
        legend.spacing.y = unit(0, "pt"),
        legend.spacing.x = unit(0, "pt"),
        legend.text=element_text(size=8),
        legend.key=element_rect(fill=NA),
        legend.margin=margin(t=1))

print(RSFplot1)

{pdf("figures/RSF-separate-repro-status-random-intercept.pdf", width=6, height=4.5)
  print(RSFplot1)
  dev.off()}
write.csv(mod910_summary, "data/RSF-separate-repro-status-random-intercept-summary.csv")



### ADDITIONAL PROCESSING, FILES FOR MAPS, AND METADATA
ntimes <- lsd %>% group_by(id) %>% summarise(n_LSD=n())
sheep_total <- data12 %>% 
  filter(id %in% info$id | id %in% info2all$id) %>%
  mutate(parturient=ifelse(id %in% info$id, 1, 0)) %>%
  select(id, year, x, y, date, year, month, day, parturient, dist, rel.angle) %>%
  filter(month==5&day>=15 | month==6 | month==7&day<=15) %>%
  group_by(id) %>%
  dplyr::summarise(start=date(min(date)), end=date(max(date)), year=max(year), n_HMM_RSF_SSF=n()) %>%
  mutate(status=ifelse(id %in% info$id, "Parturient", NA),
         status=ifelse(id %in% info2$id, "Non-Parturient", status),
         status=ifelse(is.na(status), "Uncertain", status)) %>%
  left_join(info%>%select(id, part.date, part.x, part.y), by="id") %>%
  left_join(ntimes, by="id")
write.csv(sheep_total, "sheep-sample-summary.csv")



hr <- hr_mcp(make_track(data14, x, y, parturient=parturient, id=id, crs=sp::CRS("+init=epsg:26911")), levels=0.95)
hr <- hr$mcp
hr <- hr$geometry
hr <- hr$`95%`[[1]]
write.csv(hr, "maps/95pct-hr.csv")









