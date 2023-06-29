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




## STEP SELECTION FUNCTION SEPARATE FOR EACH REPRO STATUS -----
data14 <- data12 %>% 
  filter(id %in% info$id | id %in% info2$id) %>%
  mutate(parturient=ifelse(id %in% info$id, 1, 0)) %>%
  select(id, year, x, y, date, year, month, day, parturient, dist, rel.angle) %>%
  filter(month==5&day>=15 | month==6 | month==7&day<=15)
unique(data14$id)

# Extract covariates for each point
dataSSF <- make_track(data14, x, y, date, parturient=parturient, year=year, id=id, crs=sp::CRS("+init=epsg:26911")) %>%
  nest(data=-"id") %>%
  mutate(steps=map(data, function(temp)
    temp %>% track_resample(rate=hours(2), tolerance=minutes(30)) %>% steps_by_burst %>% random_steps(n_control=15))) %>%
  select(id, steps) %>% unnest(cols=steps) %>%
  mutate(parturient=ifelse(id %in% info$id, 1, 0))
ssf <- dataSSF %>%
  extract_covariates(covariates, where="end")
ssf_unscaled <- ssf
ssf <- ssf_unscaled #Resets ssf data to before I scaled the variables. 
ssf <- ssf %>% 
  mutate(dbarren = scale(dbarren), descp = scale(descp), dforest = scale(dforest), dgrassforb = scale(dgrassforb), 
         dgrassforbshrub = scale(dgrassforbshrub), dwater = scale(dwater), 
         elevation = scale(elevation), snowdepth = scale(snowdepth), slope = scale(slope), heatload = scale(heatload), 
         vrm = scale(vrm), percentcrown = scale(percentcrown),
         cos_ta_ = cos(ta_), 
         log_sl_ = log(sl_),
         droad = scale(1 - exp(-decay.value * droad * 0.001)),
         dtrail = scale(1-exp(-decay.value * dtrail * 0.001))
         ) %>%
  mutate(case_name = paste(id, step_id_, sep="_"),
         case_name=as.numeric(factor(case_name)))

ssf_part <- ssf %>% filter(parturient==1)
ssf_nonpart <- ssf %>% filter(parturient==0)

mod7 <- mclogit(cbind(case_, case_name) ~ heatload + descp + dbarren +
                  dgrassforbshrub + droad + dtrail + dwater + elevation + 
                  percentcrown + snowdepth + vrm,
                data=ssf_part)

mod8 <- mclogit(cbind(case_, case_name) ~ heatload + descp + dbarren +
                  dgrassforbshrub + droad + dtrail + dwater + elevation + 
                  percentcrown + snowdepth + vrm,
                data=ssf_nonpart)
summary(mod7)
mod7_summary <- data.frame(summary(mod7)$coefficients) %>%
  mutate(variables=row.names(summary(mod7)$coefficients),
         AIC=AIC(mod7),
         n=summary(mod7)$N,
         CI.min=confint(mod7)[,1],
         CI.max=confint(mod7)[,2],
         parturient="parturient")
mod8_summary <- data.frame(summary(mod8)$coefficients) %>%
  mutate(variables=row.names(summary(mod8)$coefficients),
         AIC=AIC(mod8),
         n=summary(mod8)$N,
         CI.min=confint(mod8)[,1],
         CI.max=confint(mod8)[,2],
         parturient="nonparturient")
mod78_summary <- rbind(mod7_summary, mod8_summary) %>%
  left_join(covariatenames, by="variables")


SSFplot2 <- ggplot(data=mod78_summary, aes(x=variable_full, y=Estimate, fill=factor(parturient))) +
  geom_hline(yintercept=0, color="grey40") +
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
print(SSFplot2)

{pdf("figures/SSF-separate-repro-status.pdf", width=6, height=4.5)
  print(SSFplot2)
  dev.off()}
write.csv(mod78_summary, "data/SSF-separate-repro-status-summary.csv")


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







################
### OLD CODE ###
################
## LATENT SELECTION DIFFERENCE #####
#Fit the models and calculate confidence intervals. 
ids <- unique(lsd$id)
lsd_CI <- list()
lsd_summary <- list()
for(i in 1:length(ids)) {
  print(paste(i, ids[i], sep=" "))
  temp <- lsd %>% 
    filter(id == ids[i]) %>%
    select(id, withlamb, heatload, descp, dbarren, dgrassforbshrub, droad, 
           dtrail, dwater, elevation, percentcrown, snowdepth, vrm) %>%
    drop_na()
  mod <- glm(withlamb ~ heatload + descp + dbarren +
               dgrassforbshrub + droad + dtrail + dwater + elevation + 
               percentcrown + snowdepth + vrm, 
             data=temp, family=binomial)
  summary(mod)
  lsd_CI[[i]] = confint(mod, level=0.95)
  lsd_summary[[i]] = summary(mod)
}

#Store model outputs, CIs in a data frame
LSDspreadsheet <- data.frame()
for(i in 1:length(lsd_summary)) {
  int <- data.frame(lsd_CI[[i]], row.names=NULL)
  output <- data.frame(lsd_summary[[i]]$coefficients) %>%
    mutate(variable = row.names(output), 
           id = ids[i], 
           AIC = lsd_summary[[i]]$aic, 
           n = lsd_summary[[i]]$df.null + 1, 
           df.res = lsd_summary[[i]]$df.residual, 
           df.null = lsd_summary[[i]]$df.null,
           CI.min = int$X2.5..[!is.na(int$X97.5..) | !is.na(int$X2.5..)],
           CI.max = int$X97.5..[!is.na(int$X97.5..) | !is.na(int$X2.5..)]
    )
  LSDspreadsheet <- rbind(LSDspreadsheet, output)
}
row.names(LSDspreadsheet) <- NULL

# Create columns with the mean before/after values for each covariate. 
lsd_melted <- melt(lsd_unscaled, id.vars=c("id", "withlamb"), 
                   measure.vars=unique(LSDspreadsheet$variable)[2:length(unique(LSDspreadsheet$variable))])
for (i in 1:nrow(LSDspreadsheet)){
  before <- lsd_melted %>% filter(lsd_melted$variable==LSDspreadsheet$variable[i] & id==LSDspreadsheet$id[i] & withlamb==0)
  after <- lsd_melted %>% filter(lsd_melted$variable==LSDspreadsheet$variable[i] & id==LSDspreadsheet$id[i] & withlamb==1)
  LSDspreadsheet$meanbefore[i] <- mean(before$value)
  LSDspreadsheet$meanafter[i] <- mean(after$value)
}

LSDspreadsheet <- LSDspreadsheet %>%
  mutate(meanbefore = as.numeric(meanbefore), 
         meanafter = as.numeric(meanafter), 
         mean = (meanbefore + meanafter)/2)


# Figures #
# Filter out the facets I don't want to plot, tweak the id names
figuredata <- LSDspreadsheet %>% filter(variable != "(Intercept)")
variablenames_full <- c("Heat load", "Distance to escape terrain", "Distance to barren ground", "Distance to meadow", 
                        "Distance to road", "Distance to trail", "Distance to water", "Elevation", "Crown cover", 
                        "Snow depth", "Ruggedness")
variablenames <- data.frame(cbind(unique(figuredata$variable), variablenames_full))
colnames(variablenames) <- c("variable", "variablenames_full")
figuredata <- inner_join(figuredata, variablenames, by="variable")



# Give the sheep a new ID #
sheepidconversions <- cbind(sort(unique(SSF_plotdata$ordered.id)), paste("S", 1:length(unique(SSF_plotdata$id)), sep=""))
colnames(sheepidconversions) <- c("orderedID", "referencenumber")
sheepidconversions <- data.frame(sheepidconversions)
sheepidconversions$id <- substr(sheepidconversions$orderedID, 2, 10000)
write.csv(sheepidconversions, "C:/sheep_redo/sheepindexnumbers-fromoriginalid-usedtomakefigures.csv")
figuredata <- inner_join(figuredata, sheepidconversions, by="id")


# Make the plot of selection coefficients vs. id (facet=) #
LSDplot <- ggplot(figuredata%>%filter(id!="B51-2022"), aes(x=id, y=Estimate)) + 
  geom_hline(yintercept=0, color="grey") +
  geom_errorbar(color="black", aes(ymin=CI.min, ymax=CI.max), width=0) + 
  geom_point(shape=21, size=2, fill="black", color="black") + 
  xlab("Sheep ID") + ylab("Selection Coefficient") +
  facet_wrap(. ~ variablenames_full, scales="free", ncol=4) +
  theme_classic() + 
  theme(strip.background = element_blank()) +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    axis.text.x = element_text(angle = 90), 
    #axis.line = element_line(colour = 'black', size = 0.5),
    strip.text = element_text(size=8, colour="black", face=2)) 
print(LSDplot)





















# Make the plot of selection coefficients vs. availability
figuredata.2 <- LSDspreadsheet %>% filter(variable != "(Intercept)") # Filter out some facets
LSDplot2 <- ggplot(figuredata.2, aes(x=mean, y=Estimate)) + 
  geom_hline(yintercept=0, color="black") +
  #  geom_errorbar(aes(ymin=CI.min, ymax=CI.max), width=0) + 
  geom_point(shape=21, size=1.5, fill="blue") + 
  xlab("Availability") + ylab("Selection Coefficient") +
  facet_wrap(. ~ variable, scales="free", ncol=5) +
  theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x = element_text(angle = 90), 
        axis.line = element_line(colour = 'black', size = 0.5),
        strip.background = element_rect(fill="antiquewhite", colour="black", size=0.5),
        strip.text = element_text(size=8, colour="black"))
print(LSDplot2)

# Make the plot of selection coefficients vs. variables (facet=id)
LSDplot3 <- ggplot(figuredata, aes(x=variable, y=Estimate)) + 
  geom_hline(yintercept=0, color="grey") +
  #  geom_errorbar(aes(ymin=CI.min, ymax=CI.max), width=.2) + 
  geom_point(shape=21, size=1, fill="black") + 
  xlab("Variable") + ylab("Selection Coefficient") +
  facet_wrap(. ~ id, scales="free", ncol=4) +
  theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x = element_text(angle = 90), 
        axis.line = element_line(colour = 'black', size = 0.5),
        strip.background = element_rect(fill="antiquewhite", colour="black", size=0.5),
        strip.text = element_text(size=8, colour="black"))
print(LSDplot3)

{pdf("C:/_sheep/LSD-plots-one-model-per-sheep.pdf", width=10, height=6.3)
  print(LSDplot)
  print(LSDplot2)
  print(LSDplot3)
  dev.off()}
write.csv(LSDspreadsheet, "C:/_sheep/LSD-summary-scaled-one-model-per-sheep.csv")




## STEP SELECTION FUNCTION WITH PARTURITION INTERACTION TERM ----


# Interaction model (difference in selection between repro states)
mod4 <- mclogit(cbind(case_, case_name) ~ heatload*parturient + descp*parturient + dbarren*parturient +
                  dgrassforbshrub*parturient + droad*parturient + dtrail*parturient + dwater*parturient + elevation*parturient + 
                  percentcrown*parturient + snowdepth*parturient + vrm*parturient,
                data=ssf)
# Null model (no difference between repro states)
mod5 <- mclogit(cbind(case_, case_name) ~ heatload + descp + dbarren +
                  dgrassforbshrub + droad + dtrail + dwater + elevation + 
                  percentcrown + snowdepth + vrm,
                data=ssf)
# This combination of interaction terms yields the lowest AIC score. 
mod6 <- mclogit(cbind(case_, case_name) ~ heatload*parturient + descp + dbarren*parturient +
                  dgrassforbshrub + droad*parturient + dtrail*parturient + dwater + elevation*parturient + 
                  percentcrown*parturient + snowdepth*parturient + vrm*parturient,
                data=ssf)
summary(mod4)
summary(mod5)
summary(mod6)
AIC(mod4, mod5, mod6)
#     df      AIC
#mod4 22 63392.31
#mod5 11 63476.89
#mod6 19 63387.26


mod4_summary <- data.frame(summary(mod4)$coefficients) %>%
  mutate(variables=row.names(mod4_summary),
         AIC=AIC(mod4),
         n=summary(mod4)$N,
         CI.min=confint(mod4)[,1],
         CI.max=confint(mod4)[,2])

mod4_summaryX <- mod4_summary[1:11,] %>%
  mutate(parturient="nonparturient")
mod4_summaryY <- mod4_summary[12:22,] %>%
  mutate(parturient="parturient",
         Estimate=Estimate+mod4_summaryX$Estimate,
         #       Std..Error=sqrt(Std..Error^2+mod4_summaryX$Std..Error^2),
         CI.min = Estimate - qnorm(0.975)*Std..Error,
         CI.max = Estimate + qnorm(0.975)*Std..Error,
         variables=mod4_summaryX$variables)
mod4_summary_figures <- rbind(mod4_summaryX, mod4_summaryY)
rm(mod4_summaryX, mod4_summaryY)

SSFplot1 <- ggplot(data=mod4_summary_figures, aes(x=variables, y=Estimate, fill=factor(parturient))) +
  geom_hline(yintercept=0, color="grey") +
  geom_errorbar(aes(x=variables, ymin=CI.min, ymax=CI.max), width=0, position=position_dodge(width=0.5)) + 
  geom_point(shape=21, size=3, position=position_dodge(width=0.5)) +
  xlab("Variable") + ylab("Selection Coefficient (95% CI)") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90,  
                                   hjust = 1), 
        axis.line = element_line(colour = 'black', size = 0.5),
        strip.background = element_rect(fill="antiquewhite", colour="black", size=0.5),
        strip.text = element_text(size=8, colour="black")) +
  scale_fill_manual(values = c("black", "white"))
print(SSFplot1)


{pdf("SSF-all-sheep_coefficients_with-parturient-interaction.pdf", width=6, height=4)
  print(SSFplot1)
  dev.off()}
write.csv(mod4_summary, "SSF-all-sheep_summary_with-parturient-interaction.csv")



###### EVERYTHING UNDER HERE IS UNEDITED CODE ######
## RESOURCE SELECTION FUNCTION (SINGLE MODEL WITH INTERACTION) -----



#Sample random points from the pooled HR iteratively. 
#Assign indicator for parturition to available sites. 
sheeptrack <- make_track(data2, x, y, status=status, parturient=parturient, year=year, id=id, crs=sp::CRS("+init=epsg:26911"))
sheephr <- hr_mcp(sheeptrack, levels=0.99)
plot(sheephr)

rsfpoints <- data.frame()
for (i in unique(sheeptrack$id)){
  int <- sheeptrack %>% filter(id == i) %>%
    mutate(case_ = TRUE) %>%
    select(case_, x_, y_, parturient, year, id)
  int2 <- random_points(sheephr, n=10*nrow(int))
  int2 <- int2 %>%
    mutate(parturient = int$parturient[1], year = int$year[1], id = int$id[1]) %>%
    select(case_, x_, y_, parturient, year, id)
  
  rsfpoints <- rbind(rsfpoints, int, int2)
}

sheeprsf <- extract_covariates(rsfpoints, covariates)
sheeprsf_backup <- sheeprsf #Store the un-scaled version
sheeprsf <- sheeprsf_backup #Convert back to un-scaled version

{ #Scale the variables
  attach(sheeprsf)
  sheeprsf$dbarren <- scale(dbarren)
  sheeprsf$descp <- scale(descp)
  sheeprsf$dforest <- scale(dforest)
  sheeprsf$dgrassforb <- scale(dgrassforb)
  sheeprsf$dgrassforbshrub <- scale(dgrassforbshrub)
  sheeprsf$droad <- scale(droad)
  sheeprsf$dtrail <- scale(dtrail)
  sheeprsf$dwater <- scale(dwater)
  sheeprsf$elevation <- scale(elevation)
  detach(sheeprsf)
}

rsf1 <- glm(case_ ~ aspecte*parturient + aspectn*parturient + aspectw*parturient + heatload*parturient + 
              descp*parturient + dgrassforbshrub*parturient + droad*parturient + dtrail*parturient + dwater*parturient + 
              elevation*parturient + percentcrown*parturient + snowdepth*parturient + vrm*parturient, 
            data=sheeprsf, family=binomial(link='logit'))

summary(rsf1)

rsf1coefficients <- data.frame(summary(rsf1)$coefficients)
rsf1coefficients$variable <- row.names(rsf1coefficients)
rownames(rsf1coefficients) <- NULL
write.csv(rsf1coefficients, "C:/sheep_redo/RSF-allsheep-indicator-fullbreedingseason-summary.csv")

rsf1coefficients <- rsf1coefficients %>% filter(variable != c("(Intercept)", "parturient")) # might need to run this twice
rsf1plot <- data.frame(matrix(nrow=13, ncol=6))
colnames(rsf1plot) <- c("variable", "est", "se.est", "interaction", "se.interaction", "sum.est")
rsf1plot$variable <- rsf1coefficients$variable[1:13]
rsf1plot$est <- rsf1coefficients$Estimate[1:13]
rsf1plot$se.est <- rsf1coefficients$Std..Error[1:13]
rsf1plot$interaction <- rsf1coefficients$Estimate[14:26]
rsf1plot$se.interaction <- rsf1coefficients$Std..Error[14:26]
rsf1plot$sum.est <- rsf1plot$est + rsf1plot$interaction

rsfplot1 <- ggplot(rsf1plot, aes(x=variable)) + 
  geom_hline(yintercept=0, color="grey") +
  geom_point(aes(y=est), shape=21, size=3, fill="blue") +
  geom_point(aes(y=sum.est), shape=21, size=3, fill="red") +
  xlab("Variable") + ylab("Selection Coefficient") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90), 
        axis.line = element_line(colour = 'black', size = 0.5),
        strip.background = element_rect(fill="antiquewhite", colour="black", size=0.5),
        strip.text = element_text(size=8, colour="black"))

pdf("C:/sheep_redo/RSFplot-60daybreedingseasonwindow.pdf", width=5, height=3)
print(rsfplot1)
dev.off()


.














## RESOURCE SELECTION FUNCTION (SEPARATE FOR EACH INDIVIDUAL) -----

#Filter to May 15-Jul 15
data2 <- data1 %>%
  filter(id %in% info$id) %>%
  inner_join(info, by='id') %>%
  mutate(year=year.x, individual=individual.x, parturient=Parturient, 
         diff=as.numeric(date-midseasondate)/(60*60*24),
         status=ifelse(parturient==1 & diff>=0, 1, 0)) %>%
  filter(diff>=-30 & diff<=30) %>% # Filters to 60 day window around breeding season.
  select(individual, year, id, x, y, date, parturient, status, diff, dist, rel.angle)

# Extract covariates for each point
# EXPLANATION: https://conservancy.umn.edu/bitstream/handle/11299/218272/AppB_SSF_examples.html?sequence=26
RSFindiv_summary <- data.frame()

for(i in unique(data2$id)){
  print(i)
  data3 <- data2 %>% filter(id==i)
  sheeptrack <- make_track(data3, x, y, date, status=status, parturient=parturient, year=year, id=id, crs=sp::CRS("+init=epsg:26911"))
  print("extracting covariates")
  sheepRSFindiv <- sheeptrack %>% 
    random_points(n=10*nrow(sheeptrack), hr="mcp", presence=sheeptrack) %>%
    extract_covariates(covariates) %>%
    mutate(dbarren = scale(dbarren),
           descp = scale(descp),
           dforest = scale(dforest),
           dgrassforb = scale(dgrassforb),
           dgrassforbshrub = scale(dgrassforbshrub),
           droad = scale(droad),
           dtrail = scale(dtrail),
           dwater = scale(dwater),
           elevation = scale(elevation),
           slope = slope/100,
    )
  
  sheepRSFindiv$case_ <- as.logical(ifelse(sheepRSFindiv$case_==TRUE, 1, 0))
  print("Fitting Model")
  RSFindiv1 <- glm(case_ ~ aspecte + aspects + aspectw + 
                     heatload + descp + dbarren +
                     dgrassforbshrub + droad + dtrail + dwater + elevation + 
                     percentcrown + snowdepth + vrm,
                   #sl_ + log_sl_ + cos_ta_ 
                   family=binomial(link="logit"), data=sheepRSFindiv
  )
  
  LSDplot
  # THINGS TO FIX: 
  # CHECK DOCUMENTATION FOR FIT_RSF
  # MAKE SURE THE MODEL IS CORRECT. 
  # I NEED TO CHANGE THE SAMPLING REGIME TO BE NOT BASED ON STEP DISTRIBUTION ANYMORE
  # CONFIDENCE INTERVALS FOR RSF (PRESENTED DIFFERENTLY, NEED TO UPDATE CODE TO REFLECT)
  int <- summary(RSFindiv1)
  
  modeloutput <- data.frame(cbind(int$coefficients, confint(RSFindiv1)))
  modeloutput <- modeloutput %>%
    mutate(id = i, status = data3$parturient[1], year=data3$year[1],
           status = ifelse(status==1, "Parturient", "Non-Parturient"),
           variable = rownames(modeloutput))
  modeloutput$variable <- str_remove_all(modeloutput$variable, "_end")
  rownames(modeloutput) <- NULL
  RSFindiv_summary <- rbind(RSFindiv_summary, modeloutput)
}

print("hello world")


RSFindiv_plotdata <- RSFindiv_summary[order(RSFindiv_summary$status),]
RSFindiv_plotdata$index <- ifelse(RSFindiv_plotdata$status=="Non-Parturient", "a", "b")
RSFindiv_plotdata$ordered.id <- paste(RSFindiv_plotdata$index, RSFindiv_plotdata$id, sep="")
var.remove <- c("sl_", "log_sl_", "cos_ta_")
RSFindiv_plotdata <- RSFindiv_plotdata %>% filter(! variable %in% var.remove)


RSFindivplot <- ggplot(RSFindiv_plotdata, aes(x=ordered.id, y=Estimate, fill=as.factor(status))) + 
  geom_hline(yintercept=0, color="grey") + 
  geom_errorbar(aes(ymin=X2.5.., ymax=X97.5.., color=as.factor(status)), width=0) + 
  geom_point(aes(color=as.factor(status)), shape=21, size=3) +
  xlab("Sheep ID") + ylab("Selection Coefficient") +
  facet_wrap(. ~ variable, scales="free", ncol=5) +
  theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x = element_text(angle = 90, size=7), 
        axis.line = element_line(colour = 'black', size = 0.5),
        strip.background = element_rect(fill="grey90", colour="black", size=0.5),
        strip.text = element_text(size=8, colour="black")) +
  scale_x_discrete(labels=substr(unique(RSFindiv_plotdata$ordered.id), 2, 10000))

RSFindivplot

pdf("C:/sheep_redo/RSF-onindividuals-lambingseason.pdf", width=12, height=7)
print(RSFindivplot)
dev.off()
write.csv(RSFindiv_summary, "C:/sheep_redo/RSF-onindividuals-lambingseason.pdf")





## RESOURCE SELECTION FUNCTION (SEPARATE FOR REPRODUCTIVE STATUS) -----

data2 <- data1 %>%
  filter(id %in% info$id) %>%
  inner_join(info, by='id') %>%
  mutate(year=year.x, individual=individual.x, parturient=Parturient, 
         diff=as.numeric(date-midseasondate)/(60*60*24),
         status=ifelse(parturient==1 & diff>=0, 1, 0)) %>%
  filter(diff>=-30 & diff<=30) %>% # Filters to 60 day window around breeding season.
  select(individual, year, id, x, y, date, parturient, status, diff, dist, rel.angle)

## RSF ON PARTURIENT == 1
data3 <- data2 %>% filter(parturient==1)
sheeptrack <- make_track(data3, x, y, status=status, parturient=parturient, year=year, id=id, crs=sp::CRS("+init=epsg:26911"))
sheephr <- hr_mcp(sheeptrack, levels=0.99)
plot(sheephr)

int <- sheeptrack %>%
  mutate(case_ = TRUE) %>%
  select(case_, x_, y_, parturient, year, id)
int2 <- random_points(sheephr, n=10*nrow(sheeptrack)) %>%
  mutate(parturient = int$parturient[1], year = int$year[1], id = int$id[1]) %>%
  select(case_, x_, y_, parturient, year, id)
dataRSF1 <- rbind(int, int2)
RSF1points <- extract_covariates(dataRSF1, covariates)

rsf1 <- glm(case_ ~ aspecte + aspectn + aspectw + heatload + descp + dbarren +
              dgrassforbshrub + droad + dtrail + dwater + elevation + 
              percentcrown + snowdepth + vrm, 
            data=RSF1points, family=binomial)
summary(rsf1)
rsf1summary <- data.frame(summary(rsf1)$coefficients)
rsf1summary$status <- "Parturient"
rsf1confint <- confint(rsf1)
colnames(rsf1confint) <- c("CI.min", "CI.max")
rsf1summary$variable <- rownames(rsf1summary)
rsf1summary <- cbind(rsf1summary, rsf1confint)


## RSF ON PARTURIENT == 0

data3 <- data2 %>% filter(parturient==0)
sheeptrack <- make_track(data3, x, y, status=status, parturient=parturient, year=year, id=id, crs=sp::CRS("+init=epsg:26911"))
sheephr <- hr_mcp(sheeptrack, levels=0.99)
plot(sheephr)

int <- sheeptrack %>%
  mutate(case_ = TRUE) %>%
  select(case_, x_, y_, parturient, year, id)
int2 <- random_points(sheephr, n=10*nrow(data3)) %>%
  mutate(parturient = int$parturient[1], year = int$year[1], id = int$id[1]) %>%
  select(case_, x_, y_, parturient, year, id)
dataRSF2 <- rbind(int, int2)
RSF2points <- extract_covariates(dataRSF2, covariates)

rsf2 <- glm(case_ ~ aspecte + aspectn + aspectw + heatload + descp + dbarren +
              dgrassforbshrub + droad + dtrail + dwater + elevation + 
              percentcrown + snowdepth + vrm, 
            data=RSF2points, family=binomial)

summary(rsf2)
rsf2summary <- data.frame(summary(rsf2)$coefficients)
rsf2summary$status <- "Non-Parturient"
rsf2confint <- confint(rsf2)
colnames(rsf2confint) <- c("CI.min", "CI.max")
rsf2summary$variable <- rownames(rsf2summary)
rsf2summary <- cbind(rsf2summary, rsf2confint)


## COMBINE DATA, MAKE A PLOT COMPARING CIs of COEFFICIENTS

rsfsummary_both <- rbind(rsf1summary, rsf2summary)

rsfplotdata <- rsfsummary_both %>% filter(variable != "(Intercept)")
rsfplot <- ggplot(data=rsfplotdata, aes(x=variable, y=Estimate, fill=status)) + 
  geom_hline(yintercept=0, color="grey") +
  geom_errorbar(aes(ymin=CI.min, ymax=CI.max, color=factor(status)), width=0, position=position_dodge(width=0.6)) + 
  geom_point(aes(color=factor(status)), shape=21, size=3, position=position_dodge(width=0.6)) +  
  xlab("Landscape Variable") + ylab("Selection Coefficient (95% CI)") +
  theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x = element_text(angle = 90), 
        axis.line = element_line(colour = 'black', size = 0.5),
        strip.background = element_rect(fill="antiquewhite", colour="black", size=0.5),
        strip.text = element_text(size=8, colour="black"))

print(rsfplot)


{pdf("C:/sheep_redo/RSFplot-separatebyreproductivestatus.pdf", width=7, height=5)
  print(rsfplot)
  dev.off()}
write.csv(rsfsummary_both, "C:/sheep_redo/RSFsummary-separatebyreproductivestatus.pdf")




## RESOURCE SELECTION FUNCTION (OLD CODE SNIPPETS) -----

# Make a track with pooled individuals
sheeptrack <- make_track(data2, x, y, status=status, parturient=parturient, year=year, id=id, crs=sp::CRS("+init=epsg:26911"))
# Kernel density home range estimate (99%)
# sheephr <- hr_kde(sheeptrack, h=hr_kde_ref(sheeptrack), trast=make_trast(sheeptrack), levels=0.99, keep.data=T, rand_buffer=1e-05)
# Assign random points within the home range
sheephr <- hr_mcp(sheeptrack, levels=0.95)
plot(sheephr)
sheeppoints <- random_points(sheephr, n=10*nrow(sheeptrack))
# OR... SAMPLE USING THE TRACK
#sheeppoints <- random_points(sheeptrack, hr="mcp", n=10*nrow(sheeptrack))

# Code the indicators/binary variables, combine the used/available. 
sheeptrack <- sheeptrack %>% 
  mutate(case_ = TRUE) %>%
  select(case_, x_, y_, status, parturient, year, id)
sheeppoints <- sheeppoints %>% # This is where I could deal with the 0s and 1s for status
  mutate(status=NA, parturient=NA, year=NA, id=NA) %>%
  select(case_, x_, y_, status, parturient, year, id)

rsfpoints <- rbind(sheeptrack, sheeppoints)

# Extract covariates
sheeprsf <- extract_covariates(rsfpoints, covariates)

# Scale the distance variables so that the model can be fit better. 
sheeprsf_backup <- sheeprsf #Store the un-scaled version
sheeprsf <- sheeprsf_backup #Convert back to un-scaled version
{
  attach(sheeprsf)
  sheeprsf$dbarren <- scale(dbarren)
  sheeprsf$descp <- scale(descp)
  sheeprsf$dforest <- scale(dforest)
  sheeprsf$dgrassforb <- scale(dgrassforb)
  sheeprsf$dgrassforbshrub <- scale(dgrassforbshrub)
  sheeprsf$droad <- scale(droad)
  sheeprsf$dtrail <- scale(dtrail)
  sheeprsf$dwater <- scale(dwater)
  sheeprsf$elevation <- scale(elevation)
  detach(sheeprsf)
}

# Fit the model
rsf1 <- glm(case_ ~ aspecte*parturient + aspectn*parturient + aspectw*parturient + heatload*parturient + 
              descp*parturient + dgrassforbshrub*parturient + droad*parturient + dtrail*parturient + dwater*parturient + 
              elevation*parturient + percentcrown*parturient + snowdepth*parturient + vrm*parturient, 
            data=sheeprsf, family=binomial(link='logit'))

summary(rsf1)

# include random intercept for individual?
# can run separate RSF for separate behavioural states
# can run with interaction with repro state to see if there are differences.
# How do you deal with available points? Do you make part status all zeros?


















