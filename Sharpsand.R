#Sharpsand LCBH-FSG estimates

#setup, packages
.libPaths('c:/Dan/RPackages')

rm(list=ls())

require(dplyr)
require(ggplot2)

setwd('C:/Dan/_Remote_projects/ccp_2021/analysis')

#get biomass tree data
shp.biomass <- read.csv('./fuel_structure_calcs/sharpsand/sharp-biomass.csv') %>%
  mutate(cd=ht-lcbh, cr=cd/ht)

shp.bio <- filter(shp.biomass, status=='L')
shp.bioD <- filter(shp.biomass, status=='D')

#plot data for overall BA, s/ha
shp.plots <- read.csv('./fuel_structure_calcs/sharpsand/sharp-plots.csv') %>%
  mutate(ba=(dbhMID/2)^2*pi/10000,
         ba.st=ba*live.s.ha,
         trt=case_when(
           Plot %in% c(1, 10, 15, 16) ~ 'th',
           TRUE ~ 'imm' ))

plots <- shp.plots %>% group_by(Plot) %>% summarize(BA=sum(ba.st),
                                              #      DBH=Hmisc::wtd.mean(dbhMID, weights=live.s.ha),
                                                    s.ha.l=sum(live.s.ha),
                                                    s.ha.d=sum(dead.s.ha),
                                                    Trt=first(trt))
tot.meanBA <- mean(plots$BA)
imm.s.ha <- plots %>% filter(Trt=='imm') %>% summarize(s.ha=mean(s.ha.l))


#biomass height model

#params jp (random):
theta1=1.1583
delta1=0.9721
beta1=0.0401
phi1=0.2289
gamma1=0.9399
sig.sq1=1.1421

sht=9
tph=imm.s.ha %>% pull() 
ba=tot.meanBA 

jp.Sharma <- nls(ht ~ 1.3 + (theta1 + u)*sht^delta1 * (1-exp(-beta1*((tph/ba)^phi1)*dbh))^gamma1, 
               start=list(u=1), data=shp.bio)

sharmaHt.fun <- function(x) {
  predict(jp.Sharma, newdata=list(dbh=x))
}

summary(jp.Sharma)$coef

#save(sharma.jp, 'c:/Dan/_Remote_projects/ccp_2020/analysis/models_outputs/sharma_jp')
u_jp <- summary(jp.Sharma)$coef[[1]]

#Show height model

#Crown ratio & LCBH 

cor.test(shp.bio$dbh, shp.bio$cr)

cr.mod <-lm(cr~dbh, data=shp.bio)
#marginally sig. p=0.075
pred.cr <- function(dbh) {
  predict(cr.mod, newdata=list(dbh=dbh), type='response') %>% unname()
}

lcbh.fun <- function(dbh) {
  sharmaHt.fun(dbh)-(pred.cr(dbh) * sharmaHt.fun(dbh))
}

#Graph biomass data, with height and lcbh models
#lcbh.fun uses linear cr increase with dbh
ggplot(shp.bio, aes(x=dbh)) +
  geom_point(aes(y=ht), colour='red') +
  geom_point(aes(y=lcbh), colour='blue') +
  stat_function(inherit.aes=TRUE, fun=sharmaHt.fun, colour='red', n=200) +
  stat_function(inherit.aes=TRUE, fun=lcbh.fun, colour='blue', n=200) +
  labs(x='DBH', y='Height (m)', title='Sharpsand Jack Pine')+
  expand_limits(x=c(0, 12), y=c(0, 12))

##Plots
shp.plots.calc <- shp.plots %>% mutate(ht=sharmaHt.fun(dbhMID),
                                       lcbh=lcbh.fun(dbhMID))

plots2 <- shp.plots.calc %>% group_by(Plot) %>% 
  summarize(TRT=first(trt), 
            HT=Hmisc::wtd.mean(ht, weights=live.s.ha), 
            LCBH=Hmisc::wtd.mean(lcbh, weights=live.s.ha)) %>% 
  arrange(TRT)

#Check plot/burn numbers!! Not sure how to translate...
#Still need to convert:
#Imm: -1
#Imm (1981): -0.5
#Th: -0