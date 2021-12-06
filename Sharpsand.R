#Sharpsand LCBH-FSG estimates
#some confusion over Plot vs. Burn numbers
#Data for Plots 8 & 9 (thinned) are missing
#
#Use overall site mean stand height to first model tree heights
#then calculate Lorey's height for each plot
#Use this HT.L as input for more precise cohort modelling for LCBH


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

plots <- shp.plots %>% group_by(Plot) %>% 
  summarize(BA=sum(ba.st),
            DBH=Hmisc::wtd.mean(dbhMID, weights=live.s.ha),
            s.ha.l=sum(live.s.ha),
            s.ha.d=sum(dead.s.ha),
            TRT=first(trt))

tot.meanBA <- mean(plots$BA)
imm.s.ha <- plots %>% filter(TRT=='imm') %>% summarize(s.ha=mean(s.ha.l))


#biomass height model

#params jp (random):
theta1=1.1583
delta1=0.9721
beta1=0.0401
phi1=0.2289
gamma1=0.9399
sig.sq1=1.1421

st.height=9

sht=st.height
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

#Function using calibrated u and tph and ba as inputs
sharmaHt.full <- function(dbh, tph=imm.s.ha %>% pull(), ba=tot.meanBA, sht=st.height) {
  1.3 + (theta1 + u_jp)*sht^delta1 * (1-exp(-beta1*((tph/ba)^phi1)*dbh))^gamma1
}


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

lcbh.full <- function(dbh, tph=imm.s.ha %>% pull(), ba=tot.meanBA, sht=st.height) {
  sharmaHt.full(dbh, tph, ba, sht)-(pred.cr(dbh) * sharmaHt.full(dbh, tph, ba, sht))
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

##Plots - calculate overall model height and LCBH
shp.plots.calc <- shp.plots %>% mutate(ht=sharmaHt.fun(dbhMID),
                                       lcbh=lcbh.fun(dbhMID))

##Plots - calculate calibrated (mixed effect) ht, lcbh; use for imm plots; 
#First, long vectors for BA, s.ha
ba.long <- rep(plots$BA, each=8)
s.ha.long <- rep(plots$s.ha.l, each=8)


shp.plots.mixed <- mutate(shp.plots.calc, ba.plot=ba.long, 
                          s.ha.plot=s.ha.long)


plots2 <- shp.plots.mixed %>% group_by(Plot) %>% 
  summarize(TRT=first(trt), 
            HT=Hmisc::wtd.mean(ht, weights=live.s.ha), 
            HT.L=Hmisc::wtd.mean(ht, weights=ba.st),
            LCBH=Hmisc::wtd.mean(lcbh, weights=live.s.ha))
#            LCBH.mix=Hmisc::wtd.mean(lcbh.mix, weights=live.s.ha))

#plot-based Lorey's height for mean stand height for modelling
sht.L <- rep(plots2$HT.L, each=8)

#final plot-based mixed effects model lcbh using LH as stand height for each plot
shp.plots.mixed <- mutate(shp.plots.mixed, sht.plot=sht.L)
shp.plots.mixed2 <- rowwise(mutate(shp.plots.mixed, lcbh.mix=lcbh.full(
  dbh=dbhMID, tph=s.ha.plot, ba=ba.plot, sht=sht.plot))) #%>%

plots3 <- shp.plots.mixed2 %>% group_by(Plot) %>% 
  summarize(TRT=first(trt), 
            HT=Hmisc::wtd.mean(ht, weights=live.s.ha), 
            HT.L=Hmisc::wtd.mean(ht, weights=ba.st),
            LCBH=Hmisc::wtd.mean(lcbh, weights=live.s.ha),
            LCBH.mix=Hmisc::wtd.mean(lcbh.mix, weights=live.s.ha))



trt.means <- plots2 %>% group_by(TRT) %>% summarize(HT=mean(HT), LCBH=mean(LCBH))

plots.add <- data.frame(Plot=c(8,9), TRT='th', HT=trt.means$HT[trt.means$TRT=='th'], 
                        LCBH=trt.means$LCBH[trt.means$TRT=='th'])

plots3 <- left_join(plots, plots2, by=c('Plot', 'TRT')) 





plots_final <- full_join(plots3, plots.add) %>% arrange(TRT, Plot)
   
# write.csv(plots_final,  
#         'c:/Dan/_Remote_projects/ccp_2021/analysis/models_outputs/sharpsand_fsg.csv')



#Still need to convert:
#Imm: -1
#Imm (1981): -0.5
#Th: -0