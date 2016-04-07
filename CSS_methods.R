###CSS bird exclusion
###methods data
##ie vaccuum time, methods checks
library(dplyr)
library(ggplot2)

setwd("/Users/colleennell/Documents/R/CSS/data")
meth<-read.csv("CSS_methods.csv")
View(meth)

##how many replicates were vaccuummed for each plant species/treatment?
##treatment T= bird exclusion, C= control
plantcount<-meth %>%
  group_by(SPECIES,treat)%>%
  summarize(nplants=length(plant))
View(plantcount)
##looks good, generally 8 plants in exclusion, 5 control
##what is the avg vaccum time, plant wt by sp and treatment?
vactime<-meth %>%
  group_by(SPECIES,treat)%>%
  summarize(avgvac=mean(min))
View(vactime)
##Does vaccuum time differ by species, treatment?
Anova(lm(avgvac~treat+SPECIES,data=vactime))
##yes varies by species but not treatment. GOOD. 

##read in arthropod ID data
arths<-read.csv("CSS_arths.csv")
View(arths)
arthssamp<-arths%>%
  group_by(SP,Sample,ID,feed)%>%
  summarize(abun=sum(abundance))

View(arthssamp)  

cssmatrix<-dcast(arthssamp,Sample~ID,value.var="abun")
View(cssmatrix)##now in sample x species matrix for multivariate, but 1st sanity checks
cssmatrix$totalarth<-rowSums(cssmatrix[,2:23],na.rm=T)

#calculate mean arthropod density by species, treatment
plants<-read.csv("CSS_plants.csv")
View(plants)
CSS<-full_join(cssmatrix,plants[,c("Sample","species","treat","arth_abun","WT_multiplier","WT_sample")],by="Sample")
View(CSS)
CSS$WT_plant<-CSS$WT_multiplier*CSS$WT_sample


##cleaned dataframe
write.csv(CSS,"CSS_IDmatrix.csv")

clean<-read.csv("CSS_matrix_cleaned.csv")

clean$totalarths<-rowSums(clean[,2:22],na.rm=T)
clean$WT_plant<-clean$WT_multiplier*clean$WT_sample
clean$arth_dens<-clean$totalarths/clean$WT_plant
###arth dens by plant, treat
cleanmean<-clean%>%
  group_by(species,treat)%>%
  summarize(meandens=mean(totalarths))

clean$arth_dens<-clean$arth_dens*1000

###test herbivore density by species, treat
Anova(lm(clean$arth_dens~clean$species+clean$treat),type="III")

##species,treatment herbivore densities
#generate means
std <- function(x) sd(x)/sqrt(length(x))
sptreatmean<-clean%>%
  group_by(species,treat)%>%
  summarize(arthmean=mean(arth_dens),arthse=std(arth_dens),arthsd=sd(arth_dens),
            n=length(arth_dens))
View(sptreatmean)
sptreatmean$var<-((sptreatmean$arthsd)^2)/(length(arth_dens))

densrank<-reorder(sptreatmean$species,sptreatmean$arthmean)
sptreat<-ggplot(sptreatmean,aes(x=densrank,y=arthmean,group=treat,fill=treat))+geom_bar(stat="identity",position=position_dodge())
sptreat
View(sptreatmean)

#calculating pooled SE for bird effect LRR (control/exclusion)
forID<-dcast(sptreatmean,species~treat,value.var="arthmean")
forID2<-dcast(sptreatmean,species~treat,value.var="arthsd")
forIDn<-dcast(sptreatmean,species~treat,value.var="n")
forID2$Cse<-forID2$C
forID2$Tse<-forID2$T
forID<-left_join(forID,forID2,by="species")
forIDn$Cn<-forIDn$C
forIDn$Tn<-forIDn$T
forID<-left_join(forID,forIDn,by="species")
View(forID)

#bird effects
forID$ID<-log(forID$C.x/forID$T.x) ##negative value reflects arthropod removal by birds
forID$var<-(forID$Cse)^2/(forID$Cn*forID$C.x)+(forID$Tse)^2/(forID$Tn*forID$T.x)
forID$IDsd<-sqrt(forID$var)
forID$IDse<-forID$IDsd/(forID$Cn+forID$Tn)

##add to plot
IDrank<-(reorder(forID$species,forID$ID))
spID<-ggplot(forID,aes(x=IDrank,y=ID))+geom_point(size=3)+
  labs(x="Plant Species",y="Bird Effects\nlog(Control/Bird exclusion")+
  scale_x_discrete(limits = rev(levels(IDrank)))+
  geom_hline(yintercept=0,lty="dashed")+
  geom_errorbar(aes(ymin=ID-IDse,ymax=ID+IDse),width=.2)+
  theme_minimal()
spID


###randomly select T plants for arth density in exclusion (4) and calculating ID (4)
set.seed(45)

##take random smaple of plants for ID calculation
sample.df<-function(df,)
IDarca<-clean[sample((1:nrow(clean))[which("species"=="ARCA" & "treat"=="T"),],4),]
View(IDarca)

ARCAT<-clean%>%
  filter(species=="ARCA",treat=="T")
View(ARCAT)
##take 4 rows to cal ID
getgroups<-function(sp,enddf) {
  spdf<-clean%>%
    filter(species==sp,treat=="T")
  spDD<-as.data.frame(spdf[sample(nrow(spdf),4),])
  spDD$group<-"ID"
  spdf<-left_join(spdf,spDD[,c("Sample","group")],by="Sample")
  spdf$group[is.na(spdf$group)]<-"DD"
  enddf<-data.frame(spdf)
}
getgroups("ARCA",ARCATT)
View(ARCATT)##cannot get to return a dataframe to then bind 'group' variable with 'clean' df by sample...
##do long way for each species then bind together
spdf<-clean%>%
  filter(species=="ARCA",treat=="T")
spDD<-as.data.frame(spdf[sample(nrow(spdf),4),])
spDD$group<-"ID"
spdf<-left_join(spdf,spDD[,c("Sample","group")],by="Sample")
spdf$group[is.na(spdf$group)]<-"DD"
ARCAT<-data.frame(spdf)
View(ARCAT)
#ARDO
spdf<-clean%>%
  filter(species=="ARDO",treat=="T")
spDD<-as.data.frame(spdf[sample(nrow(spdf),4),])
spDD$group<-"ID"
spdf<-left_join(spdf,spDD[,c("Sample","group")],by="Sample")
spdf$group[is.na(spdf$group)]<-"DD"
ARDOT<-data.frame(spdf)
View(ARDOT)
#ENCA
spdf<-clean%>%
  filter(species=="ENCA",treat=="T")
spDD<-as.data.frame(spdf[sample(nrow(spdf),4),])
spDD$group<-"ID"
spdf<-left_join(spdf,spDD[,c("Sample","group")],by="Sample")
spdf$group[is.na(spdf$group)]<-"DD"
ENCAT<-data.frame(spdf)
#ERFA
spdf<-clean%>%
  filter(species=="ERFA",treat=="T")
spDD<-as.data.frame(spdf[sample(nrow(spdf),4),])
spDD$group<-"ID"
spdf<-left_join(spdf,spDD[,c("Sample","group")],by="Sample")
spdf$group[is.na(spdf$group)]<-"DD"
ERFAT<-data.frame(spdf)
#ERPA
spdf<-clean%>%
  filter(species=="ERPA",treat=="T")
spDD<-as.data.frame(spdf[sample(nrow(spdf),4),])
spDD$group<-"ID"
spdf<-left_join(spdf,spDD[,c("Sample","group")],by="Sample")
spdf$group[is.na(spdf$group)]<-"DD"
ERPAT<-data.frame(spdf)
#ISME
spdf<-clean%>%
  filter(species=="ISME",treat=="T")
spDD<-as.data.frame(spdf[sample(nrow(spdf),4),])
spDD$group<-"ID"
spdf<-left_join(spdf,spDD[,c("Sample","group")],by="Sample")
spdf$group[is.na(spdf$group)]<-"DD"
ISMET<-data.frame(spdf)
#LUAL
spdf<-clean%>%
  filter(species=="LUAL",treat=="T")
spDD<-as.data.frame(spdf[sample(nrow(spdf),4),])
spDD$group<-"ID"
spdf<-left_join(spdf,spDD[,c("Sample","group")],by="Sample")
spdf$group[is.na(spdf$group)]<-"DD"
LUALT<-data.frame(spdf)
#SAAP
spdf<-clean%>%
  filter(species=="SAAP",treat=="T")
spDD<-as.data.frame(spdf[sample(nrow(spdf),4),])
spDD$group<-"ID"
spdf<-left_join(spdf,spDD[,c("Sample","group")],by="Sample")
spdf$group[is.na(spdf$group)]<-"DD"
SAAPT<-data.frame(spdf)
#SAME
spdf<-clean%>%
  filter(species=="SAME",treat=="T")
spDD<-as.data.frame(spdf[sample(nrow(spdf),4),])
spDD$group<-"ID"
spdf<-left_join(spdf,spDD[,c("Sample","group")],by="Sample")
spdf$group[is.na(spdf$group)]<-"DD"
SAMET<-data.frame(spdf)
#bind these together, write csv
groupeddata<-list(ARCAT,ARDOT,ENCAT,ERFAT,ERPAT,ISMET,LUALT,SAAPT,SAMET)##create a list of the dfs
CSS<-do.call(rbind,groupeddata) ##bind rows, should be using do.call mroe often!
View(CSS)
write.csv(CSS,"CSS_data_groups")##can use this or use to bind groupings to previous df, do not rerun sample selection

##calculate species values to DD (density it treatment), density in control, bird effects(LRR)

##for calcualting ID and DD by species
ARCADDmean<-ARCADD%>%
  group_by(species)%>%
  summarize(DD=mean(arth_dens),DDsd=sd(arth_dens),DDse=std(arth_dens))
ARCAID<-ARCAT
