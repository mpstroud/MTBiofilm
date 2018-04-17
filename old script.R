

#PCA of environmental variables
#PRINCIPAL COMPINENT ANALYSIS WITH STANDARDIZATION----------------------------------------------------
enviro <- enviro_comp[, 4:24]
enviro <- na.omit(enviro)

enviro_class <- enviro[, 1:3]
enviro <- enviro[, 4:21]
enviro <- enviro[c(-1,-5,-8)]

enviro.pca <- prcomp(enviro,center = TRUE,scale. = TRUE, na.action=na.omit)
summary(enviro.pca)



library(ggplot2)
dataplot <- as.data.frame(enviro.pca$x)

ggplot(dataplot, aes(x=PC1, y=PC2, colour=enviro_class$river)) + geom_point(size=4) +
  scale_colour_manual(values=c("#CCCCCC","#000000","#336633","#99CC33","#CC6600","#FF0000","#00FF99","#0000CC","#FFCC33","#663300")) +
  xlim(-5,3) +
  ylim(-5,5) +
  theme(axis.title.x=element_blank()) + 
  theme(axis.title.y=element_blank()) +
  theme(legend.position = "top") +
  theme(legend.direction = "horizontal") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(legend.text = element_text(size=16)) +
  theme(
    panel.grid.major = element_line(colour = NA),
    panel.grid.minor = element_line(colour = NA),
    panel.background = element_rect(fill='white', colour='black'),
    axis.ticks = element_line(colour = "black"),
    axis.line = element_line(colour = "black"))


ggplot(dataplot, aes(x=PC1, y=PC2, colour=enviro_class$zone)) + geom_point(size=2)


ggplot(dataplot, aes(x=PC1, y=PC2, colour=enviro_class$habitat)) + geom_point(size=2)

ggplot(dataplot, aes(x=PC1, y=PC2, colour=enviro_class$habitat)) + geom_point(size=4) +
  scale_colour_manual(values=c("#CCCCCC","#000000","#336633","#99CC33","#CC6600","#FF0000","#00FF99","#0000CC","#FFCC33","#663300")) +
  xlim(-5,3) +
  ylim(-5,5) +
  theme(axis.title.x=element_blank()) + 
  theme(axis.title.y=element_blank()) +
  theme(legend.position = "top") +
  theme(legend.direction = "horizontal") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(legend.text = element_text(size=16)) +
  theme(
    panel.grid.major = element_line(colour = NA),
    panel.grid.minor = element_line(colour = NA),
    panel.background = element_rect(fill='white', colour='black'),
    axis.ticks = element_line(colour = "black"),
    axis.line = element_line(colour = "black"))



biplot(enviro.pca, pc.biplot = TRUE)
pca_epilithon$rotation
pca_epilithon$x


#Plotting VAriance Components
comp<-read.table("clipboard",header=TRUE, sep="\t")
ggplot(comp, aes(x=level, y=y)) + 
  geom_bar(aes(fill=var),  
           stat="identity",
           colour="black",   
           position=position_dodge()) +
theme(axis.title.x=element_blank()) + 
  theme(axis.title.y=element_blank()) +
  theme(legend.position = "top") +
  theme(legend.direction = "horizontal") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(legend.text = element_text(size=16)) +
  theme(
    panel.grid.major = element_line(colour = NA),
    panel.grid.minor = element_line(colour = NA),
    panel.background = element_rect(fill='white', colour='black'),
    axis.ticks = element_line(colour = "black"),
    axis.line = element_line(colour = "black"))
species2 <- aggregate(. ~ taxonomy, data=species, FUN=sum)
write.table(species2, "E:/A_University of Montana/PROJECTS/MARC project/Biofilm 16s/Genetic analysis 2014/species3.txt", sep="\t") 



#Principal Coordinate Analysis (PCoA)---------------------------------------------------------------------------------------------------
species <- read.table("clipboard",header=T, sep="\t")
enviro_class <- read.table("clipboard",header=T, sep="\t")


dist_biofilm <- vegdist(species)
pcoa_biofilm <- cmdscale(dist_biofilm,eig=T)
calc_eigen(pcoa_biofilm$eig)

x <- pcoa_biofilm$points[,1]
y <- pcoa_biofilm$points[,2]

library(ggplot2)
dataplot <- as.data.frame(pcoa_biofilm$points)

ggplot(dataplot, aes(x=V1, y=V2, colour=enviro_class$river)) + geom_point(size=4) +
  scale_colour_manual(values=c("#CCCCCC","#000000","#336633","#99CC33","#CC6600","#FF0000","#00FF99","#0000CC","#FFCC33","#663300")) +
  theme(axis.title.x=element_blank()) + 
  theme(axis.title.y=element_blank()) +
  theme(legend.position = "top") +
  theme(legend.direction = "horizontal") +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(legend.text = element_text(size=16)) +
  theme(
    panel.grid.major = element_line(colour = NA),
    panel.grid.minor = element_line(colour = NA),
    panel.background = element_rect(fill='white', colour='black'),
    axis.ticks = element_line(colour = "black"),
    axis.line = element_line(colour = "black"))


#NMDS

getwd()
species = read.csv("family_rarefied.csv")

species <- species[-c(55), ]
species_NMDS <- species[,6:800]
rownames(species_NMDS) <- paste(species$river, species$sampleID,sep=" ")


example_NMDS=metaMDS(species_NMDS,k=2,trymax=10000)
par(mfrow=c(1,1))
figNMDS<-ordiplot(example_NMDS,type="n",xlim = c(-0.5,0.5),ylim = c(-0.5,0.5))
text(figNMDS, "sites", col="blue", cex=0.7)

ordisurf(example_NMDS,species$np,main="",col="forestgreen")
orditorp(example_NMDS,display="sites",col="grey30",air=0.1,cex=1)
anosim(species_NMDS, species$River_1, permutations=999)

stressplot(example_NMDS)
plot(example_NMDS)

#Diversity metrics--------------------------------------------------------------------------------------
library(vegan)
library(betapart)
#shanon diversity
spec <- species[,6:819]
shannon<-diversity(spec, index = "shannon")
shannon<-as.vector(shannon)
shannon<-as.data.frame(shannon)

plot(shannon)
diversity<-cbind(shannon,species[,1:6])
plot(diversity$river,diversity$shannon)

#for betadiversity
presabs<-ifelse(spec>0,1,0)
dist<-beta.pair(presabs, index.family="sorensen")
# To get the pairwise Jaccard index turnover partition between communities, type: dist[[1]]. 
#To get nestedness partition, type: dist[[2]]. To get all beta diversity: dist[[3]].

# If we want to compare the beta diversities of communities aggregated by groups
#we can use "betadisper" analysis.
groups <- factor(species$zone)
bd<-betadisper(dist[[3]],groups)
plot(bd)
boxplot(bd)
anova(bd)


# Using bray-curtis dissimilarity and thus abudnacne data as well
#The function will calculate beta diversity (and its partitions) based on Bray-Curtis dissimilarity index.

dist<-bray.part(spec)
groups <- factor(species$zone)
bd<-betadisper(dist[[3]],groups)
plot(bd)
boxplot(bd)
anova(bd)

#LMM models--------------------------------------------------------------------------------------------------------------

library(lme4)
library(nlme)
library(piecewiseSEM)

lmm.data<-cbind(shannon,species) #creat the dataset for the model

#BIOFILM DIVERSITY-------------------------------------------------------------------------------------------
#all variables
lmm.res1 <- lmer(formula = response ~ oxygen + oxygen_sat +	temp +	dic +	doc	+ chloride +	ammonium +	
                                    srp	+ np + nitrate	+	dry_mass +	afdm +	perc_om  + (1|region/river/zone/hab), 
                                    data = lmm.data, REML = FALSE)
summary(lmm.res1)
sem.model.fits(lmm.res1)

#selecting fixed effects
lmm.res2 <- lmer(formula = response ~  ammonium + np + temp + (1|river/zone/hab), data = lmm.data, REML = FALSE)
anova(lmm.res1,lmm.res2)
summary(lmm.res2)
anova(lmm.res2)
plot(residuals(lmm.res2))
sem.model.fits(lmm.res2)
coef(lmm.res2)



#CYANOBACTERIA DIVERSITY-------------------------------------------------------------------------------------------
#all variables
lmm.res1 <- lmer(formula = response2 ~ oxygen + oxygen_sat +  temp +	dic +	doc	+ chloride +	ammonium +	
                   srp	+ nitrate	+	np + dry_mass +	afdm +	perc_om  + (1|region/river/zone/hab), 
                 data = lmm.data, REML = FALSE)
summary(lmm.res1)

#selecting fixed effects
lmm.res2 <- lmer(formula = response2 ~  doc + np + (np|river/zone/hab), 
                 data = lmm.data, REML = FALSE)
anova(lmm.res1,lmm.res2)
summary(lmm.res2)
plot(residuals(lmm.res2))
sem.model.fits(lmm.res2)


lme.nested <- lme(response2 ~ oxygen + temp +   nitrate  ,random= ~nitrate|river + ~1|hab, data=lmm.data,na.action=na.omit,method="ML")#Nested without random effects
summary(lme.nested)

