#=============================================================================================================#
# Script created by Jennifer C Molloy, jenny.molloy@cantab.net
# Script created in version R 3.1.2 
# This script is for analyzing and displaying data related to the Molloy et al. paper describing
# Wolbachia modulation of lipid metabolism in Aedes albopictus mosquito cells
# Datasets can be found at http://www.ebi.ac.uk/metabolights/MTBLS210
#=============================================================================================================#

# Set working directory
setwd("/home/jenny/Dropbox/LipidPaper/")

#clear R environment
rm(list=ls())

library("plyr") 
library("lattice")
library("ggplot2")
library("reshape2")
library("extrafont")
library("stringr")
library("dplyr")
library("cairoDevice")
library("matrixStats")
library("scales")
library("grid")
load(fonts)

#=============================================================================================================#
##################### Predefine necessary functions for plotting #####################################
#=============================================================================================================#
## Function for arranging ggplots. use png(); arrange(p1, p2, ncol=1); dev.off() to save.
### From Stephen Turner via http://www.gettinggeneticsdone.com/2010/03/arrange-multiple-ggplot2-plots-in-same.html


vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
arrange_ggplot2 <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {
  dots <- list(...)
  n <- length(dots)
  if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}
  if(is.null(nrow)) { nrow = ceiling(n/ncol)}
  if(is.null(ncol)) { ncol = ceiling(n/nrow)}
  ## NOTE see n2mfrow in grDevices for possible alternative
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )
  ii.p <- 1
  for(ii.row in seq(1, nrow)){
    ii.table.row <- ii.row	
    if(as.table) {ii.table.row <- nrow - ii.table.row + 1}
    for(ii.col in seq(1, ncol)){
      ii.table <- ii.p
      if(ii.p > n) break
      print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))
      ii.p <- ii.p + 1
    }
  }
}

# Need to display fold-change in lipid abundance with the y-intercept at 1 (no change)
# so remove -1 to 1 in the axis as these are redundant values.
# Function scaleBreaker "takes a chunk" out of the series
# written by dsparks and available from https://gist.github.com/dsparks/4086721#file-rcp_plot-r
scaleBreaker <- function(x, mn, mx){
  y <- x  #              ^ vector, low end, high end
  y[x > mn & x < mx] <- mn
  y[x >= mx] <- x[x >= mx] - mx + mn
  return(y)
}

# Make it into a scale breaking between -1 and 1, with the package scales
break_trans = function() trans_new("break",  #             min v  v max
                                   function(x) scaleBreaker(x, -1, 1),
                                   function(x) x)

#=============================================================================================================#
################# COUNT DATA FOR FOLD-CHANGE #########################
#=============================================================================================================#

# import data from file S16 Table saved as a csv
full.lcms <- read.csv("/home/jenny/Dropbox/LipidPaper/Test/S16\ Table.csv",
                      sep=",", 
                      header=TRUE,
                      stringsAsFactors = FALSE,
                      skip=17
                      )
head(full.lcms)

#Split lcms$Annotation.note into lipid class, carbon chain length and saturation using regular expressions
lcms <- subset(full.lcms, grepl("^(\\w+\\-?\\w+)\\s([1-9][0-9]):([0-9])", Annotation.note)) #select rows containing annotation
matches <- str_match(lcms$Annotation.note,"^(\\w+\\-?\\w+)\\s([1-9][0-9]):([0-9])")
lcms.annotation <- data.frame(na.omit(matches[,-1]),stringsAsFactors=FALSE) # transform to data frame
colnames(lcms.annotation) <- c("Class","Chain","Saturation") # add a header
lcms.annotation$Chain <- as.integer(lcms.annotation$Chain) # convert type of chain length from character to integer
lcms.annotation$Saturation <- as.integer(lcms.annotation$Saturation) # convert type of saturation from character to integer

#Bind lcms.annotation data to lcms
lcms <- cbind(lcms,lcms.annotation)


# Make a table with numbers of species in each class that:
# "-|-" = Decrease in both wMel and wMelPop infection with q.value < 0.05
# "+|+" = Increase in both wMel and wMelPop infection with q.value < 0.05
# "+|-" = Increase in wMel and decrease in wMelPop with q.value < 0.05 
# "-|+" = Increase in wMelPop and decrease in wMel with q.value < 0.05
# "None" = No significant changes in either infection with q.value < 0.05

ddply(lcms,.(Class),
      summarise,
      "-|-" = length(Class[q.value < 0.05 & Change.A.to.AM < 1 & Change.A.to.AP < 1]),
      "+|+" = length(Class[q.value < 0.05 & Change.A.to.AM > 1 & Change.A.to.AP > 1]),
      "+|-" = length(Class[q.value < 0.05 & Change.A.to.AM > 1 & Change.A.to.AP < 1]),
      "-|+" = length(Class[q.value < 0.05 & Change.A.to.AM < 1 & Change.A.to.AP > 1]),
      "None" = length(Class[q.value > 0.05]),
      "Total" = length(Class)
)

#=============================================================================================================#
################### MEAN FOLD-CHANGE DATA PLOTS #################################
#=============================================================================================================#

#Select only lipids manually annotated and flagged for plotting
#Load list of selected MZ values
plot.selection <- read.csv("/home/jenny/Dropbox/LipidPaper/Test/manual-plotting-selection.csv",
                           sep=",", 
                           header=TRUE,
                           stringsAsFactors = FALSE,
)
## Rename lcms$M.Z column to MZ
lcms <- rename(lcms, MZ = M.Z)
## Select only lipid species in manually curated list 
lcms.subset <- subset(lcms,  MZ %in% plot.selection$MZ)

#Convert change ratios to fold change and bind new column
lcms.subset  <- lcms.subset  %>% 
mutate(wMel.foldchange = ifelse(Change.A.to.AM<1,
                         (-(1/Change.A.to.AM)),
                         Change.A.to.AM))  %>%
mutate(wMelPop.foldchange = ifelse(Change.A.to.AP<1,
                                   (-(1/Change.A.to.AP)),
                                   Change.A.to.AP)) 

# make vectors of value and variable names
lcms.names <- (names(lcms))
lcms.valuenames <- c("wMel.foldchange","wMelPop.foldchange")
lcms.variablenames <- lcms.names[! lcms.names %in% lcms.valuenames]

#Melt dataset
lcms.full.long <- melt(lcms.subset, id=lcms.variablenames)
str(lcms.full.long)
glimpse(lcms.full.long)

#Define lipid class categories for plotting and select only significant changes (q.value < 0.05)
## Only plotting [M+Hac-H]- ions for Cer so remove all other ions
lcms.long.cer <- filter(lcms.full.long, Class == "Cer" & Ion.form == "[M+Hac-H]-" & q.value < 0.05)
## Only plotting [M-H]- ions for HexCer so remove all other ions
lcms.long.hexcer <- filter(lcms.full.long, Class == "HexCer" & Ion.form == "[M-H]-" & q.value < 0.05)
## Only plotting [M+Hac-H]- ions for LacCer so remove all other ions
lcms.long.laccer <- filter(lcms.full.long, Class == "LacCer" & Ion.form == "[M+Hac-H]-" & q.value < 0.05)
## Only plotting [M-H]- ions for PE-Cer so remove all other ions
lcms.long.pecer <- filter(lcms.full.long, Class == "PE-Cer" & Ion.form == "[M-H]-" & q.value < 0.05)
## Only plotting [M+Hac-H]- ions for SM so remove all other ions
lcms.long.sm <- filter(lcms.full.long, Class == "SM" & Ion.form == "[M+Hac-H]-" & q.value < 0.05)
## Only plotting [M+Hac-H]- ions for DG so remove all other ions
lcms.long.dg <- filter(lcms.full.long, Class == "DG" & Ion.form == "[M+Hac-H]-" & q.value < 0.05)
## Only plotting [M+Hac-H]- ions for PC so remove all other ions
lcms.long.pc <- filter(lcms.full.long, Class == "PC" & Ion.form == "[M+Hac-H]-" & q.value < 0.05)
## Only plotting [M+Hac-H]- ions for LPC so remove all other ions
lcms.long.lpc <- filter(lcms.full.long, Class == "LPC" & Ion.form == "[M+Hac-H]-" & q.value < 0.05)
## Only plotting [M-H]- ions for PE so remove all other ions
lcms.long.pe <- filter(lcms.full.long, Class == "PE" & Ion.form == "[M-H]-" & q.value < 0.05)
## Only plotting [M-H]- ions for LPE so remove all other ions
lcms.long.lpe <- filter(lcms.full.long, Class == "LPE" & Ion.form == "[M-H]-" & q.value < 0.05)
## Only plotting [M-H]- ions for PG so remove all other ions
lcms.long.pg <- filter(lcms.full.long, Class == "PG" & Ion.form == "[M-H]-" & q.value < 0.05)
## Only plotting [M-H]- ions for PI so remove all other ions
lcms.long.pi <- filter(lcms.full.long, Class == "PI" & Ion.form == "[M-H]-" & q.value < 0.05)
## Too few PS per ion form to filter on this
lcms.long.ps <- filter(lcms.full.long, Class == "PS" & q.value < 0.05)

lcms.long <- rbind(lcms.long.cer,
                   lcms.long.hexcer,
                   lcms.long.pecer,
                   lcms.long.sm,
                   lcms.long.dg,
                   lcms.long.pc,
                   lcms.long.lpc,
                   lcms.long.pe,
                   lcms.long.pg,
                   lcms.long.pi,
                   lcms.long.ps )
# LacCer and LPE not included due to lack of values (only five datapoints)

## Set lipid class groups to enable grouping of data by class in graphs
cer.class <- c("Cer","HexCer","PE-Cer","SM")
phospho.class <- c("PS","PI","PG","PE","PC","PA", "LPC") 
fa.class <-  c("DG","TG")

#Subset long dataset by class groups
lipid.class.cer <- subset(lcms.long , Class %in% cer.class)
lipid.class.phospho <- subset(lcms.long , Class %in% phospho.class)
lipid.class.fa <- subset(lcms.long , Class %in% fa.class)
glimpse(lipid.class.cer)

#=============================================================================================================#
############ MEAN FOLD-CHANGE CALCULATIONS ############################
#=============================================================================================================#

# Calculate table of mean raw ratio change in cell lines by lipid class 
ratio.change <- function(x){c(
  mean.ratio.wMel<- mean(x$Change.A.to.AM),
  mean.ratio.wMelPop <- mean(x$Change.A.to.AP))
  cbind(mean.ratio.wMel,mean.ratio.wMelPop)}
ratio.change.byclass <- ddply(lcms, c("Class"),function(df)ratio.change(df))

#=============================================================================================================#
####### SPHINGOLIPIDS FOLD-CHANGE PLOTS##############
#=============================================================================================================#

#Plot mean fold-change per sphingolipid class
## Pre-define y axis breaks and labels
cer.breaks <- myLabels <- c(-12:-2, 1, 2:8)

cer.plot <- ggplot(lipid.class.cer , aes(x = Class, 
                                         y = value, 
                                         fill=variable)) +
  scale_fill_manual(values=c("#999999", "#CCCCCC"),
                    labels=c("Aa23.wMel", "Aa23.wMelPop"))+
  geom_boxplot(outlier.size = 0) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, 
                                             jitter.height = 0, 
                                             dodge.width = 0.75), 
             alpha=0.6,
             pch=19)+
  
  # format axes
  scale_y_continuous(breaks = cer.breaks)+
  coord_trans(y = "break")+  # Apply custom scale
  geom_hline(yintercept=1, linetype="dotted", size=1)+
  xlab("Sphingolipid Subclass")+
  ylab("Mean fold-change compared to Aa23-T")+
  # format legend and theme  
  guides(fill=guide_legend(title=NULL))+
  theme_bw(base_size = 14)+
  theme(legend.justification=c(1,0), 
        legend.position=c(1,0),
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=0.9),
        legend.background = element_rect(fill="transparent"))
cer.plot

#=============================================================================================================#
####### PHOSPHOLIPIDS FOLD-CHANGE PLOTS ##############
#=============================================================================================================#

#Plot mean fold-change per phospholipid class
phospho.breaks <- myLabels <- c(-7:-2, 1, 2:5)
phospho.plot <- ggplot(lipid.class.phospho , aes(x = Class, 
                                                 y = value, 
                                                 fill=variable)) +
  scale_fill_manual(values=c("#999999", "#CCCCCC"),
                    labels=c("Aa23.wMel", "Aa23.wMelPop"))+
  geom_boxplot(outlier.size = 0) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, 
                                             jitter.height = 0, 
                                             dodge.width = 0.75),
             alpha=0.6,
             pch=19)+
  # format axes
  scale_y_continuous(breaks = cer.breaks)+
  coord_trans(y = "break")+  # Apply custom scale
  geom_hline(yintercept=0, linetype="dotted", size=1)+
  xlab("Phospholipid Subclass")+
  ylab("Mean fold-change compared to Aa23-T")+
  guides(fill=guide_legend(title=NULL))+
  theme_bw(base_size = 14)+
  theme(legend.justification=c(0,1), 
        legend.position=c(0,1),
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=0.9),
        legend.background = element_rect(fill="transparent"))
phospho.plot

#=============================================================================================================#
####### FATTY ACIDS FOLD-CHANGE PLOTS ##############
#=============================================================================================================#

#Plot mean fold-change per FA class
phospho.breaks <- myLabels <- c(-6:-2, 1, 2:4)
fa.plot <- ggplot(lipid.class.fa , aes(x = Class, 
                                       y = value, 
                                       fill=variable)) +
  scale_fill_manual(values=c("#999999", "#CCCCCC"),
                    labels=c("Aa23.wMel", "Aa23.wMelPop"))+
  geom_boxplot(outlier.size = 0) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, 
                                             jitter.height = 0, 
                                             dodge.width = 0.75),
             alpha=0.6,
             pch=19)+
  # format axes
  scale_y_continuous(breaks = cer.breaks)+
  coord_trans(y = "break")+  # Apply custom scale
  geom_hline(yintercept=0, linetype="dotted", size=1)+
  ylab("Mean fold-change compared to Aa23-T")+
  guides(fill=guide_legend(title=NULL))+
  theme_bw(base_size = 14)+
  theme(legend.justification=c(0,1), 
        legend.position=c(0,0.3),
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=0.9),
        legend.background = element_rect(fill="transparent"))
fa.plot


#=============================================================================================================#
############ DETAILED FOLD-CHANGE PLOTS BY CHAIN LENGTH AND SATURATION ######################
#=============================================================================================================#

############ Ceramides ############
cer.detailed.breaks <- myLabels <- c(-12:-2, 1, 2:4)
cer.detailed.plot <- ggplot( subset(lipid.class.cer, Class == "Cer" & Saturation < 3), 
                             #use Cer class only and ignore polyunsaturated lipids with >2 double bonds (only one species) 
                             aes(x = as.factor(Chain), #make the x axis categorical
                                 y = value, 
                                 fill=variable,
                                 width=0.75)) +
  facet_grid(~Saturation
             ,scales = "free_x"
  )+ #free scales so no gaps along x-axis
  scale_fill_manual(values=c("#4D4D4D", "#CCCCCC"),
                    labels=c("Aa23.wMel", "Aa23.wMelPop"))+
  geom_bar(stat="identity", position="dodge")+
  # format axes
  scale_y_continuous(breaks = cer.detailed.breaks)+
  coord_trans(y = "break")+  # Apply custom scale
  geom_hline(yintercept=0, size=0.5)+
  ylab("Mean fold-change compared to Aa23-T")+
  xlab("Carbon Chain Length")+
  guides(fill=guide_legend(title=NULL))+
  theme_bw(base_size = 12, base_family = "sans")+
  theme(legend.justification=c(0,1), 
        legend.position=c(0,0.2),
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=0.9),
        legend.background = element_rect(fill="transparent"),
        strip.text.x = element_text(size=12, face="bold"))
cer.detailed.plot

############ DGs ############

dg.detailed.breaks <- myLabels <- c(-6:-2, 1, 2:4)
dg.detailed.plot <- ggplot(lipid.class.fa, 
                           #use Cer class only and ignore polyunsaturated lipids with >2 double bonds (only one species) 
                           aes(x = as.factor(Chain), #make the x axis catergorical
                               y = value, 
                               fill=variable,
                               width=0.75)) +
  facet_grid(~Saturation, scales = "free_x")+ #free scales so no gaps along x-axis
  scale_fill_manual(values=c("#4D4D4D", "#CCCCCC"),
                    labels=c("Aa23.wMel", "Aa23.wMelPop"))+
  geom_bar(stat="identity", position="dodge")+
  # format axes
  scale_y_continuous(breaks = dg.detailed.breaks)+
  coord_trans(y = "break")+  # Apply custom scale
  geom_hline(yintercept=0, size=0.5)+
  ylab("Mean fold-change compared to Aa23-T")+
  xlab("Carbon Chain Length")+
  guides(fill=guide_legend(title=NULL))+
  theme_bw(base_size = 12, base_family = "sans")+
  theme(legend.justification=c(0,0.75), 
        legend.position=c(0.75,0.3),
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=0.9),
        legend.background = element_rect(fill="transparent"),
        strip.text.x = element_text(size=12, face="bold"))
dg.detailed.plot

############ PCs by saturation and chain length ############
pc.chain.detailed.breaks <- myLabels <- c(-5:-2, 1, 2:4)
pc.chain.detailed.plot <- ggplot(subset(lipid.class.phospho, Class == "PC"& Saturation < 3), 
                                 aes(x = as.factor(Chain), #make the x axis catergorical
                                     y = value, 
                                     fill=variable,
                                     width=0.75)) +
  facet_grid(~Saturation, scales = "free_x")+ #free scales so no gaps along x-axis
  scale_fill_manual(values=c("#4D4D4D", "#CCCCCC"),
                    labels=c("Aa23.wMel", "Aa23.wMelPop"))+
  geom_bar(stat="identity", position="dodge")+
  # format axes
  scale_y_continuous(breaks = pc.chain.detailed.breaks)+
  coord_trans(y = "break")+  # Apply custom scale
  geom_hline(yintercept=0, size=0.5)+
  ylab("Mean fold-change compared to Aa23-T")+
  xlab("Carbon Chain Length")+
  ggtitle("A) PC fold-change by chain length")+
  guides(fill=guide_legend(title=NULL))+
  theme_bw(base_size = 12, base_family = "sans")+
  theme(legend.justification=c(0,1), 
        legend.position=c(0.75,0.3),
        plot.title = element_text(vjust=0.9, size=14, face="bold"),
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=0.9),
        legend.background = element_rect(fill="transparent"),
        strip.text.x = element_text(size=12, face="bold"))
pc.chain.detailed.plot


pc.saturation.detailed.breaks <- myLabels <- c(-5:-2, 1, 2:4)
pc.saturation.detailed.plot <- ggplot(subset(lipid.class.phospho, Class == "PC" & Chain %in% c("36","38","40")), 
                                      aes(x = as.factor(Saturation), #make the x axis catergorical
                                          y = value, 
                                          fill=variable,
                                          width=0.75)) +
  facet_grid(~Chain, scales = "free_x")+ #free scales so no gaps along x-axis
  scale_fill_manual(values=c("#4D4D4D", "#CCCCCC"),
                    labels=c("Aa23.wMel", "Aa23.wMelPop"))+
  geom_bar(stat="identity", position="dodge")+
  # format axes
  scale_y_continuous(breaks = pc.saturation.detailed.breaks)+
  coord_trans(y = "break")+  # Apply custom scale
  geom_hline(yintercept=0, size=0.5)+
  ylab("Mean fold-change compared to Aa23-T")+
  xlab("Number of double bonds")+
  ggtitle("B) PC fold-change by saturation")+
  guides(fill=guide_legend(title=NULL))+
  theme_bw(base_size = 12, base_family = "sans")+
  scale_y_continuous(breaks=scales::pretty_breaks(n = 10))+
  theme(legend.position = "none",
        plot.title = element_text(vjust=0.9, size=14, face="bold"),
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=0.9),
        legend.background = element_rect(fill="transparent"),
        strip.text.x = element_text(size=12, face="bold"))
pc.saturation.detailed.plot

pc.detailed.plot <- arrange_ggplot2(pc.chain.detailed.plot, pc.saturation.detailed.plot, ncol=1)

############ PEs by saturation and chain length ############

pe.saturation <- c(0,1,2,3,4) #set vector of saturations levels to include
pe.saturation.detailed.breaks <- myLabels <- c(-5:-2, 1, 2:5)
pe.detailed.plot <- ggplot(subset(lipid.class.phospho, Class == "PE" & Saturation %in% pe.saturation), 
                           aes(x = as.factor(Chain), #make the x axis catergorical
                               y = value, 
                               fill=variable,
                               width=0.75)) +
  facet_grid(~Saturation, scales = "free_x")+ #free scales so no gaps along x-axis
  scale_fill_manual(values=c("#4D4D4D", "#CCCCCC"),
                    labels=c("Aa23.wMel", "Aa23.wMelPop"))+
  geom_bar(stat="identity", position="dodge")+
  # format axes
  scale_y_continuous(breaks = pe.saturation.detailed.breaks)+
  coord_trans(y = "break")+  # Apply custom scale
  geom_hline(yintercept=0,size=0.5)+
  ylab("Mean fold-change compared to Aa23-T")+
  xlab("Carbon Chain Length")+
  guides(fill=guide_legend(title=NULL))+
  theme_bw(base_size = 14, base_family = "sans")+
  theme(legend.justification=c(0,1), 
        legend.position=c(0.8,1),
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=0.9),
        legend.background = element_rect(fill="transparent"),
        strip.text.x = element_text(size=12, face="bold"))
pe.detailed.plot

#=============================================================================================================#
############# COMBINE FOLD-CHANGE AND INTENSITY DATASETS #################
#=============================================================================================================#

# Read in normalised matrix of signal intensities from LCMS 
# available from http://www.ebi.ac.uk/metabolights/MTBLS210
normmatrix <- t(read.table("/home/jenny/Dropbox/LipidPaper/Tables/TableS4.csv",
                           sep=",", 
                           check.names = TRUE))

# Add MZ column label and remove unecessary rows
normmatrix[1,1]="MZ"
normmatrix<- normmatrix[-c(2:3),]
head(normmatrix)

# Take column names from first row then delete
colnames(normmatrix) <- normmatrix[1,]
normmatrix<- (normmatrix)[-1,] 
normdata <- as.data.frame(normmatrix)

# Subset for experimental samples
normdata.Aa23 <- normdata[,1:20]
head(normdata.Aa23 )
#create column of numbers 1:number of experimental samples
id <- c(1:ncol(normdata.Aa23)) 
# Change factors to numeric in those columns
normdata.Aa23[,id] <- as.numeric(as.character(unlist(normdata.Aa23[,id])))

# Create a Peak column and add the ID number list
normdata.Aa23$Peak = c(1:nrow(normdata.Aa23))

# Prepare to merge datasets by MZ values

## Round up MZ in lcms dataset to 3 d.p for merging (otherwise unable to match MZ index columns)
lcms$MZ  = round(lcms$MZ, digits=3)
## Merge LCMS annotations and normalised intensity matrix
lcms.merged <- join(lcms, normdata.Aa23, by="MZ", type = "inner")

#=============================================================================================================#
########### Plot mean intensities per treatment per species #############
#=============================================================================================================#

## Define column names for each treatment type: Aa23T, Aa23TwMel and Aa23TwMelPop
Aa23T.names <- c("A1","A2","A3","A4","A5","A6")  
Aa23TwMel.names <- c("B1","B2","B3","B4","B5","B6")
Aa23TwMelPop.names <- c("C1","C2","C3","C4","C5","C6")

lcms.merged.selection <- subset(lcms.merged, select = c("MZ","INTENSITY",Aa23T.names,Aa23TwMel.names,Aa23TwMelPop.names))
Aa23T.selection <- subset(lcms.merged, select = Aa23T.names)
Aa23TwMel.selection <- subset(lcms.merged, select = Aa23TwMel.names)
Aa23TwMelPop.selection <- subset(lcms.merged, select = Aa23TwMelPop.names)

## Calculate mean per treatment per signal and add as additional column
lcms.merged$Aa23T.mean <- rowMeans(Aa23T.selection, na.rm = TRUE)
lcms.merged$Aa23TwMel.mean <-  rowMeans(Aa23TwMel.selection, na.rm = TRUE)
lcms.merged$Aa23TwMelPop.mean <-  rowMeans(Aa23TwMelPop.selection, na.rm = TRUE)

# Subset the four columns of interest: "MZ", "Class", "Aa23T.mean", "Aa23TwMel.mean", "Aa23TwMelPop.mean"
# from the lcms.merged dataset 
lcms.merged.means <-(subset(lcms.merged, 
                            select = c("MZ", 
                                       "Class", 
                                       "Aa23T.mean", 
                                       "Aa23TwMel.mean", 
                                       "Aa23TwMelPop.mean")))

# Melt lcms.merged.means into long format ready to feed into ggplot 
lcms.merged.means.long <- melt(lcms.merged.means, id.vars=c("MZ", "Class"))
head(lcms.merged.means.long)
# Plot the mean intensity for each signal by lipid class and treatment 
# Used a boxplot to show distribution, not much meaning added by including all points
# Tried reducing ymax but no real change to message and many outliers not displayed

intensity.boxplot <- 
  ggplot(lcms.merged.means.long, aes(x = Class, 
                                     y = value, 
                                     fill=variable))+
  geom_boxplot(pch=19)+  
  scale_fill_manual(values=c("white", "#999999", "#CCCCCC"),
                    labels=c("Aa23-T", "Aa23.wMel", "Aa23.wMelPop"))+
  xlab("Lipid Class")+
  ylab("Mean Intensity")+
  ylim(0,350000000)+
  guides(fill=guide_legend(title=NULL))+
  theme_bw(base_size = 12, 
           base_family = "sans")+
  theme(legend.justification=c(0,1), 
        legend.position=c(0,1),
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=0.9),
        legend.background = element_rect(fill="transparent"))
intensity.boxplot

#=============================================================================================================#
###### SAVE PLOTS TO FILE ################
#=============================================================================================================#

cairo_ps("/home/jenny/Dropbox/LipidPaper/fig3.eps", 
         width = 4, 
         height = 5, 
         family = "sans")

cer.plot
dev.off()

postscript("/home/jenny/Dropbox/LipidPaper/S17_fig.eps", 
           width = 8, 
           height = 5, 
           horizontal = FALSE, 
           onefile = FALSE, 
           paper = "special", 
           colormodel = "rgb", 
           family = "sans")

cer.detailed.plot
dev.off()

cairo_ps("/home/jenny/Dropbox/LipidPaper/fig5.eps", 
         width = 3, 
         height = 5, 
         family = "sans")

fa.plot
dev.off()

postscript("/home/jenny/Dropbox/LipidPaper/S18_fig.eps", 
           width = 8, 
           height = 5, 
           horizontal = FALSE, 
           onefile = FALSE, 
           paper = "special", 
           colormodel = "rgb", 
           family = "sans")

dg.detailed.plot
dev.off()


cairo_ps("/home/jenny/Dropbox/LipidPaper/fig6.eps", 
         width = 6, 
         height = 5, 
         family = "sans")

phospho.plot
dev.off()

postscript("/home/jenny/Dropbox/LipidPaper/S19_fig.eps", 
           width = 8, 
           height = 10, 
           horizontal = FALSE, 
           onefile = FALSE, 
           paper = "special", 
           colormodel = "rgb", 
           family = "sans")

arrange_ggplot2(pc.chain.detailed.plot, pc.saturation.detailed.plot, ncol=1)
dev.off()


postscript("/home/jenny/Dropbox/LipidPaper/S21_fig.eps", 
           width = 10, 
           height = 5, 
           horizontal = FALSE, 
           onefile = FALSE, 
           paper = "special", 
           colormodel = "rgb", 
           family = "sans")

intensity.boxplot
dev.off()


postscript("/home/jenny/Dropbox/LipidPaper/S20_fig.eps", 
           width = 10, 
           height = 5, 
           horizontal = FALSE, 
           onefile = FALSE, 
           paper = "special", 
           colormodel = "rgb", 
           family = "sans")

pe.detailed.plot
dev.off()
