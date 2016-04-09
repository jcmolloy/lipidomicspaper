rm(list=ls())

library("plyr") 
library("lattice")
library("ggplot2")
library("car")


#import data from plate file

Wol.Ct<- read.csv("/Users/jenny/Dropbox/DPhil/Lipids/Lipidomics/LipidomicsWolDensity/20140930_JennyLipidomicsWolDensityData.csv")


#####################    CALCULATE STANDARD CURVE    ##########################

#subset data for Standards only
HTH.standards <- subset(Wol.Ct, Type == "Standard" & Gene== "HTH" & Sample != "T")
Wsp.standards <- subset(Wol.Ct, Type == "Standard" & Gene== "Wsp" & Sample != "T")
 

#Plot standard curve
xyplot(log(Copies)~Ct, data=HTH.standards ,
       main="Copy Number Standard Curve", 
       ylab="Copies", 
       xlab="Ct Value",
       type = c("p","r"))

xyplot(log(Copies)~Ct, data=Wsp.standards ,
       main="Copy Number Standard Curve", 
       ylab="Copies", 
       xlab="Ct Value",
       type = c("p","r"))

#Create linear model to describe standard curve
lm.HTH.copies <- lm(log(Copies)~Ct, data = HTH.standards )
lm.Wsp.copies <- lm(log(Copies)~Ct, data = Wsp.standards )

#Extract equation and fit data 
cf.HTH <- coef(lm.HTH.copies) 
r2.HTH<- summary(lm.HTH.copies)$r.squared 
equation.HTH <- data.frame(intercept = cf.HTH[1], 
             slope = cf.HTH[2],  
             Rsq = r2.HTH) 

cf.Wsp <- coef(lm.Wsp.copies) 
r2.Wsp <- summary(lm.Wsp.copies)$r.squared 
equation.Wsp <- data.frame(intercept = cf.Wsp[1], 
                           slope = cf.Wsp[2],  
                           Rsq = r2.Wsp) 

#####################    APPLY LINEAR MODEL AND PREDICT SAMPLE VALUES    ##########################

#Subset data by samples only (excluding standards and controls)
HTH.samples <- subset(Wol.Ct, Type == "Sample" & Gene == "HTH" & Sample %in% c("Mel","Pop","T"))
HTH.samples.cleaned <- droplevels(HTH.samples)
Wsp.samples <- subset(Wol.Ct, Type == "Sample" & Gene== "Wsp" & Sample %in% c("Mel","Pop","T"))
Wsp.samples.cleaned <- droplevels(Wsp.samples)


# #Take mean of technical replicates
# HTH.expression <- ddply(HTH.samples, c("Sample", "BioRep"),
#                                     function(df)c(
#                                       Ct=mean(df$Ct)))
# 
# Wsp.expression <- ddply(Wsp.samples, c("Sample", "BioRep"),
#                                      function(df)c(
#                                       Ct=mean(df$Ct)))

#Predict expression values based on linear model 

HTH.samples.cleaned$predicted <- exp(predict(lm.HTH.copies, HTH.samples.cleaned))
Wsp.samples.cleaned$predicted <- exp(predict(lm.Wsp.copies, Wsp.samples.cleaned))

# Calculate ratio of Wsp to HTH

Wsp.samples.cleaned$WolDensity <- (Wsp.samples.cleaned$predicted/HTH.samples.cleaned$predicted)

# Calculate mean and sd for predicted values across biological replicates

Wsp.density.predicted <- ddply(Wsp.samples.cleaned, c("Sample"),
                                  function(df) c( mean = mean(df$WolDensity),
                                                  sd = sd(df$WolDensity))
                                  )

                                                
############################ PLOT EXPRESSION  ############################################

# Plot exponent of log transformed data

#generate bar plot
expression.plot.linear<- 
  ggplot(NULL, aes(Sample, mean))+ 
  geom_bar(aes(), data=Wsp.density.predicted , stat="identity", position = "dodge", fill="grey")+
  ggtitle(bquote(italic(Wolbachia) ~~ .("density in  Aa23 cell lines")))+
  scale_x_discrete(labels=c(expression(paste( italic(w),"Mel")), 
                             expression(paste( italic(w),"MelPop")),
                             "Uninfected"))+ 
  xlab(bquote(italic(Wolbachia) ~~ .("Strain")))+
  ylab(bquote(italic(Wolbachia) ~~ .("genomes per host genome")))+
  geom_errorbar(data=Wsp.density.predicted, aes(ymin=mean-sd, ymax=mean+sd),
                width=.1,                    # Width of the error bars
                position=position_dodge(.9))+

  theme_bw(base_size = 14, base_family = "sans")+
  theme_bw()
  
expression.plot.linear

#generate box plot
#Open postscript device, generate graphic and close
x.labels = c(expression(paste( italic(w),"Mel")), 
             expression(paste( italic(w),"MelPop")),
             "Uninfected")

postscript("/Users/jenny/Dropbox/DPhil/Lipids/Lipidomics/LipidomicsWolDensity/fig2.eps", 
           width = 3.4252, 
           height = 5, 
           horizontal = FALSE, 
           onefile = FALSE, 
           paper = "special", 
           colormodel = "rgb", 
           family = "sans",
           pointsize = 10)

plot.new() 
par(mar=c(5,4,2,2)+0.1)
plot(WolDensity~Sample, 
     data=Wsp.samples.cleaned,
     ylim=c(0,1600),
     xaxt="n",
     ylab=bquote(italic(Wolbachia) ~~ .("genomes/Host genomes")),
     xlab=bquote(italic(Wolbachia) ~~ .("Strain"))
     )
xaxt="n"
axis(1, at=c(1,2,3), labels=x.labels)
abline(v=c(1.5,2.5), col="grey40", lty="dotted")

dev.off()


########### TEST DATA FOR NORMALITY ####################################

#Check normality of dataset. Shapiro test indicates non-normality if p<0.05.
qqnorm(HTH.samples$predicted)
qqline(HTH.samples$predicted)
hist(HTH.samples$predicted)
shapiro.test(HTH.samples$predicted)


# Bartlett Test of Homogeneity of Variances. Unequal variance if p<0.05
# Bartlett Test is very sensitive to non-normality, so use Levene test for 
# non-normal data as this is more robust

bartlett.test(predicted~Sample, data=HTH.samples)
leveneTest(predicted~Sample, data=HTH.samples)

########## SIGNIFICANCE TESTING ############################################

#Log transformed data is normal and has equal variance so use ANOVA
expression.anova <- aov(predicted~Sample, data=HTH.samples)
summary(expression.anova)

#Check plots for fit
layout(matrix(c(1,2,3,4),2,2))
plot(expression.anova)

#Post-hoc testing of factors
TukeyHSD(expression.anova)