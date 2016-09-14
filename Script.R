# THis script performs all analyses and generates all figures for
# Paine, Beck and Terborgh (2016) Ecology
# Last updated 2016-08-31 by C. E. Timothy Paine
# Contact c.e.t.paine@stir.ac.uk with any comments or concerns



##################################
######## Initial Set up     ######
##################################

## Read in libraries
library(vegan)
library(lme4)
library(nlme)
library(survival)
library(doSNOW)
library(compiler)
library(broom)
library(mvtnorm)
library(RColorBrewer)
library(pbkrtest)

# Set working directory
setwd("~/Documents/professional/projects/in progress/multi-exclosures")

#Start a cluster
ncores <- parallel:::detectCores(all.tests = T)
clust <- makeCluster(ncores, type = 'SOCK')
registerDoSNOW(clust)
clusterEvalQ(clust, library(lme4, pbkrtest))
clust

# read in data
dat <- read.csv('exclo data for Ecology.csv')
str(dat)
sppdat <- read.csv('sppdat for Ecology.csv')
str(sppdat)

# modify datasets
sppdat$Sowing.date <- as.Date(sppdat$Sowing.date)
sppdat$Ending.date <- as.Date(sppdat$Ending.date)
sppdat$spcode <- toupper(paste(substr(sppdat$Genus, 0,2), substr(sppdat$species, 0,2), sep = ''))
dat$date      <- as.Date(dat$date)
dat$obs.period<- as.Date(dat$obs.period)
dat$year      <- 1+dat$obs.day%/%365.25 # what year of the experiment was each obs made?
dat$sum       <- dat$seed+dat$seedling
dat$ID        <- paste(dat$species, dat$treatment, dat$block)
dat$Sowing.group <- factor(dat$Sowing.group)
head(dat)
dat$SM <- sppdat$SM[match(dat$species, sppdat$name)]
dat$WD <- sppdat$WD[match(dat$species, sppdat$name)]
dat$TD <- sppdat$total.density[match(dat$species, sppdat$name)]

# make some additional nice useful variables
trts   <- levels(dat$treatment)
trtCOL <- brewer.pal(5, 'Set1')
times  <- as.Date(sort(unique(dat$obs.period)))
blocks <- levels(dat$block)
IDs    <- sort(unique(dat$ID))
n.IDs  <- length(IDs)
COL    <- c(Small = 'black', Medium = 'red', Large = 'blue')
n.spp  <- nrow(sppdat)

# Define the contrast matrix to test effects of mammal groups
# for small   mammals: MCH - MSH
# for medium  mammals: CCA - MSH
# for large   mammals: CS  - CCA
contrast.matrix <- matrix(c(
  # S   M   L
    0,  1, -1,  # CCA
    0,  0,  1,  # CS
    1,  0,  0,  # MCH
   -1, -1,  0   # MSH
    ), ncol = 3, byrow = T) # 1 means allowed in, -1 means excluded
colnames(contrast.matrix) <- c('S', 'M', 'L')
rownames(contrast.matrix) <- c('CCA', 'CS', 'MCH', 'MSH')
contrast.matrix


#######################################################
######## Define some helper functions             #####
#######################################################

# this takes numeric p-values and returns significance stars
pstar_calculator <- function(x){
  frac <- length(x[x < 0])/length(x)
  frac <- ifelse(frac < 0.5, frac, 1-frac)
  out  <- ifelse(frac > 0.05, 'ns', ifelse(frac < 0.0001, "***", "*"))
  return(out)
}

# A function to compute bray-curtis distance between two vectors.
bray <- function(x, y){sum(abs(x-y))/sum (x,y)}

# Calculate the effects of each mammal size-class, given four treatments: MSH, MCH, CCA and CS
#Logs are now base 10. So, -1 = log(1/10). In other words, a 10x reduction is represented by -1, and an 10x increase is indicated by +1
eff.LR.calc <- function(x, col){
  data.frame(
    small  = log(x[x$treatment == 'MCH', col]/x[x$treatment == 'MSH', col], 10),
    medium = log(x[x$treatment == 'CCA', col]/x[x$treatment == 'MSH', col], 10),
    large  = log(x[x$treatment == 'CS',  col]/x[x$treatment == 'CCA', col], 10))
}

# this and the following are taken from the help file of ?pairs
panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks
  nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "gray", border = 'darkgray', ...)
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = 'pairwise.complete.obs'))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.5/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor)
}

# Function to summarize output of bootMer objects
# accepts any arbitrary summarizing function (FUN)
# for graphing USEING log-ratios to compare effects of mammals
#Logs are now base 10. So, -1 = log(1/10). In other words, a 10x reduction is represented by -1, and an 10x increase is indicated by +1
boot.LR.calc <- function(x, nd, FUN, names = paste("SG", 1:4, sep = "_"), ...){
  out <- data.frame(rbind(
    small.509 =  apply(log(x[nd$time == 1 & nd$treatment == 'MCH'] / x[nd$time == 1 & nd$treatment == 'MSH'], 10), 2, FUN, ...),
    medium.509=  apply(log(x[nd$time == 1 & nd$treatment == 'CCA'] / x[nd$time == 1 & nd$treatment == 'MSH'], 10), 2, FUN, ...),
    large.509 =  apply(log(x[nd$time == 1 & nd$treatment == 'CS' ] / x[nd$time == 1 & nd$treatment == 'CCA'], 10), 2, FUN, ...),
    small.END =  apply(log(x[nd$time == 2 & nd$treatment == 'MCH'] / x[nd$time == 2 & nd$treatment == 'MSH'], 10), 2, FUN, ...),
    medium.END=  apply(log(x[nd$time == 2 & nd$treatment == 'CCA'] / x[nd$time == 2 & nd$treatment == 'MSH'], 10), 2, FUN, ...),
    large.END =  apply(log(x[nd$time == 2 & nd$treatment == 'CS' ] / x[nd$time == 2 & nd$treatment == 'CCA'], 10), 2, FUN, ...)))
  names(out) <- names
  return(out)
}


#######################################################
######## END Define some helper functions         #####
#######################################################




#######################################################
######## Change in seed and seedling SURVIVAL     #####
#######################################################
#summarize the data appropriately.
ddd <- numeric()
for(i in 1:length(IDs)){
  dat.i <- dat[dat$ID == IDs[i] & !is.na(dat$sum),]
  if(nrow(dat.i)>1){
    dat.i$end <- c(dat.i$obs.day[-1], NA)
    dat.i$weight.sum  <- -c(diff(dat.i$sum), NA)
    dat.i$event.sum   <- as.numeric(dat.i$weight.sum != 0)
    dat.i.sum <- dat.i[!is.na(dat.i$event.sum),]
    ddd <- rbind(ddd, dat.i.sum)
  }
}
ddd <- ddd[ddd$treatment != 'OPEN',]
ddd$treatment <- factor(ddd$treatment)
contrasts(ddd$treatment) <- contrast.matrix
str(ddd)

survreg.sum1 <- survreg(Surv(end, event.sum) ~ treatment*species, weights = weight.sum, data = ddd[ddd$weight.sum>0,])
summary(survreg.sum1) # Weibull scale parameter = 0.090. THis is < 1, meaning that risks decrease through time
hist(residuals(survreg.sum1))
plot(predict(survreg.sum1), residuals(survreg.sum1), xlab = 'fitted', ylab = 'Residuals') # residuals look OK

############
# FIGURE 1 #
############
#Summary plot to show effects of mammals on seeds as well as seeds+seedlings.
pdf('/Fig 1. summary of mammal effects on  survival times.pdf', width = 4.5, height = 7, paper = 'a4', useDingbats = F)
par(mfrow = c(1, 1), mar = c(4, 1, 1, 1), oma = c(0, 11, 0, 0), las = 1, bty = 'n', tcl = 0.2, lend = 1, cex = 0.8, xpd = NA, mgp = c(2, 0.5, 0))
nd <- expand.grid(treatment = trts[1:4], species = sppdat$name)
pred <- predict(survreg.sum1, newdata = nd, type = 'quantile', p = 0.5, se = T)
nd$predmean <- pred$fit/365.25
nd$predse   <- (pred$se.fit/365.25)^2
toplot <- data.frame(
    species = sppdat$name,
    small.m   = (nd$predmean[nd$treatment == 'MCH'] - nd$predmean[nd$treatment == 'MSH']),
    medium.m  = (nd$predmean[nd$treatment == 'CCA'] - nd$predmean[nd$treatment == 'MSH']),
    large.m   = (nd$predmean[nd$treatment == 'CS' ] - nd$predmean[nd$treatment == 'CCA']),
    small.se  = sqrt(nd$predse[nd$treatment == 'MCH'] + nd$predse[nd$treatment == 'MSH']),
    medium.se = sqrt(nd$predse[nd$treatment == 'CCA'] + nd$predse[nd$treatment == 'MSH']),
    large.se  = sqrt(nd$predse[nd$treatment == 'CS' ] + nd$predse[nd$treatment == 'CCA']))
toplot <- toplot[order(toplot$small.m),]
toplot$small.sig  <- ifelse(toplot$small.m - toplot$small.se > 0 | toplot$small.m + toplot$small.se < 0, 16, 1)
toplot$medium.sig <- ifelse(toplot$medium.m- toplot$medium.se> 0 | toplot$medium.m+ toplot$medium.se< 0, 16, 1)
toplot$large.sig  <- ifelse(toplot$large.m - toplot$large.se > 0 | toplot$large.m + toplot$large.se < 0, 16, 1)
matplot(toplot[,2:4], n.spp:1, type = 'n', axes = F, ylab = '', xlab = '', ylim = c(0, 25), xlim = c(-2, 1))
abline(h = 24.5:0.5, col = 'gray70', xpd = NA, lty = 3)
segments(0, -1, 0, 25, col = 'gray', xpd = NA)
axis(1)
axis(2, at = 24:1, labels = toplot$species, font = 3, cex.axis = 0.9, lwd = 0)
axis(2, at = 0.5:24.5, labels = NA)
for(j in 1:3){points(toplot[,j+1], n.spp:1+0.1*j-0.2, col = COL[j], pch = toplot[,j+7], cex = 1, lwd = 1)}
for(j in 1:3){segments(toplot[,j+1] + toplot[,j+4], n.spp:1+0.1*j-0.2, toplot[,j+1] - toplot[,j+4], n.spp:1+0.1*j-0.2, col = COL[j], lwd = 1)}
title(xlab = 'Change in median survival time (years)', outer = T, line = -2)
legend('bottomleft', col = COL, legend = names(COL), lty = 1, cex = 1, lwd = 1, inset = 0.05, box.col = 'black', bg = 'white', pch = 16)
mtext(adj = 0.05, family = 'Times', 'Figure 1')
dev.off()

apply(toplot, 2, range)
nrow(toplot[toplot$small.m + toplot$small.se < 0,])   # 17 negative effects of small mammals
nrow(toplot[toplot$small.m - toplot$small.se > 0,])   # 1  positive effects of small mammals
nrow(toplot[toplot$medium.m + toplot$medium.se < 0,]) # 14 negative effects of MEDIUM mammals
nrow(toplot[toplot$medium.m - toplot$medium.se > 0,]) # 2  positive effects of MEDIUM mammals
nrow(toplot[toplot$large.m + toplot$large.se < 0,])   # 4  negative effects of LARGE mammals
nrow(toplot[toplot$large.m - toplot$large.se > 0,])   # 8  positive effects of LARGE mammals





#######################################################
############ Change in evenness and diversity: H3 & H3a
#######################################################
#summarize the data appropriately.
datw <- numeric()
for(i in 1:4){
	dat.i <- dat[dat$species %in% sppdat$name[sppdat$Sowing.group==i],]
	datw.i <- reshape(dat.i, v.names = c('seed', 'seedling', 'sum'), timevar = 'species', idvar = c('bout', 'block', 'treatment'), direction = 'wide', drop = c('sort','sort2', 'sp.code','date','ID', 'year'))
	datw.i.sum  <- datw.i[,grep('sum',  names(datw.i))]
	datw.i.sum [is.na(datw.i.sum)]  <- 0
	datw.i$sum.div  <- exp(diversity(datw.i.sum))
	datw.i$sum.even <- diversity(datw.i.sum) /log(specnumber(datw.i.sum ))
	datw.i$sum.div [rowSums(datw.i.sum ) == 0] <- NA # don't calculate diversity if there are no seeds or seedlings
	datw.i <- 	datw.i[, c("obs.period", "obs.day", "Sowing.group", "block", "treatment", "sum.div",   "sum.even")]
	datw <- rbind(datw, datw.i)
}
str(datw)

div.dat <- datw[datw$treatment != 'OPEN',] # drop the open treatment
div.dat <- div.dat[!(div.dat$obs.period < '2000-07-28' | (div.dat$Sowing.group == 3 & div.dat$obs.period < '2004-07-23') )  ,]
div.dat$treatment <- factor(div.dat$treatment)
contrasts(div.dat$treatment) <- contrast.matrix
div.dat$Sowing.group <- factor(div.dat$Sowing.group)
table(div.dat$block, div.dat$treatment) # all blocks and trts are present
str(div.dat)
head(div.dat)


### Analysis of changes in diversity caused by the mammals
# Time-varying analysis of DIVERSITY
lmer.div.sum7  <- lmer(log(sum.div) ~ obs.day * treatment + (obs.day|Sowing.group)+ (1|block), data = div.dat, verbose = 2) # This allows sowing groups to have different intercepts(initial diversity), but also different rates of change through time. This focuses analysis on factor of interest: treatments
car:::Anova(lmer.div.sum7)
summary(lmer.div.sum7)
plot(lmer.div.sum7)

# Test significance of results
div.nd <- expand.grid(time = 1:2, treatment = trts[1:4])
div.nd$obs.day[div.nd$time == 1] <-  509
div.nd$obs.day[div.nd$time == 2] <- max(sppdat$obs.length)
clusterExport(clust, c('div.nd')) # export the newdata out to the cluster
boot.div <- data.frame(bootMer(lmer.div.sum7, FUN = function(x, newdata){exp(predict(x, newdata = div.nd, re.form = NA))}, nsim = 1000, parallel = 'snow', ncpus = ncores, cl = clust, verbose = T))
div.CI   <- boot.LR.calc(boot.div, div.nd, FUN = quantile, names = 'sum', c(0.025, 0.975))
row.names(div.CI) <- paste(rep(names(COL), each = 2), rep(c(509, "END"), each = 6), rep(c('lo', 'hi'), 6), sep = "_")
div.pval <- boot.LR.calc(boot.div, div.nd, FUN = pstar_calculator, names = 'sum')


# Time-varying analysis of EVENNESS
lmer.even.sum7  <- lmer(log(sum.even+0.0001) ~ obs.day * treatment  + (obs.day|Sowing.group)+ (1|block), data = div.dat, verbose = 2) # This allows sowing groups to have different intercepts(initial diversity), but also different rates of change through time. This is nice, because it focuses analysis on factor of interest: treatments
car:::Anova(lmer.even.sum7)
summary(lmer.even.sum7)
plot(lmer.even.sum7)
hist(resid(lmer.even.sum7))

# the nd is the same as for the div analysis and the even analysis.
boot.even <- data.frame(bootMer(lmer.even.sum7, FUN = function(x, newdata){exp(predict(x, newdata = div.nd, re.form = NA))}, nsim = 1000, parallel = 'snow', ncpus = ncores, cl = clust, verbose = T))
even.CI   <- boot.LR.calc(boot.even, div.nd, FUN = quantile, names = 'sum', c(0.025, 0.975), na.rm = T)
row.names(even.CI) <- paste(rep(names(COL), each = 2), rep(c(509, "END"), each = 6), rep(c('lo', 'hi'), 6), sep = "_")
even.pval <- boot.LR.calc(boot.even, div.nd, FUN = pstar_calculator, names = 'sum')



###########
# Fig 2.  #
###########
pdf('Figure 2. Predicted effects of mammals on sum diversity & evenness over time 2016 04 14.pdf', paper = 'a4', width = 3, height = 7, useDingbats = F)
par(mfcol = c(2, 1), mar = c(1, 4, 1, 0), las = 1, bty = 'n', xpd = NA, oma = c(3, 0, 2, 0), tcl = 0.2, mgp = c(2.5, 0.5, 0))
nd <- expand.grid(obs.day = 0:max(sppdat$obs.length), treatment = trts[1:4])
Xpos.pred <- outer(unique(nd$obs.day), -1:1*20, '+')
Xpos <- outer(-1:1*20, unique(div.nd$obs.day), '+')
nd$pred <- exp(predict(lmer.even.sum7, newdata = nd, re.form = NA))
toplot.pred <- eff.LR.calc(nd, 'pred')

matplot(Xpos.pred, toplot.pred, col = COL, lty = 1, type = 'l', ylab = 'Log ratio of Pielou\'s evenness (J)', xlab = '', xaxt = 'n', ylim = range(0.2, toplot.pred, even.CI), xlim = c(0, 1826), lwd = 2)
axis(1, at = 365.25*0:5, labels = NA)
abline(h= 0, col = 'gray', lty = 3, xpd = F)
segments(Xpos, even.CI$sum[grepl('_lo', rownames(even.CI))], Xpos, even.CI$sum[grepl('_hi', rownames(even.CI))], col = COL)
text(Xpos+30,  even.CI$sum[grepl('_hi', rownames(even.CI))], even.pval$sum,  col = COL, cex = 1, adj = 0, xpd =NA)
mtext('A)', adj = 0.05, xpd = NA, line = -1)
mtext('Figure 2', family = 'Times', adj = 0, xpd = NA, line = 1.8)

nd$pred <- exp(predict(lmer.div.sum7, newdata = nd, re.form = NA))
toplot.pred <- eff.LR.calc(nd, 'pred')
matplot(Xpos.pred, toplot.pred, col = COL, lty = 1, type = 'l', ylab =  expression(`Log ratio of number of effective species `(italic(e)^H)), xlab = 'Years since exposure', xaxt = 'n', ylim = range(0.05, toplot.pred, div.CI), xlim = c(0, 1826), lwd = 2)
axis(1, at = 365.25*0:5, labels = 0:5)
abline(h= 0, col = 'gray', lty = 3, xpd = F)
segments(Xpos, div.CI$sum[grepl('_lo', rownames(div.CI))], Xpos, div.CI$sum[grepl('_hi', rownames(div.CI))], col = COL)
text(Xpos+30,  div.CI$sum[grepl('_hi', rownames(div.CI))], div.pval$sum,  col = COL, cex = 1, adj = 0)
mtext('B)', adj = 0.05, xpd = NA, line = -1)
legend('bottomleft', lwd = 2, col = COL, legend = names(COL), inset = 0.05)
dev.off()



###########################################################
#### Change in Species Composition & BETA DIVERSITY    ####
###########################################################
# Through increased among-site variation in species composition.
# This would indicate that mammals affect beta diversity - site-to-site variatiuon in composition. - species turnover through space.
# Beta diversity is defined as the dissimilarty among blocks within treatments at a certain time.
# the changes in beta diversity imposed by mammals could be positive (increasing spatial heterogeneity in species composition) or nagetive (increasing homogeneity in species composition).
# How do mammals affect  dissimilarity in species composition, as compared to  control plots (from which they're excluded)?. Dissimilarity should increase through time.
# I also need to compute difference simply between treatments within blocks (and obs.periods)
# One way to do this is to analyze the difference in sp comp bwteen treatments at any point in time.
# this is bray-curtis distance between relevant treatments, WITHIN blocks. THus, different from beta-diversity calculation
###########################################################
comp.dat <- numeric()
for (i in 1:4) {
  dat.i <- dat[!(dat$obs.period < '2000-07-28' | (dat$Sowing.group == 3 & dat$obs.period < '2004-07-23') )  , ]# remove the observations before all seeds were sown. it's a problem in SG 1, not in 2, 3 or 4.
  dat.i <- dat.i[dat.i$Sowing.group == i,]
  datw.i <- reshape(dat.i, v.names = c("seed", "sum"), timevar = "species", idvar = c("bout", "block", "treatment", "Sowing.group"), direction = "wide", drop = c("sp.code", "date", "seedling", "ID", "year", "obs.day", "SM", "WD", "TD"))
  datw.i <- datw.i[order(datw.i$treatment, datw.i$block, datw.i$obs.period),]
  datw.i$obs.day <- c(datw.i$obs.period - min(datw.i$obs.period))
  datw.seed.i <- datw.i[, grep("seed|obs.period|obs.day|block|treatment", names(datw.i))]
  datw.sum.i  <- datw.i[, grep("sum|obs.period|obs.day|block|treatment" , names(datw.i))]
  ID.i <- paste(datw.sum.i$block,  datw.sum.i$obs.period)
  IDs.i <- unique(ID.i)
  for (j in 1:length(IDs.i)) {
    dat.sum.ij    <- datw.sum.i [ID.i == IDs.i[j], ]
    dat.sum.ij    <- dat.sum.ij [,!apply(dat.sum.ij, 2, function(x){all(is.na(x))})]
    sum.table.ij  <- dat.sum.ij[,grepl('sum', names(dat.sum.ij))]
    small.sum.ij   <- bray(sum.table.ij [dat.sum.ij $treatment == 'MSH',], sum.table.ij [dat.sum.ij$treatment == 'MCH',])
    medium.sum.ij  <- bray(sum.table.ij [dat.sum.ij $treatment == 'MSH',], sum.table.ij [dat.sum.ij$treatment == 'CCA',])
    large.sum.ij   <- bray(sum.table.ij [dat.sum.ij $treatment == 'CCA',], sum.table.ij [dat.sum.ij$treatment == 'CS',])
    comp.dat <- rbind(comp.dat, data.frame(block = dat.sum.ij$block[1:3], obs.period = dat.sum.ij$obs.period[1:3], obs.day = as.numeric(dat.sum.ij$obs.day[1:3]), Sowing.group = i, mammal = names(COL), sum = c(small.sum.ij, medium.sum.ij, large.sum.ij)))
  }
}
comp.dat$Sowing.group <- factor(comp.dat$Sowing.group)
table(comp.dat$block, comp.dat$mammal) # all blocks and mammals are present
comp.dat$mamm_SG <- paste(comp.dat$mammal, comp.dat$Sowing.group, sep = "_")
str(comp.dat)
head(comp.dat)

comp.dat1 <- na.omit(comp.dat[,c('sum', 'obs.day', 'mammal', 'Sowing.group', 'block')])
comp.dat1$SG_block <- factor(paste(comp.dat1$Sowing.group, comp.dat1$block, sep = "_"))

comp.dat2 <- groupedData(sum ~ obs.day | Sowing.group, data = comp.dat1)
comp.dat2$Sowing.group <- factor(comp.dat2$Sowing.group, ordered = F)
str(comp.dat2)

comp.sum.nlmer10 <- nlme(sum ~ SSasympOrig(obs.day, Asym, lrc),
                                data   = comp.dat2,
                                fixed  = list(Asym ~ mammal, lrc ~ 1),
                                random = Asym +lrc ~ 1, verbose = F)

comp.sum.nlmer13 <- update(comp.sum.nlmer10,
      fixed = list(Asym + lrc ~ mammal),
      start = c(0.6710243, 0, 0, -4.4, 0, 0))
comp.sum.nlmer14 <- update(comp.sum.nlmer13, groups = ~ block)
comp.sum.nlmer15 <- update(comp.sum.nlmer13, groups = ~ SG_block)
anova(comp.sum.nlmer10, comp.sum.nlmer13, comp.sum.nlmer14, comp.sum.nlmer15)
summary(comp.sum.nlmer15)
plot(augPred(comp.sum.nlmer15, level = 0:1))

# To simulate for CI estimates
comp.coef <- fixef(comp.sum.nlmer15)
comp.coef[2:3] <- comp.coef[2:3]+ comp.coef[1]
comp.coef[5:6] <- comp.coef[5:6]+ comp.coef[4]
comp.varcov <- data.frame(comp.sum.nlmer15$varFix)
names(comp.varcov) <- row.names(comp.varcov) <- names(comp.coef) <- outer(rev(names(COL)), c('Asym', 'lrc'), paste, sep = "_")

beta.dat <- numeric()
for (i in 1:4) {
  dat.i <- dat[dat$Sowing.group == i, ]
  dat.i <- dat.i[!(dat.i$obs.period < '2000-07-28' | (dat.i$Sowing.group == 3 & dat.i$obs.period < '2004-07-23') )  , ]# remove the observations before all seeds were sown. it's a problem in SG 1, not in 2, 3 or 4.
  datw.i <- reshape(dat.i, v.names = c("seed", "sum"), timevar = "species", idvar = c("bout", "block", "treatment", "Sowing.group"), direction = "wide", drop = c("year", 'WD', "SM", "sp.code", "date", "obs.day", "seedling", "ID"))
  datw.i <- datw.i[order(datw.i$treatment, datw.i$block, datw.i$obs.period),]
  datw.i$obs.day <- c(datw.i$obs.period - min(datw.i$obs.period)) # the c makes it not be a difftime object
  datw.sum.i  <- datw.i[, grep( "sum|obs.period|obs.day|block|treatment", names(datw.i))]
  ID.i <- paste(datw.sum.i$obs.period, datw.sum.i$treatment)
  IDs.i <- unique(ID.i)
  for (j in 1:length(IDs.i)) {
    dat.sum.ij    <- datw.sum.i [ID.i == IDs.i[j], ]
    sum.table.ij  <- dat.sum.ij [,grepl('sum', names(dat.sum.ij))]
    sum.rowSums.ij  <- rowSums(sum.table.ij)
    sum.table.ij  <- sum.table.ij [!is.na(sum.rowSums.ij ),] # need two or more plots to compute inter-plot dissimilarity
    if(nrow(sum.table.ij )>=2){dissim.sum.ij  <- vegdist(sum.table.ij,  method = 'bray')} else{dissim.sum.ij  <- NA}

    beta.dat <- rbind(beta.dat, data.frame(
      obs.period       = as.numeric(dat.sum.ij$obs.period[1]),
      obs.day          = as.numeric(dat.sum.ij$obs.day[1]),
      treatment        = dat.sum.ij$treatment [1],
      Sowing.group     = i,
      ID               = IDs.i[j],
      dissim.sum.mean  = mean(dissim.sum.ij, na.rm = T),
      dissim.sum.sd    = sd(dissim.sum.ij, na.rm = T),
      n.sum.pairs      = length(dissim.sum.ij[is.finite(dissim.sum.ij)])))
  }
}
beta.dat            <- beta.dat[beta.dat$treatment != 'OPEN',]
beta.dat$treatment  <- factor(beta.dat$treatment)
beta.dat$Sowing.group<-factor(beta.dat$Sowing.group)
beta.dat$trt_SG      <- factor(paste(beta.dat$Sowing.group, beta.dat$treatment, sep = "_"))
beta.dat$dissim.sum.sd [is.na(beta.dat$dissim.sum.sd )] <- 0 # this is a fudge to be able to plot SDs
str(beta.dat)


beta.dat1 <- na.omit(beta.dat[,c('dissim.sum.mean', 'obs.day', 'treatment', 'Sowing.group')])
beta.dat2 <- groupedData(dissim.sum.mean ~ obs.day | Sowing.group, data = beta.dat1)
beta.dat2$Sowing.group <- factor(beta.dat2$Sowing.group, ordered = F)
contrasts(beta.dat2$treatment) <- contrast.matrix
str(beta.dat2)
beta.sum.1 <- nlme(dissim.sum.mean ~ SSasympOrig(obs.day, Asym, lrc),
                   data   = beta.dat2,
                   fixed  = list(Asym ~ treatment, lrc ~ 1),
                   random = Asym +lrc ~ 1, verbose = T)
beta.sum.2 <- update(beta.sum.1,
                     fixed = list(Asym ~ treatment, lrc ~ 1),
                     start = c(0.6710243, 0,0, 0, -4.4))
beta.sum.3 <- update(beta.sum.1,
                     fixed = list(lrc ~ treatment, Asym ~ 1),
                     start = c(0.6710243, -4.4, 0, 0, 0))
beta.sum.4 <- update(beta.sum.1,
                     fixed = list(Asym + lrc ~ treatment),
                     start = c(0.6710243, 0, 0, 0, -4.4, 0, 0, 0))
anova(beta.sum.1, beta.sum.2,  beta.sum.4)
summary(beta.sum.4)
plot(augPred(beta.sum.4, level = 0:1))

# Simulate to get CI's
beta.coef <- fixef(beta.sum.4)
beta.coef[2:4] <- beta.coef[2:4]+ beta.coef[1]
beta.coef[6:8] <- beta.coef[6:8]+ beta.coef[5]
beta.varcov <- data.frame(beta.sum.4$varFix)
names(beta.varcov) <- row.names(beta.varcov) <- names(beta.coef) <- outer(trts[1:4], c('Asym', 'lrc'), paste, sep = "_")


# Plot up predicted mammal effects on beta-diversity
pred.beta <- list()
obs.day.i = 0:1612
for(i in 1:4){
    trt.i <- trts[i]
    coef.beta.i  <- rmvnorm(1000, mean = unlist(beta.coef[grep(trt.i, names(beta.coef))]), sigma = as.matrix(beta.varcov[grep(trt.i, row.names(beta.varcov)), grep(trt.i, row.names(beta.varcov))]))
    pred.beta[[i]] <- sapply(obs.day.i, SSasympOrig, Asym = coef.beta.i[,1], lrc = coef.beta.i[,2])
}
names(pred.beta)<- trts[1:4]

pred.beta.eff <- list(
    small  = pred.beta$MCH-pred.beta$MSH,
    medium = pred.beta$CCA-pred.beta$MSH,
    large  = pred.beta$CS -pred.beta$CCA)
pred.beta.eff.CI  <- lapply(pred.beta.eff,  FUN = function(x){apply(x, 2, quantile, c(0.025, 0.975))})
pred.beta.eff.mean <-sapply(pred.beta.eff,  FUN = function(x){apply(x, 2, mean)})
pred.beta.eff.pstar<-sapply(pred.beta.eff,  FUN = function(x){apply(x, 2, pstar_calculator)})




###########
# Fig 3   #
###########
pdf('Figure 3. predicted effects of mammals on species composition & beta_diversity 2016 03 18.pdf', paper = 'a4', width = 4, height = 8, useDingbats = F)
par(mfcol = c(2, 1),  mar = c(1, 5, 1, 0), las = 1, bty = 'n', xpd = NA, oma = c(3, 0, 2, 0), tcl = 0.2, mgp = c(2.5, 0.5, 0))
obs.day <- c(509, 1612)
XPOS <- outer(-1:1*10, obs.day, "+")
# Plot Seed+Seedling species composition predictions
plot(1, xlim = c(0, 1826),  ylim = c(0, 1), ylab = 'Bray-Curtis dissimilarity in\ncommunity composition ', xlab = '', type = 'n', xaxt = 'n')
axis(1, at = 365.25*0:5, labels = NA)
for(i in 1:3){
    mamm.i <- names(COL)[i]
    nd <- expand.grid(mammal = mamm.i, obs.day = 0:1612)
    nd$pred <- predict(comp.sum.nlmer15, newdata = nd, level = 0)
    lines(nd$obs.day+(i-2)*10, nd$pred, col = COL[i], lwd = 2)
    coef.i <- rmvnorm(10000, mean = unlist(comp.coef[grep(mamm.i, names(comp.coef))]), sigma = as.matrix(comp.varcov[grep(mamm.i, names(comp.varcov)),grep(mamm.i, names(comp.varcov))]))
    pred.i <- sapply(obs.day, SSasympOrig, Asym = coef.i[,1], lrc = coef.i[,2])
    comp.CI <- apply(pred.i, 2, quantile, c(0.025, 0.975))
    comp.Pstar <- apply(pred.i, 2, pstar_calculator)
    segments(XPOS[i,], comp.CI[1, ], XPOS[i,], comp.CI[2, ], col = COL[i])
    text(XPOS[i,], comp.CI[2, ], comp.Pstar, col = COL[i], pos = 4)
}
legend('bottom', legend = names(COL), col = COL, lwd = 2, ncol = 1, inset = 0.1)

obs.day.i = 0:1612
# Plot Seed+Seedling beta-diversity predictions
plot(1, xlim = c(0, 1826), type = 'n', ylim = c(-0.2, 0.7), ylab = 'Log ratio of beta-diversity', xlab = "Years since exposure", xaxt = 'n')
axis(1, at = 365.25*0:5, labels = 0:5)
abline(h = 0, col = 'gray', lty = 3, xpd = F)
for(i in 1:3){
    lines(obs.day.i+ (i-2)*10, pred.beta.eff.mean[,i], col = COL[i], lwd = 2)
    segments(c(509, 1612)+ (i-2)*10, pred.beta.eff.CI[[i]][1,obs.day.i %in% c(509, 1612)], c(509, 1612)+ (i-2)*10, pred.beta.eff.CI[[i]][2, obs.day.i %in% c(509, 1612)], col = COL[i])
    text(c(509, 1612), pred.beta.eff.CI[[i]][2, c(509, 1612)], pred.beta.eff.pstar[c(509, 1612),i], col = COL[i], pos = 4)
}
mtext('Figure 3', family = 'Times', adj = 0, line = 0, outer= T)
dev.off()



#############################################################
############ Change in FUNCTIONAL TRAITS & abundance: H3b ###
#############################################################

# Summarise exclosure-mean trait values
trt.dat <- dat
trt.dat <- trt.dat[!(trt.dat$obs.period < '2000-07-28' | (trt.dat$Sowing.group == 3 & trt.dat$obs.period < '2004-07-23') )  ,] # remove the observations before all seeds were sown. it's a problem in SG 1 & 3, not in 2 or 4.
table(trt.dat$obs.period, trt.dat$Sowing.group)
# this claculates exclosure-mean trait values, weighted by the number of SEEDS+SEEDLINGS of each species. This is more useful than SEED alone.
# Note that SM and TD are UNLOGGED logged before means are taken. They get logged in analysis.
trt.agg <- data.frame()
for(i in 1:4){
  dat.i <- trt.dat[trt.dat$Sowing.group == i,]
  periods.i <- sort(unique(dat.i$obs.period))
  for(j in 1:4){
    for(k in 1:8){
      for(m in 1:length(periods.i)){
        dat.ijkm <- dat.i[dat.i$treatment == trts[j] & dat.i$block == blocks[k] & dat.i$obs.period == periods.i[m],]
        if(sum(dat.ijkm$sum) > 0 & !is.na(sum(dat.ijkm$sum))){
          dat.ijkm4SM <- dat.ijkm[!is.na(dat.ijkm$SM),]
          dat.ijkm4WD <- dat.ijkm[!is.na(dat.ijkm$WD),]
          dat.ijkm4TD <- dat.ijkm[!is.na(dat.ijkm$TD),]
          SM.ijkm <- weighted.mean(dat.ijkm4SM$SM, w = dat.ijkm4SM$sum, na.rm = T)
          WD.ijkm <- weighted.mean(dat.ijkm4WD$WD, w = dat.ijkm4WD$sum, na.rm = T)
          TD.ijkm <- weighted.mean(dat.ijkm4TD$TD, w = dat.ijkm4TD$sum, na.rm = T)
          trt.agg <- rbind(trt.agg, data.frame(dat.ijkm[1,c('obs.period', 'obs.day', 'Sowing.group', 'block', 'treatment')], SM = SM.ijkm, WD = WD.ijkm, TD = TD.ijkm))
        }
      }
    }
  }
}
trt.agg$Sowing.group <- factor(trt.agg$Sowing.group)
trt.agg$treatment <- trt.agg$treatment2 <- factor(trt.agg$treatment)
contrasts(trt.agg$treatment) <- contrast.matrix
trt.agg.mean <- aggregate(trt.agg[,c('SM', 'WD', 'TD')], list(obs.period = trt.agg$obs.period, Sowing.group = trt.agg$Sowing.group, treatment = trt.agg$treatment), mean, na.rm = T)
trt.agg.CI   <- aggregate(trt.agg[,c('SM', 'WD', 'TD')], list(obs.period = trt.agg$obs.period, Sowing.group = trt.agg$Sowing.group, treatment = trt.agg$treatment), quantile, c(0.025, 0.975), na.rm = T)
str(trt.agg)
summary(trt.agg)



#analyse log-transformed CWM SEED MASS
lmer.SM5  <- lmer(log(SM) ~ treatment * obs.day + (1|block) + (obs.day|Sowing.group), data = trt.agg, verbose = 2)
# This allows sowing groups to have different intercepts(initial diversity), but also different rates of change through time. This is nice, because it shoudl focus analysis on factor of interest: treatments.
summary(lmer.SM5)
PBmodcomp(lmer.SM2, lmer.SM1, details = 10) # test significance of interaction
car:::Anova(lmer.SM5) # not as good as boottrapping, buut much faster
plot(lmer.SM5)
hist(residuals(lmer.SM5,type="pearson"), breaks = 60) # nicely normal


boot.lmer.SM5   <- data.frame(bootMer(lmer.SM5, FUN = function(x)fixef(x), nsim = 1000, cl = clust, ncpus = ncores, parallel = 'snow', verbose = T))
p.val.SM5 <- apply(boot.lmer.SM5, 2, FUN = function(x){
  frac <- length(x[x< 0])/length(x)
  ifelse(frac<0.5, frac, 1-frac)
})
meanSM<- apply(boot.lmer.SM5, 2, mean)
CI_SM<- apply(boot.lmer.SM5, 2, quantile, c(0.025, 0.975))
TableSM5 <- data.frame(variable = names(boot.lmer.SM5), mean = meanSM, p.val = p.val.SM5, CI_low = CI_SM[1,], CI_high = CI_SM[2,])
rownames(TableSM5) <- NULL



trt.nd <- expand.grid(treatment = trts[1:4], obs.day = c(509, max(sppdat$obs.length))) # use same new data for all three traits
trt.nd$time <- rep(1:2, each = 4)
clusterExport(clust, c('trt.nd')) # export the newdata out to the cluster
boot.SM <- data.frame(bootMer(lmer.SM5, FUN = function(x, newdata){exp(predict(x, newdata = trt.nd, re.form = NA))}, nsim = 1000, parallel = 'snow', ncpus = ncores, cl = clust, verbose = T))
SM.CI    <- boot.LR.calc(boot.SM, trt.nd, FUN= quantile, names = 'sum', c(0.025, 0.975), na.rm = T)
row.names(SM.CI) <- paste(rep(names(COL), each = 2), rep(c(509, "END"), each = 6), rep(c('lo', 'hi'), 6), sep = "_")
SM.pstar <- boot.LR.calc(boot.SM, trt.nd, FUN=pstar_calculator, names = 'sum')

#analyse mean WOOD DENSITY
lmer.WD5 <- lmer(WD ~ treatment * obs.day + (1|block) + (obs.day|Sowing.group), data = trt.agg, verbose = 2)
anova(lmer.WD1, lmer.WD2, lmer.WD3, lmer.WD4, lmer.WD5, lmer.WD6)
summary(lmer.WD5)
plot(lmer.WD5)
car:::Anova(lmer.WD5)
hist(residuals(lmer.WD5,type="pearson"), breaks = 60) # nicely normal

boot.lmer.WD5   <- data.frame(bootMer(lmer.WD5, FUN = function(x)fixef(x), nsim = 1000, cl = clust, ncpus = ncores, parallel = 'snow', verbose = T))
p.val.WD5 <- apply(boot.lmer.WD5, 2, FUN = function(x){
  frac <- length(x[x< 0])/length(x)
  ifelse(frac<0.5, frac, 1-frac)
})
meanWD<- apply(boot.lmer.WD5, 2, mean)
CI_WD<- apply(boot.lmer.WD5, 2, quantile, c(0.025, 0.975))
TableWD5 <- data.frame(variable = names(boot.lmer.WD5), mean = meanWD, p.val = p.val.WD5, CI_low = CI_WD[1,], CI_high = CI_WD[2,])
rownames(TableWD5) <- NULL


boot.WD <- data.frame(bootMer(lmer.WD5, FUN = function(x, newdata){predict(x, newdata = trt.nd, re.form = NA)}, nsim = 1000, parallel = 'snow', ncpus = ncores, cl = clust, verbose = T))
WD.CI    <- boot.LR.calc(boot.WD, trt.nd, FUN= quantile, names = 'sum', c(0.025, 0.975), na.rm = T)
row.names(WD.CI) <- paste(rep(names(COL), each = 2), rep(c(509, "END"), each = 6), rep(c('lo', 'hi'), 6), sep = "_")
WD.pstar <- boot.LR.calc(boot.WD, trt.nd, FUN=pstar_calculator, names = 'sum')

#analyse mean ADULT ABUNDANCE
lmer.TD5  <- lmer(log(TD) ~ treatment * obs.day + (1|block) + (obs.day|Sowing.group), data = trt.agg, verbose = 2)
anova(lmer.TD1, lmer.TD2, lmer.TD3, lmer.TD4, lmer.TD5, lmer.TD6) # including Sowing.group helps A LOT, and allowing
summary(lmer.TD5)
plot(lmer.TD5)
car:::Anova(lmer.TD5)
hist(residuals(lmer.TD5, type="pearson"), breaks = 60) # nicely normal
boot.TD <- data.frame(bootMer(lmer.TD5, FUN = function(x, newdata){exp(predict(x, newdata = trt.nd, re.form = NA))}, nsim = 1000, parallel = 'snow', ncpus = ncores, cl = clust, verbose = T))
TD.CI    <- boot.LR.calc(boot.TD, trt.nd, FUN= quantile, names = 'sum', c(0.025, 0.975), na.rm = T)
row.names(TD.CI) <- paste(rep(names(COL), each = 2), rep(c(509, "END"), each = 6), rep(c('lo', 'hi'), 6), sep = "_")
TD.pstar <- boot.LR.calc(boot.TD, trt.nd, FUN=pstar_calculator, names = 'sum')




###########
# Fig 4   #
###########
pdf('Figure 4. predicted effects of mammals on functional traits in each treatment one panel 2016 04 16.pdf', paper = 'a4', width = 3, height = 8, useDingbats = F)
par(mfcol = c(3, 1), mar = c(1, 4, 1, 0), las = 1, bty = 'n', xpd = NA, oma = c(3, 0, 2, 0), tcl = 0.2, mgp = c(2.5, 0.5, 0))
Xpos.pred <- outer(unique(nd$obs.day), -1:1*20, '+')
Xpos <- outer(-1:1*10, unique(trt.nd$obs.day), '+')
nd <- expand.grid(treatment = trts[1:4], obs.day = 0:max(sppdat$obs.length))

# Plot Seed Mass predictions
nd$pred <- exp(predict(lmer.SM5, newdata = nd, re.form = NA))
toplot.pred <- eff.LR.calc(nd, 'pred')
matplot(Xpos.pred, toplot.pred, col = COL, lty = 1, type = 'l', ylab = 'Log ratio seed mass', xlab = '',  ylim = c(-1, 0.1), xlim = c(0, 1826), xaxt = 'n', lwd = 2)
axis(1, at = 365.25*0:5, labels = NA)
abline(h= 0, col = 'gray', lty = 3, xpd = F)
segments(Xpos, SM.CI[grepl('_lo', rownames(SM.CI)),1], Xpos, SM.CI[grepl('_hi', rownames(SM.CI)),1], col = COL)
text(Xpos+30,  SM.CI[grepl('_hi', rownames(SM.CI)) ,1], SM.pstar$sum,col = COL, cex = 1, adj = 0)
legend('bottomleft', legend = names(COL), col = COL, lty = 1, lwd = 2, ncol = 1, inset = 0.1)
mtext('Figure 4', family = 'Times', adj = 0, xpd = NA, line = 1.8)
mtext('A)', adj = 0.05, xpd = NA, line = -1)

# Plot Wood density predictions
nd$pred <- predict(lmer.WD5, newdata = nd, re.form = NA)
toplot.pred <- eff.LR.calc(nd, 'pred')
matplot(Xpos.pred, toplot.pred, col = COL, lty = 1, type = 'l', ylab = 'Log ratio wood density', xlab = '',  ylim = c(-0.1, 0.1), xlim = c(0, 1826), xaxt = 'n', lwd = 2)
axis(1, at = 365.25*0:5, labels = NA)
abline(h= 0, col = 'gray', lty = 3, xpd = F)
segments(Xpos, WD.CI[grepl('_lo', rownames(WD.CI)),1], Xpos, WD.CI[grepl('_hi', rownames(WD.CI)),1], col = COL)
text(Xpos+30,  WD.CI[grepl('_hi', rownames(WD.CI)) ,1], WD.pstar$sum,col = COL, cex = 1, adj = 0)
mtext('B)', adj = 0.05, xpd = NA, line = -1)

# Plot Adult Density predictions
nd$pred <- exp(predict(lmer.TD5, newdata = nd, re.form = NA))
toplot.pred <- eff.LR.calc(nd, 'pred')
matplot(Xpos.pred, toplot.pred, col = COL, lty = 1, type = 'l', ylab = 'Log ratio in adult density', xlab = 'Years since exposure',  ylim = c(-2.3, 0.5), xlim = c(0, 1826), xaxt = 'n', lwd = 2)
axis(1, at = 365.25*0:5, labels = 0:5)
abline(h= 0, col = 'gray', lty = 3, xpd = F)
segments(Xpos, TD.CI[grepl('_lo', rownames(TD.CI)),1], Xpos, TD.CI[grepl('_hi', rownames(TD.CI)),1], col = COL)
text(Xpos+30,  TD.CI[grepl('_hi', rownames(TD.CI)) ,1], TD.pstar$sum, col = COL, cex = 1, adj = 0)
mtext('C)', adj = 0.05, xpd = NA, line = -1)
dev.off()



#######################
# Supplemental Fig 1  #
#######################
# Visulaize correlations among traits
pdf('Supplemental  Fig. 1. 2016 04 15.pdf', paper = 'a4', width = 5, height = 5, useDingbats = F)
par(bty = 'o', mar = c(0, 0, 0, 0), pty = 'm', las = 1)
to_plot <- sppdat[,c('SM', 'WD', 'total.density')]
to_plot$SM <- log(to_plot$SM)
to_plot$total.density <- log(to_plot$total.density)
pairs(to_plot, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, diag.panel = panel.hist, labels = c('Seed mass (mg)', 'Wood density (g·cm^-3)', 'Adult density (ind·ha^-1)'), span = 7/8)
mtext('Supplemental Figure 1', adj= 0.05, line = -1, family = 'Times', outer = T)
dev.off()



#######################
# Supplemental Fig 2  #
#######################
# Make plots of seed and seedling survival by species
TRTS <- c('Medium', 'Medium+Large', 'Small', 'CLOSED', 'ALL')
pdf('Supplemental Fig. 2. seed+seedling survival by species  2016 04 14.pdf', width = 8, height = 11, useDingbats = F, paper = 'a4')
par(mfrow = c(6, 4), bty = 'n', las = 1, mar = c(1, 1, 0, 0), oma = c(2, 3, 2, 1), xpd = F, tcl = 0.3, mgp = c(2, 0.5, 0), xps = NA)
for (i in 1:n.spp){
    dat.i <- dat[dat$species == sppdat$name[i],]
    if(i == 7){dat.i <- dat.i[dat.i$Sowing.group ==3,]} # cleans up CALOBR.
    dat.i$obs.period <- dat.i$obs.period-min(dat.i$obs.period)
    #if(sppdat$Sowing.group[i]<=2){dat.i$sum <- dat.i$sum/6*10}
    plot(log(sum+1) ~ obs.period, data = dat.i, type = 'n', axes =F, xlim = c(0, 1612), ylim = log(c(1, 11)))
    if(i %%4 ==1){axis(2, at = log(c(1, 2, 6, 11)), labels = c(0, 1, 5, 10))}else{axis(2, labels = NA, at = log(c(1, 2, 6, 11)))}
    if(i >20    ){axis(1,  at = 0:5*365.25, labels = 0:5)}else{axis(1, at = 0:5*365.25, labels = NA)}

    #axis(1, dat.i$obs.period, format(dat.i$obs.period, "%b %Y"), cex.axis = .7, col = 'white')
    mtext(sub(" ", "\n", sppdat$name[i]), adj = 0.8, font = 3, line = -3, cex = 0.7)
    for(j in 1:5){
        dat.ij <- dat.i[dat.i$treatment == trts[j],]
        all.ag.i  <- aggregate(log(dat.ij$sum+1), list(date = dat.ij$obs.period), mean, na.rm = T)
        if(any(all.ag.i$x == 0)){
            all.ag.i <- all.ag.i[1:which(all.ag.i$x == 0)[1], ]
            points(all.ag.i[nrow(all.ag.i),], col = trtCOL[j])
        }
        lines(all.ag.i, col = trtCOL[j], lwd = 1)
    }
    rug(dat.i$obs.period, ticksize = -0.04)
}
legend('right', col = trtCOL, lwd = 2, legend = TRTS, seg.len = 3, cex = 0.8, inset = 0.05)
title(ylab = "Observed number of seeds plus seedlings", xlab = "Years since exposure", outer = T, xpd = NA, cex.lab = 2, line = 1)
mtext('Supplemental Figure 2', family = 'Times', outer = T, adj = 0, xpd = NA)
dev.off()


