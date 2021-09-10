## ----setup, include=FALSE------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ------------------------------------------------------
library(dsims)
library(leaflet)
library(knitr)
library(sf)


## ----studyarea-----------------------------------------
northsea <- make.region(region.name = "minkes",
                      shape = "\\shapefiles\\Strataprj.shp",
                      strata.name = c("South", "North"),
                      units = "km")


## ----leafletmap----------------------------------------
m <- leaflet() %>% addProviderTiles(providers$Esri.OceanBasemap)
m <- m %>% 
  setView(1.4, 55.5, zoom=5)
minkes <- read_sf("\\shapefiles\\Strataprj.shp")
study.area.trans <- st_transform(minkes, '+proj=longlat +datum=WGS84')
m <- addPolygons(m, data=study.area.trans$geometry, weight=2)
m


## ----popspec-------------------------------------------
areas <- northsea@area
prop.south <- areas[1]/sum(areas)
prop.north <- areas[2]/sum(areas)
total.abundance <- 3000
abund.south <- round(total.abundance * prop.south)
abund.north <- round(total.abundance * prop.north)
constant <- make.density(region = northsea, x.space = 10, constant = 1)
minkepop <- make.population.description(region = northsea,
                                           density = constant,
                                           N = c(abund.south, abund.north))


## ----survdesign----------------------------------------
coverage.grid <- make.coverage(northsea, n.grid.points = 100)
equal.cover <- make.design(region = northsea,
                           transect.type = "line",
                           design = "systematic",
                           samplers=40, 
                           design.angle = c(50, 40),
                           truncation = 0.8,
                           coverage.grid = coverage.grid)


## ----designprop----------------------------------------
design.properties <- run.coverage(equal.cover, reps = 10, quiet=TRUE)
mine <- data.frame(Num.transects=design.properties@design.statistics$sampler.count[3,],
                   Proportion.covered=design.properties@design.statistics$p.cov.area[3,])
kable(mine)


## ----whatanalysis--------------------------------------
pooled.hn <- make.ds.analysis(dfmodel = list(~1),
                                key = "hn",
                                criteria = "AIC",
                                truncation = 0.8)
strat.specific.or.not <- make.ds.analysis(dfmodel = list(~1, ~Region.Label),
                                key = "hn",
                                criteria = "AIC",
                                truncation = 0.8)


## ----sigma---------------------------------------------
delta.multiplier <- c(seq(from=0.6, to=0.8, by=0.2), 1,
#                      seq(from=0.85, to=1.15, by=0.1),
                      seq(from=1.2, to=1.8, by=0.2))
sigma.south <- 0.5
north.sigma <- sigma.south*delta.multiplier


## ----sigmarange----------------------------------------
hn <- function(sigma, x) {return(exp(-x^2/(2*sigma^2)))}
for (i in seq_along(north.sigma)) {
  curve(hn(north.sigma[i],x),from=0,to=0.8,add=i!=1,  
        xlab="Distance", ylab="Detection probability", 
        main="Range of detection probability disparity\nSouth function in blue")
}
curve(hn(sigma.south,x),from=0,to=0.8, lwd=2, col='blue', add=TRUE)


## ----loopsigma-----------------------------------------
equalcover <- list()
whichmodel <- list()
num.sims <- 300
for (i in seq_along(delta.multiplier)) {
  sigma.strata <- c(sigma.south, sigma.south*delta.multiplier[i])
  detect <- make.detectability(key.function = "hn",
                               scale.param = sigma.strata,
                               truncation = 0.8)
  equalcover.sim <- make.simulation(reps = num.sims,
                                    design = equal.cover,
                                    population.description = minkepop,
                                    detectability = detect,
                                    ds.analysis = pooled.hn)
  whichmodel.sim <- make.simulation(reps = num.sims,
                                    design = equal.cover,
                                    population.description = minkepop,
                                    detectability = detect,
                                    ds.analysis = strat.specific.or.not)
  equalcover[[i]] <- run.simulation(equalcover.sim, run.parallel = TRUE, max.cores=6)
  whichmodel[[i]] <- run.simulation(whichmodel.sim, run.parallel = TRUE, max.cores=6)
}


## ----meaning, echo=FALSE-------------------------------
pctbias <- mat.or.vec(length(delta.multiplier),3)
ci.cover <- mat.or.vec(length(delta.multiplier),3)
num.detects <- mat.or.vec(length(delta.multiplier),3)
modelsel <- mat.or.vec(length(delta.multiplier),1)
bias.modsel <- mat.or.vec(length(delta.multiplier),3)
cover.modsel <- mat.or.vec(length(delta.multiplier),3)
for (i in seq_along(delta.multiplier)) {
  simsum <- summary(equalcover[[i]], description.summary=FALSE)
  selsum <- summary(whichmodel[[i]], description.summary=FALSE)
  pctbias[i, ] <- simsum@individuals$D[,3]
  ci.cover[i, ] <- simsum@individuals$D[,5]
  num.detects[i, ] <- simsum@individuals$summary[,3]
  bias.modsel[i, ] <- selsum@individuals$D[,3]
  cover.modsel[i, ] <- selsum@individuals$D[,5]
  modelsel[i] <- sum(whichmodel[[i]]@results$Detection[1,5, 1:num.sims]==2) / num.sims
}
pooled <- as.data.frame(cbind(pctbias, ci.cover, num.detects, modelsel))
names(pooled) <- c("Bias.N", "Bias.S", "Bias.Tot",
                "CICover.N", "CICover.S", "CICover.Tot",
                "Detects.N", "Detects.S", "Detects.Tot", "Prop.strat.model")
row.names(pooled) <- delta.multiplier
kable(pooled, digits=c(rep(1,3), rep(2,3), rep(0,3), 2), 
      caption="Results using pooled detection function")


## ----selection, echo=FALSE-----------------------------
with.selection <- as.data.frame(cbind(bias.modsel, cover.modsel))
names(with.selection) <- c("Bias.N", "Bias.S", "Bias.Tot",
                "CICover.N", "CICover.S", "CICover.Tot")
row.names(with.selection) <- delta.multiplier
kable(with.selection, digits=c(rep(1,3), rep(2,3)), 
      caption="Results using chosen model")


## ----matplot, echo=FALSE-------------------------------
matplot(delta.multiplier, pctbias, type="b", lwd=3, pch=20, xlab="Delta", 
        ylab="Percent bias", col=1, main="Bias in estimated density")
legend("bottom", lty=1:3, lwd=3, legend=c("North", "South", "Total"))
abline(h=0)

matplot(delta.multiplier, ci.cover, type="b", lwd=3, pch=20, xlab="Delta", 
        ylim=c(0.75,1),
        ylab="Confidence interval coverage", col=1, main="Confidence interval coverage")
legend("bottom", lty=1:3, lwd=3, legend=c("North", "South", "Total"))
abline(h=0.95)


## ----modsel--------------------------------------------
plot(delta.multiplier, modelsel, 
     main="Stratum-specific model chosen", type="b", pch=20,
     xlab="Delta", ylab="Stratum covariate chosen")
abline(h=0.50)


## ----poorerdesign, echo=FALSE--------------------------
poorer.design <- make.design(region = northsea,
                           transect.type = "line",
                           design = "systematic",
                           samplers=25, 
                           design.angle = c(50, 40),
                           truncation = 0.8,
                           coverage.grid = coverage.grid)
poorer.equalcover <- list()
poorer.whichmodel <- list()
num.sims <- 300
for (i in seq_along(delta.multiplier)) {
  sigma.strata <- c(sigma.south, sigma.south*delta.multiplier[i])
  detect <- make.detectability(key.function = "hn",
                               scale.param = sigma.strata,
                               truncation = 0.8)
  poorer.equalcover.sim <- make.simulation(reps = num.sims,
                                    design = poorer.design,
                                    population.description = minkepop,
                                    detectability = detect,
                                    ds.analysis = pooled.hn)
  poorer.whichmodel.sim <- make.simulation(reps = num.sims,
                                    design = poorer.design,
                                    population.description = minkepop,
                                    detectability = detect,
                                    ds.analysis = strat.specific.or.not)
  poorer.equalcover[[i]] <- run.simulation(poorer.equalcover.sim, 
                                           run.parallel = TRUE, max.cores=6)
  poorer.whichmodel[[i]] <- run.simulation(poorer.whichmodel.sim, 
                                           run.parallel = TRUE, max.cores=6)
}

ppctbias <- mat.or.vec(length(delta.multiplier),3)
pci.cover <- mat.or.vec(length(delta.multiplier),3)
pnum.detects <- mat.or.vec(length(delta.multiplier),3)
pmodelsel <- mat.or.vec(length(delta.multiplier),1)
pbias.modsel <- mat.or.vec(length(delta.multiplier),3)
pcover.modsel <- mat.or.vec(length(delta.multiplier),3)
for (i in seq_along(delta.multiplier)) {
  simsum <- summary(poorer.equalcover[[i]], description.summary=FALSE)
  selsum <- summary(poorer.whichmodel[[i]], description.summary=FALSE)
  ppctbias[i, ] <- simsum@individuals$D[,3]
  pci.cover[i, ] <- simsum@individuals$D[,5]
  pnum.detects[i, ] <- simsum@individuals$summary[,3]
  pbias.modsel[i, ] <- selsum@individuals$D[,3]
  pcover.modsel[i, ] <- selsum@individuals$D[,5]
  pmodelsel[i] <- sum(poorer.whichmodel[[i]]@results$Detection[1,5, 1:num.sims]==2) / num.sims
}
pooled <- as.data.frame(cbind(ppctbias, pci.cover, pnum.detects, pmodelsel))
names(pooled) <- c("Bias.N", "Bias.S", "Bias.Tot",
                "CICover.N", "CICover.S", "CICover.Tot",
                "Detects.N", "Detects.S", "Detects.Tot", "Prop.strat.model")
row.names(pooled) <- delta.multiplier
kable(pooled, digits=c(rep(1,3), rep(2,3), rep(0,3), 2), 
      caption="Results using pooled detection function with lower sampling intensity")


## ----poorerselection, echo=FALSE-----------------------
pwith.selection <- as.data.frame(cbind(pbias.modsel, pcover.modsel))
names(pwith.selection) <- c("Bias.N", "Bias.S", "Bias.Tot",
                "CICover.N", "CICover.S", "CICover.Tot")
row.names(pwith.selection) <- delta.multiplier
kable(pwith.selection, digits=c(rep(1,3), rep(2,3)), 
      caption="Results using chosen model with lower sampling intensity")


## ----poorerplots---------------------------------------
matplot(delta.multiplier, ppctbias, type="b", lwd=3, pch=20, xlab="Delta", 
        ylab="Percent bias", col=1, 
        main="Bias in estimated density\npoorer design")
legend("bottom", lty=1:3, lwd=3, legend=c("North", "South", "Total"))
abline(h=0)

matplot(delta.multiplier, pci.cover, type="b", lwd=3, pch=20, xlab="Delta", 
        ylim=c(0.75,1),
        ylab="Confidence interval coverage", col=1, 
        main="Confidence interval coverage\npoorer design")
legend("bottom", lty=1:3, lwd=3, legend=c("North", "South", "Total"))
abline(h=0.95)
plot(delta.multiplier, pmodelsel, 
     main="Stratum-specific model chosen\npoorer design", type="b", pch=20,
     xlab="Delta", ylab="Stratum covariate chosen")
abline(h=0.50)


## ----hissfn, echo=FALSE--------------------------------
hiss <- function(nsims, thelist, index, deltmult) {
  mm <- as.data.frame(t(thelist[[index]]@results$Detection[1,5:6,1:nsims]))
  mm$DeltaCriteria[mm$SelectedModel==1] <- mm$DeltaCriteria[mm$SelectedModel==1] * -1
  nonmissing <- sum(!is.na(mm$DeltaCriteria))
  prop.wrong <- round(sum(mm$DeltaCriteria<0,na.rm=TRUE)/nonmissing,3)
  prop.ambig <- round(sum(mm$DeltaCriteria>0 & mm$DeltaCriteria<2,na.rm=TRUE)/nonmissing,3)
  prop.right <- round(sum(mm$DeltaCriteria>2,na.rm=TRUE)/nonmissing,3)
  hist(mm$DeltaCriteria, main=paste("Delta mult on sigma =", deltmult), 
       xlab=expression(paste(Delta,"AIC")),
       sub=paste("strat=", prop.right, "ambig=", prop.ambig, "common=", prop.wrong),
       nc=15)
  abline(v=2, lwd=2); abline(v=0, lwd=2)
}


## ----hissgood, fig.cap="Distribution of deltaAIC with 40 transects in the study area.", layout="l-body-outset", fig.width=6, fig.height=6.5----
par(mfrow=c(4,2))
for (i in seq_along(delta.multiplier)) {
  hiss(num.sims, whichmodel, i, delta.multiplier[i])
}


## ----hissbad, fig.cap="Distribution of deltaAIC with 25 transects in the study area.",layout="l-body-outset", fig.width=6, fig.height=6.5----
par(mfrow=c(4,2))
for (i in seq_along(delta.multiplier)) {
  hiss(num.sims, poorer.whichmodel, i, delta.multiplier[i])
}

