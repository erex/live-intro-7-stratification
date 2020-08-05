## ---- message=FALSE--------------------------------------------------------------------
library(Distance)
library(kableExtra)
# Load data
data(minke)
head(minke, n=3)
# Specify truncation distance
minke.trunc <- 1.5


## ---- echo=T, eval=T, message=FALSE----------------------------------------------------
## Fit to each region separately - full geographical stratification
# Create data set for South
minke.S <- minke[minke$Region.Label=="South", ]
minke.df.S.strat <- ds(minke.S, truncation=minke.trunc, key="hr", adjustment=NULL)
summary(minke.df.S.strat)
# Combine selection and detection function fitting for North
minke.df.N.strat <- ds(minke[minke$Region.Label=="North", ],
                       truncation=minke.trunc, key="hr", adjustment=NULL)
summary(minke.df.N.strat)


## ---- echo=T, eval=T, message=FALSE----------------------------------------------------
minke.df.all <- ds(minke, truncation=minke.trunc, key="hr", adjustment=NULL)
summary(minke.df.all)


## ---- echo=F, eval=T-------------------------------------------------------------------
# Save AIC values
aic.all <- summary(minke.df.all$ddf)$aic
aic.S <- summary(minke.df.S.strat$ddf)$aic
aic.N <- summary(minke.df.N.strat$ddf)$aic
aic.SN <- aic.S + aic.N



## ---- echo=F---------------------------------------------------------------------------
par(mfrow=c(2,2))
plot(minke.df.S.strat, main="South")
plot(minke.df.N.strat, main="North")
plot(minke.df.all, main="All")


## ---- echo=F---------------------------------------------------------------------------
# Harvest abundance estimates
est.full <- data.frame(Label=c("North","South","Total"),Estimate=rep(NA,3))
est.full[1,2] <- minke.df.N.strat$dht$individuals$N$Estimate
est.full[2,2] <- minke.df.S.strat$dht$individuals$N$Estimate
est.full[3,2] <- est.full[1,2] + est.full[2,2]

est.dfpool <- minke.df.all$dht$individuals$N[ ,c(1,2)]


## ---- echo=FALSE-----------------------------------------------------------------------
knitr::kable(est.full, caption="Abundance estimates using full geographical stratification.", digits=0) %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE)


## ---- echo=FALSE-----------------------------------------------------------------------
knitr::kable(est.dfpool, caption="Abundance estimates calculating encounter rate by strata and a pooled detection function.", digits=0) %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE)


## ---- echo=T, eval=T, message=FALSE----------------------------------------------------
# Geographical stratification with stratum-specific detection function 
strat.specific.detfn <- ds(data=minke, truncation=minke.trunc, key="hr", 
                           adjustment=NULL, formula=~Region.Label)
abund.by.strata <- dht2(ddf=strat.specific.detfn, flatfile=minke, 
                        strat_formula=~Region.Label, stratification="geographical")
print(abund.by.strata, report="abundance")

