#-----------------------------------------------------------------
# Assis tree growth analysis
# Analysis by Sam Flake, swflake@ncsu.edu
# Abstract of methods
# License


#-----------------------------------------------------------------
# Load libraries
set.seed(129021)

library("plyr")
library("effects")
library("lme4")
library("MuMIn")
library("lmerTest")
library("ggpubr")
library("reshape2")
library("ggplot2")
library("ResourceSelection")

#-----------------------------------------------------------------
# Helper functions
is.blank <- function(x){
  #test if a cell is blank
  return(ifelse(x == "", T, F))
}

make_code <- function(species){
  #take the first two letters from genus and species and make a 4-letter code
  species <- unlist(species)
  first <- substr(species[1], start = 1, stop = 2)
  second <- substr(species[2], start = 1, stop = 2)
  return(toupper(paste0(first, second)))
}

get_heights_csas <- function(plot_no, year){
  # retrieve tree numbers, heights, and cross-sectional areas of all trees in a plot. 
  # tell it what year (2006, 2011, 2016) to retrieve data from that year. 
  # Yes this is hacky
  plot <- tree_data[tree_data$Plot == plot_no, ]
  if(year == 2006) return(plot[plot$Alive2006 == "v", c("Current.number", "H.2006", "G2006..m2.")])
  if(year == 2011) return(plot[plot$Alive2011 == "v", c("Current.number", "H.2011..m.", "G2011..m2.")])
}


vif.mer <- function (fit) {
  ## adapted from rms::vif
  ## borrowed from https://raw.githubusercontent.com/aufrank/R-hacks/master/mer-utils.R
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

addTrans <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}

logit <- function(p){
  log(p/(1-p))
}

inv.logit <- function(x){
  1/(1 + exp(-x))
}

# function to find the nearest y-value for an x data point, to add to partial residuals
# for plotting. Stolen from effects package
closest <- function(x, x0){
  apply(outer(x, x0, FUN = function(x, x0) abs(x - x0)), 1, which.min)
}

#-----------------------------------------------------------------
# Data prep
tree_data <- read.csv("./Giselda_data_per_tree_sf_edit.csv", stringsAsFactors = FALSE)
plot_data <- read.csv("./Giselda_Data_per_plot.csv")[1:30, ]

tree_data <- tree_data[!(tree_data$Plot %in% c(21:25)), ] #remove plots with incomplete data

#remove the "forest" functional type
tree_data[tree_data$C.Cerrado..G.generalist. == "F", ] <- "G"

# replace missing data with NA
# tree_data[tree_data$dbh.annual.increment.cm._10.years == "na", "dbh.annual.increment.cm._10.years"] <- NA
blanks <- apply(tree_data, c(1,2), is.blank)
tree_data[blanks] <- NA

tree_data$dbh.annual.increment.cm._10.years <- as.numeric(as.character(tree_data$dbh.annual.increment.cm._10.years))

# create a code for species names
# split up the species names
species <- strsplit(tree_data$Species..names.not.updated., split = " ")
# add 4-letter code to tree dataframe
tree_data$Code <- as.factor(unlist(lapply(species, make_code)))

#fix up variable types
tree_data$Shade.Tolerance <- as.factor(tree_data$Shade.Tolerance)
tree_data$C.Cerrado..G.generalist. <- as.factor(tree_data$C.Cerrado..G.generalist.)
tree_data$H.2006 <- as.numeric(tree_data$H.2006)
tree_data$H.2011..m. <- as.numeric(tree_data$H.2011..m.)
tree_data$G2006..m2. <- as.numeric(tree_data$G2006..m2.)
tree_data$G2011..m2. <- as.numeric(tree_data$G2011..m2.)
tree_data$dg.2006.cm <- as.numeric(tree_data$dg.2006.cm)
tree_data$dg.2016.cm <- as.numeric(tree_data$dg.2016.cm)
tree_data$dg.2011.cm <- as.numeric(tree_data$dg.2011.cm)


tree_data$bai <- tree_data$dg.2016.cm - tree_data$dg.2006.cm
# Calculate the one-sided competition index for 2006 and for 2011
# this uses the entire dataset, not removing trees with negative growth 
# and using only live trees

for(i in 1:nrow(tree_data)){
  tree_no <- tree_data[i, "Current.number"]
  tree_h <- tree_data[i, "H.2006"]
  plot_no <- tree_data[i, "Plot"]
  trees <- get_heights_csas(plot_no, 2006)
  
  #BA.above is the sum of basal area of taller trees
  tree_data$BA.above2006[i] <- sum(trees[trees$H.2006 > tree_h, "G2006..m2."], na.rm = TRUE)
}

for(i in 1:nrow(tree_data)){
  tree_no <- tree_data[i, "Current.number"]
  tree_h <- as.numeric(tree_data[i, "H.2011..m."])
  plot_no <- tree_data[i, "Plot"]
  trees <- get_heights_csas(plot_no, 2011)
  
  #BA.above is the sum of basal area of taller trees
  tree_data$BA.above2011[i] <- sum(trees[trees$H.2011..m. > tree_h, "G2011..m2."], na.rm = TRUE)
}


tree_data$Died2011 <- as.factor(ifelse(tree_data$Alive2006 == "v" & tree_data$Alive2011 == "m",
                                       "y", "n"))
tree_data$Died2016 <- as.factor(ifelse(tree_data$Alive2011 == "v" & tree_data$Alive2016 == "m",
                                       "y", "n"))


# recalculate growth rates
# calculate growth rate for ingrowths
for(i in 1:nrow(tree_data)){
  row <- tree_data[i, ]
  if(row$Added == 2006){
    tree_data[i, "dbh.annual.increment.cm._10.years"] <- (as.numeric(row$dg.2016.cm) - as.numeric(row$dg.2006.cm))/10
  }
  if(row$Added == 2011){
    tree_data[i, "dbh.annual.increment.cm._10.years"] <- (as.numeric(row$dg.2016.cm) - as.numeric(row$dg.2011.cm))/5
  }
}


t.test(tree_data[tree_data$Added == "2006", "dbh.annual.increment.cm._10.years"],
         tree_data[tree_data$Added == "2011", "dbh.annual.increment.cm._10.years"])

#remove negative growth
tree_data_pos <- subset(tree_data, dbh.annual.increment.cm._10.years > 0)

#only include trees originally in the plot
tree_data_pos <- subset(tree_data_pos, Added == 2006)

#only trees with 2006 diameter
tree_data_pos <- subset(tree_data_pos, !is.na(dg.2006.cm))
tree_data_pos <- subset(tree_data_pos, dg.2006.cm > 0)

# PCA on soils

#sand fraction was calculated weirdly
plot_data$Total_sand <- plot_data$Coarse.sand.. + plot_data$Fine.sand..
pca_data <- plot_data[, c(5,6,7,8,9,13,17,19)]

pca <- princomp(scale(pca_data))
summary(pca)
biplot(pca, var.axes = TRUE)
pca$loadings
pca$scores

plot_data$PC1 <- pca$scores[, 1]
plot_data$PC2 <- pca$scores[, 2]



#merge all data together
all_data_pos <- join(tree_data_pos, plot_data, by = c("Plot"))
all_data<- join(tree_data, plot_data, by = c("Plot"))


sp_list <- all_data[which(!duplicated(all_data$Species..names.not.updated.)), ]


#----------------------------------------------------------------------
# Model comparison for table 1
#----------------------------------------------------------------------
mod1 <- lmer(log(dbh.annual.increment.cm._10.years) ~ 
               (1|Code),
             data = all_data_pos)
  summary(mod1)
  r.squaredGLMM(mod1)
  AIC(mod1)
  
mod2 <- lmer(log(dbh.annual.increment.cm._10.years) ~ scale(log(dg.2006.cm)) + 
               (1|Code),
             data = all_data_pos)
  vif.mer(mod2)
  summary(mod2)
  r.squaredGLMM(mod2)
  AICc(mod2)
  
mod3 <- lmer(log(dbh.annual.increment.cm._10.years) ~ scale(PC1)  + scale(PC2) + 
               (1|Code),
             data = all_data_pos)
  vif.mer(mod3)
  summary(mod3)
  r.squaredGLMM(mod3)
  AICc(mod3)  
  
mod4 <- lmer(log(dbh.annual.increment.cm._10.years) ~ scale(PC1)  +
               (1|Code),
             data = all_data_pos)
vif.mer(mod4)
summary(mod4)
r.squaredGLMM(mod4)
AICc(mod4)  
  
mod5 <- lmer(log(dbh.annual.increment.cm._10.years) ~ scale(log(dg.2006.cm)) +
               scale(PC1) +  
               (1|Code),
             data = all_data_pos)
  vif.mer(mod5)
  summary(mod5)
  r.squaredGLMM(mod5)
  AICc(mod5)  
  
mod6 <- lmer(log(dbh.annual.increment.cm._10.years) ~ scale(log(dg.2006.cm)) +
               scale(PC1)  +  C.Cerrado..G.generalist. + 
               (1|Code),
             data = all_data_pos)
  vif.mer(mod6)
  summary(mod6)
  r.squaredGLMM(mod6)
  AICc(mod6)  
  
mod7 <- lmer(log(dbh.annual.increment.cm._10.years) ~ scale(log(dg.2006.cm)) +
               scale(PC1)  + C.Cerrado..G.generalist.*TBA.2006.m2ha.1 + 
               (1|Code),
             data = all_data_pos)
  vif.mer(mod7)
  summary(mod7)
  r.squaredGLMM(mod7)
  AICc(mod7)  
  
mod8 <- lmer(log(dbh.annual.increment.cm._10.years) ~ scale(log(dg.2006.cm)) +
               scale(PC1)  +  C.Cerrado..G.generalist.*Canopy.cover.. + 
               (1|Code),
             data = all_data_pos)
  vif.mer(mod8)
  summary(mod8)
  r.squaredGLMM(mod8)
  AICc(mod8)  
  
mod9 <- lmer(log(dbh.annual.increment.cm._10.years) ~ scale(log(dg.2006.cm)) +
               scale(PC1)  +  scale(BA.above2006)*C.Cerrado..G.generalist. + 
               (1|Code),
             data = all_data_pos)
  vif.mer(mod9)
  summary(mod9)
  r.squaredGLMM(mod9)
  AICc(mod9)  
  
mod10 <- lmer(log(dbh.annual.increment.cm._10.years) ~ scale(log(dg.2006.cm))*C.Cerrado..G.generalist. +
               scale(PC1)  +  scale(BA.above2006)*C.Cerrado..G.generalist. + 
               (1|Code),
             data = all_data_pos)
  vif.mer(mod10)
  summary(mod10)
  r.squaredGLMM(mod10)
  AICc(mod10)  
  
mod11 <- lmer(log(bai) ~ scale(log(dg.2006.cm)) +
                scale(PC1) + scale(BA.above2011)*C.Cerrado..G.generalist. + 
                (1|Code),
              data = all_data_pos)
  vif.mer(mod11)
  summary(mod11)
  r.squaredGLMM(mod11)
  AICc(mod11)  
  
mod12 <- lmer(log(dbh.annual.increment.cm._10.years) ~ scale(log(dg.2006.cm))*C.Cerrado..G.generalist. +
                scale(PC1)  + scale(BA.above2011)*C.Cerrado..G.generalist. + 
                (1|Code),
              data = all_data_pos)
  vif.mer(mod12)
  summary(mod12)
  r.squaredGLMM(mod12)
  AICc(mod12)  


#----------------------------------------------------------------------
# Full model: FT x CI interaction with soils

ft_ci_soils_lmm <- lmer(log(dbh.annual.increment.cm._10.years) ~ scale(log(dg.2006.cm)) +
                          scale(PC1) + scale(BA.above2011)*C.Cerrado..G.generalist. + 
                          (1|Code),
                        data = all_data_pos)

vif.mer(ft_ci_soils_lmm)
summary(ft_ci_soils_lmm)
r.squaredGLMM(ft_ci_soils_lmm)
AIC(ft_ci_soils_lmm)

plot(residuals(ft_ci_soils_lmm) ~ fitted(ft_ci_soils_lmm))
abline(h = 0)
ggqqplot(residuals(ft_ci_soils_lmm))
plot(allEffects(ft_ci_soils_lmm, partial.residuals = FALSE))
plot(Effect(focal.predictors = c("BA.above2011", "C.Cerrado..G.generalist."), 
            mod = ft_ci_soils_lmm, transformation = list(link = log, inverse = exp)))
plot(Effect(focal.predictors = c("dg.2006.cm", "C.Cerrado..G.generalist."), 
            mod = ft_ci_soils_lmm, transformation = list(link = log, inverse = exp)))
plot(Effect(focal.predictors = c("PC1"), 
            mod = ft_ci_soils_lmm, transformation = list(link = log, inverse = exp)))
plot(Effect(focal.predictors = c("PC2"), 
            mod = ft_ci_soils_lmm, transformation = list(link = log, inverse = exp)))


derivs1 <- ft_ci_soils_lmm@optinfo$derivs
sc_grad1 <- with(derivs1,solve(Hessian,gradient))
max(abs(sc_grad1))

max(pmin(abs(sc_grad1),abs(derivs1$gradient)))

relgrad <- with(ft_ci_soils_lmm@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))

library("numDeriv")
dd <- update(ft_ci_soils_lmm,devFunOnly=TRUE)
pars <- unlist(getME(ft_ci_soils_lmm,c("theta","fixef")))
grad2 <- grad(dd,pars, method = "simple")
hess2 <- hessian(dd,pars)
sc_grad2 <- solve(hess2,grad2)
max(pmin(abs(sc_grad2),abs(grad2)))


library(lattice)
qqmath(ranef(ft_ci_soils_lmm, condVar = TRUE), strip = FALSE)$Plot
qqmath(ranef(ft_ci_soils_lmm, condVar = TRUE), strip = FALSE)$Code


## Species growth rates (appendix table 1)
ranef_sp <- exp(ranef(ft_ci_soils_lmm, condVar = TRUE)$Code)
ranef_sp_df <- data.frame(rel_growth = ranef_sp$`(Intercept)`,
                          Code = rownames(ranef_sp))

ranef_sp_df <- join(ranef_sp_df, 
                    tree_data[, c("Code",
                                  "Species..names.not.updated.",
                                  "Family",
                                  "C.Cerrado..G.generalist.",
                                  "Shade.Tolerance")],
                    by = c("Code"),
                    type = "left",
                    match = "first")

ranef_sp_df <- ranef_sp_df[order(ranef_sp_df$rel_growth), ]

sp_table <- as.data.frame(table(all_data_pos$Code))
names(sp_table) <- c("Code", "Freq")
ranef_sp_df <- join(ranef_sp_df, sp_table, by = c("Code"))

write.csv(ranef_sp_df, "species_growth_rates.csv")



#------------------------------------------------------------------------------
# Subset with detailed soils data
#------------------------------------------------------------------------------
rain <- read.csv("./rain.csv")
soils <- read.csv("./soils.csv")
soils$Horizon2 <- ifelse(grepl("A", soils$Horizon), "Upper", "Lower")
soils2 <- cbind(soils[soils$Horizon2 == "Upper", ], soils[soils$Horizon2 == "Lower", ])
names(soils2)[24:46] <- paste("Lower", names(soils2)[24:46])

plot_data_sub <- plot_data[plot_data$Plot %in% soils2$Plot, ]
plot_data_sub <- join(plot_data_sub, soils2, by = "Plot")
plot_data_sub <- join(plot_data_sub, rain, by = "Plot")


pca_sub <- princomp(scale(plot_data_sub[, c(5,6,7,8,9,13,31,32,33,34,54,55,56,57)]))
summary(pca_sub)
plot_data_sub$PC1_sub <- pca_sub$scores[, 1]
plot_data_sub$PC2_sub <- pca_sub$scores[, 2]

all_data_sub <- join(plot_data_sub, tree_data_pos, by = "Plot")

soils_full_lmm <- lmer(log(dbh.annual.increment.cm._10.years) ~ scale(log(dg.2006.cm)) +
                         scale(BA.above2011)*C.Cerrado..G.generalist. + scale(PC1_sub) + scale(PC2_sub) + 
                         scale(Net.rainfall..) +
                         (1|Code),
                       data = all_data_sub)
summary(soils_full_lmm)
AICc(soils_full_lmm)
vif.mer(soils_full_lmm)

soils_full_ba_lmm <- lmer(log(dbh.annual.increment.cm._10.years) ~ scale(log(dg.2006.cm)) +
                         scale(BA.above2011)*C.Cerrado..G.generalist. + scale(PC1_sub) + scale(PC2_sub) + 
                         scale(TBA.2011.m2.ha.1) +
                         (1|Code),
                       data = all_data_sub)
summary(soils_full_ba_lmm)
AICc(soils_full_ba_lmm)
vif.mer(soils_full_ba_lmm)


plot(Effect(focal.predictors = c("Net.rainfall.."), 
            mod = soils_full_lmm, transformation = list(link = log, inverse = exp)), 
     partial.residuals = TRUE)

#exploration
summary(lm(log(all_data_sub$dbh.annual.increment.cm._10.years) ~ 
             all_data_sub$Net.rainfall..))


#------------------------------------------------------------------------------
# Mortality modeling
#------------------------------------------------------------------------------
library("pROC")


mort_ft_ci_soils_glm_2011 <- glm(Died2011 ~ scale(log(dg.2006.cm)) +
                                   scale(BA.above2006)*C.Cerrado..G.generalist.,
                                 data = all_data[all_data$dg.2006.cm > 0, ],
                                 family = binomial(link = "logit"))

mort_ft_ci_soils_glm_2016 <- glm(Died2016 ~ scale(log(dg.2011.cm)) +
                                   scale(BA.above2011)*C.Cerrado..G.generalist.,
                                 data = all_data[all_data$dg.2011.cm > 0, ],
                                 family = binomial(link = "logit"))


# vif.mer(mort_ft_ci_soils_glm_2011)

summary(mort_ft_ci_soils_glm_2011)
r.squaredGLMM(mort_ft_ci_soils_glm_2011)
AIC(mort_ft_ci_soils_glm_2011)
plot.roc(mort_ft_ci_soils_glm_2011$y, 
         fitted(mort_ft_ci_soils_glm_2011),print.auc = TRUE, 
         col = "green", lty = 2)
plot(Effect(mort_ft_ci_soils_glm_2011, 
            focal.predictors = c("BA.above2006", "C.Cerrado..G.generalist.")))


summary(mort_ft_ci_soils_glm_2016)
r.squaredGLMM(mort_ft_ci_soils_glm_2016)
AIC(mort_ft_ci_soils_glm_2016)
plot.roc(mort_ft_ci_soils_glm_2016$y, 
         fitted(mort_ft_ci_soils_glm_2016),print.auc = TRUE, 
         col = "green", lty = 2)
plot(Effect(mort_ft_ci_soils_glm_2016, 
            focal.predictors = c("BA.above2011", "C.Cerrado..G.generalist.")))



#hosmer-lemeshow test seems pretty okay!
hl1 <- hoslem.test(mort_ft_ci_soils_glm_2011$y, 
                   fitted(mort_ft_ci_soils_glm_2011), g=10)
hl2 <- hoslem.test(mort_ft_ci_soils_glm_2016$y, 
                   fitted(mort_ft_ci_soils_glm_2016), g=10)


#-----------------------------------------------------
# Effects plots
#-----------------------------------------------------

#------------------------------------------------------------------------------
# Figure 1
# Overall plot-level results
#------------------------------------------------------------------------------
tiff(filename="./plots/Figure 1 overall growth basal area.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 3, 
     height=3, 
     pointsize=7, 
     res=600)

par(mar = c(5.1, 5.1, 2, 1))

plot(plot_data$Diameter.increment.cm.yr ~ plot_data$TBA.2006.m2ha.1, 
     pch = 21,
     bg = "darkgray",
     xlab = "",
     ylab = "",
     xaxt = "n",
     yaxt = "n")

abline(coef(lm(plot_data$Diameter.increment.cm.yr ~ plot_data$TBA.2006.m2ha.1)))

mtext(side = 1, text = expression(paste("Stand basal area (m"^"2"," ha"^"-1",")")),
      line = 2.8, cex = 1.6)
mtext(side = 2, text = expression(paste("Mean diameter increment (cm year"^"-1", ")")), 
      line = 2.8, cex = 1.6)

axis(side =1, cex.axis = 1.5)
axis(side =2, cex.axis = 1.5)

dev.off()




#------------------------------------------------------------------------------
# Figure 2
# Growth ~ soil variables
#------------------------------------------------------------------------------
tiff(filename="./plots/Figure 2 soil variables effect plot with partial residuals.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 3, 
     height=3, 
     pointsize=7, 
     res=600)

par(mar = c(5.1, 5.1, 2, 1))

eff <- Effect(focal.predictors = c("PC1"), 
              mod = ft_ci_soils_lmm, xlevels = 100, transformation = list(link = log, inverse = exp),
              partial.residuals = TRUE)

y <- eff$fit
x <- eff[["x"]]$PC1
x.fit <- eff[["x.all"]]$PC1

fitted <- y[closest(x.fit, x)]
resids <- eff$residuals + fitted
exp_resids <- exp(resids)

resid_df <- data.frame(x.fit,
                       resids)

e1_r_agg <- aggregate(resid_df, by = list(resid_df$x.fit), FUN = mean)[, c(2,3)]
e1_r_agg[, 3] <- aggregate(resid_df, by = list(resid_df$x.fit), FUN = sd)[, 3]
names(e1_r_agg) <- c("x.fit", "mean_part_resids", "sd_resids")
e1_r_agg$low <- e1_r_agg$mean_part_resids - e1_r_agg$sd_resids
e1_r_agg$high <- e1_r_agg$mean_part_resids + e1_r_agg$sd_resids
e1_r_agg$exp_low <- exp(e1_r_agg$mean_part_resids - e1_r_agg$sd_resids)
e1_r_agg$exp_high <- exp(e1_r_agg$mean_part_resids + e1_r_agg$sd_resids)
e1_r_agg$exp_mean_resids <- exp(e1_r_agg[, 2])


test <- plot(eff)


eff2 <- Effect(focal.predictors = c("PC1"), 
               mod = mod3, xlevels = 100, transformation = list(link = log, inverse = exp))

dat <- data.frame(y = exp(eff$fit),
                  lower = exp(eff$lower),
                  upper = exp(eff$upper),
                  soil = eff[["x"]][["PC1"]])

dat2 <- data.frame(y = exp(eff2$fit),
                   lower = exp(eff2$lower),
                   upper = exp(eff2$upper),
                   soil = eff2[["x"]][["PC1"]])

plot(x = NA, 
     xlim = c(-3,3.4),
     ylim = c(0,0.6),
     xlab = "",
     ylab = "",
     xaxt = "n", 
     yaxt = "n")

points(e1_r_agg$exp_mean_resids ~ e1_r_agg$x.fit, 
       pch = 21, 
       bg= "darkblue",
       col = "darkblue")
arrows(y0 = e1_r_agg$exp_low, x0 = e1_r_agg$x.fit, 
       y1 = e1_r_agg$exp_high, x1 = e1_r_agg$x.fit,
       code = 3,
       angle = 90,
       length = 0.02,
       lwd = 0.5)

lines(dat$y ~ dat$soil, col = "darkblue", lwd = 2)
polygon(c(dat$soil, rev(dat$soil)), c(dat$upper, rev(dat$lower)),
        col = addTrans("darkblue",30), border = NA)
# points(exp_resids ~ x.fit)


lines(dat2$y ~ dat2$soil, col = "magenta", lwd = 2)
polygon(c(dat2$soil, rev(dat2$soil)), c(dat2$upper, rev(dat2$lower)),
        col = addTrans("magenta",30), border = NA)

mtext(side = 1, text = "Soil PC1", line = 2.7, cex = 1.8)
mtext(side = 2, text = expression(paste("Diameter increment (cm year"^"-1", ")")), 
      line = 2.7, cex = 1.8)

axis(side =1, cex.axis = 1.5)
axis(side =2, cex.axis = 1.5)

legend("topleft", legend = c("Soils Alone", "Full Model"), col = c("magenta", "darkblue"),
       lty = 1)

dev.off()


tiff(filename="./plots/Figure 2 soil variables effect plot.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 3, 
     height=3, 
     pointsize=7, 
     res=600)

par(mar = c(5.1, 5.1, 2, 1))

eff <- Effect(focal.predictors = c("PC1"), 
              mod = ft_ci_soils_lmm, xlevels = 100, transformation = list(link = log, inverse = exp),
              partial.residuals = TRUE)

eff2 <- Effect(focal.predictors = c("PC1"), 
               mod = mod3, xlevels = 100, transformation = list(link = log, inverse = exp))

dat <- data.frame(y = exp(eff$fit),
                  lower = exp(eff$lower),
                  upper = exp(eff$upper),
                  soil = eff[["x"]][["PC1"]])

dat2 <- data.frame(y = exp(eff2$fit),
                   lower = exp(eff2$lower),
                   upper = exp(eff2$upper),
                   soil = eff2[["x"]][["PC1"]])

plot(x = NA, 
     xlim = c(-3,3.4),
     ylim = c(0.08, 0.25),
     xlab = "",
     ylab = "",
     xaxt = "n", 
     yaxt = "n")


lines(dat$y ~ dat$soil, col = "darkblue", lwd = 2)
polygon(c(dat$soil, rev(dat$soil)), c(dat$upper, rev(dat$lower)),
        col = addTrans("darkblue",30), border = NA)
# points(exp_resids ~ x.fit)


lines(dat2$y ~ dat2$soil, col = "magenta", lwd = 2)
polygon(c(dat2$soil, rev(dat2$soil)), c(dat2$upper, rev(dat2$lower)),
        col = addTrans("magenta",30), border = NA)

mtext(side = 1, text = "Soil PC1", line = 2.7, cex = 1.8)
mtext(side = 2, text = expression(paste("Diameter increment (cm year"^"-1", ")")), 
      line = 2.7, cex = 1.8)

axis(side =1, cex.axis = 1.5)
axis(side =2, cex.axis = 1.5)

legend("topleft", legend = c("Soils Alone", "Full Model"), col = c("magenta", "darkblue"),
       lty = 1)

dev.off()

#------------------------------------------------------------------------------
# Figure 3
# Growth ~ diam by FT
#------------------------------------------------------------------------------

eff <- Effect(focal.predictors = c("dg.2006.cm", "C.Cerrado..G.generalist."), 
              mod = ft_ci_soils_lmm, xlevels = 100, transformation = list(link = log, inverse = exp),
              partial.residuals = TRUE)

y <- eff$fit
x <- eff[["x"]]$"dg.2006.cm"
x.fit <- eff[["x.all"]]$"dg.2006.cm"
fg <- eff[["x.all"]]$"C.Cerrado..G.generalist."

data_orig <- data.frame(x.fit, 
                        fg)

#make sure we match the correct line
fitted <- ifelse(data_orig$fg == "C", y[1:100][closest(x.fit, x[1:100])],
                 y[101:200][closest(x.fit, x[101:200])])
  
resids <- eff$residuals + fitted
exp_resids <- exp(resids)



tiff(filename="./plots/Figure 3 increment_diam_by_fg_with_partial_residuals.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 5, 
     height=5, 
     pointsize=7, 
     res=600)

dat <- data.frame(y = exp(eff$fit),
                  lower = exp(eff$lower),
                  upper = exp(eff$upper),
                  diam = eff[["x"]][["dg.2006.cm"]],
                  fg = eff[["x"]][["C.Cerrado..G.generalist."]])

par(mar = c(5.1, 5.1, 2, 1))

plot(x = NA, 
     xlim = c(5, 60),
     ylim = c(0.001, 1.2),
     xlab = "",
     ylab = "",
     xaxt = "n", 
     yaxt = "n",
     log = "xy")

points(exp_resids ~ x.fit, 
       col = ifelse(fg == "G", addTrans("#1b9e77",50), addTrans("#d95f02",50)),
       bg = ifelse(fg == "G", addTrans("#1b9e77",50), addTrans("#d95f02",50)),
       pch = 21,
       cex = 0.5)

dat_sub <- dat[dat$fg == "G", ]
lines(dat_sub$y ~ dat_sub$diam, col = "#1b9e77", lwd = 2)
polygon(c(dat_sub$diam, rev(dat_sub$diam)), c(dat_sub$upper, rev(dat_sub$lower)),
        col = addTrans("#1b9e77",30), border = NA)

dat_sub <- dat[dat$fg == "C", ]
lines(dat_sub$y ~ dat_sub$diam, col = "#d95f02", lwd = 2)
polygon(c(dat_sub$diam, rev(dat_sub$diam)), c(dat_sub$upper, rev(dat_sub$lower)),
        col = addTrans("#d95f02",30), border = NA)

mtext(side = 1, text = "Diameter (cm)", line = 2.7, cex = 1.8)
mtext(side = 2, text = expression(paste("Diameter increment (cm year"^"-1", ")")), 
      line = 2.7, cex = 1.8)

axis(side =1, cex.axis = 1.5)
axis(side =2, cex.axis = 1.5)

legend("topright", legend = c("Forest species", "Savanna species"), col = c("#1b9e77", "#d95f02"),
       lty = 1)

dev.off()

tiff(filename="./plots/Figure 3 increment_diam_by_fg.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 3, 
     height=3, 
     pointsize=7, 
     res=600)

dat <- data.frame(y = exp(eff$fit),
                  lower = exp(eff$lower),
                  upper = exp(eff$upper),
                  diam = eff[["x"]][["dg.2006.cm"]],
                  fg = eff[["x"]][["C.Cerrado..G.generalist."]])

par(mar = c(5.1, 5.1, 2, 1))

plot(x = NA, 
     xlim = c(5, 60),
     ylim = c(0.05, 0.25),
     xlab = "",
     ylab = "",
     xaxt = "n", 
     yaxt = "n")

dat_sub <- dat[dat$fg == "G", ]
lines(dat_sub$y ~ dat_sub$diam, col = "#1b9e77", lwd = 2)
polygon(c(dat_sub$diam, rev(dat_sub$diam)), c(dat_sub$upper, rev(dat_sub$lower)),
        col = addTrans("#1b9e77",30), border = NA)

dat_sub <- dat[dat$fg == "C", ]
lines(dat_sub$y ~ dat_sub$diam, col = "#d95f02", lwd = 2)
polygon(c(dat_sub$diam, rev(dat_sub$diam)), c(dat_sub$upper, rev(dat_sub$lower)),
        col = addTrans("#d95f02",30), border = NA)

mtext(side = 1, text = "Diameter (cm)", line = 2.7, cex = 1.8)
mtext(side = 2, text = expression(paste("Diameter increment (cm year"^"-1", ")")), 
      line = 2.7, cex = 1.8)

axis(side =1, cex.axis = 1.5)
axis(side =2, cex.axis = 1.5)

legend("topright", legend = c("Forest species", "Savanna species"), col = c("#1b9e77", "#d95f02"),
       lty = 1)

dev.off()

# pdf(file = "growth_figure_for_video.pdf", width = 3, height =3, pointsize = 7)
# 
# dat <- data.frame(y = exp(eff$fit),
#                   lower = exp(eff$lower),
#                   upper = exp(eff$upper),
#                   diam = eff[["x"]][["dg.2006.cm"]],
#                   fg = eff[["x"]][["C.Cerrado..G.generalist."]])
# 
# par(mar = c(5.1, 5.1, 2, 1))
# 
# plot(x = NA, 
#      xlim = c(5, 60),
#      ylim = c(0.05, .23),
#      xlab = "",
#      ylab = "",
#      xaxt = "n", 
#      yaxt = "n")
# 
# dat_sub <- dat[dat$fg == "G", ]
# lines(dat_sub$y ~ dat_sub$diam, col = "#1b9e77", lwd = 2)
# polygon(c(dat_sub$diam, rev(dat_sub$diam)), c(dat_sub$upper, rev(dat_sub$lower)),
#         col = addTrans("#1b9e77",30), border = NA)
# 
# dat_sub <- dat[dat$fg == "C", ]
# lines(dat_sub$y ~ dat_sub$diam, col = "#d95f02", lwd = 2)
# polygon(c(dat_sub$diam, rev(dat_sub$diam)), c(dat_sub$upper, rev(dat_sub$lower)),
#         col = addTrans("#d95f02",30), border = NA)
# 
# mtext(side = 1, text = "Diameter (cm)", line = 2.7, cex = 1.8)
# mtext(side = 2, text = expression(paste("Diameter increment (cm year"^"-1", ")")), 
#       line = 2.7, cex = 1.8)
# 
# axis(side =1, cex.axis = 1.5)
# axis(side =2, cex.axis = 1.5)
# 
# legend("topright", legend = c("Forest species", "Savanna species"), col = c("#1b9e77", "#d95f02"),
#        lty = 1)
# 
# dev.off()


#------------------------------------------------------------------------------
# Figure 4
# Growth ~ BA above by FT
#------------------------------------------------------------------------------

eff <- Effect(focal.predictors = c("BA.above2011", "C.Cerrado..G.generalist."), 
              mod = ft_ci_soils_lmm, xlevels = 100, 
              transformation = list(link = log, inverse = exp),
              partial.residuals = TRUE)

tiff(filename="./plots/Figure 4 increment_BA_above_by_fg.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 3, 
     height=3, 
     pointsize=7, 
     res=600)

dat <- data.frame(y = exp(eff$fit),
                  lower = exp(eff$lower),
                  upper = exp(eff$upper),
                  ba = eff[["x"]][["BA.above2011"]],
                  fg = eff[["x"]][["C.Cerrado..G.generalist."]])

y <- eff$fit
x <- eff[["x"]]$"BA.above2011"
x.fit <- eff[["x.all"]]$"BA.above2011"
fg <- eff[["x.all"]]$"C.Cerrado..G.generalist."

data_orig <- data.frame(x.fit, 
                        fg)

#make sure we match the correct line
fitted <- ifelse(data_orig$fg == "C", y[1:100][closest(x.fit, x[1:100])],
                 y[101:200][closest(x.fit, x[101:200])])

resids <- eff$residuals + fitted
exp_resids <- exp(resids)



par(mar = c(5.1, 5.1, 2, 1))

plot(x = NA, 
     xlim = c(0.01, 30),
     ylim = c(0.01, 2),
     xlab = "",
     ylab = "",
     xaxt = "n", 
     yaxt = "n",
     log = "y")

points(exp_resids ~ I(x.fit*10), 
       col = ifelse(fg == "G", addTrans("#1b9e77",50), addTrans("#d95f02",50)),
       bg = ifelse(fg == "G", addTrans("#1b9e77",50), addTrans("#d95f02",50)),
       pch = 21,
       cex = 1)


dat_sub <- dat[dat$fg == "G", ]
lines(dat_sub$y ~ I(dat_sub$ba * 10), col = "#1b9e77", lwd = 2)
polygon(c(I(dat_sub$ba*10), rev(I(dat_sub$ba*10))), c(dat_sub$upper, rev(dat_sub$lower)),
        col = addTrans("#1b9e77",30), border = NA)

dat_sub <- dat[dat$fg == "C", ]
lines(dat_sub$y ~ I(dat_sub$ba * 10), col = "#d95f02", lwd = 2)
polygon(c(I(dat_sub$ba*10), rev(I(dat_sub$ba*10))), c(dat_sub$upper, rev(dat_sub$lower)),
        col = addTrans("#d95f02",30), border = NA)

mtext(side = 1, text = expression(paste("BA above (m"^"2"," ha"^"-1",")")), line = 3.1, cex = 1.8)
mtext(side = 2, text = expression(paste("Diameter increment (cm year"^"-1", ")")), 
      line = 2.7, cex = 1.8)

axis(side =1, cex.axis = 1.5)
axis(side =2, cex.axis = 1.5)

legend("topright", legend = c("Forest species", "Savanna species"), col = c("#1b9e77", "#d95f02"),
       lty = 1)

dev.off()


#------------------------------------------------------------------------------
# Figure 6
# Mortality ~ BA by FT
#------------------------------------------------------------------------------
eff <- Effect(focal.predictors = c("BA.above2006", "C.Cerrado..G.generalist."), 
              mod = mort_ft_ci_soils_glm_2011, xlevels = 100)
eff2 <- Effect(focal.predictors = c("BA.above2011", "C.Cerrado..G.generalist."),
mod = mort_ft_ci_soils_glm_2016, xlevels = 100)


tiff(filename="./plots/Figure 5 mortality_ba_above_by_fg.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 3, 
     height=3, 
     pointsize=7, 
     res=600)

dat1 <- data.frame(y = inv.logit(eff$fit)/5,
                   lower = inv.logit(eff$lower)/5,
                   upper = inv.logit(eff$upper)/5,
                   ba = eff[["x"]][["BA.above2006"]]*10,
                   fg = eff[["x"]][["C.Cerrado..G.generalist."]])

dat2 <- data.frame(y = inv.logit(eff2$fit)/5,
                   lower = inv.logit(eff2$lower)/5,
                   upper = inv.logit(eff2$upper)/5,
                   ba = eff2[["x"]][["BA.above2011"]]*10,
                   fg = eff2[["x"]][["C.Cerrado..G.generalist."]])

par(mar = c(5.1, 5.1, 2, 1))

plot(x = NA, 
     xlim = c(0, 30),
     ylim = c(0, 0.12),
     xlab = "",
     ylab = "",
     xaxt = "n", 
     yaxt = "n")

dat_sub <- dat1[dat1$fg == "G", ]
lines(dat_sub$y ~ dat_sub$ba, col = "#1b9e77", lwd = 2)
polygon(c(dat_sub$ba, rev(dat_sub$ba)), c(dat_sub$upper, rev(dat_sub$lower)),
        col = addTrans("#1b9e77",30), border = NA)

dat_sub <- dat1[dat1$fg == "C", ]
lines(dat_sub$y ~ dat_sub$ba, col = "#d95f02", lwd = 2)
polygon(c(dat_sub$ba, rev(dat_sub$ba)), c(dat_sub$upper, rev(dat_sub$lower)),
        col = addTrans("#d95f02",30), border = NA)

dat_sub <- dat2[dat2$fg == "G", ]
lines(dat_sub$y ~ dat_sub$ba, col = "#1b9e77", lwd = 2, lty = 2)
polygon(c(dat_sub$ba, rev(dat_sub$ba)), c(dat_sub$upper, rev(dat_sub$lower)),
        col = addTrans("#1b9e77",30), border = NA)

dat_sub <- dat2[dat2$fg == "C", ]
lines(dat_sub$y ~ dat_sub$ba, col = "#d95f02", lwd = 2, lty = 2)
polygon(c(dat_sub$ba, rev(dat_sub$ba)), c(dat_sub$upper, rev(dat_sub$lower)),
        col = addTrans("#d95f02",30), border = NA)

mtext(side = 1, text = expression(paste("BA above (m"^"2"," ha"^"-1",")")), line = 3.1, cex = 1.8)
mtext(side = 2, text = expression(paste("p(mortality)")), 
      line = 2.7, cex = 1.8)

axis(side =1, cex.axis = 1.5)
axis(side =2, cex.axis = 1.5)

legend("topleft", legend = c("Forest species, 2006", "Forest species, 2011", "Savanna species, 2006", "Savanna species, 2011"), col = c("#1b9e77", "#1b9e77", "#d95f02", "#d95f02"),
       lty = c(1, 2, 1, 2))


dev.off()


tiff(filename="./plots/Figure 6 mortality_ba_above_by_fg_for_talk.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 3, 
     height=3, 
     pointsize=7, 
     res=600)

dat1 <- data.frame(y = inv.logit(eff$fit)/5,
                   lower = inv.logit(eff$lower)/5,
                   upper = inv.logit(eff$upper)/5,
                   ba = eff[["x"]][["BA.above2006"]]*10,
                   fg = eff[["x"]][["C.Cerrado..G.generalist."]])

# dat2 <- data.frame(y = inv.logit(eff2$fit),C
#                    lower = inv.logit(eff2$lower),
#                    upper = inv.logit(eff2$upper),
#                    ba = eff2[["x"]][["BA.above2011"]],
#                    fg = eff2[["x"]][["C.Cerrado..G.generalist."]])

par(mar = c(5.1, 5.1, 2, 1))

plot(x = NA, 
     xlim = c(0, 30),
     ylim = c(0, 0.12),
     xlab = "",
     ylab = "",
     xaxt = "n", 
     yaxt = "n")

dat_sub <- dat1[dat1$fg == "G", ]
lines(dat_sub$y ~ dat_sub$ba, col = "#1b9e77", lwd = 2)
polygon(c(dat_sub$ba, rev(dat_sub$ba)), c(dat_sub$upper, rev(dat_sub$lower)),
        col = addTrans("#1b9e77",30), border = NA)

dat_sub <- dat1[dat1$fg == "C", ]
lines(dat_sub$y ~ dat_sub$ba, col = "#d95f02", lwd = 2)
polygon(c(dat_sub$ba, rev(dat_sub$ba)), c(dat_sub$upper, rev(dat_sub$lower)),
        col = addTrans("#d95f02",30), border = NA)

# dat_sub <- dat2[dat2$fg == "G", ]
# lines(dat_sub$y ~ dat_sub$ba, col = "#1b9e77", lwd = 2, lty = 2)
# polygon(c(dat_sub$ba, rev(dat_sub$ba)), c(dat_sub$upper, rev(dat_sub$lower)),
#         col = addTrans("#1b9e77",30), border = NA)
# 
# dat_sub <- dat2[dat2$fg == "C", ]
# lines(dat_sub$y ~ dat_sub$ba, col = "#d95f02", lwd = 2, lty = 2)
# polygon(c(dat_sub$ba, rev(dat_sub$ba)), c(dat_sub$upper, rev(dat_sub$lower)),
#         col = addTrans("#d95f02",30), border = NA)

mtext(side = 1, text = expression(paste("BA above (m"^"2"," ha"^"-1",")")), line = 3.1, cex = 1.8)
mtext(side = 2, text = expression(paste("p(mortality)")), 
      line = 2.7, cex = 1.8)

axis(side =1, cex.axis = 1.5)
axis(side =2, cex.axis = 1.5)

# legend("topleft", legend = c("Forest species, 2006", "Forest species, 2011", "Savanna species, 2006", "Savanna species, 2011"), col = c("#1b9e77", "#1b9e77", "#d95f02", "#d95f02"),
#        lty = c(1, 2, 1, 2))


dev.off()

#------------------------------------------------------------------------------
# Figure 5
# Growth ~ rain throughfall
#------------------------------------------------------------------------------
tiff(filename="./plots/Figure 5 rainfall effect plot.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 3, 
     height=3, 
     pointsize=7, 
     res=600)

par(mar = c(5.1, 5.1, 2, 1))

eff <- Effect(focal.predictors = c("Net.rainfall.."), 
              mod = soils_full_lmm, xlevels = 100, transformation = list(link = log, inverse = exp))

dat <- data.frame(y = exp(eff$fit),
                  lower = exp(eff$lower),
                  upper = exp(eff$upper),
                  rain = eff[["x"]][["Net.rainfall.."]])

plot(x = NA, 
     xlim = c(75, 100),
     ylim = c(0.1, 0.25),
     xlab = "",
     ylab = "",
     xaxt = "n", 
     yaxt = "n")

lines(dat$y ~ dat$rain, col = "darkblue", lwd = 2)
polygon(c(dat$rain, rev(dat$rain)), c(dat$upper, rev(dat$lower)),
        col = addTrans("darkblue",30), border = NA)

mtext(side = 1, text = "Percent rain throughfall", line = 2.7, cex = 1.8)
mtext(side = 2, text = expression(paste("Diameter increment (cm year"^"-1", ")")), 
      line = 2.7, cex = 1.8)

axis(side =1, cex.axis = 1.5)
axis(side =2, cex.axis = 1.5)

dev.off()



#------------------------------------------------------------------------------
# Supplemental Figure S2
# Size-growth relationship without accounting for covariates
#------------------------------------------------------------------------------
tiff(filename="./plots/Figure S2 diameter effect plot.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 3, 
     height=3, 
     pointsize=7, 
     res=600)

par(mar = c(5.1, 5.1, 2, 1))

eff <- Effect(focal.predictors = c("dg.2006.cm"), 
              mod = mod7, xlevels = 100, 
              transformation = list(link = log, inverse = exp))
dat <- data.frame(y = exp(eff$fit),
                  lower = exp(eff$lower),
                  upper = exp(eff$upper),
                  diam = eff[["x"]][["dg.2006.cm"]])

plot(x = NA, 
     xlim = c(0, 55),
     ylim = c(0.1, 0.35),
     xlab = "",
     ylab = "",
     xaxt = "n", 
     yaxt = "n")

lines(dat$y ~ dat$diam, col = "darkblue", lwd = 2)
polygon(c(dat$diam, rev(dat$diam)), c(dat$upper, rev(dat$lower)),
        col = addTrans("darkblue",30), border = NA)

mtext(side = 1, text = "Diameter", line = 2.7, cex = 1.8)
mtext(side = 2, text = expression(paste("Diameter increment (cm year"^"-1", ")")), 
      line = 2.7, cex = 1.8)

axis(side =1, cex.axis = 1.5)
axis(side =2, cex.axis = 1.5)

dev.off()


