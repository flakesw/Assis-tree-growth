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

se <- function(x){sd(x)/sqrt(length(x))}

# function to find the nearest y-value for an x data point, to add to partial residuals
# for plotting. Stolen from effects package
closest <- function(x, x0){
  apply(outer(x, x0, FUN = function(x, x0) abs(x - x0)), 1, which.min)
}

#-----------------------------------------------------------------
# Data prep
tree_data <- read.csv("./Giselda_data_per_tree_sf_edit.csv", stringsAsFactors = FALSE)
plot_data <- read.csv("./Giselda_Data_per_plot.csv")[1:30, ]

tree_data <- (tree_data[!(tree_data$Plot %in% c(21:25)), ]) #remove plots with incomplete data

#remove the "forest" functional type and relabel as generalist
tree_data[tree_data$C.Cerrado..G.generalist. == "F", ]$C.Cerrado..G.generalist. <- "G"

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
tree_data$BA.above2006 <- NA
tree_data$BA.above2011 <- NA

for(i in 1:nrow(tree_data)){
  tree_no <- tree_data[i, "Current.number"]
  tree_h <- tree_data[i, "H.2006"]
  plot_no <- tree_data[i, "Plot"]
  trees <- get_heights_csas(plot_no, 2006)
  trees <- trees[trees$Current.number != tree_no, ]
  #BA.above is the sum of basal area of taller trees
  if(!is.na(tree_h)){
    tree_data$BA.above2006[i] <- sum(trees[trees$H.2006 >= tree_h, "G2006..m2."], na.rm = TRUE)
  }
}

for(i in 1:nrow(tree_data)){
  tree_no <- tree_data[i, "Current.number"]
  tree_h <- as.numeric(tree_data[i, "H.2011..m."])
  plot_no <- tree_data[i, "Plot"]
  trees <- get_heights_csas(plot_no, 2011)
  trees <- trees[trees$Current.number != tree_no, ]
  #BA.above is the sum of basal area of taller trees
  if(!is.na(tree_h)){
    tree_data$BA.above2011[i] <- sum(trees[trees$H.2011..m. >= tree_h, "G2011..m2."], na.rm = TRUE)
  }
}


tree_data$Died2011 <- as.factor(ifelse(tree_data$Alive2006 == "v" & tree_data$Alive2011 == "m",
                                       "y", "n"))
tree_data$Died2016 <- as.factor(ifelse(tree_data$Alive2011 == "v" & tree_data$Alive2016 == "m",
                                       "y", "n"))
tree_data$Died_all <- as.factor(ifelse(tree_data$Died2011 == "y" | tree_data$Died2016 == "y", "y", "n"))


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

#only use trees which were present in 2006
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

pca_scores_std <- scale(pca$scores)

plot_data$PC1 <- pca$scores[, 1]
plot_data$PC2 <- pca$scores[, 2]
plot_data$PC1_std <- pca_scores_std[, 1]
plot_data$PC2_std <- pca_scores_std[, 2]

#merge all data together
all_data_pos <- join(tree_data_pos, plot_data, by = c("Plot"))
all_data<- join(tree_data, plot_data, by = c("Plot"))


sp_list <- all_data[which(!duplicated(all_data$Species..names.not.updated.)), ]

#how does QMD change over the gradient?
qmd <- aggregate(as.numeric(tree_data$dg.2006.cm), by = list(tree_data$Plot), FUN = function(x){sqrt(mean(x^2, na.rm = TRUE))})
plot(qmd$x ~ plot_data$TBA.2006.m2ha.1)

#do generalists grow faster?
t.test(dbh.annual.increment.cm._10.years ~ C.Cerrado..G.generalist., data = all_data_pos)

#----------------------------------------------------------------------
# Model comparison for table 1
#----------------------------------------------------------------------
mod1 <- expression(lmer(log(bai) ~ (1|Code), data = all_data_pos))
  
mod2 <- expression(lmer(log(bai) ~ scale(log(dg.2006.cm)) + (1|Code), data = all_data_pos))
  
mod3 <- expression(lmer(log(bai) ~ scale(PC1)  + scale(PC2) + 
               (1|Code),
             data = all_data_pos))
  
mod4 <- expression(lmer(log(bai) ~ scale(PC1)  +
               (1|Code),
             data = all_data_pos))

mod5 <- expression(lmer(log(bai) ~ poly(PC1,2)  +
               (1|Code),
             data = all_data_pos))
  
mod6 <- expression(lmer(log(bai) ~ scale(log(dg.2006.cm)) +
               poly(PC1, 2) +  
               (1|Code),
             data = all_data_pos))
  
mod7 <- expression(lmer(log(bai) ~ scale(log(dg.2006.cm)) +
               poly(PC1, 2)  +  C.Cerrado..G.generalist. + 
               (1|Code),
             data = all_data_pos))
  
mod8 <- expression(lmer(log(bai) ~ scale(log(dg.2006.cm)) +
               poly(PC1, 2)  + C.Cerrado..G.generalist.*TBA.2006.m2ha.1 + 
               (1|Code),
             data = all_data_pos))
  
mod9 <- expression(lmer(log(bai) ~ scale(log(dg.2006.cm)) +
               poly(PC1, 2)  +  C.Cerrado..G.generalist.*Canopy.cover.. + 
               (1|Code),
             data = all_data_pos))
  
mod10 <- expression(lmer(log(bai) ~ scale(log(dg.2006.cm)) +
               poly(PC1, 2)  +  scale(BA.above2006)*C.Cerrado..G.generalist. + 
               (1|Code),
             data = all_data_pos)) 
  
mod11 <- expression(lmer(log(bai) ~ scale(log(dg.2006.cm))*C.Cerrado..G.generalist. +
                poly(PC1, 2)  +  scale(BA.above2006)*C.Cerrado..G.generalist. + 
               (1|Code),
             data = all_data_pos))
  
mod12 <- expression(lmer(log(bai) ~ scale(log(dg.2006.cm)) +
                  scale(PC1) + scale(BA.above2011)*C.Cerrado..G.generalist. + 
                  (1|Code),
                data = all_data_pos))
  
mod13 <- expression(lmer(log(bai) ~ scale(log(dg.2006.cm)) +
                poly(PC1, 2) + scale(BA.above2011)*C.Cerrado..G.generalist. + 
                (1|Code),
              data = all_data_pos))
  
mod14 <- expression(lmer(log(bai) ~ scale(log(dg.2006.cm))*C.Cerrado..G.generalist. +
                poly(PC1, 2)  + scale(BA.above2011)*C.Cerrado..G.generalist. + 
                (1|Code),
              data = all_data_pos))

mod15 <- expression(lmer(log(bai) ~ scale(log(dg.2006.cm)) +
                           scale(PC1)  + scale(BA.above2011)*C.Cerrado..G.generalist. + 
                           (1|Code),
                         data = all_data_pos))


model_list <- c(paste0("mod", seq(1:15)))

model_summary <- list(model_names = model_list,
                            models_fit = lapply(model_list, FUN = function(x){eval(eval(str2expression(x)))}))

model_summary$n_params <- lapply(model_summary$models_fit, function(x){attributes(logLik(x))$df})
model_summary$AICc <- lapply(model_summary$models_fit, function(x){AICc(x)})
model_summary$R2m <- lapply(model_summary$models_fit, function(x){r.squaredGLMM(x)[1,1]})
model_summary$R2c <- lapply(model_summary$models_fit, function(x){r.squaredGLMM(x)[1,2]})

model_table <- data.frame(model_number = model_list,
                          desc = NA,
                          n_params = unlist(model_summary$n_params),
                          AICc = unlist(model_summary$AICc),
                          R2m = unlist(model_summary$R2m),
                          R2c = unlist(model_summary$R2c))

write.csv(model_table, "./model outputs/model_table.csv")

#----------------------------------------------------------------------
# Full model: FT x CI interaction with soils
ft_ci_soils_lmm <- lmer(log(bai) ~ scale(log(dg.2006.cm)) +
                               poly(PC1, 2, raw = FALSE) + scale(BA.above2011)*C.Cerrado..G.generalist. + 
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
# qqmath(ranef(ft_ci_soils_lmm, condVar = TRUE), strip = FALSE)$Plot
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

mort_ft_ci_soils_glm_died_all<- glm(Died_all ~ scale(log(dg.2006.cm)) +
                                   scale(BA.above2006)*C.Cerrado..G.generalist.,
                                 data = all_data[all_data$dg.2006.cm > 0, ],
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


summary(mort_ft_ci_soils_glm_died_all)
r.squaredGLMM(mort_ft_ci_soils_glm_died_all)
AIC(mort_ft_ci_soils_glm_died_all)
plot.roc(mort_ft_ci_soils_glm_died_all$y, 
         fitted(mort_ft_ci_soils_glm_died_all),print.auc = TRUE, 
         col = "green", lty = 2)
plot(Effect(mort_ft_ci_soils_glm_died_all, 
            focal.predictors = c("BA.above2006", "C.Cerrado..G.generalist.")))

#hosmer-lemeshow test seems pretty okay!
hl1 <- hoslem.test(mort_ft_ci_soils_glm_2011$y, 
                   fitted(mort_ft_ci_soils_glm_2011), g=10)
hl2 <- hoslem.test(mort_ft_ci_soils_glm_2016$y, 
                   fitted(mort_ft_ci_soils_glm_2016), g=10)
hl3 <- hoslem.test(mort_ft_ci_soils_glm_died_all$y, 
                   fitted(mort_ft_ci_soils_glm_died_all), g=10)

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
# Figure 3
# Growth ~ soil variables
#------------------------------------------------------------------------------
tiff(filename="./plots/Figure 3 soil variables effect plot with partial residuals.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 3, 
     height=3, 
     pointsize=7, 
     res=600)

par(mar = c(7.4, 5.1, 2, 1),
    oma = c(0,1,1,0))

eff <- Effect(focal.predictors = c("PC1"), 
              mod = ft_ci_soils_lmm, xlevels = 100, transformation = list(link = log, inverse = exp),
              partial.residuals = TRUE,
              se = TRUE)

y <- eff$fit
x <- eff[["x"]]$PC1
x.fit <- eff[["x.all"]]$PC1

fitted <- y[closest(x.fit, x)]
resids <- eff$residuals + fitted
exp_resids <- exp(resids)

resid_df <- data.frame(x.fit,
                       resids)

e1_r_agg <- aggregate(resid_df, by = list(resid_df$x.fit), FUN = mean)[, c(2,3)]
e1_r_agg[, 3] <- aggregate(resid_df, by = list(resid_df$x.fit), FUN = se)[, 3]
names(e1_r_agg) <- c("x.fit", "mean_part_resids", "se_resids")
e1_r_agg$low <- e1_r_agg$mean_part_resids - e1_r_agg$se_resids
e1_r_agg$high <- e1_r_agg$mean_part_resids + e1_r_agg$se_resids
e1_r_agg$exp_low <- exp(e1_r_agg$mean_part_resids - e1_r_agg$se_resids)
e1_r_agg$exp_high <- exp(e1_r_agg$mean_part_resids + e1_r_agg$se_resids)
e1_r_agg$exp_mean_resids <- exp(e1_r_agg[, 2])


eff2 <- Effect(focal.predictors = c("PC1"), 
               mod = eval(eval(mod3)), xlevels = 100, transformation = list(link = log, inverse = exp))

dat <- data.frame(y = exp(eff$fit),
                  lower = exp(eff$lower),
                  upper = exp(eff$upper),
                  soil = eff[["x"]][["PC1"]])

dat2 <- data.frame(y = exp(eff2$fit),
                   lower = exp(eff2$lower),
                   upper = exp(eff2$upper),
                   soil = eff2[["x"]][["PC1"]])

plot(x = NA, 
     xlim = c(-3.5,3.5),
     ylim = c(0.08, 0.25),
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

lines(dat2$y ~ dat2$soil, col = "magenta", lwd = 2)
polygon(c(dat2$soil, rev(dat2$soil)), c(dat2$upper, rev(dat2$lower)),
        col = addTrans("magenta",30), border = NA)

mtext(side = 1, text = "Soil PC1", line = 2.7, cex = 1.8)
mtext(side = 2, text = expression(paste("Diameter increment (cm year"^"-1", ")")), 
      line = 2.7, cex = 1.8)

axis(side =1, cex.axis = 1.5)
axis(side =2, cex.axis = 1.5)

legend("topright", legend = c("Soils Alone", "Full Model"), col = c("magenta", "darkblue"),
       lty = 1)

mtext(side = 1, line = 5, text = "Sand                                         Fines", cex = 1.4)

arrows(x0 = -2.5, x1 = 2.5, y0 = 0.007, y1 = 0.007, xpd = TRUE, code = 3, length = 0.125)

dev.off()
# 
# 
# tiff(filename="./plots/Figure 3 soil variables effect plot.tiff", 
#      type = "cairo",
#      antialias = "gray",
#      compression = "lzw",
#      units="in", 
#      width = 3, 
#      height=3, 
#      pointsize=7, 
#      res=600)
# 
# par(mar = c(5.1, 5.1, 2, 1))
# 
# eff <- Effect(focal.predictors = c("PC1"), 
#               mod = ft_ci_soils_lmm, xlevels = 100, transformation = list(link = log, inverse = exp),
#               partial.residuals = TRUE)
# 
# eff2 <- Effect(focal.predictors = c("PC1"), 
#                mod = mod3, xlevels = 100, transformation = list(link = log, inverse = exp))
# 
# dat <- data.frame(y = exp(eff$fit),
#                   lower = exp(eff$lower),
#                   upper = exp(eff$upper),
#                   soil = eff[["x"]][["PC1"]])
# 
# dat2 <- data.frame(y = exp(eff2$fit),
#                    lower = exp(eff2$lower),
#                    upper = exp(eff2$upper),
#                    soil = eff2[["x"]][["PC1"]])
# 
# plot(x = NA, 
#      xlim = c(-3,3.4),
#      ylim = c(0.08, 0.25),
#      xlab = "",
#      ylab = "",
#      xaxt = "n", 
#      yaxt = "n")
# 
# 
# lines(dat$y ~ dat$soil, col = "darkblue", lwd = 2)
# polygon(c(dat$soil, rev(dat$soil)), c(dat$upper, rev(dat$lower)),
#         col = addTrans("darkblue",30), border = NA)
# # points(exp_resids ~ x.fit)
# 
# 
# lines(dat2$y ~ dat2$soil, col = "magenta", lwd = 2)
# polygon(c(dat2$soil, rev(dat2$soil)), c(dat2$upper, rev(dat2$lower)),
#         col = addTrans("magenta",30), border = NA)
# 
# mtext(side = 1, text = "Soil PC1", line = 2.7, cex = 1.8)
# mtext(side = 2, text = expression(paste("Diameter increment (cm year"^"-1", ")")), 
#       line = 2.7, cex = 1.8)
# 
# axis(side =1, cex.axis = 1.5)
# axis(side =2, cex.axis = 1.5)
# 
# legend("topleft", legend = c("Soils Alone", "Full Model"), col = c("magenta", "darkblue"),
#        lty = 1)
# 
# dev.off()

#------------------------------------------------------------------------------
# Figure 4 
# Growth ~ diam by FT
#------------------------------------------------------------------------------
tiff(filename="./plots/Figure 4 differences in growth by fg.tiff",
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in",
     width = 6,
     height=3,
     pointsize=7,
     res=600)

par(mar = c(5.1, 5.1, 2, 1),
    mfrow = c(1,2))


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

dat <- data.frame(y = exp(eff$fit),
                  lower = exp(eff$lower),
                  upper = exp(eff$upper),
                  diam = eff[["x"]][["dg.2006.cm"]],
                  fg = eff[["x"]][["C.Cerrado..G.generalist."]])

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

text(x = 8, y = 0.242, labels = "(a)", cex = 1.5)



# growth by BA above

eff <- Effect(focal.predictors = c("BA.above2011", "C.Cerrado..G.generalist."),
              mod = ft_ci_soils_lmm, xlevels = 100,
              transformation = list(link = log, inverse = exp),
              partial.residuals = TRUE)
dat <- data.frame(y = exp(eff$fit),
                  lower = exp(eff$lower),
                  upper = exp(eff$upper),
                  ba = eff[["x"]][["BA.above2011"]]*10,
                  fg = eff[["x"]][["C.Cerrado..G.generalist."]])

plot(x = NA,
     xlim = c(0, 25),
     ylim = c(0, 0.6),
     xlab = "",
     ylab = "",
     xaxt = "n",
     yaxt = "n")

dat_sub <- dat[dat$fg == "G", ]
lines(dat_sub$y ~ dat_sub$ba, col = "#1b9e77", lwd = 2)
polygon(c(dat_sub$ba, rev(dat_sub$ba)), c(dat_sub$upper, rev(dat_sub$lower)),
        col = addTrans("#1b9e77",30), border = NA)

dat_sub <- dat[dat$fg == "C", ]
lines(dat_sub$y ~ dat_sub$ba, col = "#d95f02", lwd = 2)
polygon(c(dat_sub$ba, rev(dat_sub$ba)), c(dat_sub$upper, rev(dat_sub$lower)),
        col = addTrans("#d95f02",30), border = NA)

mtext(side = 1, text = expression(paste("BA above (m"^"2"," ha"^"-1",")")), line = 3.1, cex = 1.8)
mtext(side = 2, text = expression(paste("Diameter increment (cm year"^"-1", ")")),
      line = 2.7, cex = 1.8)

axis(side =1, cex.axis = 1.5)
axis(side =2, cex.axis = 1.5)


legend("topright", legend = c("Generalist species", "Savanna species"), col = c("#1b9e77", "#d95f02"),
       lty = 1)


text(x = 1.7, y = 0.58, labels = "(b)", cex = 1.5)

dev.off()

#------------------------------------------------------------------------------
# Figure 6
# Mortality ~ BA by FT
#------------------------------------------------------------------------------
eff <- Effect(focal.predictors = c("BA.above2006", "C.Cerrado..G.generalist."), 
              mod = mort_ft_ci_soils_glm_2011, xlevels = 100)
eff2 <- Effect(focal.predictors = c("BA.above2011", "C.Cerrado..G.generalist."),
              mod = mort_ft_ci_soils_glm_2016, xlevels = 100)


tiff(filename="./plots/Figure 6 mortality_ba_above_by_fg.tiff", 
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

legend("topleft", legend = c("Generalist species, 2006", "Generalist species, 2011", "Savanna species, 2006", "Savanna species, 2011"), col = c("#1b9e77", "#1b9e77", "#d95f02", "#d95f02"),
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

# legend("topleft", legend = c("Generalist species, 2006", "Generalist species, 2011", "Savanna species, 2006", "Savanna species, 2011"), col = c("#1b9e77", "#1b9e77", "#d95f02", "#d95f02"),
#        lty = c(1, 2, 1, 2))


dev.off()

#------------------------------------------------------------------------------
# Version with just one model
#------------------------------------------------------------------------------
eff <- Effect(focal.predictors = c("BA.above2006", "C.Cerrado..G.generalist."), 
              mod = mort_ft_ci_soils_glm_died_all, xlevels = 100)


tiff(filename="./plots/Figure 6 mortality_ba_above_by_fg.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 3, 
     height=3, 
     pointsize=7,
     res=600)

dat1 <- data.frame(y = 1 - (1 - inv.logit(eff$fit))^(1/10),
                   lower = 1 - (1-inv.logit(eff$lower))^(1/10),
                   upper = 1 - (1-inv.logit(eff$upper))^(1/10),
                   ba = eff[["x"]][["BA.above2006"]]*10,
                   fg = eff[["x"]][["C.Cerrado..G.generalist."]])

par(mar = c(1,7,1,1),
    oma = c(5,0,0,0),
    ps = 9,
    cex = 1,
    yaxs = "i",
    xaxs = "i")

plot(x = NA, 
     xlim = c(0, 30),
     ylim = c(-0.03, 0.18),
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


mtext(side = 1, text = expression(paste("BA above (m"^"2"," ha"^"-1",")")), line = 3.1, cex = 1.8)
mtext(side = 2, text = expression(paste("p(mortality)")), 
      line = 2.7, cex = 1.8)

axis(side =1, cex.axis = 1.5)
axis(side =2, cex.axis = 1.5, at = c(0, 0.05, .1, .15))

abline(h = 0)
abline(h = 0.15)

legend(x = 0, y = 0.15, legend = c("Generalist species", "Savanna species"), 
       col = c("#1b9e77", "#d95f02"),
       lty = c(1, 1))

h <- density(all_data[all_data$Died_all == "y" & all_data$C.Cerrado..G.generalist. == "C" & 
                        !is.na(all_data$BA.above2006), "BA.above2006"]*10, na.rm = TRUE)
h$y <- 0.18 - h$y/4
lines(h, col = c("#d95f02"))

h <- density(all_data[all_data$Died_all == "n" & all_data$C.Cerrado..G.generalist. == "C" & 
                        !is.na(all_data$BA.above2006), "BA.above2006"]*10, na.rm = TRUE)
h$y <- -0.03 + h$y/4
lines(h, col = c("#d95f02"))

h <- density(all_data[all_data$Died_all == "y" & all_data$C.Cerrado..G.generalist. == "G" & 
                        !is.na(all_data$BA.above2006), "BA.above2006"]*10, na.rm = TRUE)
h$y <- 0.18 - h$y/4
lines(h, col = c("#1b9e77"))

h <- density(all_data[all_data$Died_all == "n" & all_data$C.Cerrado..G.generalist. == "G" & 
                        !is.na(all_data$BA.above2006), "BA.above2006"]*10, na.rm = TRUE)
h$y <- -0.03 + h$y/4
lines(h, col = c("#1b9e77"))


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
              mod = eval(eval(mod7)), xlevels = 100, 
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

mtext(side = 1, text = "Original diameter (cm)", line = 2.7, cex = 1.8)
mtext(side = 2, text = expression(paste("Diameter increment (cm year"^"-1", ")")), 
      line = 2.7, cex = 1.8)

axis(side =1, cex.axis = 1.5)
axis(side =2, cex.axis = 1.5)

dev.off()



#------------------------------------------------------------------------------
# Basal area increment by functional type
#------------------------------------------------------------------------------
change_2006 <- data.frame(Plot = unique(tree_data$Plot),
                     BA_s = numeric(length(unique(tree_data$Plot))),
                     BA_g = numeric(length(unique(tree_data$Plot))),
                     BA_total = numeric(length(unique(tree_data$Plot))),
                     BA_s_died = numeric(length(unique(tree_data$Plot))),
                     BA_g_died = numeric(length(unique(tree_data$Plot)))
                     )
change_2011 <- change_2006
change_2016 <- change_2006


#sum the basal area by each functional type
for(plot in unique(tree_data$Plot)){
  temp <- tree_data[tree_data$Plot == plot, ]
  change_2006[change_2006$Plot == plot, ]$BA_s <- sum((((temp[temp$C.Cerrado..G.generalist. == "C", "dg.2006.cm"]/2)^2)*pi), na.rm = TRUE)/10000
  change_2006[change_2006$Plot == plot, ]$BA_g <- sum((((temp[temp$C.Cerrado..G.generalist. == "G", "dg.2006.cm"]/2)^2)*pi), na.rm = TRUE)/10000
  change_2006[change_2006$Plot == plot, ]$BA_total <- change_2006[change_2006$Plot == plot, ]$BA_s + change_2006[change_2006$Plot == plot, ]$BA_g
  change_2006[change_2006$Plot == plot, ]$BA_s_died <- NA
  change_2006[change_2006$Plot == plot, ]$BA_g_died <- NA
  
  
  change_2011[change_2011$Plot == plot, ]$BA_s <- sum((((temp[temp$C.Cerrado..G.generalist. == "C", "dg.2011.cm"]/2)^2)*pi), na.rm = TRUE)/10000
  change_2011[change_2011$Plot == plot, ]$BA_g <- sum((((temp[temp$C.Cerrado..G.generalist. == "G", "dg.2011.cm"]/2)^2)*pi), na.rm = TRUE)/10000
  change_2011[change_2011$Plot == plot, ]$BA_total <- change_2011[change_2011$Plot == plot, ]$BA_s + change_2011[change_2011$Plot == plot, ]$BA_g
  change_2011[change_2011$Plot == plot, ]$BA_s_died <- sum((((temp[temp$C.Cerrado..G.generalist. == "C" & temp$Died2011 == "y", "dg.2006.cm"]/2)^2)*pi), na.rm = TRUE)/10000
  change_2011[change_2011$Plot == plot, ]$BA_g_died <- sum((((temp[temp$C.Cerrado..G.generalist. == "G"& temp$Died2011 == "y", "dg.2006.cm"]/2)^2)*pi), na.rm = TRUE)/10000
  
  
  change_2016[change_2016$Plot == plot, ]$BA_s <- sum((((temp[temp$C.Cerrado..G.generalist. == "C", "dg.2016.cm"]/2)^2)*pi), na.rm = TRUE)/10000
  change_2016[change_2016$Plot == plot, ]$BA_g <- sum((((temp[temp$C.Cerrado..G.generalist. == "G", "dg.2016.cm"]/2)^2)*pi), na.rm = TRUE)/10000
  change_2016[change_2016$Plot == plot, ]$BA_total <- change_2016[change_2016$Plot == plot, ]$BA_s + change_2016[change_2016$Plot == plot, ]$BA_g
  change_2016[change_2016$Plot == plot, ]$BA_s_died <- sum((((temp[temp$C.Cerrado..G.generalist. == "C" & temp$Died2016 == "y", "dg.2011.cm"]/2)^2)*pi), na.rm = TRUE)/10000
  change_2016[change_2016$Plot == plot, ]$BA_g_died <- sum((((temp[temp$C.Cerrado..G.generalist. == "G"& temp$Died2016 == "y", "dg.2011.cm"]/2)^2)*pi), na.rm = TRUE)/10000
  
}

change_ba <- data.frame(Plot = change_2006$Plot,
                        ba_tot_2006 = change_2006$BA_total,
                        ba_inc_s = change_2016$BA_s - change_2006$BA_s,
                        ba_inc_g = change_2016$BA_g - change_2006$BA_g,
                        ba_died_s = change_2011$BA_s_died + change_2016$BA_s_died,
                        ba_died_g = change_2011$BA_g_died + change_2016$BA_g_died,
                        ba_inc_s_w_dead = change_2016$BA_s - change_2006$BA_s + change_2011$BA_s_died + change_2016$BA_s_died,
                        ba_inc_g_w_dead = change_2016$BA_g - change_2006$BA_g + change_2011$BA_g_died + change_2016$BA_g_died
                        )

change_ba_prop <- change_ba
change_ba_prop[, c(3, 5, 7)] <- change_ba_prop[, c(3,5,7)]/change_2006$BA_s
change_ba_prop[, c(4, 6, 8)] <- change_ba_prop[, c(4,6,8)]/change_2006$BA_g

#make three-panel figure

tiff(filename="./plots/Figure 1 change in basal area proportional.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="mm", 
     width = 90, 
     height=200, 
     pointsize=7, 
     res=600)


par(mar = c(1,7,1,1),
    oma = c(5,0,0,0),
    mfrow = c(3,1),
    ps = 9,
    cex = 1,
    yaxs = "i")

new <- data.frame(ba_tot_2006 = seq(0.5, 2.5, length.out = 100))

plot(NA,
     ylim = c(-0.7, 1.8),
     xlim = c(0.4, 2.6),
     xlab = expression(paste("Original basal area (m"^2, " ha"^-1, ")")),
     ylab = expression(paste("Net change in basal area (%)")),
     xaxt = "n",
     yaxt = "n",
     cex.lab = 1.2)

abline(h = 0)

points(ba_inc_s ~ ba_tot_2006, data = change_ba_prop, 
       pch = 21,
       col = "#d95f02",
       bg = "#d95f02")

mod <- lm(log(ba_inc_s + 1) ~ log(ba_tot_2006), data = change_ba_prop)
pred <- predict(mod, newdata = new)
lines(exp(pred) - 1 ~ new$ba_tot_2006,
      col = "#d95f02",
      lty = 1)


points(ba_inc_g ~ ba_tot_2006, data = change_ba_prop, 
       pch = 22,
       col = "#1b9e77",
       bg = "#1b9e77")
mod <- lm(log(ba_inc_g + 1) ~ log(ba_tot_2006), data = change_ba_prop)
pred <- predict(mod, newdata = new)
lines(exp(pred) - 1 ~ new$ba_tot_2006,
      col = "#1b9e77",
      lty = 1)

axis(side = 2, at = pretty(c(-0.7, 1.8)), labels = pretty(c(-0.7, 1.8))*100)
axis(side = 1, at = pretty(change_ba$ba_tot_2006), labels = NA)

legend(x = 0.5, y = -0.2, legend = c("Savanna species", "Generalist species"),
       pch = c(16,15),
       col =  c("#d95f02","#1b9e77"),
       cex = 0.9)

text(x = 0.5, y = 1.67, labels = "(a)", cex = 1.2)

#just losses due to mortality

plot(NA,
     ylim = c(0, .8),
     xlim = c(0.4, 2.6),
     xaxt = "n",
     yaxt = "n",
     xlab = "",
     ylab = expression(paste("Mortality basal area (%)")),
     cex.lab = 1.2)

abline(h = 0)

points(ba_died_s ~ ba_tot_2006, data = change_ba_prop, 
       pch = 21,
       col = "#d95f02",
       bg = "#d95f02")

mod <- lm(log(ba_died_s + 1) ~ log(ba_tot_2006), data = change_ba_prop)
pred <- predict(mod, newdata = new)
lines(exp(pred) - 1 ~ new$ba_tot_2006,
      col = "#d95f02",
      lty = 1)


points(ba_died_g ~ ba_tot_2006, data = change_ba_prop, 
       pch = 22,
       col = "#1b9e77",
       bg = "#1b9e77")
mod <- lm(log(ba_died_g + 1) ~ log(ba_tot_2006), data = change_ba_prop)
pred <- predict(mod, newdata = new)
lines(exp(pred) - 1 ~ new$ba_tot_2006,
      col = "#1b9e77",
      lty = 1)

axis(side = 2, at = pretty(c(0,0.8)), labels = pretty(c(0,0.8))*100)
axis(side = 1, at = pretty(change_ba$ba_tot_2006), labels = NA)


text(x = 0.5, y = 0.75, labels = "(b)", cex = 1.2)

#with resurrected trees
plot(NA,
     ylim = c(-0.7, 1.8),
     xlim = c(0.4, 2.6),
     xaxt = "n",
     yaxt = "n",
     xlab = "",
     ylab = expression(paste("Change in basal area \n without mortality (%)")),
     cex.lab = 1.2)

abline(h = 0)

points(ba_inc_s_w_dead ~ ba_tot_2006, data = change_ba_prop,
       pch = 21,
       col = "#d95f02",
       bg = "#d95f02")
mod <- lm(log(ba_inc_s_w_dead + 1) ~ log(ba_tot_2006), data = change_ba_prop)
pred <- predict(mod, newdata = new)
lines(exp(pred) - 1 ~ new$ba_tot_2006,
      col = "#d95f02",
      lty = 1)

points(ba_inc_g_w_dead ~ ba_tot_2006, data = change_ba_prop,
       pch = 22,
       col = "#1b9e77",
       bg = "#1b9e77")
mod <- lm(log(ba_inc_g_w_dead + 1) ~ log(ba_tot_2006), data = change_ba_prop)
pred <- predict(mod, newdata = new)
lines(exp(pred) - 1 ~ new$ba_tot_2006,
      col = "#1b9e77",
      lty = 1)

axis(side = 1, at = pretty(change_ba$ba_tot_2006), labels = pretty(change_ba$ba_tot_2006)*10)
axis(side = 2, at = pretty(c(-0.7, 1.8)), labels = pretty(c(-0.7, 1.8))*100)

mtext(side = 1, line = 2.9, text = expression(paste("Original basal area (m"^2, " ha"^-1, ")")),
      cex = 1.3)


text(x = 0.5, y = 1.67, labels = "(c)", cex = 1.2)

dev.off()


#  
# #combined into one figure
# png(filename="./plots/Figure 1 change in basal area.png", 
#     type = "cairo",
#     antialias = "gray",
#     # compression = "lzw",
#     units="in", 
#     width = 4, 
#     height=3, 
#     pointsize=7, 
#     res=600)
# 
# par(mar = c(5,5,1,1),
#     oma = c(1,1,0,0))
# 
# new <- data.frame(BA_total = seq(0.5, 2.5, length.out = 100))
# 
# plot(NA,
#      ylim = c(-0.1, 1.2),
#      xlim = c(0.4, 2.6),
#      xaxt = "n",
#      yaxt = "n",
#      xlab = expression(paste("Original basal area (m"^2, " ha"^-1, ")")),
#      ylab = expression(paste("Change in basal area (m"^2, " ha"^-1, ")")),
#      cex.lab = 1.2)
# 
# abline(h = 0)
# 
# points(ba_inc_s ~ ba_tot_2006, data = change_ba, 
#        pch = 21,
#        col = "#d95f02",
#        bg = "#d95f02")
# points(ba_inc_s_w_dead ~ ba_tot_2006, data = change_ba,
#        pch = 21,
#        col = "#d95f02")
# segments(x0 = change_ba$ba_tot_2006, x1 = change_ba$ba_tot_2006,
#          y0 = change_ba$ba_inc_s,
#          y1 = change_ba$ba_inc_s_w_dead,
#          col = "#d95f02")
# abline(coef(lm(ba_inc_s ~ ba_tot_2006, data = change_ba)),
#        col = "#d95f02",
#        lty = 1)
# abline(coef(lm(ba_inc_s_w_dead ~ ba_tot_2006, data = change_ba)),
#        col = "#d95f02",
#        lty = 2)
# 
# points(ba_inc_g ~ ba_tot_2006, data = change_ba, 
#        pch = 22,
#        col = "#1b9e77",
#        bg = "#1b9e77")
# points(ba_inc_g_w_dead ~ ba_tot_2006, data = change_ba,
#        pch = 22,
#        col = "#1b9e77")
# segments(x0 = change_ba$ba_tot_2006, x1 = change_ba$ba_tot_2006,
#          y0 = change_ba$ba_inc_g,
#          y1 = change_ba$ba_inc_g_w_dead,
#          col = "#1b9e77")
# 
# abline(coef(lm(ba_inc_g ~ ba_tot_2006, data = change_ba)),
#        col = "#1b9e77",
#        lty = 1)
# abline(coef(lm(ba_inc_g_w_dead ~ ba_tot_2006, data = change_ba)),
#        col = "#1b9e77",
#        lty = 2)
# 
# legend(x = .5, y = 1.2, legend = c("Savanna", "Savanna + Dead", "Forest", "Forest + Dead"),
#        pch = c(16, 21, 15, 22),
#        col =  c("#d95f02","#d95f02","#1b9e77","#1b9e77"),
#        cex = 0.9)
# 
# axis(side = 1, at = pretty(change_ba$ba_tot_2006), labels = pretty(change_ba$ba_tot_2006)*10)
# axis(side = 2, at = pretty(c(-0.1, 1.2)), labels = pretty(c(-0.1, 1.2))*10)
# 
# dev.off()
# 
# curve <- curvefit(change_2006$BA_total, I((change_2016$BA_g - change_2006$BA_g)/change_2006$BA_g), y.max = NULL, extrapol = NULL, 
#                   plot.curves = FALSE, print.results = FALSE)

# #------------------------------------------------------------------------------
# # Leaf area increment by FT
# #------------------------------------------------------------------------------
# 
# la_mod <- readRDS("./diam_fg_lm.RDS")
# summary(la_mod)
# 
# all_data$la_2006 <- predict(la_mod, newdata = list(FG = all_data$FG, D30_eff = all_data$dg.2006.cm))
# all_data$la_2011 <- predict(la_mod, newdata = list(FG = all_data$FG, D30_eff = all_data$dg.2011.cm))
# all_data$la_2016 <- predict(la_mod, newdata = list(FG = all_data$FG, D30_eff = all_data$dg.2016.cm))
# all_data$delta_la <- exp(all_data$la_2016) - exp(all_data$la_2006)
# all_data$delta_la <- all_data$delta_la*10
# 
# delta_la_plot <- aggregate(all_data$delta_la, by = list(all_data$FG, all_data$Plot), FUN = function(x){sum(x, na.rm = TRUE)})[-61, ]
# names(delta_la_plot)[2] <- "Plot"
# delta_la_plot <- join(delta_la_plot, plot_data[, c("Plot", "TBA.2011.m2.ha.1")], by = c("Plot"))
# 
# # plot(NA, 
# #      xlim = c(5, 27),
# #      ylim = c(0, 800),
# #      ylab = "LA Increment",
# #      xlab = "Stand BA")
# # points(x ~ TBA.2011.m2.ha.1, data = delta_la_plot[delta_la_plot$Group.1 == "F", ], col = "#1b9e77")
# # points(x ~ TBA.2011.m2.ha.1, data = delta_la_plot[delta_la_plot$Group.1 == "S", ], col = "#d95f02")
# 
# ggplot(data = delta_la_plot, aes(x = TBA.2011.m2.ha.1, y = x)) +
#   geom_hline(yintercept=0) + 
#   geom_point(aes(color = Group.1)) +
#   # geom_line(aes(color = Group.1)) + 
#   geom_smooth(aes(color = Group.1)) +
#   scale_color_manual(values=c("#1b9e77", "#d95f02")) + 
#   xlab(expression(paste("Stand basal area (m"^2, " ha"^-1, ")"))) + 
#   ylab(expression(paste("Change in leaf area (m"^2, " ha"^-1, ")"))) + 
#   theme(legend.position = "none")
# 
