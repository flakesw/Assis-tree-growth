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
library("maptools")
library("vegan")
library("dplyr")
library("lattice")
library("flora")

options(na.action = "na.omit")

#------------------------------------------------------------------------------
# Set some parameters

min.growth <- -0.2 #omit trees which grow more slowly than this value

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
plot_data <- read.csv("./Giselda_Data_per_plot.csv", stringsAsFactors = FALSE)[1:30, ]
plot_data  <- plot_data[, c(1:3, 14, 18)]

soil_data1 <- read.csv("./soil_nutrients_0-20.csv", header = TRUE, 
                       fileEncoding= "UTF-8-BOM", stringsAsFactors = FALSE)[-c(1,2), ]
soil_data1[, -c(1,2)] <- apply(soil_data1[, -c(1,2)], 2, FUN = as.numeric)
soil_data1$Fines <- soil_data1$Argila.Silte
soil_data1$Sand <- 100 - soil_data1$Fines
names(soil_data1)[-c(1,2)] <- paste0(names(soil_data1)[-c(1,2)], "_20")
soil_data1 <- soil_data1[, c(1, 5, 6 ,9:25, 27)]

soil_data2 <- read.csv("./soil_nutrients_60-80.csv", header = TRUE, 
                       fileEncoding= "UTF-8-BOM", stringsAsFactors = FALSE)[-c(1,2), -24]
soil_data2[, -c(1,2)] <- apply(soil_data2[, -c(1,2)], 2, FUN = as.numeric)
soil_data2$Fines <- soil_data2$Argila.Silte
soil_data2$Sand <- 100 - soil_data2$Fines
names(soil_data2)[-c(1,2)] <- paste0(names(soil_data2)[-c(1,2)], "_60")
soil_data2 <- soil_data2[, c(1, 3, 4, 7:23, 25)]

plot_data  <- plyr::join(plot_data, soil_data1, by = "Plot")
plot_data  <- plyr::join(plot_data, soil_data2, by = "Plot")

tree_data <- tree_data[!(tree_data$Plot %in% c(21:25)), ] #remove plots with incomplete data

#remove the "forest" functional type and relabel as generalist
tree_data[tree_data$C.Cerrado..G.generalist. == "F", ]$C.Cerrado..G.generalist. <- "G"

# replace missing data with NA

blanks <- apply(tree_data, c(1,2), is.blank)
tree_data[blanks] <- NA

# tree_data$dbh.annual.increment.cm._10.years <- as.numeric(as.character(tree_data$dbh.annual.increment.cm._10.years))

# create a code for species names
# split up the species names
species <- strsplit(tree_data$Species..names.not.updated., split = " ")
# add 4-letter code to tree dataframe
tree_data$Code <- as.factor(unlist(lapply(species, make_code)))

tree_data[is.na(tree_data$Exclude), "Exclude"] <- "n"
tree_data <- tree_data[tree_data$Exclude != "y", ]

#fix up variable types
tree_data$Shade.Tolerance <- as.factor(tree_data$Shade.Tolerance)
tree_data$C.Cerrado..G.generalist. <- as.factor(tree_data$C.Cerrado..G.generalist.)
tree_data$H.2006 <- as.numeric(tree_data$H.2006)
tree_data$H.2011..m. <- as.numeric(tree_data$H.2011..m.)
tree_data$H.2016.m <- as.numeric(tree_data$H.2016.m)
tree_data$G2006..m2. <- as.numeric(tree_data$G2006..m2.)
tree_data$G2011..m2. <- as.numeric(tree_data$G2011..m2.)
tree_data$dg.2006.cm <- as.numeric(tree_data$dg.2006.cm)
tree_data$dg.2016.cm <- as.numeric(tree_data$dg.2016.cm)
tree_data$dg.2011.cm <- as.numeric(tree_data$dg.2011.cm)

# tree_data$bai <- tree_data$dg.2016.cm - tree_data$dg.2006.cm
# plot_data$bai <- aggregate(tree_data$bai, by = list(tree_data$Plot),
#                            FUN = function(x){sum(x, na.rm = TRUE)})


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


# clean up data to estimate growth

# calculate diameter growth separately for each time interval
tree_data$inc_06_11 <- NA
tree_data$inc_11_16 <- NA
tree_data$inc_06_16 <- NA
for(i in 1:nrow(tree_data)){
  row <- tree_data[i, ]
  if(!is.na(row$Alive2006) & !is.na(row$Alive2011) & row$Alive2006 == "v" & row$Alive2011 == "v" ){
    tree_data[i, "inc_06_11"] <- (as.numeric(row$dg.2011.cm) - as.numeric(row$dg.2006.cm))/5
  }
  if(!is.na(row$Alive2011) & !is.na(row$Alive2016) & row$Alive2011 == "v" & row$Alive2016 == "v" ){
    tree_data[i, "inc_11_16"] <- (as.numeric(row$dg.2016.cm) - as.numeric(row$dg.2011.cm))/5
  }
  if(!is.na(row$Alive2006) & !is.na(row$Alive2016) & row$Alive2006 == "v" & row$Alive2016 == "v" ){
    tree_data[i, "inc_06_16"] <- (as.numeric(row$dg.2016.cm) - as.numeric(row$dg.2006.cm))/10
  }
}

#calculate height growth for each interval
tree_data$inc_06_11_ht <- NA
tree_data$inc_11_16_ht <- NA
tree_data$inc_06_16_ht <- NA
for(i in 1:nrow(tree_data)){
  row <- tree_data[i, ]
  if(!is.na(row$Alive2006) & !is.na(row$Alive2011) & row$Alive2006 == "v" & row$Alive2011 == "v" ){
    tree_data[i, "inc_06_11_ht"] <- (as.numeric(row$H.2011..m.) - as.numeric(row$H.2006))/5
  }
  if(!is.na(row$Alive2011) & !is.na(row$Alive2016) & row$Alive2011 == "v" & row$Alive2016 == "v" ){
    tree_data[i, "inc_11_16_ht"] <- (as.numeric(row$H.2016.m) - as.numeric(row$H.2011..m.))/5
  }
  if(!is.na(row$Alive2006) & !is.na(row$Alive2016) & row$Alive2006 == "v" & row$Alive2016 == "v" ){
    tree_data[i, "inc_06_16_ht"] <- (as.numeric(row$H.2016.m) - as.numeric(row$H.2006))/10
  }
}

# first interval growth is higher than second interval or overall growth
t.test(tree_data$inc_06_11_ht, tree_data$inc_11_16_ht)
t.test(tree_data$inc_06_11_ht, tree_data$inc_06_16_ht)
t.test(tree_data$inc_11_16_ht, tree_data$inc_06_16_ht)

# make a new dataframe with clean data
tree_data_clean <- tree_data[, c(2, 3, 5, 7, 8, 9, 10, 12, 49, 13, 14, 35, 25, 15, 45, 46, 47, 50:60)]

# PCA on soils

#just data from 60 cm
pca_data <- plot_data[, c(26:45)]

#remove some constant terms
pca_data <- pca_data[, !(colnames(pca_data) %in% c("Mg_60", "Zn_60", "H.Al_60"))]
pca_data <- scale(pca_data)
clean_pca_names <- c("Clay", "Silt", "P", "Org. matter", "pH", 
                     "K", "Ca", "Al", "Sum of bases", "C.E.C.", "Base sat.",
                     "Al. sat", "Cu", "Mn", "Fe", "B", "Sand")


pca <- prcomp(pca_data, scale = TRUE)

# code stolen from someone on CrossValidated; I forget who. 
# Modified from biplot
# to scale scores and loadings to same scale
choices = 1L:3L
scale = 0.5
pc.biplot = FALSE
scores<-pca$x
lam <- pca$sdev[choices]
n <- NROW(scores)
lam <- lam * sqrt(n)
lam <- lam^scale
yy<-t(t(pca$rotation[, choices]) * lam)
xx<-t(t(scores[, choices])/lam)
biplot(xx,yy)

plot_scores <- xx
# vegan::scores(pca, choices = c(1:3), display = c("sites"))

soil_scores <- yy
# vegan::scores(pca, choices = c(1:3), display = c("species"))

# biplot_dat <- stats::biplot(pca, choices = c(1,2), scale = .5)


plot_data$PC1 <- plot_scores[, 1]
plot_data$PC2 <- plot_scores[, 2]

#merge all data together
all_data <- join(tree_data_clean, plot_data, by = c("Plot"))

sp_list <- all_data[which(!duplicated(all_data$Species..names.not.updated.)), ]

#self-thinning
summary(lm(log(density[plot_data$TBA.2006.m2ha.1 > 18, ]$x) ~ log(qmd[plot_data$TBA.2006.m2ha.1 > 18, ]$x)))
abline(coef(lm(log(density[plot_data$TBA.2006.m2ha.1 > 18, ]$x) ~ log(qmd[plot_data$TBA.2006.m2ha.1 > 18, ]$x))))

#do generalists grow faster?
t.test(inc_06_11_ht ~ C.Cerrado..G.generalist., data = all_data) 
t.test(inc_11_16_ht ~ C.Cerrado..G.generalist., data = all_data)


#----------------------------------------------------------------------
# Model comparison for table 1
# height
#----------------------------------------------------------------------
# just species
mod1 <- expression(lmer(inc_06_11_ht ~ (1|Code) + (1|Plot), data = subset(all_data, inc_06_11_ht >= min.growth & dg.2006.cm > 0))) # just species and diameter

mod2 <- expression(lmer(inc_06_11_ht ~ log(dg.2006.cm) + (1|Code) + (1|Plot), data = subset(all_data, inc_06_11_ht >= min.growth & dg.2006.cm > 0)))
  
# linear soils
mod3 <- expression(lmer(inc_06_11_ht ~ scale(PC1) + scale(PC2) + 
               (1|Code) + (1|Plot),
             data = subset(all_data, inc_06_11_ht >= min.growth & dg.2006.cm > 0)))
  
# quadratic soils
mod4 <- expression(lmer(inc_06_11_ht ~ poly(PC1,2)  +
                          poly(PC2,2) +  (1|Code) + (1|Plot),
             data = subset(all_data, inc_06_11_ht >= min.growth & dg.2006.cm > 0)))

# just size and competition
mod5 <- expression((lmer(inc_06_11_ht ~ 
                           scale(log(dg.2006.cm)) +
                           scale(BA.above2006) +
                           (1|Code) + (1|Plot),
                         data = subset(all_data, inc_06_11_ht >= min.growth & dg.2006.cm > 0))))

# let FTs differ in competition susceptibility
mod6 <- expression((lmer(inc_06_11_ht ~ 
                           scale(log(dg.2006.cm)) +
                           scale(BA.above2006) * C.Cerrado..G.generalist. +
                           (1|Code) + (1|Plot),
                         data = subset(all_data, inc_06_11_ht >= min.growth & dg.2006.cm > 0))))

# let FTs also differ in diameter-growth relationships
# doesn't help AIC
mod7 <- expression((lmer(inc_06_11_ht ~ 
                           scale(log(dg.2006.cm))* C.Cerrado..G.generalist. +
                           scale(BA.above2006) * C.Cerrado..G.generalist. +
                           (1|Code) + (1|Plot),
                         data = subset(all_data, inc_06_11_ht >= min.growth & dg.2006.cm > 0))))
  
# add in PC1
mod8 <- expression((lmer(inc_06_11_ht ~ 
                           scale(log(dg.2006.cm)) +
                           scale(BA.above2006) * C.Cerrado..G.generalist. +
                           poly(PC1, 2) + 
                           (1|Code) + (1|Plot),
                         data = subset(all_data, inc_06_11_ht >= min.growth & dg.2006.cm > 0))))
# add in linear PC2  
mod9 <-  expression((lmer(inc_06_11_ht ~ 
                 scale(log(dg.2006.cm)) +
                 scale(BA.above2006) * C.Cerrado..G.generalist. +
                 poly(PC1, 2) + PC2 +
                 (1|Code) + (1|Plot),
               data = subset(all_data, inc_06_11_ht >= min.growth & dg.2006.cm > 0))))

#best model
# quadratic PC2
mod10 <- expression((lmer(inc_06_11_ht ~ 
                            scale(log(dg.2006.cm)) +
                            scale(BA.above2006) * C.Cerrado..G.generalist. +
                            poly(PC1, 2) + poly(PC2, 2) +
                            (1|Code) + (1|Plot),
                          data = subset(all_data, inc_06_11_ht >= min.growth & dg.2006.cm > 0))))

mod11 <- expression((lmer(inc_06_11_ht ~ 
                            scale(log(dg.2006.cm)) +
                            scale(TBA.2006.m2ha.1) * C.Cerrado..G.generalist. +
                            poly(PC1, 2) + poly(PC2, 2) +
                            (1|Code) + (1|Plot),
                          data = subset(all_data, inc_06_11_ht >= min.growth & dg.2006.cm > 0))))

model_list <- c(paste0("mod", seq(1:11)))

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

write.csv(model_table, "./model outputs/model_table_06_11.csv")

#------------------------------------------------------------------------------
# redo for second time period

# just species
mod1 <- expression(lmer(inc_11_16_ht ~ (1|Code) + (1|Plot), data = subset(all_data, inc_11_16_ht >= min.growth & dg.2011.cm > 0))) # just species and diameter

mod2 <- expression(lmer(inc_11_16_ht ~ log(dg.2011.cm) + (1|Code) + (1|Plot), data = subset(all_data, inc_11_16_ht >= min.growth & dg.2011.cm > 0)))

# linear soils
mod3 <- expression(lmer(inc_11_16_ht ~ scale(PC1) + scale(PC2) + 
                          (1|Code) + (1|Plot),
                        data = subset(all_data, inc_11_16_ht >= min.growth & dg.2011.cm > 0)))

# quadratic soils
mod4 <- expression(lmer(inc_11_16_ht ~ poly(PC1,2)  +
                          poly(PC2,2) +  (1|Code) + (1|Plot),
                        data = subset(all_data, inc_11_16_ht >= min.growth & dg.2011.cm > 0)))

# just size and competition
mod5 <- expression((lmer(inc_11_16_ht ~ 
                           scale(log(dg.2011.cm)) +
                           scale(BA.above2011) +
                           (1|Code) + (1|Plot),
                         data = subset(all_data, inc_11_16_ht >= min.growth & dg.2011.cm > 0))))

# let FTs differ in competition susceptibility
mod6 <- expression((lmer(inc_11_16_ht ~ 
                           scale(log(dg.2011.cm)) +
                           scale(BA.above2011) * C.Cerrado..G.generalist. +
                           (1|Code) + (1|Plot),
                         data = subset(all_data, inc_11_16_ht >= min.growth & dg.2011.cm > 0))))

# let FTs also differ in diameter-growth relationships
# doesn't help AIC
mod7 <- expression((lmer(inc_11_16_ht ~ 
                           scale(log(dg.2011.cm))* C.Cerrado..G.generalist. +
                           scale(BA.above2011) * C.Cerrado..G.generalist. +
                           (1|Code) + (1|Plot),
                         data = subset(all_data, inc_11_16_ht >= min.growth & dg.2011.cm > 0))))

# add in PC1
mod8 <- expression((lmer(inc_11_16_ht ~ 
                           scale(log(dg.2011.cm)) +
                           scale(BA.above2011) * C.Cerrado..G.generalist. +
                           poly(PC1, 2) + 
                           (1|Code) + (1|Plot),
                         data = subset(all_data, inc_11_16_ht >= min.growth & dg.2011.cm > 0))))
# add in linear PC2  
mod9 <-  expression((lmer(inc_11_16_ht ~ 
                            scale(log(dg.2011.cm)) +
                            scale(BA.above2011) * C.Cerrado..G.generalist. +
                            poly(PC1, 2) + PC2 +
                            (1|Code) + (1|Plot),
                          data = subset(all_data, inc_11_16_ht >= min.growth & dg.2011.cm > 0))))

#best model
# quadratic PC2
mod10 <- expression((lmer(inc_11_16_ht ~ 
                            scale(log(dg.2011.cm)) +
                            scale(BA.above2011) * C.Cerrado..G.generalist. +
                            poly(PC1, 2) + poly(PC2, 2) +
                            (1|Code) + (1|Plot),
                          data = subset(all_data, inc_11_16_ht >= min.growth & dg.2011.cm > 0))))

mod11 <- expression((lmer(inc_11_16_ht ~ 
                            scale(log(dg.2011.cm)) +
                            scale(TBA.2011.m2.ha.1) * C.Cerrado..G.generalist. +
                            poly(PC1, 2) + poly(PC2, 2) +
                            (1|Code) + (1|Plot),
                          data = subset(all_data, inc_11_16_ht >= min.growth & dg.2011.cm > 0))))

model_list <- c(paste0("mod", seq(1:11)))

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

write.csv(model_table, "./model outputs/model_table_11_16.csv")

#----------------------------------------------------------------------
# Full model: FT x CI interaction with soils
# = model 10 from above
# all_data$inc_06_11_ht_rel <- all_data$inc_06_11_ht/all_data$dg.2006.cm
# all_data$inc_11_16_ht_rel <- all_data$inc_11_16_ht/all_data$dg.2011.cm
options(na.action = na.omit)
data_06 <- subset(all_data, inc_06_11_ht >= min.growth & dg.2006.cm > 0)
data_11 <- subset(all_data, inc_11_16_ht >= min.growth & dg.2011.cm > 0)

ft_ci_soils_06_ht <- lmer(inc_06_11_ht ~ 
                         scale(H.2006) +
                         scale(BA.above2006) * C.Cerrado..G.generalist. +
                         # poly(PC1, 2) + poly(PC2, 2) +
                           PC1 + PC2 + 
                         (1|Code) + (1|Plot),
                       data = data_06)

ft_ci_soils_11_ht <- lmer(inc_11_16_ht ~ 
                         scale(H.2011..m.) +
                         scale(BA.above2011) * C.Cerrado..G.generalist. +
                         # poly(PC1, 2) + poly(PC2, 2) +
                         (1|Code) + (1|Plot),
                       data = data_11)

summary(ft_ci_soils_06_ht)
summary(ft_ci_soils_11_ht)

write.csv(summary(ft_ci_soils_06_ht)$coef, "./model outputs/model_coefs_06_ht.csv")
write.csv(summary(ft_ci_soils_11_ht)$coef, "./model outputs/model_coefs_11_ht.csv")


plot(residuals(ft_ci_soils_06_ht) ~ fitted(ft_ci_soils_06_ht))
abline(h = 0)
ggqqplot(residuals(ft_ci_soils_06_ht))
plot(allEffects(ft_ci_soils_06_ht, partial.residuals = FALSE))
plot(Effect(focal.predictors = c("BA.above2006", "C.Cerrado..G.generalist."),
            mod = ft_ci_soils_06_ht))
plot(Effect(focal.predictors = c("dg.2006.cm", "C.Cerrado..G.generalist."),
            mod = ft_ci_soils_0_ht))
plot(Effect(focal.predictors = c("PC1"),
            mod = ft_ci_soils_06_ht))

plot(residuals(ft_ci_soils_11_ht) ~ fitted(ft_ci_soils_11_ht))
abline(h = 0)
ggqqplot(residuals(ft_ci_soils_11_ht))
plot(allEffects(ft_ci_soils_11_ht, partial.residuals = FALSE))
plot(Effect(focal.predictors = c("BA.above2006", "C.Cerrado..G.generalist."), 
            mod = ft_ci_soils_11_ht, transformation = list(link = log, inverse = exp)))
plot(Effect(focal.predictors = c("dg.2006.cm", "C.Cerrado..G.generalist."), 
            mod = ft_ci_soils_11_ht, transformation = list(link = log, inverse = exp)))
plot(Effect(focal.predictors = c("PC1"), 
            mod = ft_ci_soils_11_ht, transformation = list(link = log, inverse = exp)))

#making sure everything converged okay; I think I stole this from Ben Bolker
# derivs1 <- ft_ci_soils_06@optinfo$derivs
# sc_grad1 <- with(derivs1,solve(Hessian,gradient))
# max(abs(sc_grad1))
# max(pmin(abs(sc_grad1),abs(derivs1$gradient)))
# 
# relgrad <- with(ft_ci_soils_06@optinfo$derivs,solve(Hessian,gradient))
# max(abs(relgrad))

# qqmath(ranef(ft_ci_soils_06, condVar = TRUE), strip = FALSE)$Code

## Species growth rates (appendix table 1)
ranef_sp_06 <- exp(ranef(ft_ci_soils_06_ht, condVar = TRUE)$Code)
ranef_sp_11 <- exp(ranef(ft_ci_soils_11_ht, condVar = TRUE)$Code)
ranef_sp_df_06 <- data.frame(rel_growth_06_ht = ranef_sp_06_ht$`(Intercept)`,
                          Code = rownames(ranef_sp_06_ht))
ranef_sp_df_11 <- data.frame(rel_growth_11_ht = ranef_sp_11_ht$`(Intercept)`,
                             Code = rownames(ranef_sp_11_ht))

ranef_sp_df_ht <- join(ranef_sp_df_06_ht, ranef_sp_df_11_ht, by = "Code", type = "full")

ranef_sp_df_ht <- join(ranef_sp_df_ht, 
                    tree_data[, c("Code",
                                  "Species..names.not.updated.",
                                  "Family",
                                  "C.Cerrado..G.generalist.",
                                  "Shade.Tolerance")],
                    by = c("Code"),
                    type = "left",
                    match = "first")

ranef_sp_df_ht <- ranef_sp_df_ht[order(ranef_sp_df_ht$rel_growth_06), ]

sp_table <- as.data.frame(table(all_data$Code))
names(sp_table) <- c("Code", "Freq")
ranef_sp_df_ht <- join(ranef_sp_df_ht, sp_table, by = c("Code"))

#update names using Flora package
old_names <- sapply(ranef_sp_df_ht$Species..names.not.updated., FUN = remove.authors)
new_names <- sapply(old_names, FUN = suggest.names)
new_names <- get.taxa(new_names)$`scientific.name`

ranef_sp_df_ht$Species..names.not.updated. <- new_names

write.csv(ranef_sp_df_ht, "./model outputs/species_growth_rates_height.csv")


#-------------------------------------------------------------------------------
# height growth figure
#-------------------------------------------------------------------------------
data_11 <- subset(all_data, inc_11_16_ht >= min.growth & dg.2011.cm > 0) %>%
  dplyr::rename(height = H.2011..m., diam = dg.2011.cm)

ft_ci_soils_11_ht <- lmer(inc_11_16_ht ~ 
                         scale(height) +
                         scale(BA.above2011) * C.Cerrado..G.generalist. +
                         poly(PC1, 2) + poly(PC2, 2) +
                         (1|Code) + (1|Plot),
                       data = data_11)

ft_ci_soils_11_diam <- lmer(inc_11_16 ~ 
                         scale(log(diam)) +
                         scale(BA.above2011) * C.Cerrado..G.generalist. +
                         poly(PC1, 2) + poly(PC2, 2) +
                         (1|Code),
                       data = data_11)

nyear <- 20
c_growth <- data.frame(year = numeric(nyear),
                       height = numeric(nyear))

  height <- mean(all_data[all_data$C.Cerrado..G.generalist. == "C" & all_data$Added == 2011, ]$H.2011..m., na.rm = TRUE)
  ft <- "C"
  PC1 <- mean(plot_data$PC1)
  PC2 <- mean(plot_data$PC2)
  BA.above2011 <- 0
for(year in 1:nyear){
  newdat <- data.frame(height = height,
                       C.Cerrado..G.generalist. = ft,
                       PC1 = PC1,
                       PC2 = PC2,
                       BA.above2011 = BA.above2011)
  growth <- predict(ft_ci_soils_11_ht, newdata = newdat, re.form = NA)
  height <- height + growth
  c_growth$year[year] <- year
  c_growth$height[year] <- height
  
}

g_growth <- data.frame(year = numeric(nyear),
                       height = numeric(nyear))

  height <- mean(all_data[all_data$C.Cerrado..G.generalist. == "G" & all_data$Added == 2011, ]$H.2011..m., na.rm = TRUE)
  ft <- "G"
  PC1 <- mean(plot_data$PC1)
  PC2 <- mean(plot_data$PC2)
  BA.above2011 <- 0
for(year in 1:nyear){
  newdat <- data.frame(height = height,
                       C.Cerrado..G.generalist. = ft,
                       PC1 = PC1,
                       PC2 = PC2,
                       BA.above2011 = BA.above2011)
  growth <- predict(ft_ci_soils_11_ht, newdata = newdat, re.form = NA)
  height <- height + growth
  g_growth$year[year] <- year
  g_growth$height[year] <- height
}

plot(c_growth$height ~ c_growth$year, type = "l", ylim = c(0,15))
lines(g_growth$height ~ g_growth$year, lty = 2)

#-------------------------------------------------------------------------------
# simulating a whole stand

recruits <- data_11 %>%
  dplyr::filter(Added == 2011 & !is.na(height) & !is.na(diam))

new_stand <- data.frame(Current.number = character(),
                        Year = numeric(),
                        height = numeric(),
                        C.Cerrado..G.generalist. = character(),
                        diam = numeric(),
                        PC1 = numeric(),
                        PC2 = numeric(),
                        BA.above2011 = numeric())
BAtot <- 0

while(BAtot < 5){
  c_tree <- recruits[recruits$C.Cerrado..G.generalist. == "C", ] %>%
    '['(sample(nrow(.), 1), ) %>% #this might be clunky, but it selects one row at random
    dplyr::select(Current.number, height, C.Cerrado..G.generalist., diam) %>%
    dplyr::mutate(PC1 = mean(plot_data$PC1),
                  PC2 = mean(plot_data$PC2),
                  BA.above2011 = 0,
                  Year = 0) %>%
    dplyr::mutate(height = height + rnorm(n = nrow(c_tree), mean = 0, sd = 0.02))
  g_tree <- recruits[recruits$C.Cerrado..G.generalist. == "G", ] %>%
    '['(sample(nrow(.), 1), ) %>%
    dplyr::select(Current.number, height, C.Cerrado..G.generalist., diam) %>%
    dplyr::mutate(PC1 = mean(plot_data$PC1),
                  PC2 = mean(plot_data$PC2),
                  BA.above2011 = 0,
                  Year = 0) %>%
    dplyr::mutate(height = height + rnorm(n = nrow(g_tree), mean = 0, sd = 0.02))
  
  new_stand <- dplyr::bind_rows(new_stand, c_tree, g_tree)
  
  BAtot <- sum((new_stand$diam / 2)^2 * pi) / 10000 #convert to m^2 
  BAtot <- BAtot * 10 #convert to equivalent BA for a 0.1-ha plot
  
}

for(i in 1:nrow(new_stand)){
  new_stand$BA.above2011[i] <- sum(((new_stand[new_stand$height >= new_stand$height[i], "diam"]/2)^2 * pi), na.rm = TRUE) / 10000
}

stand_data <- new_stand
nyear <- 20
for(year in 1:nyear){
  height_inc <- predict(ft_ci_soils_11_ht, newdata = stand_data[stand_data$Year == (year - 1), ], re.form = NA)
  diam_inc <- predict(ft_ci_soils_11_diam, newdata = stand_data[stand_data$Year == (year - 1), ], re.form = NA)

  update <- stand_data[stand_data$Year == (year - 1), ] %>%
    mutate(Year = year,
           height = height + height_inc,
           diam = diam + diam_inc)
  
  for(i in 1:nrow(update)){
    update$BA.above2011[i] <- sum((update[update$height >= update$height[i], "diam"]/2)^2 * pi, na.rm = TRUE) / 10000
  }
  
  stand_data <- dplyr::bind_rows(stand_data, update)
}

ggplot(stand_data, aes(x = Year, y = height)) + 
  geom_jitter(aes(col = C.Cerrado..G.generalist.)) +
  stat_smooth(aes(col = C.Cerrado..G.generalist.), method = "loess")

ggplot(stand_data, aes(x = Year, y = diam)) + 
  geom_jitter(aes(col = C.Cerrado..G.generalist.)) +
  stat_smooth(aes(col = C.Cerrado..G.generalist.), method = "loess")

ggplot(stand_data, aes(x = Year, y = BA.above2011)) + 
  # geom_jitter(aes(col = C.Cerrado..G.generalist.)) +
  geom_line(aes(group = Current.number, col = C.Cerrado..G.generalist.))+ 
  stat_smooth(aes(col = C.Cerrado..G.generalist.), method = "loess")



#test to check if a repeat measures model would work
all_data$Unique_id <- paste0(all_data$Plot, ".", all_data$Current.number, ".", 1:nrow(all_data))

#there absolutely has to be a way to pivot_longer this, but I couldn't figure it out
test <- data.frame(Plot = character(),
                   Unique_id = character(),
                   FT = character(),
                   Species = character(),
                   Year = integer(),
                   Height = numeric(),
                   Diam = numeric(),
                   BA_above = numeric(),
                   PC1 = numeric(),
                   PC2 = numeric(),
                   diam_inc = numeric(),
                   ht_inc = numeric())
for(i in 1:length(unique(all_data$Unique_id))){
  id <- all_data$Unique_id[i]
  test[(2*i-1):(2*i), "Plot"] <- subset(all_data, Unique_id == id)$"Plot"
  test[(2*i-1):(2*i), "Unique_id"] <- id
  test[(2*i-1):(2*i), "FT"] <- subset(all_data, Unique_id == id)$"C.Cerrado..G.generalist."
  test[(2*i-1):(2*i), "Species"] <- subset(all_data, Unique_id == id)$Code
  test[(2*i-1):(2*i), "Year"] <- c(2006, 2011)
  test[(2*i-1):(2*i), "PC1"] <- subset(all_data, Unique_id == id)$PC1
  test[(2*i-1):(2*i), "PC2"] <- subset(all_data, Unique_id == id)$PC2
  test[(2*i-1), "Height"] <- subset(all_data, Unique_id == id)$H.2006
  test[(2*i), "Height"] <- subset(all_data, Unique_id == id)$H.2011..m.
  test[(2*i-1), "Diam"] <- subset(all_data, Unique_id == id)$dg.2006.cm
  test[(2*i), "Diam"] <- subset(all_data, Unique_id == id)$dg.2011.cm
  test[(2*i-1), "BA_above"] <- subset(all_data, Unique_id == id)$BA.above2006
  test[(2*i), "BA_above"] <- subset(all_data, Unique_id == id)$BA.above2011
  test[(2*i-1), "diam_inc"] <- subset(all_data, Unique_id == id)$inc_06_11
  test[(2*i), "diam_inc"] <- subset(all_data, Unique_id == id)$inc_11_16
  test[(2*i-1), "ht_inc"] <- subset(all_data, Unique_id == id)$inc_06_11_ht
  test[(2*i), "ht_inc"] <- subset(all_data, Unique_id == id)$inc_11_16_ht
}

model <- lmer(ht_inc ~ Height + BA_above*FT + PC1 + PC2 +
                (1|Plot) + (1|Species) + (1|Year),
              data = test)

summary(model)

