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
  if(year == 2016) return(plot[plot$Alive2016 == "v", c("Current.number", "H.2016.m", "G2016..m2.")])
  
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
tree_data$G2016..m2. <- as.numeric(tree_data$G2016..m2.)
tree_data$dg.2006.cm <- as.numeric(tree_data$dg.2006.cm)
tree_data$dg.2016.cm <- as.numeric(tree_data$dg.2016.cm)
tree_data$dg.2011.cm <- as.numeric(tree_data$dg.2011.cm)

# Calculate the one-sided competition index for 2006 and for 2011
# this uses the entire dataset, not removing trees with negative growth 
# and using only live trees

tree_data$BA.above2006 <- NA
tree_data$BA.above2006_prop <- NA #finish this!
tree_data$BA.above2011 <- NA
tree_data$BA.above2011_prop <- NA
tree_data$BA.above2016 <- NA
tree_data$BA.above2016_prop <- NA

for(i in 1:nrow(tree_data)){
  tree_no <- tree_data[i, "Current.number"]
  tree_h <- tree_data[i, "H.2006"]
  plot_no <- tree_data[i, "Plot"]
  trees <- get_heights_csas(plot_no, 2006)
  trees <- trees[trees$Current.number != tree_no, ]
  #BA.above is the sum of basal area of taller trees
  if(!is.na(tree_h)){
    tree_data$BA.above2006[i] <- sum(trees[trees$H.2006 >= tree_h, "G2006..m2."], na.rm = TRUE)
    tree_data$BA.above2006_prop[i] <- tree_data$BA.above2006[i] / sum(trees$G2006..m2., na.rm = TRUE)
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
    tree_data$BA.above2011_prop[i] <- tree_data$BA.above2011[i] / sum(trees$G2011..m2., na.rm = TRUE)
  }
}


for(i in 1:nrow(tree_data)){
  tree_no <- tree_data[i, "Current.number"]
  tree_h <- as.numeric(tree_data[i, "H.2016.m"])
  plot_no <- tree_data[i, "Plot"]
  trees <- get_heights_csas(plot_no, 2016)
  trees <- trees[trees$Current.number != tree_no, ]
  #BA.above is the sum of basal area of taller trees
  if(!is.na(tree_h)){
    tree_data$BA.above2016[i] <- sum(trees[trees$H.2016.m >= tree_h, "G2016..m2."], na.rm = TRUE)
    tree_data$BA.above2016_prop[i] <- tree_data$BA.above2016[i] / sum(trees$G2016..m2., na.rm = TRUE)
  }
}


tree_data$Died2011 <- as.factor(ifelse(tree_data$Alive2006 == "v" & tree_data$Alive2011 == "m",
                                       "y", "n"))
tree_data$Died2016 <- as.factor(ifelse(tree_data$Alive2011 == "v" & tree_data$Alive2016 == "m",
                                       "y", "n"))
tree_data$Died_all <- as.factor(ifelse(tree_data$Died2011 == "y" | tree_data$Died2016 == "y", "y", "n"))


# clean up data to estimate growth

# calculate growth separately for each time interval
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
t.test(tree_data$inc_06_11, tree_data$inc_11_16)
t.test(tree_data$inc_06_11, tree_data$inc_06_16)
t.test(tree_data$inc_11_16, tree_data$inc_06_16)

# make a new dataframe with clean data
tree_data_clean <- tree_data[, c(2, 3, 5, 7, 8, 9, 10, 12, 49, 13, 14, 35, 25, 15, 45, 46, 47, 50:64)]

#-------------------------------------------------------------------------------
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
# plot_data$PC1_std <- pca_scores_std[, 1]
# plot_data$PC2_std <- pca_scores_std[, 2]


#merge all data together
all_data <- join(tree_data_clean, plot_data, by = c("Plot"))

sp_list <- all_data[which(!duplicated(all_data$Species..names.not.updated.)), ]

#how does QMD change over the gradient?
qmd <- aggregate(as.numeric(tree_data$dg.2006.cm), by = list(tree_data$Plot), FUN = function(x){sqrt(mean(x^2, na.rm = TRUE))})
plot(qmd$x ~ plot_data$TBA.2006.m2ha.1)

density <- aggregate(tree_data$dg.2006.cm, by = list(tree_data$Plot), FUN = function(x){sum(x, na.rm = TRUE)})

plot(log(density$x) ~ log(qmd$x))

#self-thinning
summary(lm(log(density[plot_data$TBA.2006.m2ha.1 > 18, ]$x) ~ log(qmd[plot_data$TBA.2006.m2ha.1 > 18, ]$x)))
abline(coef(lm(log(density[plot_data$TBA.2006.m2ha.1 > 18, ]$x) ~ log(qmd[plot_data$TBA.2006.m2ha.1 > 18, ]$x))))

#do generalists grow faster?
t.test(inc_06_11 ~ C.Cerrado..G.generalist., data = all_data) 
t.test(inc_11_16 ~ C.Cerrado..G.generalist., data = all_data)

#------------------------------------------------------------------------------
# Make a figure of the PCA
#------------------------------------------------------------------------------
library("viridis")

tiff(filename="./plots/Figure S2 soil_PCA_axes_1_2.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 4, 
     height= 4, 
     pointsize=9, 
     res=600)

par(mar = c(3.9,3.9,1,0),
    oma = c(0,0,0,0))

col <- viridis(30)[order(plot_data$TBA.2006.m2ha.1, decreasing = TRUE)]

layout(mat = matrix(data = c(1,2), nrow = 1, ncol = 2),
       widths = c(.85, .15))

plot(NA,
     xlab = "", ylab = "",
     xaxt = "n",
     yaxt = "n",
     xlim = c( (min(c(plot_scores[, 1], soil_scores[, 1]) - 0.1)), 
               (max(c(plot_scores[, 1], soil_scores[, 1])) + 0.1)),
     ylim = c( (min(c(plot_scores[, 2], soil_scores[, 2]) - 0.1)), 
               (max(c(plot_scores[, 2], soil_scores[, 2])) + 0.1)))
  
  abline(h = 0)
  abline(v = 0)
  axis(1)
  mtext(side = 1, line = 2, text = "PC1 (43%)")
  axis(2)
  mtext(side = 2, line = 2, text = "PC2 (31%)")
  # x<-pointLabel(x = plot_scores[, 1], y = plot_scores[, 2], labels = as.character(plot_data$Plot), cex = 0.7)
  
  points(x = plot_scores[, 1], y = plot_scores[, 2], cex = 1,
         col = col,
         bg = col,
         pch = 21)
  
  
  segments(x0 = 0, y0 = 0, x1 = soil_scores[, 1]*.9, y1 = soil_scores[, 2]*.9)
  yoffsets <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0.05,0,0,0)
  xoffsets <- c(0,0,0,0,0,0,0,0,-0.2,0,0.2,0,0,0,0,0,0)
  
  text(x = soil_scores[, 1]+xoffsets, y = soil_scores[, 2] + yoffsets, 
       labels = clean_pca_names, cex = 1)
  
  par(mar = c(0,0,0,0))
  
  plot.new()
  plot.window(ylim = c(0,10), xlim = c(0,10))
  
  
  lgd_ = rep(NA, 30)
  lgd_[c(1,30)] = c(round(max(plot_data$TBA.2006.m2ha.1),1), 
                    round(min(plot_data$TBA.2006.m2ha.1),1))
  legend(x = 0, y = 7,
         legend = lgd_,
         fill = viridis(30),
         border = NA,
         bty = "n",
         y.intersp = 0.2,
         x.intersp = .5,
         cex = 1, text.font = 2)
  text(x = 0, y = 8, labels = "Basal\narea", 
       pos = 4, cex = 1)
  text(x = 0, y = 7.4, labels = expression(paste("(m"^"2"," ha"^"-1",")")), 
       pos = 4, cex = 1)

dev.off()


#----------------------------------------------------------------------
# Model comparison for table S2
#----------------------------------------------------------------------

data_06 <- subset(all_data, inc_06_11 >= min.growth & inc_06_11 <= 3 &  dg.2006.cm > 5)

data_11 <- subset(all_data, inc_11_16 >= min.growth & inc_06_11 <= 3 & dg.2011.cm > 5 &
                    !(Current.number %in% c(2775, 2030)))


#----------------------------------------------------------------------
# Model comparison for table 1
#----------------------------------------------------------------------
# just species
mod1 <- expression(lmer(inc_06_11 ~ (1|Code) + (1|Plot), data = data_06)) # just species and diameter

mod2 <- expression(lmer(inc_06_11 ~ log(dg.2006.cm) + (1|Code) + (1|Plot), data = data_06))

# linear soils
mod3 <- expression(lmer(inc_06_11 ~ scale(PC1) + scale(PC2) + 
                          (1|Code) + (1|Plot),
                        data = data_06))

# quadratic soils
mod4 <- expression(lmer(inc_06_11 ~ poly(PC1,2)  +
                          poly(PC2,2) +  (1|Code) + (1|Plot),
                        data = data_06))

# just size and competition
mod5 <- expression((lmer(inc_06_11 ~ 
                           scale(log(dg.2006.cm)) +
                           scale(BA.above2006) +
                           (1|Code) + (1|Plot),
                         data = data_06)))

# let FTs differ in competition susceptibility
mod6 <- expression((lmer(inc_06_11 ~ 
                           scale(log(dg.2006.cm)) +
                           scale(BA.above2006) * C.Cerrado..G.generalist. +
                           (1|Code) + (1|Plot),
                         data = data_06)))

# let FTs also differ in diameter-growth relationships
# doesn't help AIC
mod7 <- expression((lmer(inc_06_11 ~ 
                           scale(log(dg.2006.cm))* C.Cerrado..G.generalist. +
                           scale(BA.above2006) * C.Cerrado..G.generalist. +
                           (1|Code) + (1|Plot),
                         data = data_06)))

# add in PC1
mod8 <- expression((lmer(inc_06_11 ~ 
                           scale(log(dg.2006.cm)) +
                           scale(BA.above2006) * C.Cerrado..G.generalist. +
                           poly(PC1, 2) + 
                           (1|Code) + (1|Plot),
                         data = data_06)))
# add in linear PC2  
mod9 <-  expression((lmer(inc_06_11 ~ 
                            scale(log(dg.2006.cm)) +
                            scale(BA.above2006) * C.Cerrado..G.generalist. +
                            poly(PC1, 2) + PC2 +
                            (1|Code) + (1|Plot),
                          data = data_06)))

#best model
# quadratic PC2
mod10 <- expression((lmer(inc_06_11 ~ 
                            scale(log(dg.2006.cm)) +
                            scale(BA.above2006) * C.Cerrado..G.generalist. +
                            poly(PC1, 2) + poly(PC2, 2) +
                            (1|Code) + (1|Plot),
                          data = data_06)))

mod11 <- expression((lmer(inc_06_11 ~ 
                            scale(log(dg.2006.cm)) +
                            scale(TBA.2006.m2ha.1) * C.Cerrado..G.generalist. +
                            poly(PC1, 2) + poly(PC2, 2) +
                            (1|Code) + (1|Plot),
                          data = data_06)))

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
mod1 <- expression(lmer(inc_11_16 ~ (1|Code) + (1|Plot), data = data_11)) # just species and diameter

mod2 <- expression(lmer(inc_11_16 ~ log(dg.2011.cm) + (1|Code) + (1|Plot), data = data_11))

# linear soils
mod3 <- expression(lmer(inc_11_16 ~ scale(PC1) + scale(PC2) + 
                          (1|Code) + (1|Plot),
                        data = data_11))

# quadratic soils
mod4 <- expression(lmer(inc_11_16 ~ poly(PC1,2)  +
                          poly(PC2,2) +  (1|Code) + (1|Plot),
                        data = data_11))

# just size and competition
mod5 <- expression((lmer(inc_11_16 ~ 
                           scale(log(dg.2011.cm)) +
                           scale(BA.above2011) +
                           (1|Code) + (1|Plot),
                         data = data_11)))

# let FTs differ in competition susceptibility
mod6 <- expression((lmer(inc_11_16 ~ 
                           scale(log(dg.2011.cm)) +
                           scale(BA.above2011) * C.Cerrado..G.generalist. +
                           (1|Code) + (1|Plot),
                         data = data_11)))

# let FTs also differ in diameter-growth relationships
# doesn't help AIC
mod7 <- expression((lmer(inc_11_16 ~ 
                           scale(log(dg.2011.cm))* C.Cerrado..G.generalist. +
                           scale(BA.above2011) * C.Cerrado..G.generalist. +
                           (1|Code) + (1|Plot),
                         data = data_11)))

# add in PC1
mod8 <- expression((lmer(inc_11_16 ~ 
                           scale(log(dg.2011.cm)) +
                           scale(BA.above2011) * C.Cerrado..G.generalist. +
                           poly(PC1, 2) + 
                           (1|Code) + (1|Plot),
                         data = data_11)))
# add in linear PC2  
mod9 <-  expression((lmer(inc_11_16 ~ 
                            scale(log(dg.2011.cm)) +
                            scale(BA.above2011) * C.Cerrado..G.generalist. +
                            poly(PC1, 2) + PC2 +
                            (1|Code) + (1|Plot),
                          data = data_11)))

#best model
# quadratic PC2
mod10 <- expression((lmer(inc_11_16 ~ 
                            scale(log(dg.2011.cm)) +
                            scale(BA.above2011) * C.Cerrado..G.generalist. +
                            poly(PC1, 2) + poly(PC2, 2) +
                            (1|Code) + (1|Plot),
                          data = data_11)))

mod11 <- expression((lmer(inc_11_16 ~ 
                            scale(log(dg.2011.cm)) +
                            scale(TBA.2011.m2.ha.1) * C.Cerrado..G.generalist. +
                            poly(PC1, 2) + poly(PC2, 2) +
                            (1|Code) + (1|Plot),
                          data = data_11)))

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
# all_data$inc_06_11_rel <- all_data$inc_06_11/all_data$dg.2006.cm
# all_data$inc_11_16_rel <- all_data$inc_11_16/all_data$dg.2011.cm
options(na.action = na.omit)
min.growth <- -0.2

data_06 <- subset(all_data, inc_06_11 >= min.growth & inc_06_11 <= 3 &  dg.2006.cm > 5)

data_11 <- subset(all_data, inc_11_16 >= min.growth & inc_06_11 <= 3 & dg.2011.cm > 5 &
                    !(Current.number %in% c(2775, 2030)))


ft_ci_soils_06 <- lmer(inc_06_11 ~ 
                         log(dg.2006.cm) +
                         scale(BA.above2006) * C.Cerrado..G.generalist. +
                         poly(PC1, 2) + poly(PC2, 2) +
                         (1|Code) + (1|Plot),
                       data = data_06)


ft_ci_soils_11 <- lmer(inc_11_16~ 
                         scale(log(dg.2011.cm)) +
                         scale(BA.above2011) * C.Cerrado..G.generalist. +
                         poly(PC1, 2) + poly(PC2, 2) +
                         (1|Code) + (1|Plot),
                       data = data_11)

summary(ft_ci_soils_06, ddf = "Satterthwaite")
summary(ft_ci_soils_11)

write.csv(summary(ft_ci_soils_06)$coef, "./model outputs/model_coefs_06_diam.csv")
write.csv(summary(ft_ci_soils_11)$coef, "./model outputs/model_coefs_11_diam.csv")

plot(residuals(ft_ci_soils_06) ~ fitted(ft_ci_soils_06))
# plot(residuals(ft_ci_soils_06) ~ ft_ci_soils_06@frame$`scale(BA.above2006)`)
abline(h = 0)

ggqqplot(residuals(ft_ci_soils_06))
plot(allEffects(ft_ci_soils_06, partial.residuals = FALSE))
plot(Effect(focal.predictors = c("BA.above2006", "C.Cerrado..G.generalist."), 
            mod = ft_ci_soils_06))
plot(Effect(focal.predictors = c("dg.2006.cm", "C.Cerrado..G.generalist."), 
            mod = ft_ci_soils_06))
plot(Effect(focal.predictors = c("PC1"), 
            mod = ft_ci_soils_06))

plot(residuals(ft_ci_soils_11) ~ fitted(ft_ci_soils_11))
abline(h = 0)
ggqqplot(residuals(ft_ci_soils_11))
plot(allEffects(ft_ci_soils_11, partial.residuals = FALSE))
plot(Effect(focal.predictors = c("BA.above2011", "C.Cerrado..G.generalist."), 
            mod = ft_ci_soils_11, transformation = list(link = log, inverse = exp)))
plot(Effect(focal.predictors = c("dg.2011.cm", "C.Cerrado..G.generalist."), 
            mod = ft_ci_soils_11, transformation = list(link = log, inverse = exp)))
plot(Effect(focal.predictors = c("PC1"), 
            mod = ft_ci_soils_11, transformation = list(link = log, inverse = exp)))

#making sure everything converged okay; I think I stole this from Ben Bolker
# derivs1 <- ft_ci_soils_11@optinfo$derivs
# sc_grad1 <- with(derivs1,solve(Hessian,gradient))
# max(abs(sc_grad1))
# max(pmin(abs(sc_grad1),abs(derivs1$gradient)))
# 
# relgrad <- with(ft_ci_soils_11@optinfo$derivs,solve(Hessian,gradient))
# max(abs(relgrad))

# qqmath(ranef(ft_ci_soils_11, condVar = TRUE), strip = FALSE)$Code

## Species growth rates (appendix table 1)
ranef_sp_06 <- exp(ranef(ft_ci_soils_06, condVar = TRUE)$Code)
ranef_sp_11 <- exp(ranef(ft_ci_soils_11, condVar = TRUE)$Code)
ranef_sp_df_06 <- data.frame(rel_growth_06 = ranef_sp_06$`(Intercept)`,
                          Code = rownames(ranef_sp_06))
ranef_sp_df_11 <- data.frame(rel_growth_11 = ranef_sp_11$`(Intercept)`,
                             Code = rownames(ranef_sp_11))

ranef_sp_df <- join(ranef_sp_df_06, ranef_sp_df_11, by = "Code", type = "full")

ranef_sp_df <- join(ranef_sp_df, 
                    tree_data[, c("Code",
                                  "Species..names.not.updated.",
                                  "Family",
                                  "C.Cerrado..G.generalist.",
                                  "Shade.Tolerance")],
                    by = c("Code"),
                    type = "left",
                    match = "first")

ranef_sp_df <- ranef_sp_df[order(ranef_sp_df$rel_growth_06), ]

sp_table <- as.data.frame(table(all_data$Code))
names(sp_table) <- c("Code", "Freq")
ranef_sp_df <- join(ranef_sp_df, sp_table, by = c("Code"))

#update names using Flora package
old_names <- sapply(ranef_sp_df$Species..names.not.updated., FUN = remove.authors)
new_names <- sapply(old_names, FUN = suggest.names)
new_names <- get.taxa(new_names)$`scientific.name`

ranef_sp_df$Species..names.not.updated. <- new_names

write.csv(ranef_sp_df, "./model outputs/species_growth_rates.csv")



#-----
# bootstrap standard errors
# code from http://svmiller.com/blog/2020/03/bootstrap-standard-errors-in-r/
# 
# detach("package:lmerTest", unload = TRUE)
# 
# n_samples <- 5000
# 
# data_06 %>%
#   modelr::bootstrap(n_samples) %>%
#   dplyr::mutate(lm = purrr::map(strap, ~lmer(inc_06_11 ~ 
#                                  log(dg.2006.cm) +
#                                  scale(BA.above2006) * C.Cerrado..G.generalist. +
#                                  poly(PC1, 2) + poly(PC2, 2) +
#                                  (1|Code) + (1|Plot),
#                                data = .)),
#          tidy = purrr::map(lm, broom.mixed::tidy)) -> bootCrime
# 
# bootCrime %>%
#   pull(tidy) %>%
#   purrr::map2_df(., # map to return a data frame
#           seq(1, n_samples), # make sure to get this seq right. We did this 1000 times.
#           ~mutate(.x, resample = .y)) -> tidybootCrime
# 
# tidybootCrime %>%
#   # group by term, naturally
#   dplyr::group_by(term) %>%
#   # This is the actual bootstrapped standard error you want
#   dplyr::summarize(bse = sd(estimate)) -> bseM1
# 
# 
# 
# data_11 %>%
#   modelr::bootstrap(n_samples) %>%
#   mutate(lm = purrr::map(strap, ~lmer(inc_11_16~ 
#                                       scale(log(dg.2011.cm)) +
#                                       scale(BA.above2011) * C.Cerrado..G.generalist. +
#                                       poly(PC1, 2) + poly(PC2, 2) +
#                                       (1|Code) + (1|Plot),
#                                       data = .)),
#          tidy = purrr::map(lm, broom.mixed::tidy)) -> bootCrime
# 
# bootCrime %>%
#   pull(tidy) %>%
#   purrr::map2_df(., # map to return a data frame
#                  seq(1, n_samples), # make sure to get this seq right. We did this 1000 times.
#                  ~mutate(.x, resample = .y)) -> tidybootCrime
# 
# tidybootCrime %>%
#   # group by term, naturally
#   dplyr::group_by(term) %>%
#   # This is the actual bootstrapped standard error you want
#   dplyr::summarize(bse = sd(estimate)) -> bseM2


################################################################################
# Model height growth
################################################################################

#----------------------------------------------------------------------
# Model comparison for table 1
# height
#----------------------------------------------------------------------
min.growth <- -0.5

# just species
mod1 <- expression(lmer(inc_06_11_ht ~ (1|Code) + (1|Plot), data = subset(all_data, inc_06_11_ht >= min.growth & dg.2006.cm > 5))) # just species and diameter

mod2 <- expression(lmer(inc_06_11_ht ~ log(dg.2006.cm) + (1|Code) + (1|Plot), data = subset(all_data, inc_06_11_ht >= min.growth & dg.2006.cm > 5)))

# linear soils
mod3 <- expression(lmer(inc_06_11_ht ~ scale(PC1) + scale(PC2) + 
                          (1|Code) + (1|Plot),
                        data = subset(all_data, inc_06_11_ht >= min.growth & dg.2006.cm > 5)))

# quadratic soils
mod4 <- expression(lmer(inc_06_11_ht ~ poly(PC1,2)  +
                          poly(PC2,2) +  (1|Code) + (1|Plot),
                        data = subset(all_data, inc_06_11_ht >= min.growth & dg.2006.cm > 5)))

# just size and competition
mod5 <- expression((lmer(inc_06_11_ht ~ 
                           scale(log(dg.2006.cm)) +
                           scale(BA.above2006) +
                           (1|Code) + (1|Plot),
                         data = subset(all_data, inc_06_11_ht >= min.growth & dg.2006.cm > 5))))

# let FTs differ in competition susceptibility
mod6 <- expression((lmer(inc_06_11_ht ~ 
                           scale(log(dg.2006.cm)) +
                           scale(BA.above2006) * C.Cerrado..G.generalist. +
                           (1|Code) + (1|Plot),
                         data = subset(all_data, inc_06_11_ht >= min.growth & dg.2006.cm > 5))))

# let FTs also differ in diameter-growth relationships
# doesn't help AIC
mod7 <- expression((lmer(inc_06_11_ht ~ 
                           scale(log(dg.2006.cm))* C.Cerrado..G.generalist. +
                           scale(BA.above2006) * C.Cerrado..G.generalist. +
                           (1|Code) + (1|Plot),
                         data = subset(all_data, inc_06_11_ht >= min.growth & dg.2006.cm > 5))))

# add in PC1
mod8 <- expression((lmer(inc_06_11_ht ~ 
                           scale(log(dg.2006.cm)) +
                           scale(BA.above2006) * C.Cerrado..G.generalist. +
                           poly(PC1, 2) + 
                           (1|Code) + (1|Plot),
                         data = subset(all_data, inc_06_11_ht >= min.growth & dg.2006.cm > 5))))
# add in linear PC2  
mod9 <-  expression((lmer(inc_06_11_ht ~ 
                            scale(log(dg.2006.cm)) +
                            scale(BA.above2006) * C.Cerrado..G.generalist. +
                            poly(PC1, 2) + PC2 +
                            (1|Code) + (1|Plot),
                          data = subset(all_data, inc_06_11_ht >= min.growth & dg.2006.cm > 5))))

#best model
# quadratic PC2
mod10 <- expression((lmer(inc_06_11_ht ~ 
                            scale(log(dg.2006.cm)) +
                            scale(BA.above2006) * C.Cerrado..G.generalist. +
                            poly(PC1, 2) + poly(PC2, 2) +
                            (1|Code) + (1|Plot),
                          data = subset(all_data, inc_06_11_ht >= min.growth & dg.2006.cm > 5))))

mod11 <- expression((lmer(inc_06_11_ht ~ 
                            scale(log(dg.2006.cm)) +
                            scale(TBA.2006.m2ha.1) * C.Cerrado..G.generalist. +
                            poly(PC1, 2) + poly(PC2, 2) +
                            (1|Code) + (1|Plot),
                          data = subset(all_data, inc_06_11_ht >= min.growth & dg.2006.cm > 5))))

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

write.csv(model_table, "./model outputs/model_table_06_11_height.csv")

#------------------------------------------------------------------------------
# redo for second time period

# just species
mod1 <- expression(lmer(inc_11_16_ht ~ (1|Code) + (1|Plot), data = subset(all_data, inc_11_16_ht >= min.growth & dg.2011.cm > 5))) # just species and diameter

mod2 <- expression(lmer(inc_11_16_ht ~ log(dg.2011.cm) + (1|Code) + (1|Plot), data = subset(all_data, inc_11_16_ht >= min.growth & dg.2011.cm > 5)))

# linear soils
mod3 <- expression(lmer(inc_11_16_ht ~ scale(PC1) + scale(PC2) + 
                          (1|Code) + (1|Plot),
                        data = subset(all_data, inc_11_16_ht >= min.growth & dg.2011.cm > 5)))

# quadratic soils
mod4 <- expression(lmer(inc_11_16_ht ~ poly(PC1,2)  +
                          poly(PC2,2) +  (1|Code) + (1|Plot),
                        data = subset(all_data, inc_11_16_ht >= min.growth & dg.2011.cm > 5)))

# just size and competition
mod5 <- expression((lmer(inc_11_16_ht ~ 
                           scale(log(dg.2011.cm)) +
                           scale(BA.above2011) +
                           (1|Code) + (1|Plot),
                         data = subset(all_data, inc_11_16_ht >= min.growth & dg.2011.cm > 5))))

# let FTs differ in competition susceptibility
mod6 <- expression((lmer(inc_11_16_ht ~ 
                           scale(log(dg.2011.cm)) +
                           scale(BA.above2011) * C.Cerrado..G.generalist. +
                           (1|Code) + (1|Plot),
                         data = subset(all_data, inc_11_16_ht >= min.growth & dg.2011.cm > 5))))

# let FTs also differ in diameter-growth relationships
# doesn't help AIC
mod7 <- expression((lmer(inc_11_16_ht ~ 
                           scale(log(dg.2011.cm))* C.Cerrado..G.generalist. +
                           scale(BA.above2011) * C.Cerrado..G.generalist. +
                           (1|Code) + (1|Plot),
                         data = subset(all_data, inc_11_16_ht >= min.growth & dg.2011.cm > 5))))

# add in PC1
mod8 <- expression((lmer(inc_11_16_ht ~ 
                           scale(log(dg.2011.cm)) +
                           scale(BA.above2011) * C.Cerrado..G.generalist. +
                           poly(PC1, 2) + 
                           (1|Code) + (1|Plot),
                         data = subset(all_data, inc_11_16_ht >= min.growth & dg.2011.cm > 5))))
# add in linear PC2  
mod9 <-  expression((lmer(inc_11_16_ht ~ 
                            scale(log(dg.2011.cm)) +
                            scale(BA.above2011) * C.Cerrado..G.generalist. +
                            poly(PC1, 2) + PC2 +
                            (1|Code) + (1|Plot),
                          data = subset(all_data, inc_11_16_ht >= min.growth & dg.2011.cm > 5))))

#best model
# quadratic PC2
mod10 <- expression((lmer(inc_11_16_ht ~ 
                            scale(log(dg.2011.cm)) +
                            scale(BA.above2011) * C.Cerrado..G.generalist. +
                            poly(PC1, 2) + poly(PC2, 2) +
                            (1|Code) + (1|Plot),
                          data = subset(all_data, inc_11_16_ht >= min.growth & dg.2011.cm > 5))))

mod11 <- expression((lmer(inc_11_16_ht ~ 
                            scale(log(dg.2011.cm)) +
                            scale(TBA.2011.m2.ha.1) * C.Cerrado..G.generalist. +
                            poly(PC1, 2) + poly(PC2, 2) +
                            (1|Code) + (1|Plot),
                          data = subset(all_data, inc_11_16_ht >= min.growth & dg.2011.cm > 5))))

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

write.csv(model_table, "./model outputs/model_table_11_16_height.csv")

#----------------------------------------------------------------------
# Full model: FT x CI interaction with soils
# = model 10 from above

min.growth <- -0.5
max.growth <- 1.2

options(na.action = na.omit)
data_06 <- subset(all_data, inc_06_11_ht >= min.growth & inc_06_11_ht <= max.growth & 
                    dg.2006.cm > 5)
data_11 <- subset(all_data, inc_11_16_ht >= min.growth & inc_11_16_ht <= max.growth & 
                    dg.2011.cm > 5)

ft_ci_soils_06_ht <- lmer(inc_06_11_ht ~ 
                            scale(H.2006) +
                            scale(BA.above2006) * C.Cerrado..G.generalist. +
                            poly(PC1, 2) + poly(PC2, 2) +
                            (1|Code) + (1|Plot),
                          data = data_06)

ft_ci_soils_11_ht <- lmer(inc_11_16_ht ~ 
                            scale(H.2011..m.) +
                            scale(BA.above2011) * C.Cerrado..G.generalist. +
                            poly(PC1, 2) + poly(PC2, 2) +
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
derivs1 <- ft_ci_soils_06@optinfo$derivs
sc_grad1 <- with(derivs1,solve(Hessian,gradient))
max(abs(sc_grad1))
max(pmin(abs(sc_grad1),abs(derivs1$gradient)))

relgrad <- with(ft_ci_soils_06@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))

qqmath(ranef(ft_ci_soils_06, condVar = TRUE), strip = FALSE)$Code

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

#-----
# bootstrap height models
# code from http://svmiller.com/blog/2020/03/bootstrap-standard-errors-in-r/
# 
# detach("package:lmerTest", unload = TRUE)
# 
# n_samples <- 5000
# 
# data_06 %>%
#   modelr::bootstrap(n_samples) %>%
#   mutate(lm = purrr::map(strap, ~lmer(inc_06_11_ht ~ 
#                                         scale(H.2006) +
#                                         scale(BA.above2006) * C.Cerrado..G.generalist. +
#                                         poly(PC1, 2) + poly(PC2, 2) +
#                                         (1|Code) + (1|Plot),
#                                       data = .)),
#          tidy = purrr::map(lm, broom.mixed::tidy)) -> bootCrime
# 
# bootCrime %>%
#   pull(tidy) %>%
#   purrr::map2_df(., # map to return a data frame
#                  seq(1, n_samples), # make sure to get this seq right. We did this 1000 times.
#                  ~mutate(.x, resample = .y)) -> tidybootCrime
# 
# tidybootCrime %>%
#   # group by term, naturally
#   dplyr::group_by(term) %>%
#   # This is the actual bootstrapped standard error you want
#   dplyr::summarize(bse = sd(estimate)) -> bseM3
# 
# 
# 
# data_11 %>%
#   modelr::bootstrap(n_samples) %>%
#   mutate(lm = purrr::map(strap, ~lmer(inc_11_16_ht~ 
#                                         scale(H.2011..m.) +
#                                         scale(BA.above2011) * C.Cerrado..G.generalist. +
#                                         poly(PC1, 2) + poly(PC2, 2) +
#                                         (1|Code) + (1|Plot),
#                                       data = .)),
#          tidy = purrr::map(lm, broom.mixed::tidy)) -> bootCrime
# 
# bootCrime %>%
#   pull(tidy) %>%
#   purrr::map2_df(., # map to return a data frame
#                  seq(1, n_samples), # make sure to get this seq right. We did this 1000 times.
#                  ~mutate(.x, resample = .y)) -> tidybootCrime
# 
# tidybootCrime %>%
#   # group by term, naturally
#   dplyr::group_by(term) %>%
#   # This is the actual bootstrapped standard error you want
#   dplyr::summarize(bse = sd(estimate)) -> bseM4

#------------------------------------------------------------------------------
# Mortality modeling
#------------------------------------------------------------------------------
library("pROC")


#------------------------------------------------------------------------------
# Model selection for mortality
#------------------------------------------------------------------------------
mod1 <- expression(glm(Died2011 ~ scale(log(dg.2006.cm)),
                            data = subset(all_data, dg.2006.cm > 5),
                            family = binomial(link = "logit")))

# linear soils
mod2 <- expression(glm(Died2011 ~ scale(PC1) + scale(PC2),
                       data = subset(all_data, dg.2006.cm > 5),
                       family = binomial(link = "logit")))

# quadratic soils
mod3 <- expression(glm(Died2011 ~ poly(PC1, 2) + poly(PC2, 2),
                       data = subset(all_data, dg.2006.cm > 5),
                       family = binomial(link = "logit")))

# just size and competition
mod4 <- expression(glm(Died2011 ~ scale(log(dg.2006.cm)) +
                         scale(BA.above2006),
                       data = subset(all_data, dg.2006.cm > 5),
                       family = binomial(link = "logit")))

# let FTs differ in competition susceptibility
mod5 <- expression(glm(Died2011 ~ scale(log(dg.2006.cm)) +
                         scale(BA.above2006) * C.Cerrado..G.generalist.,
                       data = subset(all_data, dg.2006.cm > 5),
                       family = binomial(link = "logit")))

# let FTs also differ in diameter-growth relationships
# doesn't help AIC
mod6 <- expression(glm(Died2011 ~ scale(log(dg.2006.cm))* C.Cerrado..G.generalist. +
                         scale(BA.above2006) * C.Cerrado..G.generalist.,
                       data = subset(all_data, dg.2006.cm > 5),
                       family = binomial(link = "logit")))

#best model
# add in PC1
mod7 <- expression(glm(Died2011 ~ scale(log(dg.2006.cm)) +
                         scale(BA.above2006) * C.Cerrado..G.generalist. +
                         poly(as.numeric(PC1), 2),
                       data = subset(all_data, dg.2006.cm > 5),
                       family = binomial(link = "logit")))

# add in linear PC2  
mod8 <-  expression(glm(Died2011 ~  scale(log(dg.2006.cm)) +
                          scale(BA.above2006) * C.Cerrado..G.generalist. +
                          poly(PC1, 2) + PC2,
                        data = subset(all_data, dg.2006.cm > 5),
                        family = binomial(link = "logit")))

# quadratic PC2
mod9 <- expression(glm(Died2011 ~  scale(log(dg.2006.cm)) +
                          scale(BA.above2006) * C.Cerrado..G.generalist. +
                          poly(PC1, 2) + poly(PC2, 2),
                        data = subset(all_data, dg.2006.cm > 5),
                        family = binomial(link = "logit")))

#only include PC1 because PC2 increases the AIC
mod10 <- expression(glm(Died2011 ~  scale(log(dg.2006.cm)) +
                          scale(TBA.2006.m2ha.1) * C.Cerrado..G.generalist. +
                          poly(PC1, 2),
                        data = subset(all_data, dg.2006.cm > 5),
                        family = binomial(link = "logit")))


model_list <- c(paste0("mod", seq(1:10)))

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

write.csv(model_table, "./model outputs/model_table_06_11_mortality.csv")

#------------------------------------------------------------------------------
# redo for second time period

mod1 <- expression(glm(Died2016 ~ scale(log(dg.2011.cm)),
                       data = subset(all_data, dg.2011.cm > 5),
                       family = binomial(link = "logit")))

# linear soils
mod2 <- expression(glm(Died2016 ~ scale(PC1) + scale(PC2),
                       data = subset(all_data, dg.2011.cm > 5),
                       family = binomial(link = "logit")))

# quadratic soils
mod3 <- expression(glm(Died2016 ~ poly(PC1, 2) + poly(PC2, 2),
                       data = subset(all_data, dg.2011.cm > 5),
                       family = binomial(link = "logit")))

# just size and competition
mod4 <- expression(glm(Died2016 ~ scale(log(dg.2011.cm)) +
                         scale(BA.above2011),
                       data = subset(all_data, dg.2011.cm > 5),
                       family = binomial(link = "logit")))

# let FTs differ in competition susceptibility
mod5 <- expression(glm(Died2016 ~ scale(log(dg.2011.cm)) +
                         scale(BA.above2011) * C.Cerrado..G.generalist.,
                       data = subset(all_data, dg.2011.cm > 5),
                       family = binomial(link = "logit")))

# let FTs also differ in diameter-growth relationships
# doesn't help AIC
mod6 <- expression(glm(Died2016 ~ scale(log(dg.2011.cm))* C.Cerrado..G.generalist. +
                         scale(BA.above2011) * C.Cerrado..G.generalist.,
                       data = subset(all_data, dg.2011.cm > 5),
                       family = binomial(link = "logit")))

# add in PC1
mod7 <- expression(glm(Died2016 ~ scale(log(dg.2011.cm)) +
                         scale(BA.above2011) * C.Cerrado..G.generalist. +
                         poly(as.numeric(PC1), 2),
                       data = subset(all_data, dg.2011.cm > 5),
                       family = binomial(link = "logit")))

# add in linear PC2  
mod8 <-  expression(glm(Died2016 ~  scale(log(dg.2011.cm)) +
                          scale(BA.above2011) * C.Cerrado..G.generalist. +
                          poly(PC1, 2) + PC2,
                        data = subset(all_data, dg.2011.cm > 5),
                        family = binomial(link = "logit")))
#best model
# quadratic PC2
mod9 <- expression(glm(Died2016 ~  scale(log(dg.2011.cm)) +
                         scale(BA.above2011) * C.Cerrado..G.generalist. +
                         poly(PC1, 2) + poly(PC2, 2),
                       data = subset(all_data, dg.2011.cm > 5),
                       family = binomial(link = "logit")))

#only include PC1 because PC2 increases the AIC
mod10 <- expression(glm(Died2016 ~  scale(log(dg.2011.cm)) +
                          scale(TBA.2011.m2.ha.1) * C.Cerrado..G.generalist. +
                          poly(PC1, 2),
                        data = subset(all_data, dg.2011.cm > 5),
                        family = binomial(link = "logit")))


model_list <- c(paste0("mod", seq(1:10)))

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

write.csv(model_table, "./model outputs/model_table_11_16_mortality.csv")

#------------------------------------------------------------------------------
#best model from selection below
mort_ft_ci_soils_glm_2011 <- glm(Died2011 ~ log(dg.2006.cm) +
                                   BA.above2006 * C.Cerrado..G.generalist. +
                                   poly(PC1, 2),
                                 data = subset(all_data, dg.2006.cm > 5),
                                 family = binomial(link = "logit"))

summary(mort_ft_ci_soils_glm_2011)
r.squaredGLMM(mort_ft_ci_soils_glm_2011)
AIC(mort_ft_ci_soils_glm_2011)
plot.roc(mort_ft_ci_soils_glm_2011$y, 
         fitted(mort_ft_ci_soils_glm_2011),print.auc = TRUE, 
         col = "green", lty = 2)
plot(Effect(mort_ft_ci_soils_glm_2011, 
            focal.predictors = c("BA.above2006", "C.Cerrado..G.generalist.")))

mort_ft_ci_soils_glm_2016 <- glm(Died2016 ~  log(dg.2011.cm) +
                                   BA.above2011 * C.Cerrado..G.generalist. +
                                   poly(PC1, 2) + poly(PC2, 2),
                                 data = subset(all_data, dg.2011.cm > 5),
                                 family = binomial(link = "logit"))

summary(mort_ft_ci_soils_glm_2016)
r.squaredGLMM(mort_ft_ci_soils_glm_2016)
AIC(mort_ft_ci_soils_glm_2016)
plot.roc(mort_ft_ci_soils_glm_2016$y, 
         fitted(mort_ft_ci_soils_glm_2016),print.auc = TRUE, 
         col = "green", lty = 2)
plot(Effect(mort_ft_ci_soils_glm_2016, 
            focal.predictors = c("BA.above2011", "C.Cerrado..G.generalist.")))

write.csv(summary(mort_ft_ci_soils_glm_2011)$coef, "./model outputs/model_coefs_06_mort.csv")
write.csv(summary(mort_ft_ci_soils_glm_2016)$coef, "./model outputs/model_coefs_11_mort.csv")


#hosmer-lemeshow test seems pretty okay!
hl1 <- hoslem.test(mort_ft_ci_soils_glm_2011$y, 
                   fitted(mort_ft_ci_soils_glm_2011), g=12)
hl2 <- hoslem.test(mort_ft_ci_soils_glm_2016$y, 
                   fitted(mort_ft_ci_soils_glm_2016), g=12)

###############################################################################
# Effects plots
#-----------------------------------------------------

#------------------------------------------------------------------------------
# Figure S1
# Overall plot-level results
#------------------------------------------------------------------------------

svg(filename="./plots/Figure S1 overall growth basal area.svg", 
     width = 3, 
     height=3, 
     pointsize=7)

par(mar = c(5.1, 5.1, 2, 1))

plot(plot_data$Diameter.increment.cm.yr ~ plot_data$TBA.2006.m2ha.1, 
     pch = 21,
     bg = "darkgray",
     xlab = "",
     ylab = "",
     axes = FALSE)

abline(coef(lm(plot_data$Diameter.increment.cm.yr ~ plot_data$TBA.2006.m2ha.1)))

mtext(side = 1, text = expression(paste("Stand basal area (m"^"2"," ha"^"-1",")")),
      line = 2.8, cex = 1.6)
mtext(side = 2, text = expression(paste("Mean diameter increment (cm year"^"-1", ")")), 
      line = 2.8, cex = 1.6)

axis(side =1, cex.axis = 1.5)
axis(side =2, cex.axis = 1.5)

dev.off()




#------------------------------------------------------------------------------
# Figure S3
# Growth and mortality ~ soil variables
#------------------------------------------------------------------------------

svg(filename="./plots/Figure S3 soil variables.svg", 
     width = 6, 
     height=9, 
     pointsize=8)

par(mar = c(2, 2, 2, 1),
    oma = c(6,3.5,1,0),
    mfcol = c(3,2))

opar <- par() #TODO add the arrows

eff <- Effect(focal.predictors = c("PC1"), 
              mod = ft_ci_soils_06, xlevels = 100,
              partial.residuals = TRUE,
              se = TRUE)

eff2 <- Effect(focal.predictors = c("PC1"), 
               mod = ft_ci_soils_11, xlevels = 100,
               partial.residuals = FALSE,
               se = TRUE)

dat <- data.frame(y = eff$fit,
                  lower = eff$lower,
                  upper = eff$upper,
                  soil = eff[["x"]][["PC1"]])

dat2 <- data.frame(y = eff2$fit,
                   lower = eff2$lower,
                   upper = eff2$upper,
                   soil = eff2[["x"]][["PC1"]])

plot(x = NA, 
     xlim = c(-1, 1.5),
     ylim = c(0, 0.35),
     xlab = "",
     ylab = "",
     xaxt = "n", 
     yaxt = "n")

lines(dat$y ~ dat$soil, col = "darkblue", lwd = 2, lty = 2)
polygon(c(dat$soil, rev(dat$soil)), c(dat$upper, rev(dat$lower)),
        col = addTrans("darkblue",30), border = NA)

lines(dat2$y ~ dat2$soil, col = "magenta", lwd = 2, lty = 2)
polygon(c(dat2$soil, rev(dat2$soil)), c(dat2$upper, rev(dat2$lower)),
        col = addTrans("magenta",30), border = NA)

mtext(side = 2, text = expression(paste("Diameter increment (cm year"^"-1", ")")), 
      line = 2.7, cex = 1.8)

axis(side =1, cex.axis = 1.5,at = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2),
     labels = c(NA, -1, NA, 0, NA, 1, NA, 2))
axis(side =2, cex.axis = 1.5)

legend("topright", legend = c("2006-2011", "2011-2016"), col = c("darkblue", "magenta"),
       lty = 1)

text(x = -1, y = 0.35, labels = "(a)", cex = 1.4)

#-------------------------------------------------
# Panel for height growth, PC1


eff <- Effect(focal.predictors = c("PC1"), 
              mod = ft_ci_soils_06_ht, xlevels = 100,
              partial.residuals = TRUE,
              se = TRUE)

eff2 <- Effect(focal.predictors = c("PC1"), 
               mod = ft_ci_soils_11_ht, xlevels = 100,
               partial.residuals = FALSE,
               se = TRUE)

dat <- data.frame(y = eff$fit,
                  lower = eff$lower,
                  upper = eff$upper,
                  soil = eff[["x"]][["PC1"]])

dat2 <- data.frame(y = eff2$fit,
                   lower = eff2$lower,
                   upper = eff2$upper,
                   soil = eff2[["x"]][["PC1"]])

plot(x = NA, 
     xlim = c(-1, 1.5),
     ylim = c(0, 0.35),
     xlab = "",
     ylab = "",
     xaxt = "n", 
     yaxt = "n")

lines(dat$y ~ dat$soil, col = "darkblue", lwd = 2)
polygon(c(dat$soil, rev(dat$soil)), c(dat$upper, rev(dat$lower)),
        col = addTrans("darkblue",30), border = NA)

lines(dat2$y ~ dat2$soil, col = "magenta", lwd = 2)
polygon(c(dat2$soil, rev(dat2$soil)), c(dat2$upper, rev(dat2$lower)),
        col = addTrans("magenta",30), border = NA)

mtext(side = 2, text = expression(paste("Height increment (m year"^"-1", ")")), 
      line = 2.7, cex = 1.8)

axis(side =1, cex.axis = 1.5,at = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2),
     labels = c(NA, -1, NA, 0, NA, 1, NA, 2))
axis(side =2, cex.axis = 1.5)

legend("topright", legend = c("2006-2011", "2011-2016"), col = c("darkblue", "magenta"),
       lty = 1)

text(x = -1, y = 0.35, labels = "(c)", cex = 1.4)

#----------------
# panel for mortality, PC1


eff <- Effect(focal.predictors = c("PC1"), 
              mod = mort_ft_ci_soils_glm_2011, xlevels = 100,
              partial.residuals = TRUE,
              se = TRUE)

eff2 <- Effect(focal.predictors = c("PC1"), 
               mod = mort_ft_ci_soils_glm_2016, xlevels = 100,
               partial.residuals = TRUE,
               se = TRUE)

dat <- data.frame(y = inv.logit(eff$fit)/5,
                  lower = inv.logit(eff$lower)/5,
                  upper = inv.logit(eff$upper)/5,
                  soil = eff[["x"]][["PC1"]])

dat2 <- data.frame(y = inv.logit(eff2$fit)/5,
                   lower = inv.logit(eff2$lower)/5,
                   upper = inv.logit(eff2$upper)/5,
                   soil = eff2[["x"]][["PC1"]])

plot(x = NA, 
     xlim = c(-1, 1.5),
     ylim = c(0, .05),
     xlab = "",
     ylab = "",
     xaxt = "n", 
     yaxt = "n")


lines(dat$y ~ dat$soil, col = "darkblue", lwd = 2)
polygon(c(dat$soil, rev(dat$soil)), c(dat$upper, rev(dat$lower)),
        col = addTrans("darkblue",30), border = NA)

lines(dat2$y ~ dat2$soil, col = "magenta", lwd = 2)
polygon(c(dat2$soil, rev(dat2$soil)), c(dat2$upper, rev(dat2$lower)),
        col = addTrans("magenta",30), border = NA)

mtext(side = 1, text = "Soil PC1", line = 2.7, cex = 1.8)
mtext(side = 2, text = expression(paste("P(Mortality)")), 
      line = 2.7, cex = 1.8)

axis(side =1, cex.axis = 1.5,at = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2),
     labels = c(NA, -1, NA, 0, NA, 1, NA, 2))
axis(side =2, cex.axis = 1.5)

legend("topright", legend = c("2006-2011", "2011-2016"), col = c("darkblue", "magenta"),
       lty = 1)

text(x = -1, y = 0.05, labels = "(e)", cex = 1.4)

mtext(side = 1, line = 5, text = "pH                                      C.E.C.", cex = 1.4)

# arrows(x0 = -0.8, x1 = 1.2, y0 = -0.008, y1 = -0.008, xpd = TRUE, code = 3, length = 0.125)

#------------------------------
# panel for growth, soil pc2

eff <- Effect(focal.predictors = c("PC2"), 
              mod = ft_ci_soils_06, xlevels = 100, transformation = list(link = log, inverse = exp),
              partial.residuals = TRUE,
              se = TRUE)

eff2 <- Effect(focal.predictors = c("PC2"), 
               mod = ft_ci_soils_11, xlevels = 100, transformation = list(link = log, inverse = exp),
               partial.residuals = TRUE,
               se = TRUE)

dat <- data.frame(y = eff$fit,
                  lower = eff$lower,
                  upper = eff$upper,
                  soil = eff[["x"]][["PC2"]])

dat2 <- data.frame(y = eff2$fit,
                   lower = eff2$lower,
                   upper = eff2$upper,
                   soil = eff2[["x"]][["PC2"]])

plot(x = NA, 
     xlim = c(-1.5,2),
     ylim = c(0, 0.35),
     xlab = "",
     ylab = "",
     xaxt = "n", 
     yaxt = "n")

# points(e1_r_agg$exp_mean_resids ~ e1_r_agg$x.fit, 
#        pch = 21, 
#        bg= "darkblue",
#        col = "darkblue")
# arrows(y0 = e1_r_agg$exp_low, x0 = e1_r_agg$x.fit, 
#        y1 = e1_r_agg$exp_high, x1 = e1_r_agg$x.fit,
#        code = 3,
#        angle = 90,
#        length = 0.02,
#        lwd = 0.5)

lines(dat$y ~ dat$soil, col = "darkblue", lwd = 2)
polygon(c(dat$soil, rev(dat$soil)), c(dat$upper, rev(dat$lower)),
        col = addTrans("darkblue",30), border = NA)

lines(dat2$y ~ dat2$soil, col = "magenta", lwd = 2, lty = 2)
polygon(c(dat2$soil, rev(dat2$soil)), c(dat2$upper, rev(dat2$lower)),
        col = addTrans("magenta",30), border = NA)

axis(side =1, cex.axis = 1.5, at = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2),
     labels = c(NA, -1, NA, 0, NA, 1, NA, 2))
axis(side =2, cex.axis = 1.5)

legend("topright", legend = c("2006-2011", "2011-2016"), col = c("darkblue", "magenta"),
       lty = 1)

text(x = -1.4, y = 0.35, labels = "(b)", cex = 1.4)
#-------------------------------------------------
# Panel for height growth, PC1

eff <- Effect(focal.predictors = c("PC2"), 
              mod = ft_ci_soils_06_ht, xlevels = 100,
              partial.residuals = TRUE,
              se = TRUE)

eff2 <- Effect(focal.predictors = c("PC2"), 
               mod = ft_ci_soils_11_ht, xlevels = 100,
               partial.residuals = FALSE,
               se = TRUE)

dat <- data.frame(y = eff$fit,
                  lower = eff$lower,
                  upper = eff$upper,
                  soil = eff[["x"]][["PC2"]])

dat2 <- data.frame(y = eff2$fit,
                   lower = eff2$lower,
                   upper = eff2$upper,
                   soil = eff2[["x"]][["PC2"]])

plot(x = NA, 
     xlim = c(-1.5, 2),
     ylim = c(0, 0.35),
     xlab = "",
     ylab = "",
     xaxt = "n", 
     yaxt = "n")

lines(dat$y ~ dat$soil, col = "darkblue", lwd = 2)
polygon(c(dat$soil, rev(dat$soil)), c(dat$upper, rev(dat$lower)),
        col = addTrans("darkblue",30), border = NA)

lines(dat2$y ~ dat2$soil, col = "magenta", lwd = 2)
polygon(c(dat2$soil, rev(dat2$soil)), c(dat2$upper, rev(dat2$lower)),
        col = addTrans("magenta",30), border = NA)

axis(side =1, cex.axis = 1.5,at = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2),
     labels = c(NA, -1, NA, 0, NA, 1, NA, 2))
axis(side =2, cex.axis = 1.5)

legend("topright", legend = c("2006-2011", "2011-2016"), col = c("darkblue", "magenta"),
       lty = 1)

text(x = -1, y = 0.35, labels = "(d)", cex = 1.4)

#----------------
# panel for mortality, PC1

eff2 <- Effect(focal.predictors = c("PC2"), 
               mod = mort_ft_ci_soils_glm_2016, xlevels = 100,
               partial.residuals = TRUE,
               se = TRUE)

dat2 <- data.frame(y = inv.logit(eff2$fit)/5,
                   lower = inv.logit(eff2$lower)/5,
                   upper = inv.logit(eff2$upper)/5,
                   soil = eff2[["x"]][["PC2"]])

plot(x = NA, 
     xlim = c(-1.5, 2),
     ylim = c(0, .05),
     xlab = "",
     ylab = "",
     xaxt = "n", 
     yaxt = "n")

lines(dat2$y ~ dat2$soil, col = "magenta", lwd = 2, lty = 1)
polygon(c(dat2$soil, rev(dat2$soil)), c(dat2$upper, rev(dat2$lower)),
        col = addTrans("magenta",30), border = NA)

mtext(side = 1, text = "Soil PC2", line = 2.7, cex = 1.8)

axis(side =1, cex.axis = 1.5,at = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2),
     labels = c(NA, -1, NA, 0, NA, 1, NA, 2))
axis(side =2, cex.axis = 1.5)

legend("topright", legend = c("2011-2016"), col = c("magenta"),
       lty = 1)

text(x = -1.4, y = 0.05, labels = "(f)", cex = 1.4)



mtext(side = 1, line = 5, text = "Fewer Bases                 More Bases", cex = 1.4)

# arrows(x0 = -1.5, x1 = 2, y0 = -0.015, y1 = -0.02, xpd = TRUE, code = 3, length = 0.125)

dev.off()


#------------------------------------------------------------------------------
# figure S4 soil effects for the appendix
#------------------------------------------------------------------------------

svg(filename="./plots/Figure S4 soil variables no covariates.svg", 
    width = 6, 
    height=9, 
    pointsize=8)

par(mar = c(2, 2, 2, 1),
    oma = c(6,3.5,1,0),
    mfcol = c(3,2))

opar <- par()

model1 <- lmer(inc_06_11 ~ poly(PC1,2)  +
                poly(PC2,2) +  (1|Code) + (1|Plot),
              data = data_06)

eff <- Effect(focal.predictors = c("PC1"), 
              mod = model1, xlevels = 100,
              partial.residuals = TRUE,
              se = TRUE)

model2 <- lmer(inc_11_16 ~ poly(PC1,2)  +
                 poly(PC2,2) +  (1|Code) + (1|Plot),
               data = subset(all_data, inc_11_16 >= min.growth & dg.2006.cm > 5))

eff2 <- Effect(focal.predictors = c("PC1"), 
               mod = model2, xlevels = 100,
               partial.residuals = TRUE,
               se = TRUE)

dat <- data.frame(y = eff$fit,
                  lower = eff$lower,
                  upper = eff$upper,
                  soil = eff[["x"]][["PC1"]])

dat2 <- data.frame(y = eff2$fit,
                   lower = eff2$lower,
                   upper = eff2$upper,
                   soil = eff2[["x"]][["PC1"]])

plot(x = NA, 
     xlim = c(-1, 1.5),
     ylim = c(0, 0.5),
     xlab = "",
     ylab = "",
     xaxt = "n", 
     yaxt = "n")

lines(dat$y ~ dat$soil, col = "darkblue", lwd = 2)
polygon(c(dat$soil, rev(dat$soil)), c(dat$upper, rev(dat$lower)),
        col = addTrans("darkblue",30), border = NA)

lines(dat2$y ~ dat2$soil, col = "magenta", lwd = 2)
polygon(c(dat2$soil, rev(dat2$soil)), c(dat2$upper, rev(dat2$lower)),
        col = addTrans("magenta",30), border = NA)

mtext(side = 2, text = expression(paste("Diameter increment (cm year"^"-1", ")")), 
      line = 2.7, cex = 1.8)

axis(side =1, cex.axis = 1.5,at = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2),
     labels = c(NA, -1, NA, 0, NA, 1, NA, 2))
axis(side =2, cex.axis = 1.5)

legend("topright", legend = c("2006-2011", "2011-2016"), col = c("darkblue", "magenta"),
       lty = 1)

text(x = -1, y = 0.5, labels = "(a)", cex = 1.4)

#------------------------------------------------------------------------------
# panel for height, PC1

model1 <- lmer(inc_06_11_ht ~ poly(PC1,2)  +
                 poly(PC2,2) +  (1|Code) + (1|Plot),
               data = subset(all_data, inc_06_11_ht >= min.growth & dg.2006.cm > 5))

eff <- Effect(focal.predictors = c("PC1"), 
              mod = model1, xlevels = 100,
              partial.residuals = TRUE,
              se = TRUE)

model2 <- lmer(inc_11_16_ht ~ poly(PC1,2)  +
                 poly(PC2,2) +  (1|Code) + (1|Plot),
               data = subset(all_data, inc_11_16_ht >= min.growth & dg.2006.cm > 5))

eff2 <- Effect(focal.predictors = c("PC1"), 
               mod = model2, xlevels = 100,
               partial.residuals = TRUE,
               se = TRUE)

dat <- data.frame(y = eff$fit,
                  lower = eff$lower,
                  upper = eff$upper,
                  soil = eff[["x"]][["PC1"]])

dat2 <- data.frame(y = eff2$fit,
                   lower = eff2$lower,
                   upper = eff2$upper,
                   soil = eff2[["x"]][["PC1"]])

plot(x = NA, 
     xlim = c(-1, 1.5),
     ylim = c(0, 0.5),
     xlab = "",
     ylab = "",
     xaxt = "n", 
     yaxt = "n")

lines(dat$y ~ dat$soil, col = "darkblue", lwd = 2, lty = 2)
polygon(c(dat$soil, rev(dat$soil)), c(dat$upper, rev(dat$lower)),
        col = addTrans("darkblue",30), border = NA)

lines(dat2$y ~ dat2$soil, col = "magenta", lwd = 2, lty = 2)
polygon(c(dat2$soil, rev(dat2$soil)), c(dat2$upper, rev(dat2$lower)),
        col = addTrans("magenta",30), border = NA)

mtext(side = 2, text = expression(paste("Height increment (m year"^"-1", ")")), 
      line = 2.7, cex = 1.8)

axis(side =1, cex.axis = 1.5,at = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2),
     labels = c(NA, -1, NA, 0, NA, 1, NA, 2))
axis(side =2, cex.axis = 1.5)

legend("topright", legend = c("2006-2011", "2011-2016"), col = c("darkblue", "magenta"),
       lty = 1)

text(x = -1, y = 0.5, labels = "(c)", cex = 1.4)

#----------------
# panel for mortality, PC1

model_mort1 <- glm(formula = Died2011 ~ poly(PC1, 2) + poly(PC2, 2), 
                   family = binomial(link = "logit"), 
                   data = subset(all_data, dg.2006.cm > 5))

eff <- Effect(focal.predictors = c("PC1"), 
              mod = model_mort1, xlevels = 100,
              partial.residuals = TRUE,
              se = TRUE)

model_mort2 <- glm(formula = Died2016 ~ poly(PC1, 2) + poly(PC2, 2), 
                   family = binomial(link = "logit"), 
                   data = subset(all_data, dg.2011.cm > 5))

eff2 <- Effect(focal.predictors = c("PC1"), 
               mod = model_mort2, xlevels = 100,
               partial.residuals = TRUE,
               se = TRUE)

dat <- data.frame(y = inv.logit(eff$fit)/5,
                  lower = inv.logit(eff$lower)/5,
                  upper = inv.logit(eff$upper)/5,
                  soil = eff[["x"]][["PC1"]])

dat2 <- data.frame(y = inv.logit(eff2$fit)/5,
                   lower = inv.logit(eff2$lower)/5,
                   upper = inv.logit(eff2$upper)/5,
                   soil = eff2[["x"]][["PC1"]])

plot(x = NA, 
     xlim = c(-1, 1.5),
     ylim = c(0, .05),
     xlab = "",
     ylab = "",
     xaxt = "n", 
     yaxt = "n")

# points(e1_r_agg$mean_part_resids ~ e1_r_agg$x.fit, 
#        pch = 21, 
#        bg= "darkblue",
#        col = "darkblue")
# arrows(y0 = e1_r_agg$low, x0 = e1_r_agg$x.fit, 
#        y1 = e1_r_agg$high, x1 = e1_r_agg$x.fit,
#        code = 3,
#        angle = 90,
#        length = 0.02,
#        lwd = 0.5)

lines(dat$y ~ dat$soil, col = "darkblue", lwd = 2)
polygon(c(dat$soil, rev(dat$soil)), c(dat$upper, rev(dat$lower)),
        col = addTrans("darkblue",30), border = NA)

lines(dat2$y ~ dat2$soil, col = "magenta", lwd = 2)
polygon(c(dat2$soil, rev(dat2$soil)), c(dat2$upper, rev(dat2$lower)),
        col = addTrans("magenta",30), border = NA)

mtext(side = 1, text = "Soil PC1", line = 2.7, cex = 1.8)
mtext(side = 2, text = expression(paste("P(Mortality)")), 
      line = 2.7, cex = 1.8)

axis(side =1, cex.axis = 1.5,at = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2),
     labels = c(NA, -1, NA, 0, NA, 1, NA, 2))
axis(side =2, cex.axis = 1.5)

legend("topright", legend = c("2006-2011", "2011-2016"), col = c("darkblue", "magenta"),
       lty = 1)

text(x = -1, y = 0.05, labels = "(e)", cex = 1.4)

mtext(side = 1, line = 5, text = "pH                                      C.E.C.", cex = 1.4)

# arrows(x0 = -0.8, x1 = 1.2, y0 = -0.008, y1 = -0.008, xpd = TRUE, code = 3, length = 0.125)

#------------------------------
# panel for diameter growth, soil pc2

eff <- Effect(focal.predictors = c("PC2"), 
              mod = model1, xlevels = 100, transformation = list(link = log, inverse = exp),
              partial.residuals = TRUE,
              se = TRUE)

eff2 <- Effect(focal.predictors = c("PC2"), 
               mod = model2, xlevels = 100, transformation = list(link = log, inverse = exp),
               partial.residuals = TRUE,
               se = TRUE)

dat <- data.frame(y = eff$fit,
                  lower = eff$lower,
                  upper = eff$upper,
                  soil = eff[["x"]][["PC2"]])

dat2 <- data.frame(y = eff2$fit,
                   lower = eff2$lower,
                   upper = eff2$upper,
                   soil = eff2[["x"]][["PC2"]])

plot(x = NA, 
     xlim = c(-1.5,2),
     ylim = c(0, 0.5),
     xlab = "",
     ylab = "",
     xaxt = "n", 
     yaxt = "n")

# points(e1_r_agg$exp_mean_resids ~ e1_r_agg$x.fit, 
#        pch = 21, 
#        bg= "darkblue",
#        col = "darkblue")
# arrows(y0 = e1_r_agg$exp_low, x0 = e1_r_agg$x.fit, 
#        y1 = e1_r_agg$exp_high, x1 = e1_r_agg$x.fit,
#        code = 3,
#        angle = 90,
#        length = 0.02,
#        lwd = 0.5)

lines(dat$y ~ dat$soil, col = "darkblue", lwd = 2)
polygon(c(dat$soil, rev(dat$soil)), c(dat$upper, rev(dat$lower)),
        col = addTrans("darkblue",30), border = NA)

lines(dat2$y ~ dat2$soil, col = "magenta", lwd = 2)
polygon(c(dat2$soil, rev(dat2$soil)), c(dat2$upper, rev(dat2$lower)),
        col = addTrans("magenta",30), border = NA)

axis(side =1, cex.axis = 1.5, at = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2),
     labels = c(NA, -1, NA, 0, NA, 1, NA, 2))
axis(side =2, cex.axis = 1.5)

legend("topright", legend = c("2006-2011", "2011-2016"), col = c("darkblue", "magenta"),
       lty = 1)

text(x = -1.4, y = 0.5, labels = "(b)", cex = 1.4)

#------------------------------------------------------------------------------
# panel for height, PC2

model1 <- lmer(inc_06_11_ht ~ poly(PC1,2)  +
                 poly(PC2,2) +  (1|Code) + (1|Plot),
               data = data_06)

eff <- Effect(focal.predictors = c("PC2"), 
              mod = model1, xlevels = 100,
              partial.residuals = TRUE,
              se = TRUE)

model2 <- lmer(inc_11_16_ht ~ poly(PC1,2)  +
                 poly(PC2,2) +  (1|Code) + (1|Plot),
               data = data_11)

eff2 <- Effect(focal.predictors = c("PC2"), 
               mod = model2, xlevels = 100,
               partial.residuals = TRUE,
               se = TRUE)

dat <- data.frame(y = eff$fit,
                  lower = eff$lower,
                  upper = eff$upper,
                  soil = eff[["x"]][["PC2"]])

dat2 <- data.frame(y = eff2$fit,
                   lower = eff2$lower,
                   upper = eff2$upper,
                   soil = eff2[["x"]][["PC2"]])

plot(x = NA, 
     xlim = c(-1.5, 2),
     ylim = c(0, 0.5),
     xlab = "",
     ylab = "",
     xaxt = "n", 
     yaxt = "n")

lines(dat$y ~ dat$soil, col = "darkblue", lwd = 2, lty = 2)
polygon(c(dat$soil, rev(dat$soil)), c(dat$upper, rev(dat$lower)),
        col = addTrans("darkblue",30), border = NA)

lines(dat2$y ~ dat2$soil, col = "magenta", lwd = 2, lty = 2)
polygon(c(dat2$soil, rev(dat2$soil)), c(dat2$upper, rev(dat2$lower)),
        col = addTrans("magenta",30), border = NA)

axis(side =1, cex.axis = 1.5,at = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2),
     labels = c(NA, -1, NA, 0, NA, 1, NA, 2))
axis(side =2, cex.axis = 1.5)

legend("topright", legend = c("2006-2011", "2011-2016"), col = c("darkblue", "magenta"),
       lty = 1)

text(x = -1.4, y = 0.5, labels = "(d)", cex = 1.4)


#----------------
# panel for mortality, PC2

eff <- Effect(focal.predictors = c("PC2"), 
              mod = model_mort1, xlevels = 100,
              partial.residuals = TRUE,
              se = TRUE)

eff2 <- Effect(focal.predictors = c("PC2"), 
               mod = model_mort2, xlevels = 100,
               partial.residuals = TRUE,
               se = TRUE)

dat <- data.frame(y = inv.logit(eff$fit)/5,
                  lower = inv.logit(eff$lower)/5,
                  upper = inv.logit(eff$upper)/5,
                  soil = eff[["x"]][["PC2"]])

dat2 <- data.frame(y = inv.logit(eff2$fit)/5,
                   lower = inv.logit(eff2$lower)/5,
                   upper = inv.logit(eff2$upper)/5,
                   soil = eff2[["x"]][["PC2"]])

plot(x = NA, 
     xlim = c(-1.5, 2),
     ylim = c(0, .05),
     xlab = "",
     ylab = "",
     xaxt = "n", 
     yaxt = "n")

lines(dat$y ~ dat$soil, col = "darkblue", lwd = 2)
polygon(c(dat$soil, rev(dat$soil)), c(dat$upper, rev(dat$lower)),
        col = addTrans("darkblue",30), border = NA)


lines(dat2$y ~ dat2$soil, col = "magenta", lwd = 2)
polygon(c(dat2$soil, rev(dat2$soil)), c(dat2$upper, rev(dat2$lower)),
        col = addTrans("magenta",30), border = NA)

mtext(side = 1, text = "Soil PC2", line = 2.7, cex = 1.8)

axis(side =1, cex.axis = 1.5,at = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2),
     labels = c(NA, -1, NA, 0, NA, 1, NA, 2))
axis(side =2, cex.axis = 1.5)

legend("topright", legend = c("2006-2011", "2011-2016"), col = c("darkblue", "magenta"),
       lty = 1)

text(x = -1.4, y = 0.05, labels = "(f)", cex = 1.4)



mtext(side = 1, line = 5, text = "Fewer Bases                 More Bases", cex = 1.4)

# arrows(x0 = -1.5, x1 = 2, y0 = -0.015, y1 = -0.02, xpd = TRUE, code = 3, length = 0.125)

dev.off()


#------------------------------------------------------------------------------
# Figure 4
# Mortality ~ BA by FT
#------------------------------------------------------------------------------
eff <- Effect(focal.predictors = c("BA.above2006", "C.Cerrado..G.generalist."), 
              mod = mort_ft_ci_soils_glm_2011, xlevels = 100)
eff2 <- Effect(focal.predictors = c("BA.above2011", "C.Cerrado..G.generalist."),
              mod = mort_ft_ci_soils_glm_2016, xlevels = 100)


tiff(filename="./plots/Figure 4 mortality_ba_above_by_fg.tiff", 
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

#------------------------------------------------------------------------------
# Figure 2
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
change_live <- data.frame(Plot = unique(tree_data$Plot),
                          s_2016 = numeric(length(unique(tree_data$Plot))),
                          s_2006 = numeric(length(unique(tree_data$Plot))),
                          g_2016 = numeric(length(unique(tree_data$Plot))),
                          g_2006 = numeric(length(unique(tree_data$Plot)))
                          )


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
  
  change_live[change_live$Plot == plot, ]$s_2016 <- sum((((temp[temp$C.Cerrado..G.generalist. == "C" & temp$Died2016 != "y" & temp$Died2011 != "y", 
                                     "dg.2016.cm"]/2)^2)*pi), na.rm = TRUE)/10000
  change_live[change_live$Plot == plot, ]$g_2016 <- sum((((temp[temp$C.Cerrado..G.generalist. == "G" & temp$Died2016 != "y" & temp$Died2011 != "y", 
                                     "dg.2016.cm"]/2)^2)*pi), na.rm = TRUE)/10000
  change_live[change_live$Plot == plot, ]$s_2006 <- sum((((temp[temp$C.Cerrado..G.generalist. == "C" & temp$Died2016 != "y" & temp$Died2011 != "y", 
                                          "dg.2006.cm"]/2)^2)*pi), na.rm = TRUE)/10000
  change_live[change_live$Plot == plot, ]$g_2006 <- sum((((temp[temp$C.Cerrado..G.generalist. == "G" & temp$Died2016 != "y" & temp$Died2011 != "y", 
                                         "dg.2006.cm"]/2)^2)*pi), na.rm = TRUE)/10000
  
}

change_ba <- data.frame(Plot = change_2006$Plot,
                        ba_tot_2006 = change_2006$BA_total,
                        ba_inc_s = change_2016$BA_s - change_2006$BA_s,
                        ba_inc_g = change_2016$BA_g - change_2006$BA_g,
                        ba_died_s = -(change_2011$BA_s_died + change_2016$BA_s_died),
                        ba_died_g = -(change_2011$BA_g_died + change_2016$BA_g_died),
                        ba_inc_s_w_dead = change_2016$BA_s - change_2006$BA_s + change_2011$BA_s_died + change_2016$BA_s_died,
                        ba_inc_g_w_dead = change_2016$BA_g - change_2006$BA_g + change_2011$BA_g_died + change_2016$BA_g_died,
                        ba_live_s = change_live$s_2016 - change_live$s_2006,
                        ba_live_g = change_live$g_2016 - change_live$g_2006
                        )

change_ba_prop <- change_ba
change_ba_prop[, c(3, 5, 7, 9)] <- change_ba_prop[, c(3,5,7,9)]/change_2006$BA_s
change_ba_prop[, c(4, 6, 8, 10)] <- change_ba_prop[, c(4,6,8,10)]/change_2006$BA_g

#make three-panel figure

svg(filename="./plots/Figure 2 change in basal area proportional.svg",
     width = 3.6, 
     height=8, 
     pointsize=7) 


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

text(x = 0.4, y = 1.67, labels = "(a)", cex = 1.2)

#just losses due to mortality

plot(NA,
     ylim = c(-.8, 0),
     xlim = c(0.4, 2.6),
     xaxt = "n",
     yaxt = "n",
     xlab = "",
     ylab = expression(paste("Basal area lost\nto mortality (%)")),
     cex.lab = 1.2)

abline(h = 0)

points(ba_died_s ~ ba_tot_2006, data = change_ba_prop, 
       pch = 21,
       col = "#d95f02",
       bg = "#d95f02")

mod <- lm(log(-ba_died_s + 1) ~ log(ba_tot_2006), data = change_ba_prop)
pred <- predict(mod, newdata = new)
lines(-(exp(pred) - 1) ~ new$ba_tot_2006,
      col = "#d95f02",
      lty = 1)


points(ba_died_g ~ ba_tot_2006, data = change_ba_prop, 
       pch = 22,
       col = "#1b9e77",
       bg = "#1b9e77")
mod <- lm(log(-ba_died_g + 1) ~ log(ba_tot_2006), data = change_ba_prop)
pred <- predict(mod, newdata = new)
lines(-(exp(pred) - 1) ~ new$ba_tot_2006,
      col = "#1b9e77",
      lty = 1)

axis(side = 2, at = pretty(c(-.8, 0)), labels = pretty(c(-.8, 0))*100)
axis(side = 1, at = pretty(change_ba$ba_tot_2006), labels = NA)


text(x = 0.4, y = -.05, labels = "(b)", cex = 1.2)

#without dead trees
plot(NA,
     ylim = c(-0.7, 1.8),
     xlim = c(0.4, 2.6),
     xaxt = "n",
     yaxt = "n",
     xlab = "",
     ylab = expression(paste("Change in basal area \n without mortality (%)")),
     cex.lab = 1.2)

abline(h = 0)

points(ba_live_s ~ ba_tot_2006, data = change_ba_prop,
       pch = 21,
       col = "#d95f02",
       bg = "#d95f02")
mod <- lm(log(ba_live_s + 1) ~ log(ba_tot_2006), data = change_ba_prop)
pred <- predict(mod, newdata = new)
lines(exp(pred) - 1 ~ new$ba_tot_2006,
      col = "#d95f02",
      lty = 1)

points(ba_live_g ~ ba_tot_2006, data = change_ba_prop,
       pch = 22,
       col = "#1b9e77",
       bg = "#1b9e77")
mod <- lm(log(ba_live_g + 1) ~ log(ba_tot_2006), data = change_ba_prop)
pred <- predict(mod, newdata = new)
lines(exp(pred) - 1 ~ new$ba_tot_2006,
      col = "#1b9e77",
      lty = 1)

axis(side = 1, at = pretty(change_ba$ba_tot_2006), labels = pretty(change_ba$ba_tot_2006)*10)
axis(side = 2, at = pretty(c(-0.7, 1.8)), labels = pretty(c(-0.7, 1.8))*100)

mtext(side = 1, line = 2.9, text = expression(paste("Original basal area (m"^2, " ha"^-1, ")")),
      cex = 1.3)


text(x = 0.4, y = 1.67, labels = "(c)", cex = 1.2)

dev.off()

#------------------------------------------------------------------------------
# Figure 3
# Growth ~ diam by FT
#------------------------------------------------------------------------------
tiff(filename="./plots/Figure 3 differences in growth by fg.tiff",
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
#-------------------------------------------------------------------------------
# diameter growth by BA above

get_predictions_ba <- function(model, ba_name){
  eff <- Effect(focal.predictors = c(ba_name, "C.Cerrado..G.generalist."),
                mod = model, xlevels = 100,
                partial.residuals = TRUE)
  
  y <- eff$fit
  x <- eff[["x"]][[ba_name]]
  x.fit <- eff[["x.all"]][[ba_name]]
  fg <- eff[["x.all"]]$"C.Cerrado..G.generalist."
  
  data_orig <- data.frame(x.fit,
                          fg)
  
  #make sure we match the correct line
  fitted <- ifelse(data_orig$fg == "C", y[1:100][closest(x.fit, x[1:100])],
                   y[101:200][closest(x.fit, x[101:200])])
  
  # resids <- eff$residuals + fitted
  # exp_resids <- exp(resids)
  
  dat <- data.frame(y = eff$fit,
                    lower = eff$lower,
                    upper = eff$upper,
                    ba = eff[["x"]][[ba_name]] * 10,
                    fg = eff[["x"]][["C.Cerrado..G.generalist."]])
  return(dat)
}

dat1 <- get_predictions_ba(ft_ci_soils_06, "BA.above2006")
dat2 <- get_predictions_ba(ft_ci_soils_11, "BA.above2011")

plot(x = NA,
     xlim = c(0, 28),
     ylim = c(0, 0.6),
     xlab = "",
     ylab = "",
     xaxt = "n",
     yaxt = "n")

dat_sub <- dat1[dat1$fg == "G", ]
lines(dat_sub$y ~ dat_sub$ba, col = "#1b9e77", lwd = 2)
polygon(c(dat_sub$ba, rev(dat_sub$ba)), c(dat_sub$upper, rev(dat_sub$lower)),
        col = addTrans("#1b9e77",20), border = NA)

dat_sub <- dat1[dat1$fg == "C", ]
lines(dat_sub$y ~ dat_sub$ba, col = "#d95f02", lwd = 2)
polygon(c(dat_sub$ba, rev(dat_sub$ba)), c(dat_sub$upper, rev(dat_sub$lower)),
        col = addTrans("#d95f02",20), border = NA)

dat_sub <- dat2[dat2$fg == "G", ]
lines(dat_sub$y ~ dat_sub$ba, col = "#1b9e77", lwd = 2, lty = 2)
polygon(c(dat_sub$ba, rev(dat_sub$ba)), c(dat_sub$upper, rev(dat_sub$lower)),
        col = addTrans("#1b9e77",20), border = NA)

dat_sub <- dat2[dat2$fg == "C", ]
lines(dat_sub$y ~ dat_sub$ba, col = "#d95f02", lwd = 2, lty = 2)
polygon(c(dat_sub$ba, rev(dat_sub$ba)), c(dat_sub$upper, rev(dat_sub$lower)),
        col = addTrans("#d95f02",20), border = NA)

mtext(side = 1, text = expression(paste("BA above (m"^"2"," ha"^"-1",")")), line = 3.1, cex = 1.8)
mtext(side = 2, text = expression(paste("Diameter increment (cm year"^"-1", ")")),
      line = 2.7, cex = 1.8)

axis(side =1, cex.axis = 1.5)
axis(side =2, cex.axis = 1.5)

text(x = 1.7, y = 0.58, labels = "(a)", cex = 1.5)

#-------------------------------------------------------------------------------
# do the same for height
#-------------------------------------------------------------------------------
dat1 <- get_predictions_ba(ft_ci_soils_06_ht, "BA.above2006")
dat2 <- get_predictions_ba(ft_ci_soils_11_ht, "BA.above2011")

plot(x = NA,
     xlim = c(0, 28),
     ylim = c(0, 0.4),
     xlab = "",
     ylab = "",
     xaxt = "n",
     yaxt = "n")

dat_sub <- dat1[dat1$fg == "G", ]
lines(dat_sub$y ~ dat_sub$ba, col = "#1b9e77", lwd = 2)
polygon(c(dat_sub$ba, rev(dat_sub$ba)), c(dat_sub$upper, rev(dat_sub$lower)),
        col = addTrans("#1b9e77",20), border = NA)

dat_sub <- dat1[dat1$fg == "C", ]
lines(dat_sub$y ~ dat_sub$ba, col = "#d95f02", lwd = 2)
polygon(c(dat_sub$ba, rev(dat_sub$ba)), c(dat_sub$upper, rev(dat_sub$lower)),
        col = addTrans("#d95f02",20), border = NA)

dat_sub <- dat2[dat2$fg == "G", ]
lines(dat_sub$y ~ dat_sub$ba, col = "#1b9e77", lwd = 2, lty = 2)
polygon(c(dat_sub$ba, rev(dat_sub$ba)), c(dat_sub$upper, rev(dat_sub$lower)),
        col = addTrans("#1b9e77",20), border = NA)

dat_sub <- dat2[dat2$fg == "C", ]
lines(dat_sub$y ~ dat_sub$ba, col = "#d95f02", lwd = 2, lty = 2)
polygon(c(dat_sub$ba, rev(dat_sub$ba)), c(dat_sub$upper, rev(dat_sub$lower)),
        col = addTrans("#d95f02",20), border = NA)

mtext(side = 1, text = expression(paste("BA above (m"^"2"," ha"^"-1",")")), line = 3.1, cex = 1.8)
mtext(side = 2, text = expression(paste("Height increment (m year"^"-1", ")")),
      line = 2.7, cex = 1.8)

axis(side =1, cex.axis = 1.5)
axis(side =2, cex.axis = 1.5)


legend("topright", legend = c("Generalist species, 2006", "Generalist species, 2011", "Savanna species, 2006", "Savanna species, 2011"), col = c("#1b9e77", "#1b9e77", "#d95f02", "#d95f02"),
       lty = c(1, 2, 1, 2))


text(x = 1.7, y = 0.38, labels = "(b)", cex = 1.5)


dev.off()



#-------------------------------------------------------------------------------
# Figure 5 change in BA_above for trees in the field

svg(filename="./plots/Figure 5 change in ba_above.svg", 
    width = 6, 
    height=4, 
    pointsize=9)

  plot(xlim = c(2006, 2016),
       ylim = c(0, 1.6),
       x = NA,
       y = NA,
       xlab = "",
       ylab = "",
       axes = FALSE)
  par(mar = c(5.1, 5.1, 2.1, 2.1))
  
  #individual lines
  for(i in 1:length(unique(trees_to_use1))){
    lines(BA_above ~ Year, data = test[test$Unique_id == trees_to_use1[i], ],
          col = ifelse(test[test$Unique_id == trees_to_use1[i], ]$FT == "C",
                       addTrans("#d95f02", 50),
                       addTrans("#1b9e77", 50)))
  }

  mod1 <- lm(BA_above ~ Year*FT, data = test) 
  mod2 <- lm(BA_above ~ Year*FT, data = test2) 
  segments(2006, predict(mod1, newdata = list(Year = 2006, FT = "C")),
                    2016, predict(mod1, newdata = list(Year = 2016, FT = "C")),
           col = "#d95f02",
           lwd = 3)
  segments(2006, predict(mod1, newdata = list(Year = 2006, FT = "G")),
           2016, predict(mod1, newdata = list(Year = 2016, FT = "G")),
           col = "#1b9e77",
           lwd = 3)
  segments(2011, predict(mod2, newdata = list(Year = 2011, FT = "C")),
           2016, predict(mod2, newdata = list(Year = 2016, FT = "C")),
           col = "#d95f02",
           lwd = 3,
           lty = 2)
  segments(2011, predict(mod2, newdata = list(Year = 2011, FT = "G")),
           2016, predict(mod2, newdata = list(Year = 2016, FT = "G")),
           col = "#1b9e77",
           lwd = 3,
           lty = 2)
  axis(1, at = c(2006, 2011, 2016))
  axis(2, at = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5),
       labels = c(0, 2.5, 5, 7.5, 10, 12.5, 15))
  mtext(side = 1, text = "Year", line = 3.1, cex = 1.2)
  mtext(side = 2, text = expression(paste("BA above (m"^"2"," ha"^"-1",")")), line = 3.1, cex = 1.2)
  
  legend(x = 2006, y = 1.5, legend = c("Savanna species", "Generalist species", "Added 2011", "Added 2016"),
         lty = c(1,1,1,2),
         lwd = c(2,2,2,2),
         col = c("#d95f02", "#1b9e77", "black", "black"))

dev.off()


#-------------------------------------------------------------------------------
# Table S1 correlation matrix for predictors
#-------------------------------------------------------------------------------

#2006 data
cormat06 <- cor(data_06[, c("dg.2006.cm", "H.2006", "BA.above2006", "PC1", "PC2")], use = "pairwise.complete.obs")
cormat11 <- cor(data_11[, c("dg.2011.cm", "H.2011..m.", "BA.above2011", "PC1", "PC2")], use = "pairwise.complete.obs")
cormat11[upper.tri(cormat11)] <- cormat06[upper.tri(cormat06)]

diag(cormat11) <- c(cor(all_data$dg.2006.cm, all_data$dg.2011.cm, use = "pairwise.complete.obs"),
                              cor(all_data$H.2006, all_data$H.2011..m., use = "pairwise.complete.obs"),
                              cor(all_data$BA.above2006, all_data$BA.above2011, use = "pairwise.complete.obs"),
                              1,
                              1)
colnames(cormat11) <- colnames(cormat06)

cormat11

#-------------------------------------------------------------------------------
# Figure 6 height:diameter ratios
#-------------------------------------------------------------------------------

plot(log(H.2011..m.) ~ log(dg.2011.cm), data = data_11[data_11$C.Cerrado..G.generalist. == "C", ])
plot(log(H.2011..m.) ~ log(dg.2011.cm), data = data_11[data_11$C.Cerrado..G.generalist. == "G", ])


data_11$HD_ratio <- data_11$H.2011..m/data_11$dg.2011.cm
data_11$Savanna <- as.factor(ifelse(data_11$Plot %in% c(2,3,4,5,11:17,19), TRUE, FALSE))
data_11$small_diam <- ifelse(data_11$dg.2011.cm >= 5 & data_11$dg.2011.cm <= 7, TRUE, FALSE)

data_11_backup <- data_11

data_11 <- data_11[data_11$small_diam, ]


svg(filename="./plots/Figure 6 height diam allometry.svg", 
    width = 3, 
    height=3, 
    pointsize=9)

par(mar = c(5.1, 5.1, 2, 1))

boxplot(data_11[data_11$Savanna == TRUE &
                  data_11$C.Cerrado..G.generalist. == "C",
                "HD_ratio"],
        data_11[data_11$Savanna == TRUE &
                  data_11$C.Cerrado..G.generalist. == "G",
                "HD_ratio"],
        data_11[data_11$Savanna == FALSE &
                  data_11$C.Cerrado..G.generalist. == "C",
                "HD_ratio"],
        data_11[data_11$Savanna == FALSE &
                  data_11$C.Cerrado..G.generalist. == "G",
                "HD_ratio"],
        col = c("#d95f02", "#1b9e77"),
        xlab = "",
        ylab = "",
        ylim = c(0, 2.4))

data_11$Savanna <- as.factor(data_11$Savanna)
data_11$C.Cerrado..G.generalist. <- as.factor(data_11$C.Cerrado..G.generalist.)
mod <- aov(lm(HD_ratio ~ Savanna*C.Cerrado..G.generalist., data = data_11))
TukeyHSD(mod)

# text(x =  c(1,2,3,4), y = 2.3, c("a", "b", "c", "d"), )
axis(side = 1, at = c(1,2,3,4), labels = c("S", "G", "S", "G"))
mtext(text = "Savanna", side = 1, line= 3, at = 1.5)
mtext(text = "Forest", side = 1, line= 3, at = 3.5)
title(ylab = expression(paste("Height:diameter ratio (m cm"^"-1", ")")))

dev.off()


#-------------------------------------------------------------------------------
# Figure S5 plots of residuals ~ fitted

svg(filename="./plots/Figure S5 residuals vs fitted.svg", 
    width = 6, 
    height=6, 
    pointsize=9)
par(mfrow = c(2,2))

plot(residuals(ft_ci_soils_06) ~ fitted(ft_ci_soils_06),
     xlab = "Fitted",
     ylab = "Residuals",
     main = "Diameter growth, 2006-2011",
     pch = 21, 
     bg = addTrans("gray", 50),
     col = NA)
abline(h = 0)

plot(residuals(ft_ci_soils_11) ~ fitted(ft_ci_soils_11),
     xlab = "Fitted",
     ylab = "Residuals",
     main = "Diameter growth, 2011-2016",
     pch = 21, 
     bg = addTrans("gray", 50),
     col = NA)
abline(h = 0)

plot(residuals(ft_ci_soils_06_ht) ~ fitted(ft_ci_soils_06_ht),
     xlab = "Fitted",
     ylab = "Residuals",
     main = "Height growth, 2006-2011",
     pch = 21, 
     bg = addTrans("gray", 50),
     col = NA)
abline(h = 0)

plot(residuals(ft_ci_soils_11_ht) ~ fitted(ft_ci_soils_11_ht),
     xlab = "Fitted",
     ylab = "Residuals",
     main = "Height growth, 2011-2016",
     pch = 21, 
     bg = addTrans("gray", 50),
     col = NA)
abline(h = 0)

dev.off()

