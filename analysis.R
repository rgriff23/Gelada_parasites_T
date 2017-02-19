#################
# PREPARATIONS #
###############

# Load packages
library(lme4) # GLMMs
library(MuMIn) # model comparison

# Read data
data <- read.csv("~/Desktop/GitHub/Gelada_parasites_T/data.csv")

#####################
# MODIFY VARIABLES #
###################

# Cut non adults
#data <- data[data$Status != "NonA",]

# Center age (mean age is 10.6 years)
#data$Age <- data$Age - mean(data$Age)

# Center temperature (midpoint is 13 degrees C)
#data$Temp_midpoint <- data$Temp_midpoint - mean(data$Temp_midpoint)

# New leadership variable
#data$Leader <- ifelse(data$Status=="L", 1, 0)

# Export
#write.csv(data, file="~/Desktop/GitHub/Gelada_parasites_T/data.csv", row.names=F)

###########
# MODELS #
#########

# Subset data and remove NAs
data_sub <- data[,c("Z_t", "Cyst", "Age", "Leader", "Z_gc", "Temp_midpoint", "Rain_median_30days", "Name")]
data_sub <- data_sub[complete.cases(data_sub),]

# GLMM
mod <- lmer(Z_t ~ (Cyst + Age + Leader + Z_gc + Temp_midpoint + Rain_median_30days)^2 + (1|Name), data=data_sub, na.action=na.fail)
summary(mod)
dredge1 <- dredge(mod)

# view top models (within 2 AICc of top model)
get.models(dredge1, subset=delta<2)

# model averaging
mod.average <- model.avg(get.models(dredge1, subset=delta<2))
summary(mod.average)

########
# END #
######

# Look at males with cysts, sort by name and date
data$Date <- as.Date(data$Date, format="%m/%d/%y")
cystmales <- data[data$Cyst==1,c("Name", "Date", "Status")]
cystmales <- cystmales[order(cystmales$Name, cystmales$Date),]
cystmales

