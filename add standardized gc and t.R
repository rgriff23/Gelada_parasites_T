# Read data
data <- read.csv("~/Desktop/GitHub/Gelada_parasites_T/data.csv")

# Add column for standardized log GC 
gc.princeton.mean <- mean(log(data[data$GC_loc=="Princeton","GC"]))
gc.princeton.sd <- sd(log(data[data$GC_loc=="Princeton","GC"]))
gc.um.mean <- mean(log(data[data$GC_loc=="UM","GC"]))
gc.um.sd <- sd(log(data[data$GC_loc=="UM","GC"]))
data$Z_gc <- NA
data$Z_gc[data$GC_loc=="Princeton"] <- (log(data$GC[data$GC_loc=="Princeton"]) - gc.princeton.mean)/gc.princeton.sd
data$Z_gc[data$GC_loc=="UM"] <- (log(data$GC[data$GC_loc=="UM"]) - gc.um.mean)/gc.um.sd
sum(is.na(data$Z_gc))

# visual confirmation that GC standardization worked (new data has mean 0 and sd=1)
layout(matrix(1:2,1,2))
boxplot(log(GC) ~ GC_loc, data=data, main="unstandardized")
boxplot(Z_gc ~ GC_loc, data=data, main="standardized")

# Add column for standardized log T
t.pantex.mean <- mean(log(data[data$T_antibody =="Pantex","T"]))
t.pantex.sd <- sd(log(data[data$T_antibody =="Pantex","T"]))
t.dsl.mean <- mean(log(data[data$T_antibody =="DSL","T"]))
t.dsl.sd <- sd(log(data[data$T_antibody =="DSL","T"]))
data$Z_t <- NA
data$Z_t[data$T_antibody=="Pantex"] <- (log(data$T[data$T_antibody =="Pantex"]) - t.pantex.mean)/t.pantex.sd
data$Z_t[data$T_antibody=="DSL"] <- (log(data$T[data$T_antibody =="DSL"]) - t.dsl.mean)/t.dsl.sd
sum(is.na(data$Z_t))

# visual confirmation that T standardization worked (new data has mean 0 and sd=1)
layout(matrix(1:2,1,2))
boxplot(log(T) ~ T_antibody, data=data, main="unstandardized")
boxplot(Z_t ~ T_antibody, data=data, main="standardized")

# export
write.csv(data, file="~/Desktop/GitHub/Gelada_parasites_T/data.csv", row.names=F)
