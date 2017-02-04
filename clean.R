################
# PREPARATIONS #
################

# load package plyr for ddply function
library(plyr)

# read data
data <- read.csv("~/Desktop/GitHub/Gelada_parasites_hormones/data.csv")

# check dimensions
dim(data) # 2726 rows, 23 columns

# view first 6 rows
head(data)

# looks like there are some dates, let's format them
# note that the 'format' argument tells R where the month, day, and year are
# check the documentation for 'as.POSIXct' to see more on how to specify this
data$Date <- as.POSIXct(data$Date, format="%m/%d/%y")
data$DOB <- as.POSIXct(data$DOB, format="%m/%d/%y")

######################
# GO THROUGH COLUMNS #
######################

### 1. Name column

length(unique(data$Name)) # 150 unique individuals

# Rows for each individual should have the same birthdate, sex, and cyst status
# the following code counts the number of unique "DOB", "Sex" and "Cyst" exist for
# each Name
ddply(data, .(Name), function(x) {
  bday <- length(unique(x$DOB))
  sex <- length(unique(x$Sex))
  cyst <- length(unique(x$Cyst))
  data.frame(bday, sex, cyst)
})
# Oh, it appears that DOB and Sex are fine (there is one unique value per individual),
# but some individuals apparently change their cyst status over the study period.
# I guess that makes sense, and it shouldn't be a problem (more just wanted to demonstrate
# how this function works and can reveal how many different values a variable takes over
# repeated measurements from individuals). 

### 2. Date of sample column

# start by looking at range of dates
range(data$Date) # 2006-01-06 to 2014-02-20 KST... seems right?

# missing values?
sum(is.na(data$Date)) # nope... good, each sample should have a date

# histogram shows how sampling varies over time (breaks=365 means 365 bins for histogram)
hist(data$Date, breaks=365) # there is a lot of variation over time, and a big gap
# luckily, we can account for this gap (it was intentionally removed due to assay issues)
# everything seems to make sense

### 3. Sex

table(data$Sex) # everyone is male! makes sense since this is the T dataset 
sum(is.na(data$Sex)) # there are 0 NAs, good
# since Sex is invariant, we can drop it from the data 
# I'll drop variables at the end to avoid messing up the column numbers as I go

### 4. DOB

range(data$DOB) # 1993-01-04 to 2014-11-12... seems reasonable, indicates we have some young individuals

# now I'm going to do a quick sanity check... the sample collection date should always be LATER
# than the birthdate (i.e., an individual can't be sampled before they are born)
sum(data$DOB > data$Date) # OMG!!! Looks like we have one error!!! Where is it???
data[data$DOB > data$Date,] # It is "NON"... NON was apparently sampled 8 years before he was born
# Let's look at other rows for NON
data[data$Name=="NON",] # Looks like there is only one row for this guy
# Someone familiar with the monkeys will have to determine what is going on with this row
# For the record, I ALWAYS do these kinds of checks with dates... if you know some dates must
# logically come before others, then check to make sure that is true!

# Considering we have Age in the data, DOB is redundant and could be dropped

### 5. Age.at.smp

range(data$Age.at.smp) # range is 0.51 to 16.6
# That seems reasonable, except for the fact that we know this doesn't quite match up
# With the fact that if Age.at.smp was computed by subtracting the sample date from the
# birthdate, the minimum of the range should be a negative number. So how was this computed?
# Let's compute Age.at.smp using the dates in the dataframe
Age.at.smp2 <- difftime(data$Date, data$DOB, units="days")/365.25
# now let's make a scatterplot.. these should fall on a straight line
plot(Age.at.smp2 ~ data$Age.at.smp)
# Looks like they do, except for one weird value... that would have to be NON
# So the Age.at.smp variable seems to check out in general, and the only problem is NON
# I am going to rename Age.at.smp to simple "Age"... easier to type
names(data)[5] <- "Age"
names(data)

### 6. GC

# for a continuous variable, start with a simple histogram
hist(data$GC) # Seems reasonable.. looks like log transform would make more normal
hist(log(data$GC)) # Yup, seems like this one should be logged for analysis

# does GC depend on the location of the sample?
boxplot(log(GC) ~ GC.Loc, data=data) # looks like maybe it does... let's do a model to confirm
# we can just do a regular linear model on log transformed GC
summary(lm(log(GC) ~ GC.Loc, data=data)) # yup, there is a significant effect of location

### 7. GC.Loc

# let's use 'table' to see how many there are of each value
table(data$GC.Loc) # 163 Princeton, 2563 UM

# missing values?
sum(is.na(data$GC.Loc)) # nope, good

### 8. GC.Loc1

# This should just be a binary version of GC.Loc (why is this useful? i don't know)
# Let's look at a cross tabulation just to be sure
table(data$GC.Loc, data$Gc.Loc1) # Yup, looks fine
# Note the sloppy naming of variable (lowercase 'C' in GC)
# I don't see why this is useful, so will plan to delete

### 9. T.Loc

# let's use 'table' to see how many there are of each value
table(data$T.Loc) # 1462 UM, 1264 UN

# missing values?
sum(is.na(data$T.Loc)) # nope, good

### 10. T.Loc1

# Similar to GC.Loc1, this should just be a binary version of T.Loc
# Let's look at a cross tabulation just to be sure
table(data$T.Loc, data$T.Loc1) # Yup, looks fine
# I don't see why this is useful, so will plan to delete

### 11. T

# for a continuous variable, start with a simple histogram
hist(data$T) # Seems reasonable.. looks like log transform would make more normal
hist(log(data$T)) # Yup, seems like this one should be logged for analysis

# does T depend on the location of the sample?
boxplot(log(T) ~ T.Loc, data=data) # looks like maybe it does... let's do a model to confirm
# we can just do a regular linear model on log transformed T
summary(lm(log(T) ~ T.Loc, data=data)) # yup, there is a significant effect of location

### 12. T.antibody

# I don't know what this is... I suspect it maps onto T location... let's find out
table(data$T.Loc, data$T.antibody) # looks like it almost does, but there is one weird value
# Where a Pantex assay was performed at UN (all other Pantex were done at UM, and all other
# UN samples were DSL rather than Pantex)
# Let's find that value in the data
data[data$T.Loc=="UN" & data$T.antibody=="Pantex",] # huh, looks like its the last row (2726)
# Seems likely that this is an error

### 13. avgrain

range(data$avgrain) # hmm, got an error... looks like there are some NAs, let's investigate
sum(is.na(data$avgrain)) # there are 27 NAs... let's pull them up
data[is.na(data$avgrain),] # hmm, not sure why these are here

# plot histogram for the non-NA avgrain values
hist(data$avgrain) # seems reasonable

### 14. MaxT

range(data$MaxT) # 9.9 - 27.2
hist(data$MaxT) # seems reasonable

### 15. MinT

range(data$MinT) # 4.7 - 71
hist(data$MinT) # seems reasonable

### 16. Rainfall

range(data$Rainfall) # looks like there are NAs again
sum(is.na(data$Rainfall)) # there are 11 NAs... let's pull them up
data[is.na(data$Rainfall),] # hmm, not sure why these are here

# plot histogram for the non-NA Rainfall values
hist(data$Rainfall) # seems reasonable

### 17. RMedian

range(data$RMedian) # seems fine

# plot histogram for RMedian
hist(data$RMedian) # hmm, seems like a very jumpy distribution... why?

### 18. wetvsdry

# NOTE: this is a bad variable name because it is ambiguous what 0 and 1 mean
# if 1 means "wet", then the variable name should just be "wet"
# I will rename it now
names(data)[18] <- "Wet"

table(data$wetvsdry) # 1284=0 wet, 1415=1
sum(is.na(data$wetvsdry)) # there are 27 NAs, must be because of the avgrain NAs

# Let's check that it was done right... these statements should evaluate TRUE
all(data[data$avgrain < data$RMedian,"Wet"] == 0)
all(data[data$avgrain > data$RMedian,"Wet"] == 1)
# looks right!

### 19. MaxTMedian

range(data$MaxTMedian) # 15.8-19.2 seems fine

### 20. MinTMedian

range(data$MinTMedian) # 7.1 - 8.8 seems fine

### 21. Cyst

table(data$Cyst) # 2497=0, 229=1
sum(is.na(data$Cyst)) # no NAs, good

### 22. Status

# Should be B=Bachelor, L=Leader, F=Follower, N=Natal, SA=subadult
table(data$Status) # Hmm, seems like there are only 4 categories (B,L,F, and NonA)
# Is 'NonA' combining natal and subadult?

### 23. Temp

table(data$Temp) # 924=0, 897=1, 905=2
sum(is.na(data$Temp)) # no NAs, good

# Let's check that it was done right... these statements should all evaluate TRUE
all(data[data$MinT < data$MinTMedian & data$MaxT < data$MaxTMedian,"Temp"] == 0)
all(data[data$MinT < data$MinTMedian & data$MaxT > data$MaxTMedian,"Temp"] == 1)
all(data[data$MinT > data$MinTMedian & data$MaxT > data$MaxTMedian,"Temp"] == 2)
# looks right!

###########
# SUMMARY #
###########

# Here are some potential issues/questions found from my checking of all the columns:

# 1. Problem with either NON's birthdate or sample date
# 2. Possible error in either T.Loc or T.antibody for the 2726th row
# 3. There are 27 NAs for avgrain and wetvsdry: data[is.na(data$avgrain),]
# 4. There are 11 NAs for Rainfall
# 5. What is up with the very jumpy distribution for median rainfall (RMedian)?
# 6. What is "NonA" under "Status"? The description lead me to expect N and SA, not NonA

# Several columns may as well be dropped: Sex, DOB, Gc.Loc1, T.Loc1
drop <- c(3, 4, 8, 10)

# Additional, climate variables can be dropped once we are happy with "Wet" and "Temp"
drop <- c(drop, 13, 14, 15, 16, 17, 19, 20)

# I also suspect that T.antibody is not useful, but need to figure out whether that weird
# mismatch between T.antibody and T.Loc is an error, and if so, which is wrong, T.antibody
# or T.Loc? Once this is determined, drop T.antibody too
drop <- c(drop, 12)

# This should bring us down to 11 columns
# But I won't drop the columns and export the new data file until we have addressed
# the 6 issues/questions raised above

#######
# END #
#######




