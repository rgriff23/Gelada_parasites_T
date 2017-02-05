
# read data
data <- read.csv("~/Desktop/GitHub/Gelada_parasites_T/data.csv")
weather <- read.csv("~/Dropbox/Geladas/Cortisol Analysis/weather_091516.csv")

# format dates
data$Date <- (data$Date, format="%m/%d/%y")
weather$Date <- as.POSIXct(weather$Date, format="%m/%d/%y")

# empty column
data$Rain_median_30days
data$Rain_count

# compute median rainfall from past 30 days for each row
for (i in 1:nrow(data)) {
	end <- data[i,"Date"]
	start <- end - as.difftime(30, unit="days")
	rain <- weather[weather$Date > start & weather$Date < end,"Rainfall"]
	data$Rain_count[i] <- 30 - sum(is.na(rain))
	if (data$Rain_count[i] > 9) {
		data$Rain_median_30days[i] <- median(rain, na.rm=T)
	} else {data$Rain_median_30days[i] = NA}
}

# export
write.csv(data, file="~/Desktop/GitHub/Gelada_parasites_T/data.csv", row.names=FALSE)

