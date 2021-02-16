## code to prepare `wildfires` dataset goes here

fires <- read.csv(file = "data-raw/Incendios.csv", sep = ",", dec = ".", na.strings="", header = TRUE,
                      colClasses = c("integer","character","character","character","character"))
fires <- na.omit(fires)
fires$Date <- as.character(fires$Date)
fires$Month <- substr(x = fires$Date, start = 4, stop = 5)
fires$Year_Month <- paste0(fires$Year,"-",fires$Month)

wildfires_new <- data.frame(Year = c(2007L, 2008L, 2008L, 2010L, 2010L),
                        Date = c(NA, NA, NA, NA, NA),
                        Code = c(NA, NA, NA, NA, NA),
                Municipality = c(NA, NA, NA, NA, NA),
                  Total.Area = c(0, 0, 0, 0, 0),
                       Month = c("08", "06", "11", "10", "12"),
                  Year_Month = c("2007-08","2008-06","2008-11","2010-10","2010-12"))

fires <- rbind(fires, wildfires_new)
fires <- fires[with(fires, order(Year, Month)), ]
fires$GIF <- ifelse(fires$Total.Area > 500, 1, 0)

wildfires <-aggregate(GIF ~ Year_Month, data = fires, sum)
colnames(wildfires) <- c("Date","GIF")
wildfires$GIF <- as.integer(wildfires$GIF)
