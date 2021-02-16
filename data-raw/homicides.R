## code to prepare `homicides` dataset goes here

# Transform data base
homicides <- read.csv("data-raw/Homicidios.csv", sep = ";", dec = ",")
homicides$Rate <- homicides$Murder / homicides$Population * 100000
homicides <- homicides[14:72,]

