#### Intro to R with the Cal Fire Data ####
# Simran Bawa
# 10/11/2020

#### Load Packages ####

#install.packages("tidyverse")
library(tidyverse)
#install.packages("dplyr")
library(dplyr)
#install.packages("readxl")
library(readxl)
#install.packages("praise")
library(praise)

#### Load Dataset ####

excel_sheets("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx")

metadata <- read_excel("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx", sheet = 1)

data <- read_excel("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx", sheet = 2)

#### Initial Data Exploration ####

names(data)
dim(data)
class(data)
head(data)
tail(data)
str(data)
glimpse(data)
typeof(data$Total_Acres_Burned)

max(data$Total_Acres_Burned)
max(data$Structures_Destroyed)
max(data$Structures_Destroyed, na.rm = T)

summary(data)

#### Basic data wrangling (dplyr functions) ####

df1 <- select(data, County_Unit:Controlled_Date, Total_Acres_Burned, Cause:Structures_Damaged)
unique(df1$County_Unit)

df2 <- filter(df1, County_Unit %in% c("SANTA BARBARA", "VENTURA", "LOS ANGELES", "SAN DIEGO", "ORANGE", "VENTURA/SANTA BARBARA") & Total_Acres_Burned >= 500)
unique(df2$Total_Acres_Burned)

# | = "or"
# == = "equals"/"matches" , %in% c()
# & = "and"

df3 <- arrange(df2, desc(Total_Acres_Burned))

df4 <- mutate_at(df3, vars("Structures_Destroyed", "Structures_Damaged"), replace_na, 0)

df5 <- mutate(df4, struc_impact = Structures_Damaged + Structures_Destroyed)

#mess with time
library(lubridate)

df6 <- mutate(df5, interv = interval(Start_Date, Controlled_Date), 
               dur = as.duration(interv),
               days = as.numeric(dur, "days"))

#### Introduction to Piping ####

#%>% # cmd + shift + m 
socal.fires <- data %>% 
  select(County_Unit:Controlled_Date, Total_Acres_Burned, Cause:Structures_Damaged) %>% 
  filter(County_Unit %in% c("SANTA BARBARA", "VENTURA", "LOS ANGELES", "SAN DIEGO", "ORANGE", "VENTURA/SANTA BARBARA") & Total_Acres_Burned >= 500) %>% 
  arrange(desc(Total_Acres_Burned)) %>% 
  mutate_at(vars("Structures_Destroyed", "Structures_Damaged"), replace_na, 0) %>% 
  mutate(struc_impact = Structures_Damaged + Structures_Destroyed,
         interv = interval(Start_Date, Controlled_Date), 
         dur = as.duration(interv),
         days = as.numeric(dur, "days"))
  
#### Our first graph in ggplot ####

ggplot(socal.fires, aes(x = Start_Date, y= Total_Acres_Burned)) +
  geom_point(aes(color = County_Unit)) +
  ggtitle("CA South Coast Major Fires \n2014 - 2018") +
  labs(x = "", y = "Total Acres Burned", Color = "County") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_grid(rows = "County_Unit", scales = "free")
  
incidents <- socal.fires %>% 
  rename(county = County_Unit,
         acres = Total_Acres_Burned,
         start = Start_Date,
         end = Controlled_Date) %>% 
  mutate(year = year(start),
         county = ifelse(county == "VENTURA/SANTA BARBARA", "VENTURA", county))
 
incidents2 <- incidents %>%  
  group_by(county, year) %>% 
  tally() %>% 
  ungroup()
  
incidents2.plot <- incidents2 %>% 
  ggplot(aes(x = year, y = n)) +
  geom_point(aes(color = county)) +
  geom_point() +
  geom_line() +
  labs(title = "CA South Coast Major Fire Incidents \n 2013 - 2018", x = "", y = "Incidents", color = "County") +
  theme_bw() + 
  facet_grid(rows = "county", scales = "free")
  # guides(color = F)
  
  all_incidents <- incidents2.plot %>% 
  group_by(year) %>% 
  tally() %>% 
  ungroup()
  
all_incidents.plot <- incidents %>% 
  ggplot(aes(x = year, y = n)) +
  geom_point(aes(color = "blue")) +
  geom_line(color = "blue") +
  labs(title = "CA South Coast Major Fire Incidents \n 2013 - 2018", x = "", y = "Incidents") +
  theme_bw()

#### Save Data and Plots ####

saveRDS(socal.fires, file = "~/Desktop/GitHub/144l_students/Output_Data/socal_fires_data.rds")












