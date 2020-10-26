Thomas Fire
================
Simran Bawa
10/22/2020

``` r
library(readxl)
thomas.fire <- read_excel("~/Desktop/EEMB 144L/144l_students/Input_Data/week1/Thomas_Fire_Progression.xlsx")
View(thomas.fire)

acres.burned <- thomas.fire$Acres_Burned
containment <- thomas.fire$Containment
plot(acres.burned, containment,
     xlab = "Acres Burned",
     ylab = "Containment",
     main = "Containment of Fire as Acres Burned")
```

![](ThomasFire_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

To represent the data in the given data set, I chose to use a scatter
plot with the independent variable being the amount of acres burned in
the Thomas Fire and the dependent variable being the amount of the fire
that was contained. This shows how the fire was most contained once it
has burned about 260,000 acres. The fire had spread rapidly in the very
beginning and this is why it grew out of control. Once more than 40% of
the fire was contained, there was a sharp increase of the percent of the
fire that was contained and a low increase in the number of acres
burned.
