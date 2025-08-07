# #!/bin/bash

# datadir="/mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/81499/released/2022-06-07/data"

# head -n 1 $datadir/data.51913.csv | tr ',' '\n' > ukb.datafields

# less ukb.datafields 

# grep -n "^\"129-" ukb.datafields | cut -d ":" -f 1 
# grep -n "^\"130-" ukb.datafields | cut -d ":" -f 1 


# cut -d "," -f 1,435-440 $datadir/data.51913.csv > ukb.coords
# cut -d "," -f 10005-10044 $datadir/data.51913.csv > ukb.pcs






library(dplyr)
library(data.table)
library(ggplot2)

coords <- fread("ukb.coords") %>% filter(!is.na(`130-0.0`))
names(coords)[2] <- "northing"
names(coords)[5] <- "easting"

pcs <- fread("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/principal_components/data.pca1-40.plink.txt")

linker <- fread("/mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/81499/released/2022-06-07/data/linker.81499.csv")
pcs <- inner_join(pcs, linker, by=c("V1"="ieu"))
pcs
dat <- inner_join(pcs, coords, by=c("app"="eid"))



generate_variant_based_on_location <- function(coord1, coord2, centroid, sharpness, maxfreq=0.1)  {
    # Calculate the distance from the centroid
    distance_from_centroid <- sqrt((coord1 - centroid[1])^2 + (coord2 - centroid[2])^2)
    print(summary(distance_from_centroid))
    # Frequency of the variant reduces with distance according to sharpness
    frequency <- distance_from_centroid - min(distance_from_centroid, na.rm=T)
    frequency <- 1 - frequency / max(frequency, na.rm=T)
    print(summary(frequency))
    hist(frequency, breaks=100)

    variant <- rbinom(length(coord1), size = 2, prob = frequency)
    return(variant)
}

i <- sample(1:nrow(dat), 1)
variant <- generate_variant_based_on_location(dat$northing, dat$easting, c(dat$northing[i], dat$easting[i]), 1, 0.1)
tibble(variant, northing=dat$northing, easting=dat$easting) %>%
ggplot(., aes(x=easting, y=northing)) +
stat_summary_hex(aes(z=variant), fun = mean) +
geom_point(data=dat[i,], aes(x=easting, y=northing), colour="red")





summary(lm(variant ~ as.matrix(dat[,3:42])))

phen <- 


tibble(variant, northing=dat$northing, easting=dat$easting) %>%
ggplot(., aes(x=easting, y=northing)) +
geom_hex()




