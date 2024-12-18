
library(ncdf4) # package for netcdf manipulation
library(geodata)
library(rgdal) # read shape files
library(raster)

shape <- readOGR("data/shapefiles/67100000.shp")
Files <- dir("data/landuse/netcdf")
Border <- world(path=tempdir())
Border <- world(path=tempdir())


roundUp <- function(x, to = 10){
  to*(x%/%to + as.logical(x%%to))
}

ex <- extent(shape)
ex[1] <- -64
ex[2] <- -52
ex[3] <- -23
ex[4] <- -14

MaximumLandUse <- 70
cuts <- seq(0,MaximumLandUse,10)
rgb.palette <- colorRampPalette(c(rgb(0.8,0.8,0.8,0.8),"seagreen","orange","firebrick"),
                                space = "rgb")
Time <- list()
Average_LandUse <- empty_list <- vector(mode='list', length=length(Files))
Max <- c(10,10,40,40,30,10,10)
for(i in 1:length(Files)){
  
  Title <- strsplit(strsplit(Files[i], "LU")[[1]][2], "19")[[1]][1]
  if(is.na(Title)) next
  
  nc_data <- nc_open(paste0("data/landuse/netcdf/",Files[i]))
  
  #{
  #  sink('data/landuse/netcdf/nc_metadata.txt')
  #  print(nc_data)
  #  sink()
  #}
  
  t <- ncvar_get(nc_data, "time")
  t2 <- seq_along(t)-1
  Dates <- as.Date(paste0(t,"-01-01"))
  
  Time[[i]] <- t
  
  lonIdx <- which(nc_data$dim$lon$vals > ex[1] & nc_data$dim$lon$vals < ex[2]) ##
  latIdx <- which(nc_data$dim$lat$vals > ex[3] & nc_data$dim$lat$vals < ex[4]) ##
  
  # https://rstudio-pubs-static.s3.amazonaws.com/538346_6c663fca47c24288b56d4294618fcf7d.html
  landuse <- ncvar_get(nc_data,
                       names(nc_data$var),
                       start=c(min(lonIdx),
                               min(latIdx),
                               min(t2)+1),
                       count=c((max(lonIdx) - I(min(lonIdx)+1)),
                               (max(latIdx) - I(min(latIdx)+1)),
                               (max(t2)-min(t2)+1)))
  
  
  lon <- ncvar_get(nc_data, "lon", start=c(min(lonIdx)),
                   count=c((max(lonIdx) - I(min(lonIdx)+1))))
  lat <- ncvar_get(nc_data, "lat", start=c(min(latIdx)),
                   count=c((max(latIdx) - I(min(latIdx)+1))))
  
  fillvalue <- ncatt_get(nc_data, names(nc_data$var), "_FillValue")
  landuse[landuse == fillvalue$value] <- NA
  
  Sub_Average_LandUse <- array(NA, length(t))
  
  
  for(j in 1:length(t)){
    print(paste("Map",t[j], "File =",Title))
    landuse_slice <- landuse[, , j] 
    
    r <- raster(t(landuse_slice),
                xmn=min(lon),
                xmx=max(lon),
                ymn=min(lat),
                ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
    
    rr <- flip(r)
    r_pantanal <- crop(rr, extent(shape))
    LU_paraguai <- raster::mask(r_pantanal, shape)
    med_LU_paraguai <- mean(getValues(LU_paraguai), na.rm = TRUE)
    Sub_Average_LandUse[j] <- med_LU_paraguai
    
    png(paste0("plots/landuse/",Title,"/",Title,"_",t[j],".png"),
        height = 800*5, width = 700*5, res = 75 * 9)
    
    par(mfrow = c(2,1), mar = c(4,4,2,1))
    
    plot(rr, main = paste0("Land use area for the variable ",Title, " (%)"),
         bty = "n", 
         col = rgb.palette(10),
         xlim = c(-65,-52), ylim = c(-23,-12), zlim = c(0,MaximumLandUse))
    
    plot(shape, add = TRUE, lwd = 3)
    plot(Border, add=TRUE, lwd=1, border = "darkgrey")
    
    mtext("Longitude", side = 1, line = 2.3)
    mtext("Latitude", side = 2, line = 2.3)
    
    plot(t, Sub_Average_LandUse, type = "l", bty = "n",
         ylim = c(0, Max[i]),
         xlim = c(min(t, na.rm = TRUE), 2020),
         ylab = "", xlab = "", col = "seagreen", lwd = 3)
    
    mtext("Year", side = 1, line = 2.3)
    mtext("Land use occupation (%)", side = 2, line = 2.3)
    
    dev.off()
    
  }
  
  Average_LandUse[[i]] <- Sub_Average_LandUse
  
}

Max <- sapply(Average_LandUse, max, na.rm = TRUE)
Max <- ifelse(Max == "-Inf", 8, Max)
Max <- roundUp(Max)

names(Average_LandUse) <- Files

for(i in 1:7){
  Average_LandUse[[i]] <- data.frame("Time" = Time[[i]],
                                     "LU" = Average_LandUse[[i]])
}
i <- 1
for(i in 1:7){
  Title <- strsplit(strsplit(Files[i], "LU")[[1]][2], "19")[[1]][1]
  write_csv2(Average_LandUse[[i]], paste0("tables/",Title,"_Dias.csv"))
}

