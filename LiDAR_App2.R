library(lidR)
library(shiny)
library(ggplot2)
library(rlas)
library(raster)
library(rgdal)


las <- readLAS("D:/LiDARD_Testing/0609/0609_segmented.las", select = "xyzc") # input las at X date during the growing season
las2 <- readLAS("D:/LiDARD_Testing/0401/0401_clipped.las", select = "xyzc")  # input las of bare soil before growing season 
resolution =   #define resolution for raters
k <- .45       #extinction coeficient

----------------------------------------------------------------------------------------------------------

mycsf <- csf(sloop_smooth = FALSE, class_threshold = .1, cloth_resolution = 7, time_step = .65) #classify ground in the bar soil las
las_clas1 <- classify_ground(las2, mycsf)

dtm_tin <- grid_terrain(las_clas1, res = 2, algorithm = tin()) #creates DTM

dsm <- grid_canopy(las, res = 2, algorithm = p2r()) #creates dsm

lidar_chm <- dsm - dtm_tin #creates CHM

#plot the resulting chm
plot(lidar_chm,
     main = "Lidar Canopy Height Model (CHM)")

avg_height <- cellStats(lidar_chm, 'mean') #retreive average field height
avg_height

#parameterize teh Grid Resolution (GR) based on the crop height  

if(avg_height < 0.4){gr <- 0.1
} else if (avg_height < 0.6){gr <- 3
} else if (avg_height < 0.8){gr <- 4
} else if (avg_height < 1){gr <- 6}

gr


# ground point segmentation
mycsf <- csf(sloop_smooth = FALSE, class_threshold = .1, cloth_resolution = gr, time_step = .65)
las_clas2 <- classify_ground(las, mycsf) 


#plot segmentation to verify correct classification
plot_crossection <- function(las_clas2,
                             p1 = c(min(las_clas2@data$X), mean(las_clas2@data$Y)),
                             p2 = c(max(las_clas2@data$X), mean(las_clas2@data$Y)),
                             width =1, colour_by = NULL)
{
  colour_by <- enquo(colour_by)
  data_clip <- clip_transect(las_clas2, p1, p2, width)
  p <- ggplot(data_clip@data, aes(X,Z)) + geom_point(size = 1) + coord_equal() + theme_minimal() + ylim(102.5,103.5)
  
  if (!is.null(colour_by))
    p <- p + aes(color = !!colour_by) + labs(color = "")
  
  return(p)
}

#plot crossection of ground to other point segmentation
plot_crossection(las_clas2, p1 = p1, p2 = p2, colour_by = factor(Classification))

#filter pounts to only ground points
gnd <- filter_ground(las_clas2)
#plot(gnd, size = 3, bg = "white") 

#grid ground returns
grnd_returns <- grid_metrics(gnd, ~length(Z), 1) # calculate density
#plot(grnd_returns, col = gray.colors(50,0,1)) # some plotting
x <- reclassify(grnd_returns, cbind(-Inf, NA, 1), right=TRUE)
plot(x, col = gray.colors(50,0,1)) # some plotting

#grid all returns
all_returns <- grid_metrics(las, ~length(Z), 1) # calculate density
plot(all_returns, col = gray.colors(50,0,1)) # some plotting
#grnd_returns_add <- x + 1

#calculate Gap Fraction
GF <- x/all_returns
plot(GF , col = gray.colors(50,0,1)) # some plotting

#process for extracting average scan angle
angl <- readLAS("D:/LiDARD_Testing/0609/0609_segmented.las", filter = "a") 
angl_mean <- grid_metrics(angl, mean(ScanAngleRank), 1) #Grid the mean scan angle
abs_angl <- abs(angl_mean) #convert the absolute values
layer_angl_mean <- cellStats(abs_angl, stat='mean', na.rm=TRUE) #find mean scan angle of all cells 
layer_angl_mean
angl_rad <- layer_angl_mean * pi/180 # convert to radians
angl_rad
plot(abs_angl, col = gray.colors(50,0,1))

#calulate PAI
PAI_cal= -((cos(angl_rad)*log(GF))/k) #modified beer labert equation for calulating LAI
plot(PAI_cal,
     main = "Lidar PAI")
