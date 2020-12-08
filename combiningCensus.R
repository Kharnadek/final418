#These steps will help you combine the outputs 
#from your spatial interpolation with your income data.
# Convert your interpolation into a raster and map it:
r <- raster(P.idw)
r.m <- mask(r, income.tracts)
surfaceMap <-tm_shape(r.m) + 
  tm_raster(n=5, palette="Blues",  
            title="IDW \nPredicted PM2.5 \n(in ppm)") +
  tm_shape(pm2.5) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)
surfaceMap
#If you have too many cells, 
#you can reduce the number by aggregating values
#agg <- aggregate(yourRasterFromKriging, fact=??, fun=mean)

#Extract average pm2.5 for each polygon
income.tracts$Pm2.5 <- round(extract(r, income.tracts, fun = mean)[,1], 5)

map_pm25 <- tm_shape(income.tracts) +
  tm_polygons(col = "Pm2.5",
              title = "Mean PM2.5 Per Tract",
              style = "jenks",
              palette = "seq", n = 6) +
  tm_legend(legend.position = c("LEFT", "BOTTOM")) +
  tm_layout(aes.palette = list(seq = "-RdBu"))

map_pm25

