## Extract the future climate data from WorldClim into individual bioclimatic variable tif files

library(raster)

r <- stack("wc2.1_5m_bioc_MIROC6_ssp585_2081-2100.tif")
nlayers(r)
for(i in 1:nlayers(r)){
  band <- r[[i]]
  writeRaster(band,paste('bio',i,'.tif', sep=''))
}
