#### Project Solaris ####
## Initial setup and libraries

# Clean enviromment
rm(list=ls(all=TRUE)) 
cat("\014")  
options(scipen = 999)             # Modify global options in R
options(sf_max.plot=1)
# dev.off()

# Load libraries
library(tidyverse)
library(sf)
library(osmdata)
library(here)
library(ncdf4) #library to read and process netcdf data
library(chron) #deal with chronological objects
library(lattice)
library(RColorBrewer)
library(tmap)
library(readxl)
library(writexl)
library(raster)
library(gstat)
library(terra)
library(geodata)
library(chron)
library(esquisse)
library("ggspatial")

# Write functions
'%ni%' <- Negate('%in%')

## read Indonesia boundaries
boundary = st_read(here("input","indonesia_geo.gpkg"))

#### land use analysis ####
land_types = read_excel(here("input","land_cover","lccs_types.xlsx"),sheet = "lccs")

land <- nc_open(here("input","land_cover", "C3S-LC-L4-LCCS-Map-300m-P1Y-2022-v2.1.1.area-subset.6.95.-12.180.nc"))
#land

#get Dimension :lon, lat
lon = ncvar_get(land, "lon")
lat = ncvar_get(land, "lat")

land_array = ncvar_get(land,"lccs_class") #get the land cover class

laname = ncatt_get(land,"lccs_class","long_name")
launits = ncatt_get(land,"lccs_class","units")
lafillvalue = ncatt_get(land,"lccs_class","_FillValue")

lonlat = as.matrix( (expand.grid(lon, lat))) #lon and lat are what we extracted in step 2.
land_vec = as.vector(land_array)

land_df = data.frame( cbind( lonlat,land_vec  ))
colnames(land_df) = c("lon", "lat", "lccs")

rm(lafillvalue, laname, land, land_array, launits, lonlat)
rm(land_vec)

# take away stuff that is not indonesia
tmp = st_bbox(boundary)
xmin = tmp[1]
ymin = tmp[2]
xmax = tmp[3]
ymax = tmp[4]

land_df = land_df %>% dplyr::filter(lat < ymax)
land_df = land_df %>% dplyr::filter(lat > ymin)
land_df = land_df %>% dplyr::filter(lon < xmax)
land_df = land_df %>% dplyr::filter(lon > xmin)

# reduce dataframe to be able to process
#land_df = land_df[seq(1, nrow(land_df), 100), ]
land_df = sample_n(land_df, 5000000)

# replace lccs with suitable yes/no
tmp = land_types %>% dplyr::select(lccs, suitable2) 
land_df =left_join(land_df, tmp, by = "lccs")
land_df = land_df %>% dplyr::select(-(lccs)) %>% rename(ok = suitable2)
land_sf = st_as_sf(land_df, coords = c("lon","lat"), crs = 4326) 
#st_write(land_sf, here("output","land_cover.gpkg"),layer = "lccs", append = FALSE)
land_sfU = st_transform(land_sf, 23879)

tmp = st_bbox(land_sfU)
xmin = tmp[1]
ymin = tmp[2]
xmax = tmp[3]
ymax = tmp[4]
raster_template = rast(resolution = 20000 , xmin= xmin, ymin=ymin,
                       xmax=xmax, ymax=ymax,
                       crs = st_crs(land_sfU)$wkt)

rasterized <- rasterize(land_sfU, raster_template, field = "ok")
boundaryU = st_transform(boundary,  crs = 23879)
r_mask <- mask(rasterized, boundaryU)
plot(r_mask)

# raster::writeRaster(r_mask, 
#                     filename=here("tmp", "landuse_class.tif"), 
#                     overwrite=TRUE)

rm(land_df, land_sf, raster_template)

#### solar radiance analysis with year-round values ####

# read one file for the 10th of Jan, April, July, October 2023
# time of day: 0800, 1200, 1600, 2000
### define geo area for copernicus data
tmp = st_bbox(boundary)
xmin = tmp[1] # LON / East / X => west
ymin = tmp[2] # LAT / North / Y => south
xmax = tmp[3] # east
ymax = tmp[4] # north

xmin
# west = 95
ymin
# south = -11
xmax
# east = 142
ymax
# north = 6

tera01 = rast(here("input","solar","23-01.10-data_0.nc"))
tera02 = rast(here("input","solar","23-04-10-data_0.nc"))
tera03 = rast(here("input","solar","23-07-10-data_0.nc"))
tera04 = rast(here("input","solar","23-10-10-data_0.nc"))

# merge into one file
tera_all = c(tera01, tera02, tera03, tera04)

# calculate the mean ssrd and check results
tera_mean <- mean(tera_all, na.rm = TRUE)
dim(tera_mean)
freq(tera_mean)
plot(tera_mean)

## make sf object from spatraster 
# step 3 turn into sf object
tmp2 = crds(tera_mean, df = TRUE, na.all = TRUE)
tmp3 = as.data.frame(tera_mean)
ssrd_df = cbind(tmp2, tmp3)

colnames(ssrd_df) = c("lon", "lat", "ssrd")
ssrd_df_value = na.omit (ssrd_df)
head(ssrd_df_value, 3) 

ssrd_sf = st_as_sf(ssrd_df_value, coords = c("lon","lat"), crs = 4326)
plot(ssrd_sf)

# Radiation to power function
# an example of a 1m2 (A) solar panel - kwh per square meter per hour

radiation_to_power <- function(G, A=1, r=0.175, p=0.6, hours=1){
  kWh <- G * A * r * p * (hours/3600) / 1000
  return(kWh)
}

ssrd_kwh <- as.data.frame (radiation_to_power (ssrd_df_value))
ssrd_df_value <- cbind(ssrd_df_value,ssrd_kwh$ssrd)
colnames(ssrd_df_value) [4] <- 'ssrd_kwh'
ssrd_sf$ssrd_kwh = ssrd_kwh$ssrd

ssrd_sf <- ssrd_sf[, -which(names(ssrd_sf) == "ssrd")]
plot(ssrd_sf)
#st_write(ssrd_sf, here("output","ssrd_kwh.gpkg"),layer = "kwh", append = FALSE)
rm(tera_all, tera_mean, tera01, tera02, tera03, tera04)

## IDW infill ##
# Q2: How to decide the distance decay parameters for spatial interpolation model based on IDW? 
#ssrd_sf = read_sf(here("output","ssrd_kwh.gpkg"))

# step 0: convert to UTM, EPSG:23879 for Indonesia
ssrd_sfU = st_transform(ssrd_sf, crs = 23879)
ssrd_sfU

# step 1: IDW_extract point cooridnates

coor = as.data.frame(st_coordinates(ssrd_sfU))
ssrd_sfU$x = coor$X
ssrd_sfU$y = coor$Y
ssrd_sf_nogeomU = st_drop_geometry(ssrd_sfU)
#pts_nogeom is a data frame with its original attributes and also contains 
# coordinates columns x and y

# step 2: IDW_create gstat object
#data should be in data frame format
#idp is the distance power 1/d^k, so idp here indicate k value
gs <- gstat(formula=ssrd_kwh~1, locations=~x+y, data=ssrd_sf_nogeomU, nmax=Inf,set=list(idp=2)) 
gs

#step 3: IDW_spatial interpolation
#step 3.1 IDW_create raster template

#tmp = st_bbox(ssrd_sfU)
tmp = st_bbox(land_sfU)
xmin = tmp[1]
ymin = tmp[2]
xmax = tmp[3]
ymax = tmp[4]
raster_template = rast(resolution = 20000 , xmin= xmin, ymin=ymin,
                       xmax=xmax, ymax=ymax,
                       crs = st_crs(ssrd_sfU)$wkt)
raster_template

# step 3.2 IDW_interpolation
#interpolate is the function comes with terra
idw <- interpolate(raster_template, gs, debug.level=0) 
#plot(idw$var1.pred)

# step 3.3 mask and write to file + cleanup
boundaryU = st_transform(boundary,  crs = 23879)
idw_mask <- mask(idw, boundaryU)

names(idw_mask)
plot(idw_mask$var1.pred)

idw_mask <- subset(idw_mask, "var1.var", negate=TRUE)

raster::writeRaster(idw_mask, 
                    filename=here("tmp", "idw2_class_mean.tif"), 
                    overwrite=TRUE)

rm(coor, gs, ssrd_sf, raster_template, ssrd_sfU, ssrd_sf_nogeomU)
rm(lat,lon, xmax, xmin, ymax, ymin)

#### grid analysis and buffering ####
## read grid data and clean up
grid = st_read(here("input","grid.geojson")) 
grid = grid %>% filter(is.na(disused)) # remove disused cables

grid = grid %>% dplyr::select(voltage)
st_write(grid, here("input","grid.gpkg"),layer = "grid", append = FALSE)
gridU = st_transform(grid, crs = 23879)

## calculate grid line length
grid_m <- st_length(gridU)
grid_km = as.numeric(sum(grid_m) / 1000)
print(grid_km)

grid10km = st_buffer(gridU, dist = 10000) %>% st_union()
#plot(grid500)
grid25km = st_buffer(gridU, dist = 25000) %>% st_union()

tmp = st_difference(grid25km, grid10km)

st_write(grid10km, here("tmp","grid10km.gpkg"),layer = "a", append = FALSE)
st_write(grid25km, here("tmp","grid25km.gpkg"),layer = "a", append = FALSE)
st_write(tmp, here("tmp","grid25km_net.gpkg"),layer = "a", append = FALSE)

rm(tmp, grid, grid10km, grid25km)

#### buffer roads ####
road = read_sf(here("input","IDN_roads","IDN_roads.shp"))
road = st_transform(road, crs = 23879)
road %>% count(RTT_DESCRI)
plot(road)
roadA = st_buffer(road, dist = 10000) %>% st_union()
roadB = st_buffer(road, dist = 25000) %>% st_union()
tmp = st_difference(roadB, roadA)
st_write(roadA, here("tmp","road10.gpkg"), append = FALSE)
st_write(tmp, here("tmp","road25N.gpkg"), append = FALSE)
rm(road, roadA, roadB, tmp)

#### populated areas raster ####

populated = raster(here("input","populated_areas","idn_bsgme_v0a_100m_2020.tif"))
populated = rast(populated)
plot(populated)
tmp = project(populated, crs("EPSG:23879"))

raster::writeRaster(tmp, 
                    filename=here("tmp", "popU.tif"), 
                    overwrite=TRUE)


popr = rast(here("tmp", "popu.tif"))
dim(tmp)
dim(populated)
# Crop the first raster to match the extent of the second raster
popr = crop(popr, kwh)

# align the resolution (rows and columns) to match r2
popr <- resample(popr, kwh)
plot(popr)
freq(popr)

#write to file
raster::writeRaster(popr, 
                    filename=here("tmp", "pop_resampled.tif"), 
                    overwrite=TRUE)

#### moutain areas ####
## read elevation map and remove NAs
indo_ele = read_delim(here("input","INDONESIA_altitude_data.csv"), delim = ",", locale = locale(decimal_mark = "."))
tmp = indo_ele
tmp = tmp %>% filter(!is.na(alt))

# set mountain limit
m = 1500
tmp = tmp %>% filter(alt > m)

# set as SF and local UTM 
tmp = st_as_sf(tmp,coords = c("x", "y"), crs = 4326)
tmp = st_transform(tmp, crs = 23879) 

# make polygons for mountains, unite and write to file
dist=500 # buffer each point 500m
buffer = st_buffer(tmp, dist = dist) # buffre
buffer = st_union(buffer) %>% st_cast("POLYGON") # merge
buffer = st_transform(buffer, crs = 4326) # transform to WGS84 for QGIS visuals
st_write(buffer, here("tmp","mountains.gpkg"),layer = "mountain", append = FALSE)
rm(tmp, buffer)

#### nature reserves ####
layer <- st_layers(here("input","land_cover", "nature4.gpkg"))
layer = layer %>% dplyr::select(name)
layer[1,]

nature1 = st_read(here("input","land_cover", "nature4.gpkg"), layer = layer[1,])
nature2 = st_read(here("input","land_cover", "nature4.gpkg"), layer = layer[2,])
nature3 = st_read(here("input","land_cover", "nature4.gpkg"), layer = layer[3,])

nature = rbind(nature1, nature2, nature3)

st_write(nature, here("output","nature.gpkg"), append = FALSE)

rm(layer, nature1, nature2, nature3)

#### AHP evaluation ####
# Clean enviromment
rm(list=ls(all=TRUE)) 
cat("\014")  

#### read data

# read Indonesia boundaries
boundary = st_read(here("input","indonesia_geo.gpkg"))
boundaryU = st_transform(boundary,  crs = 23879) #local UTM

# land cover
landuse = rast(here("tmp","landuse_class.tif"))
names(landuse) = "landuse"

# mountain
mountain = st_read(here("tmp","mountains.gpkg"))
mountainU = st_transform(mountain,  crs = 23879)

# nature
nature = read_sf(here("output", "nature.gpkg"))
natureU = st_transform(nature,  crs = 23879)

## populated area
# for points loop
#pop = read_sf(here("tmp", "pop_b500.gpkg"))
#popU = st_transform(pop,  crs = 23879)

# for raster
popU = rast(here("tmp","pop_resampled.tif"))
dim(popU)

## electric grid
#grid25 = read_sf(here("tmp", "grid25km_net.gpkg"), layer = "a")
#grid10 = read_sf(here("tmp", "grid10km.gpkg"), layer = "a")
grid = read_sf(here("input","grid.gpkg"))
grid = st_transform(grid, crs = 23879)

# kwh solar power
kwh = rast(here("tmp","idw2_class_mean.tif"))
names(kwh) = "kwhm2"
dim(kwh)

#### Filter exclusion areas

## test raster evaluation
plot(landuse)

## raster with landuse = 1
luone = landuse
luone[luone != 1] <- NA
kwh_lu = mask(kwh, luone)
#plot(kwh)
#plot(kwh_lu)

## mask with mountains
kwhlumo = mask(kwh_lu, vect(mountainU), inverse = TRUE)
#plot(kwhlumo)

## mask with nature
kwh_lumona = mask(kwhlumo, vect(natureU), inverse = TRUE)
#plot(kwh_lumona)

## mask with populated areas
popU
#plot(popU)

# Use `ifel` to set <0,5 (not populated) values to NA
popUpt5 <- ifel(popU >= 0.5, 1, NA)
#plot(popUpt5)

kwh_lumonapo = mask(kwh_lumona, popUpt5, inverse = TRUE)
#plot(kwh_lumona)
#plot(kwh_lumonapo)

writeRaster(kwh_lumonapo, filename=here("tmp", "kwh_all.tif"), overwrite=TRUE)

## cell raster
# Get the resolution (cell size) of the raster
res_raster <- res(kwh)

# Calculate the area of a raster cell
# In projected CRS, area = cell width * cell height
cell_area <- res_raster[1] * res_raster[2]
cell_area

#### calculate distance to nearest populated area

r1 = kwh_lumonapo

# want to be close to poplated areas to reduce transmission costs
popUpt8 <- ifel(popU >= 0.5, 1, NA) 
#target = popUpt5
target = popUpt8
target[!is.na(values(target))] <- 1  # Ensure target cells are 1
target[is.na(values(target))] <- NA # Non-target cells as NA
dist_to_target <- distance(target)
dist_kwh_to_pop <- mask(dist_to_target, r1)
#plot(dist_kwh_to_pop, main="Distance from kwh to pop")
#writeRaster(dist_kwh_to_pop, filename=here("tmp", "distance.tif"), overwrite=TRUE)

#### calculate distance to grid
# load the kwh
r = kwh_lumonapo

# Rasterize the line (cells touched by the line will have a value of 1, others NA)
line_raster <- rasterize(grid, r, field=1)

# Calculate distance from each cell to the nearest cell with the line
distance_raster <- distance(line_raster)
plot(distance_raster, main="Distance to Nearest Grid")
dist_kwh_to_grid <- mask(distance_raster, r)

#### make final matrix and write to file

# Convert kwh raster to polygons
polygon_kwh = as.polygons(kwh_lumonapo, dissolve=FALSE)  # Keep individual cells

# Convert to an sf object
sf_kwh = st_as_sf(polygon_kwh)
sf_kwh = sf_kwh %>% rename(kwh = kwhm2)

# Convert distance to population raster to sf object
polygon_pop = as.polygons(dist_kwh_to_pop, dissolve = FALSE)
sf_pop = st_as_sf(polygon_pop)
sf_pop = sf_pop %>% rename(pop_dist = idn_bsgme_v0a_100m_2020)

# Convert distance to grid raster to sf object
polygon_grid = as.polygons(dist_kwh_to_grid, dissolve = FALSE)
sf_grid = st_as_sf(polygon_grid)
sf_grid = sf_grid %>% rename(grid_dist = layer)

# merge into one sf object
matrix1 = cbind(sf_kwh, sf_grid, sf_pop) %>% 
  dplyr::select(-c(geometry.1, geometry.2)) %>%
  mutate(row_id = row_number())

st_write(matrix1, here("tmp","matrix1.gpkg"),append = FALSE)
matrix1 = read_sf(here("tmp","matrix1.gpkg"))

matrixX = matrix1

## apply quintiles for matrix evaluation

# kwh quintile
q_kwh = quantile(matrixX$kwh, probs = c(0.2, 0.4, 0.6, 0.8, 1))
q_kwh
q1 = q_kwh[1]
q2 = q_kwh[2]
q3 = q_kwh[3]
q4 = q_kwh[4]
q5 = q_kwh[5]

matrixX = matrixX %>%
  mutate(g_kwh = case_when(
    kwh <= q1 ~ 1,
    kwh <= q2 ~ 2,
    kwh <= q3 ~ 3,
    kwh <= q4 ~ 4,
    kwh <= q5 ~ 5,
    TRUE ~ 99 ))

# grid quintile
q_grid = quantile(matrixX$grid_dist, probs = c(0.2, 0.4, 0.6, 0.8, 1))
q_grid
q1 = q_grid[1]
q2 = q_grid[2]
q3 = q_grid[3]
q4 = q_grid[4]
q5 = q_grid[5]

matrixX = matrixX %>%
  mutate(g_grid = case_when(
    grid_dist <= q1 ~ 5,
    grid_dist <= q2 ~ 4,
    grid_dist <= q3 ~ 3,
    grid_dist <= q4 ~ 2,
    grid_dist <= q5 ~ 1,
    TRUE ~ 99 ))

# population quintiles
q_pop = quantile(matrixX$pop_dist, probs = c(0.2, 0.4, 0.6, 0.8, 1))
q_pop
q1 = q_pop[1]
q2 = q_pop[2]
q3 = q_pop[3]
q4 = q_pop[4]
q5 = q_pop[5]

matrixX = matrixX %>%
  mutate(g_pop = case_when(
    pop_dist <= q1 ~ 5,
    pop_dist <= q2 ~ 4,
    pop_dist <= q3 ~ 3,
    pop_dist <= q4 ~ 2,
    pop_dist <= q5 ~ 1,
    TRUE ~ 99 ))

# apply weights
w_kwh = 0.65
w_grid = 0.23
w_pop = 0.12

# calc_weights and sum
matrixX = matrixX %>% rowwise() %>%
  mutate(we_kwh = g_kwh * w_kwh,
         we_grid = g_grid * w_grid,
         we_pop = g_pop * w_pop,
         we_sum = sum(we_kwh, we_grid, we_pop))

# calculate twh per cell / site
sites_per_cell = 1
effective_hours = 24 * 0.24
m2_site = 12000000
days_yr = 365
kwh_to_twh = 1000000000

matrixX = matrixX %>% mutate(
  twh_yr = kwh * sites_per_cell * effective_hours * m2_site * days_yr / kwh_to_twh)

## sort the matrix ala suitability criteria
matrixX = matrixX %>% arrange(desc(we_sum))
st_write(matrixX, here("tmp","matrixX.gpkg"), append = FALSE)

## slice the matrix to find needed twh production + write to file
tmp = matrixX$twh_yr
tmp2 = as.data.frame(cumsum(tmp)) %>% rename(twh_cum = "cumsum(tmp)" )
matrixX = cbind(matrixX, tmp2)

twh_needed = 127
tmp = matrixX %>% filter(twh_cum <= twh_needed)
a = nrow(tmp)
matrix2 = slice(matrixX, 0:a+1)
st_write(matrix2, here("tmp","matrix2.gpkg"), append = FALSE)

## plot 
ggplot() +
  geom_sf(data = boundaryU) +
  geom_sf(data = matrix2, aes(fill = "kwh" )) +
  coord_sf(crs = 4326) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)

ggsave(here("tmp","test.png"), width = 6, height = 6, dpi = "screen")


kwh <- as.data.frame(kwh_lumonapo, xy = TRUE)  # Include x and y coordinates

ggplot() +
  geom_sf(data = boundaryU) +
  geom_tile(data = kwh, aes(x = x, y = y, fill = kwhm2)) +
  #geom_sf(data = grid) +
  scale_fill_viridis_c() +  # Optional: Viridis color scale for better visualization
  #coord_equal() +           # Ensures square raster cells
  labs(title = "kwH values outside exclusion areas", fill = "Value") +
  #coord_sf(crs = 4326) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)


#### NPV analysis ####
# Clean enviromment
rm(list=ls(all=TRUE)) 
cat("\014")  

## read data
matrix2 = read_sf(here("tmp","matrix2.gpkg"))

#### set parameters

## revenue
# Use a value of 7.67p/kWh (£0.0767/kWh) 
revenue_kwh = 7.67 / 100

## solar system capacity cost
# $30 to $80 per m² is a reasonable estimate for CAPEX for large-scale solar projects.
CpEx_capacity_mw = 1160000 * 0.8 #usd 50 converted to GBP

## network connection cost
# $30,000 to $50,000 per km for high-voltage lines
# $1 million to $5 million for building a new substation, 
# depending on the required voltage and capacity.

CpEx_nc_km = 590 * 0.8 #connection per meter converted to GBP
#nc_fixed = 2500000 * 0.8 # fixed connection costs per site

## grid distance
# include new grid distance, set to 5000 meters if original value = 0
matrix3 = matrix2 %>% 
  mutate(grid_dist2 = ifelse(grid_dist == 0, 5000, grid_dist))

## capacity factor
cap_factor = 0.24

## system lifetime
lifetime_yrs = 25

rm(matrix2)

#### calculate capex (initial costs)

## calculate mw capacity per site and mw_km per site
matrix3 = matrix3 %>% 
  mutate(mw_capacity = twh_yr * 1000 / cap_factor / 365 / 24 * 1000,
         mw_km = grid_dist2 * mw_capacity / 1000)

## calculate capacity capex
matrix3 = matrix3 %>% 
  mutate(capex = (mw_capacity * CpEx_capacity_mw) + (mw_km * CpEx_nc_km))

#### calculate yearly revenues
matrix3 = matrix3 %>%
  mutate(rev_year = twh_yr * 10^9 * revenue_kwh)

tmp = matrix3 %>% st_drop_geometry()
write_xlsx(tmp, here("tmp","matrix3.xlsx"))

#### Net Present Value
#in this project, we do not consider opex but only capex
calc_NPV <- function(annual_revenue, i=0.05, lifetime_yrs, CAPEX, OPEX=0){
  revenue <- rep(annual_revenue, lifetime_yrs) 
  t <- seq(1, lifetime_yrs, 1) #output: 1, 2, 3, ...25
  
  NPV <- sum( (revenue )/(1 + i)**t ) - CAPEX
  return(round(NPV, 0))
}

annual_revenue = sum(matrix3$rev_year)
CAPEX = sum(matrix3$capex)

## calculate the NPV

npv=calc_NPV(annual_revenue = annual_revenue,
             lifetime_yrs=lifetime_yrs,
             CAPEX=CAPEX)
ifelse(npv>0, "Support","obeject" )

## Levelized cost of electricity (LCOE)

kwh_year = sum(matrix3$twh_yr * 1000000000)

Life_span_generation_kWH <- function (yearly_generation_kWH, discount = 0.08, lifetime_yrs = 25){
  t<- seq(1, lifetime_yrs, 1)
  L_S_G <- sum(yearly_generation_kWH/(1+discount)**t)
  return (round(L_S_G,0))
}

L_S_G = Life_span_generation_kWH(yearly_generation_kWH = kwh_year)

# NPV of cost. if you don't consider the operational cost, 
# you can just use CAPEX as proxy for NPV of cost

NPV_cost = sum(matrix3$capex)

LCOE <- function(NPV_cost,Life_span_generation){
  lcoe <- NPV_cost/Life_span_generation
  return(round(lcoe,4))
}

lcoe = LCOE(NPV_cost, L_S_G)
print(lcoe)


