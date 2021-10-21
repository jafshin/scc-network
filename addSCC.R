# Add Strategic Cycling Corridor links to matsim network

# Approach - 
# 1 Filter network links to exclude pt and non-cyclable paths, and filter network
#   nodes to those connected by the cyclable links
# 2 For each cycle path in Strategic Cycle network, buffer cycle paths to 250m,
#   creating constrained network within which links can be located
# 3 Locate nearest nodes to start and end points of cycle path
# 4 Within the constrained network, find shortest network path connecting nodes,
#   in each direction
# 5 If a route can't be found within the 250m buffer, progressively increase buffer
#   up to maximum of 1.5km
# 6 Write the SCC id (row number) and type (C1 or C2) to the links file
# 7 Will not find routes where:
#   (a) there are no network links at all within 250m (eg cycle paths outside study region),
#   (b) route cannot be found within 1.5km buffer (eg disconnected network), or
#   (c) cycle path is so short that start and end points are located near same network node.


# set inputs
# networkFile <- "./data/MATSimMelbNetwork_SRL.sqlite"  # old network with SRL
# networkFile <- "./data/roadDataAllNoPT_sep24.sqlite"  # new unsimplified network (without PT or SRL)
networkFile <- "./data/network.sqlite"  # new simplified network
linkLayer <- "links"
nodeLayer <- "nodes"
# SCCFile <- "./data/2020_Strategic_Cycling_Corridors_(SCC)/2020_Strategic_Cycling_Corridors_(SCC).shp"
SCCzip <- "./data/2020_Strategic_Cycling_Corridors_(SCC).zip"

# set up environment
library(dplyr)
library(sf)
library(igraph)
library(lwgeom)
library(stringr)
library(nngeo)  # for nn (nearest)

# constants
# MinLength <- 500

source("./functions/getPathLinks.R")

# 1. Load inputs - links, nodes and SCC
# -----------------------------------------------------------------------------
# links and nodes
links <- st_read(networkFile, layer = linkLayer) %>%
  st_make_valid() %>%  # note - without this there are some invalid links (duplicated nodes)
  mutate(link_id = row_number())
nodes <- st_read(networkFile, layer = nodeLayer)

# SCC <- st_read(SCCFile) %>%
# # set CRS to be same as links
# st_transform(st_crs(links)) %>%
#   mutate(scc_id = row_number())

# SCC  https://stackoverflow.com/questions/59740419/unzipping-and-reading-shape-file-in-r-without-rgdal-installed
temp <- tempfile()
unzip(zipfile = SCCzip, exdir = temp)
SCCshp1 <- list.files(temp, pattern = ".shp", full.names = TRUE)
SCCshp2 <- list.files(temp, pattern = ".xml", full.names = TRUE)
SCCshp <- setdiff(SCCshp1, SCCshp2)  # file that contains .shp but not .xml (https://stackoverflow.com/questions/31590730/list-files-in-r-that-do-not-match-a-pattern)
SCC <- st_read(SCCshp) %>%
  # set CRS to be same as links
  st_transform(st_crs(links)) %>%
  mutate(scc_id = row_number()) %>%
  st_snap_to_grid(1)


# 2. Set up cyclable links and nodes; divide SCC into workable segments
# -----------------------------------------------------------------------------
# filter links - remove PT and motorway; add link_id
# note - can't limit to all where is_cycle==1 - some 'non-cyclable' links on main roads required for connectivity 
cyclable.links <- links %>%
  filter(!highway %in% c("pt", "motorway", "motorway_link")) %>%
  # create weight (for shortest path), with cyclepaths weighted @ [80%] and non-cyclable @150%
  # for 'simplified'
  mutate(weight = ifelse(!is.na(cycleway), length*0.85,
                         ifelse(is_cycle == 0, length*1.5,
                                length)))
  # for 'unsimplfied'
  # mutate(weight = ifelse(!is.na(cycleway), road_length_meter*0.9, 
  #                        ifelse(!str_detect(modes, "bike"), road_length_meter*1.5,
  #                              road_length_meter)))




# filter nodes to those used in links, to remove any disconnected
cyclable.nodes <- nodes %>%
  filter(id %in% cyclable.links$from_id | id %in% cyclable.links$to_id)

# ### trying to fix SCC
# ## stuff that doesn't work
# test <- st_linesubstring(SCC, 0, 500) %>%
#   mutate(testlength = st_length(geometry))
# 
# test <- st_segmentize(SCC, 500) %>%
#   st_cast(to="LINESTRING")
# 
# SCC.splits <- SCC %>%
#   st_line_sample(density = units::set_units(500, m)) %>%
#   st_sf() %>%
#   st_snap_to_grid(1)
# 
# split.SCC <- st_split(SCC, SCC.splits)
# 
# st_write(SCC.splits, "testsccsplits.sqlite", delete_dsn = TRUE)
# st_write(split.SCC, "testsplitscc.sqlite", delete_dsn = TRUE)


## intro to  '500m' approach
# shortSCC <- SCC %>%
#   filter(as.numeric(st_length(geometry)) <= MinLength)
# longSCC <- SCC %>%
#   filter(as.numeric(st_length(geometry)) > MinLength)
# 
# newSCC <- shortSCC

## approach that shortens some links, but still leaves some long
# for (i in 1:nrow(longSCC)) {
#   # split into 500m segments - doesn't work yet - still has some segments > 500m
#   # also, takes a couple of minutes
#   path <- longSCC[i,] %>%
#     # st_snap_to_grid(1) %>%
#     st_segmentize(., MinLength)
#   splitpoints1 <- path %>%
#     st_line_sample(density = units::set_units(MinLength, m)) %>%
#     st_sf() 
#   splitpoints2 <- splitpoints1 %>%
#     # st_snap_to_grid(1) %>%
#     st_snap(., path, 50) %>% # doesn't succeed in snapping every point.
#     st_cast() 
#   splitpath <- st_split(path, splitpoints2) %>%
#     # st_sf() %>%
#     st_cast() #%>%
#     # mutate(xlength = as.numeric(st_length(geometry)))
#   
#   newSCC <- rbind(newSCC, splitpath)
# }

## maybe...(splitting at intersections?)
#
newSCC <- NULL
for (i in 1:nrow(SCC)) {
  
  path <- SCC[i,]
  cat("Splitting path #", i, "of", nrow(SCC), "at intersections\n")
  intersections <- st_intersection(path, SCC)
  if (nrow(intersections) > 1) {
    splitpoints <- st_collection_extract(intersections, type = "POINT", warn = FALSE)
    path <- st_split(path, splitpoints) %>%
      st_cast()
  }

  # mutate(xlength = as.numeric(st_length(geometry)))
  
  newSCC <- rbind(newSCC, path)
}

## testing - easy load of newSCC
newSCC <- st_read("./testnewscc.sqlite")

# # testing outputs...
# newSCC <- newSCC %>%
#   mutate(xlength = as.numeric(st_length(geometry)))
# 
# st_write(splitpath, "testsplitpath.sqlite", delete_dsn = TRUE)
# st_write(splitpoints, "testsplitpoints.sqlite", delete_dsn = TRUE)
# st_write(splitpoints2, "testsplitpoints2.sqlite", delete_dsn = TRUE)
# st_write(newSCC, "testnewSCC.sqlite", delete_dsn = TRUE)





# 4. Run function to find links corresponding to each SCC path
# -----------------------------------------------------------------------------
## clearing 'links' for testing
links <- links %>%
  mutate(scc_id = NULL, scc_type = NULL)

## repeated here for easy testing
source("./functions/getPathLinks.R")

## error reporting (didn't work - tried to write to it - can't write from function)
# errorFile <- data.frame(comment = character())

for (i in 1:nrow(SCC)) {
# for (i in c(71:80)) {
  paths <- newSCC %>%  ### NEED A BETTER NAME FOR newSCC
    filter(scc_id == i)
  startpoint <- lwgeom::st_startpoint(SCC[i,])  ## used to try to paths
  cat("Finding links for cycle path #", i, "of", nrow(SCC), "(", nrow(paths), "segments )\n")
  # nextVertices <- c(0, 0)
  fromPoints <- c("0")
  toPoints <- c("0")
  
  for (j in 1:nrow(paths)) {
  # for (j in c(1)) {
    path <- paths[j,]
    bufferDistance <- 250

    # get the links for each direction
    # path.links <- getPathLinks(path, bufferDistance, nextVertices)
    path.link.outputs <- getPathLinks(path, bufferDistance, startpoint, fromPoints, toPoints)
    # path.links <- getPathLinks(path, bufferDistance)
    
    # # if there are no links in either direction, try increased bufferDistance, up to 2km
    # while (length(path.links[[1]]) == 0 | length(path.links[[2]]) == 0) {
    #   bufferDistance <- bufferDistance + 50
    #   if (bufferDistance > 2000) break
    #   path.links <- getLinks(path, bufferDistance)
    # }
    
    # combine the paths in each direction, and remove duplicates
    path.links <- c(path.link.outputs[[1]], path.link.outputs[[2]]) %>% 
      unique()
    
    # add SCC cycle path id and type to links
    for (k in 1:length(path.links)) {
      links[links$link_id == path.links[k], "scc_id"] <- st_drop_geometry(SCC[i, "scc_id"])
      links[links$link_id == path.links[k], "scc_type"] <- st_drop_geometry(SCC[i, "TYPE"])
    }
    
    # replace fromPoints and toPoints ready for next loop
    fromPoints <- path.link.outputs[[3]]
    toPoints <- path.link.outputs[[4]]
    
  }
  
}

# 5. Write outputs
# -----------------------------------------------------------------------------

st_write(links, "outputlinks.sqlite", delete_dsn = TRUE)
st_write(SCC, "outputscc.sqlite", delete_dsn = TRUE)

## testing
filtered.links <- links %>%
  filter(!is.na(scc_id))
st_write(filtered.links, "outputlinksa.sqlite", delete_dsn = TRUE)

# 6. Testing/display tools
# -----------------------------------------------------------------------------
library(ggplot2)
library(ggspatial)

ggplot() + geom_sf(data = cyclable.links)
ggplot() + annotation_map_tile(type="osm",zoom=11, alpha=0.6) +
  geom_sf(data = path)
ggplot() + annotation_map_tile(type="osm",zoom=11, alpha=0.6) +
  geom_sf(data = path.buffered)
ggplot() + annotation_map_tile(type="osm",zoom=11, alpha=0.6) +
  geom_sf(data = local.links)
ggplot() + annotation_map_tile(type="osm",zoom=11, alpha=0.6) +
  geom_sf(data = local.links) +
  geom_sf(data = point1, colour = "red") +
  geom_sf(data = point2, colour = "blue")
ggplot() + annotation_map_tile(type="osm",zoom=11, alpha=0.6) +
  geom_sf(data = local.links) +
  geom_sf(data = point1, colour = "red") +
  geom_sf(data = point2, colour = "blue") + 
  geom_sf(data = reachable.area, alpha = 0.3) +
  geom_sf(data = path.intersecs, colour = "green")
ggplot() + annotation_map_tile(type="osm",zoom=9, alpha=0.6) +
  geom_sf(data = SCC[SCC$scc_id == 8,])
ggplot() + annotation_map_tile(type="osm", zoom=11, alpha=0.6) +
  geom_sf(data = reachable.area)

st_write(local.links, "testoutputlinks.sqlite", delete_dsn = TRUE)
