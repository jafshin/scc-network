# Add Strategic Cycling Corridor links to matsim network

# Approach - 
# 1 Filter network links to exclude pt and freeways paths, and filter network
#   nodes to those connected by the cyclable links
# 2 Split each Strategic Cycle Corridors (SCC) path into segments at intersections 
#   (where SCC paths intersect – not road intersections)
# 3 For each SC segment, buffer cycle paths to 250m [parameter], 
#   creating constrained network within which links can be located
# 3 Locate nearest nodes to start and end points of cycle path (use 3 nodes - closest may be unsuitable,
#   for example because it connects only to one-way links in the wrong direction)
# 4 Within the constrained network, find shortest network path connecting nodes,
#   in each direction [getPathLinks function]
# 5 If a route can't be found within the 250m buffer, then find as much (if any) of a route 
#   as possible from/to each end (in the correct direction) [getPartialPathLinks function]
# 6 Where there is a remaining gap to be filled, create a new link to join the 
#   nodes at each end of the gap [bridgeGaps function]
# 7 If the gap exceeds 250m [parameter], break the unmatched part of the cycle path into 250m sections, and
#   create new links matching the 250m sections (each connecting the nearest nodes to the path) 
# 8 Where 2nd or 3rd closest node was used as start or end point of route, see whether a direct link
#   exists (in the correct direction) to a closer node
# 9 For each link corresponding to a cycle path (whether matched or new), add attributes for
#   the SCC id (row number) and type (C1 or C2)
# 10 In addition, if the link is new, then complete the SCC_new field with '1', and complete
#   its length field (corresponding to length of the equivalent path section)


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

# parameters
bufferDistance <- 100  # distance to buffer paths, creating constrained network to search for matching links
gapSegmentDistance <- 250  # distance (m) to break up unmatched cycle path segments to create new links
cyclableBonus <- 0.85  # input to cost function for shortest path weight: bonus for cyclable routes
nonCyclablePenalty <- 1.5  # equivalent penalty for non-cyclable (sometimes needed for connectivity)

source("./functions/getPathLinks.R")
source("./functions/getStartEndLinks.R")


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
  # cost function to create weight (for shortest path)
  # for 'simplified'
  mutate(weight = ifelse(!is.na(cycleway), length * cyclableBonus,
                         ifelse(is_cycle == 0, length * nonCyclablePenalty,
                                length)))
  # for 'unsimplfied'
  # mutate(weight = ifelse(!is.na(cycleway), road_length_meter * cyclableBonus, 
  #                        ifelse(!str_detect(modes, "bike"), road_length_meter * nonCyclablePenalty,
  #                              road_length_meter)))


# filter nodes to those used in links, to remove any disconnected
cyclable.nodes <- nodes %>%
  filter(id %in% cyclable.links$from_id | id %in% cyclable.links$to_id)

# split SCC at intersections
splitSCC <- NULL
for (i in 1:nrow(SCC)) {
  path <- SCC[i,]
  cat("Splitting path #", i, "of", nrow(SCC), "at intersections\n")
  intersections <- st_intersection(path, SCC)
  if (nrow(intersections) > 1) {
    splitpoints <- st_collection_extract(intersections, type = "POINT", warn = FALSE)
    path <- st_split(path, splitpoints) %>%
      st_cast()
  }
  splitSCC <- rbind(splitSCC, path)
}

# ## testing - easy load of splitSCC
# splitSCC <- st_read("./outputscc.sqlite")

# confine SCC to intersection with links bounding box
studyAreaSCC <- st_intersection(splitSCC, 
                                st_as_sfc(st_bbox(cyclable.links)))


# 4. Run function to find links corresponding to each SCC path
# -----------------------------------------------------------------------------
# note - 'links' needs to be re-loaded before running this section
# ## testing - clearing links, and re-loading function
# links <- st_read(networkFile, layer = linkLayer) %>%
#   st_make_valid() %>%  # note - without this there are some invalid links (duplicated nodes)
#   mutate(link_id = row_number())
# # # 
# # # links <- links %>%
# # #   mutate(scc_id = NULL, scc_type = NULL)
# # # 
# source("./functions/getPathLinks.R")
# source("./functions/getStartEndLinks.R")

# create element to hold new links returned from each path
accumulated.new.links <- c()

for (i in 1:nrow(SCC)) {
# for (i in c(1:10)) {
  paths <- studyAreaSCC %>%
    filter(scc_id == i)
  if (nrow(paths) > 0) {
    startpoint <- lwgeom::st_startpoint(SCC[i,])  # used to align segments of path end-to-end, and for start/end links
    endpoint <- lwgeom::st_endpoint(SCC[i,])  # used for start/end links
    cat("Finding links for cycle path #", i, "of", nrow(SCC), "(", nrow(paths), "segments )\n")
    firstDirFromPoints <- c("0")
    firstDirToPoints <- c("0")
    secondDirFromPoints <- c("0")
    secondDirToPoints <- c("0")
    
    for (j in 1:nrow(paths)) {
      # for (j in c(1)) {
      path <- paths[j,]
      
      # get the links for each direction: returns (1 & 2) matched links in each direction,
      # (3) new links, (4 - 5) start and end points (used to align path segments)
      path.link.outputs <- getPathLinks(path, startpoint, 
                                        firstDirFromPoints, firstDirToPoints, 
                                        secondDirFromPoints, secondDirToPoints)
      
      # combine the links in each direction, and remove duplicates
      path.links <- c(path.link.outputs[[1]], path.link.outputs[[2]]) %>% 
        unique()
      
      # add SCC cycle path id and type to links
      for (k in 1:length(path.links)) {
        links[links$link_id == path.links[k], "scc_id"] <- st_drop_geometry(SCC[i, "scc_id"])
        links[links$link_id == path.links[k], "scc_type"] <- st_drop_geometry(SCC[i, "TYPE"])
      }
      
      # add new links to accumulated new links
      new.links <- path.link.outputs[[3]]
      if (length(new.links) > 0) {
        accumulated.new.links <- bind_rows(accumulated.new.links, new.links)
      }
      
      # replace fromPoints and toPoints ready for next loop
      firstDirFromPoints <- path.link.outputs[[4]]
      firstDirToPoints <- path.link.outputs[[5]]
      secondDirFromPoints <- path.link.outputs[[6]]
      secondDirToPoints <- path.link.outputs[[7]]
    }
    
    # see whether further links can be added at start and end (selected links will start and
    # end at one of closest 3 nodes - if 2nd or 3rd has been selected, this may be able to link to a closer)
    start.end.links <- c(0)
    if (length(firstDirFromPoints) > 1 | length(firstDirToPoints) > 1) {
      first.dir.start.end.links <- getStartEndLinks(startpoint, endpoint,
                                                    firstDirFromPoints, firstDirToPoints)
    }
    if (length(secondDirFromPoints) > 1 | length(secondDirToPoints) > 1) {
      second.dir.start.end.links <- getStartEndLinks(endpoint, startpoint,
                                                     secondDirFromPoints, secondDirToPoints)
    }

    start.end.links <- c(start.end.links, first.dir.start.end.links, second.dir.start.end.links) %>%
      unique()  # remove duplicates (could arise because in both directions)
    for (j in 1:length(start.end.links)) {  
      links[links$link_id == start.end.links[j], "scc_id"] <- st_drop_geometry(SCC[i, "scc_id"])
      links[links$link_id == start.end.links[j], "scc_type"] <- st_drop_geometry(SCC[i, "TYPE"])
    }
  
  }
}

# add all accumulated new links to links
links <- bind_rows(links, accumulated.new.links)


# 5. Write outputs
# -----------------------------------------------------------------------------
st_write(links, "outputlinks100_250se.sqlite", delete_dsn = TRUE)

st_write(splitSCC, "outputscc.sqlite", delete_dsn = TRUE)

# ## testing
# filtered.links <- links %>%
#   filter(!is.na(scc_id))
#   # filter(scc_id == 4)
# st_write(filtered.links, "outputlinksa.sqlite", delete_dsn = TRUE)
# # st_write(SCC, "outputscc.sqlite", delete_dsn = TRUE)
# st_write(accumulated.new.links, "testaccumulated.sqlite", delete_dsn = TRUE)
# 
# ## testing - 'lost' sections of new links
# View(accumulated.new.links)
# lost <- accumulated.new.links %>% filter(from_id == to_id)
# View(lost)
# unique(lost$scc_id)
# sum(lost$length)


# 6. Testing/display tools [not needed to run script]
# -----------------------------------------------------------------------------
library(ggplot2)
library(ggspatial)

ggplot() + geom_sf(data = cyclable.links)
ggplot() + annotation_map_tile(type="osm",zoom=11, alpha=0.6) +
  geom_sf(data = path)
ggplot() + annotation_map_tile(type="osm",zoom=11, alpha=0.6) +
  geom_sf(data = paths) + 
  geom_sf(data = startpoint, colour = "red") +
  geom_sf(data = endpoint, colour = "blue")
ggplot() + annotation_map_tile(type="osm",zoom=11, alpha=0.6) +
  geom_sf(data = path.buffered, alpha=0.6) + 
  geom_sf(data = path, colour = "blue") +
  geom_sf(data = local.links)
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
  geom_sf(data = path.intersecs, colour = "green") +
  geom_sf(data = path, colour = "green")
ggplot() + annotation_map_tile(type="osm",zoom=11, alpha=0.6) +
  geom_sf(data = local.links) +
  geom_sf(data = point1, colour = "red") +
  geom_sf(data = point2, colour = "blue") + 
  geom_sf(data = reachable.nodes)
ggplot() + annotation_map_tile(type="osm",zoom=9, alpha=0.6) +
  geom_sf(data = SCC[SCC$scc_id == 8,])
ggplot() + annotation_map_tile(type="osm", zoom=11, alpha=0.6) +
  geom_sf(data = reachable.area)
ggplot() + annotation_map_tile(type="osm", zoom = 11, alpha=0.6) +
  geom_sf(data = pathPoints) +
  # geom_sf(data = pathSectionToBridge, colour = "red") + 
  geom_sf(data = nodeA, colour = "blue") +
  geom_sf(data = nodeB, colour = "blue")
ggplot() + annotation_map_tile(type="osm", zoom = 11, alpha=0.6) +
  geom_sf(data = pathPoints) +
  geom_sf(data = pathSectionToBridge, colour = "red") +
  geom_sf(data = newlink, colour = "blue")
ggplot() + geom_sf(data = splitSCC, colour = "blue") +
  geom_sf(data = studyAreaSCC, colour = "red")

st_write(local.links, "testoutputlinks.sqlite", delete_dsn = TRUE)
st_write(pathPoints, "testpathpoints.sqlite", delete_dsn = TRUE)
st_write(pathSectionToBridge, "testoutputbridge.sqlite", delete_dsn = TRUE)
