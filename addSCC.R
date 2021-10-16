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
networkFile <- "./data/MATSimMelbNetwork_SRL.sqlite"
linkLayer <- "links"
nodeLayer <- "nodes"
# SCCFile <- "./data/2020_Strategic_Cycling_Corridors_(SCC)/2020_Strategic_Cycling_Corridors_(SCC).shp"
SCCzip <- "./data/2020_Strategic_Cycling_Corridors_(SCC).zip"

# set up environment
library(dplyr)
library(sf)
library(igraph)
library(lwgeom)


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
  mutate(scc_id = row_number())


# 2. Set up cyclable links and nodes
# -----------------------------------------------------------------------------
# filter links - remove non-cyclable and PT; add link_id
cyclable.links <- links %>%
  filter(is_cycle == 1 & highway != "pt")

# filter nodes to those used in links, to remove any disconnected
cyclable.nodes <- nodes %>%
  filter(id %in% cyclable.links$from_id | id %in% cyclable.links$to_id)


# 3. Main function - to find network links comprising shortest path for given cycle path
# -----------------------------------------------------------------------------
getLinks <- function(path, bufferDistance) {
  # create empty vectors for outputs to be returned
  path.links1 <- c()
  path.links2 <- c()
  
  # buffer cycle path to given distance (constrains to network links near cycle path route)
  path.buffered <- st_buffer(path, bufferDistance)
  
  # identify network links within the buffer
  local.links <- cyclable.links %>% 
    filter(st_intersects(GEOMETRY, path.buffered, sparse = FALSE))
  
  # if there are local links, proceed; if not (ie paths outside study region) then return empty vectors 
  if (length(local.links) > 0) {
    # create graph to find shortest paths in network links for the given cycle path
    # for two-way links, create an extra set of links with 'from' and 'to' reversed
    links2way <- local.links %>%
      st_drop_geometry() %>%
      filter(is_oneway == 0) %>%
      mutate(from = to_id, to = from_id) %>%  # reverse order
      dplyr::select(from, to, weight = length, link_id)
    
    # combine links with the extra set of reversed links for two-way links
    directedlinks <- local.links %>%
      st_drop_geometry() %>%
      dplyr::select(from = from_id, to = to_id, weight = length, link_id) %>%
      rbind(., links2way)
    
    # create directed graph
    graph <- graph_from_data_frame(directedlinks, directed = T, vertices = cyclable.nodes)
    
    
    # get id's of nodes nearest to cycle path start and end points
    point1 <- lwgeom::st_startpoint(path)
    point2 <- lwgeom::st_endpoint(path)
    node1 <- cyclable.nodes[st_nearest_feature(point1, cyclable.nodes), "id"][[1]]
    node2 <- cyclable.nodes[st_nearest_feature(point2, cyclable.nodes), "id"][[1]]
    
    # if nodes are different, proceed; if not (ie start and end are same) then return empty vectors
    if (node1 != node2) {
      # find shortest path links for cycle path - first direction
      # find shortest path between nodes, and return edges ('epath')
      shortest1 <- shortest_paths(graph, 
                                  from = as.character(node1), 
                                  to = as.character(node2), 
                                  output = c("both"))
      
      # make vector of edges in the epath
      if (length(shortest1$epath[[1]] > 0)) {  # length will be zero if path can't be found
        for (j in 1:length(shortest1$epath[[1]])) {
          # get the link_id for the link
          link <- directedlinks %>%
            filter(link_id == shortest1$epath[[1]][j][[1]]$link_id)  # this is required to return the link_id!!
          link_id <- link$link_id
          # add the link to the path.links vector
          path.links1 <- c(path.links1, link_id)
        }
        # remove duplicates (which arise from duplication for two-way links)
        path.links1 <- unique(path.links1)
      }
      
      
      # find shortest path links for cycle path - second direction
      # get id's of nodes nearest to cycle path start and end points
      
      # find paths between nodes, and return edges ('epath')
      shortest2 <- shortest_paths(graph, 
                                  from = as.character(node2), 
                                  to = as.character(node1), 
                                  output = c("both"))
      
      # make vector of edges in the epath
      if (length(shortest2$epath[[1]] > 0)) {  # length will be zero if path can't be found
        for (j in 1:length(shortest2$epath[[1]])) {
          # get the link_id for the link
          link <- directedlinks %>%
            filter(link_id == shortest2$epath[[1]][j][[1]]$link_id)  # this is required to return the link_id!!
          link_id <- link$link_id
          # add the link to the path.links vector
          path.links2 <- c(path.links2, link_id)
        }
        # remove duplicates (which arise from duplication for two-way links)
        path.links2 <- unique(path.links2)
      }
    }
  }

    # return the output as a list
  return(list(path.links1, path.links2))
}


# 4. Run function to find links corresponding to each SCC path
# -----------------------------------------------------------------------------
for (i in 1:nrow(SCC)) {
  path <- SCC[i,]
  bufferDistance <- 250
  cat("Finding links for cycle path #", i, "of", nrow(SCC), "\n")
  
  # get the links for each direction
  path.links <- getLinks(path, bufferDistance)
  
  # if there are no links in either direction, try increased bufferDistance, up to 2km
  while (length(path.links[[1]]) == 0 | length(path.links[[2]]) == 0) {
    bufferDistance <- bufferDistance + 50
    if (bufferDistance > 2000) break
    path.links <- getLinks(path, bufferDistance)
  }
  
  # combine the paths in each direction, and remove duplicates
  path.links <- c(path.links[[1]], path.links[[2]]) %>% 
    unique()
  
  # add SCC cycle path id and type to links
  for (k in 1:length(path.links)) {
    links[links$link_id == path.links[k], "scc_id"] <- st_drop_geometry(SCC[i, "scc_id"])
    links[links$link_id == path.links[k], "scc_type"] <- st_drop_geometry(SCC[i, "TYPE"])
  }
}

# 5. Write outputs
# -----------------------------------------------------------------------------

st_write(links, "outputlinks.sqlite", delete_dsn = TRUE)
st_write(SCC, "outputscc.sqlite", delete_dsn = TRUE)


# 6. Testing/display tools
# -----------------------------------------------------------------------------
library(ggplot2)
library(ggspatial)

ggplot() + geom_sf(data = cyclable.links)
ggplot() + annotation_map_tile(type="osm",zoom=11, alpha=0.6) +
  geom_sf(data = path.buffered)
ggplot() + annotation_map_tile(type="osm",zoom=11, alpha=0.6) +
  geom_sf(data = local.links)
ggplot() + annotation_map_tile(type="osm",zoom=9, alpha=0.6) +
  geom_sf(data = SCC[SCC$scc_id == 8,])

