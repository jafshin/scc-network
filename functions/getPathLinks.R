# Function - to find network links comprising shortest path for given cycle path segment

# getPathLinks <- function(path, bufferDistance, nextVertices) {
getPathLinks <- function(path, startpoint, fromPoints, toPoints) {
  
  source("./functions/getPartialPathLinks.R")
  source("./functions/bridgeGaps.R")
  
  cat("Finding links for cycle path #", path$scc_id, "segment", j, "\n")
  
  # create empty vectors for outputs to be returned
  path.links1 <- c()  # for links in one direction
  path.links2 <- c()  # for links in the other direction
  new.links <- c()  # for new links added to bridge gaps
  
  # buffer cycle path to given distance (constrains to network links near cycle path route)
  path.buffered <- st_buffer(path, bufferDistance)
  
  # identify network links within the buffer
  local.links <- cyclable.links %>% 
    filter(st_intersects(GEOMETRY, path.buffered, sparse = FALSE))
  
  # if there are local links, proceed; if not (ie paths outside study region) then return empty vectors 
  if (nrow(local.links) > 0) {
    
    # for simplified (where there is 'is_oneway' column) -
    # create graph to find shortest paths in network links for the given cycle path
    # for two-way links, create an extra set of links with 'from' and 'to' reversed
    links2way <- local.links %>%
      st_drop_geometry() %>%
      filter(is_oneway == 0) %>%
      mutate(from = to_id, to = from_id) %>%  # reverse order
      dplyr::select(from, to, weight, link_id)  # 'weight' is from 'length' using cost function in 'addSCC.R'
    
    # combine links with the extra set of reversed links for two-way links
    directedlinks <- local.links %>%
      st_drop_geometry() %>%
      dplyr::select(from = from_id, to = to_id, weight, link_id) %>%
      rbind(., links2way)
    
    # ## for unsimplfied (where all links are already one way)
    # directedlinks <- local.links %>%
    #   st_drop_geometry() %>%
    #   dplyr::select(from = from_id, to = to_id, weight, link_id)
    
    # create directed graph
    graph <- graph_from_data_frame(directedlinks, directed = T, vertices = cyclable.nodes)
    
    
    # find cycle path start and end points, and order on basis of nearest to start 
    # of overall path, so segments are laid out with end of one matching start of next
    # helps matching start and end point (but wouldn't give good results for a U-shaped path)
    pathstart <- lwgeom::st_startpoint(path)
    pathend <- lwgeom::st_endpoint(path)
    
    if (st_distance(startpoint, pathstart) < st_distance(startpoint, pathend)) {
      point1 <- pathstart
      point2 <- pathend
    } else {
      point1 <- pathend
      point2 <- pathstart
    }

    # find 3 nearest nodes to each of start and end points, because sometimes nearest is
    # unsuitable - on one way road; note the igraph functions expect character input
    ## however, sometimes it would be better to use nearest - for future enhancement, add
    ## function to add back in link(s) from #1 to #2 and/or #3 where nearest not selected 
    ## and there is a direct connection
    nodes1 <- as.character(cyclable.nodes[st_nn(point1, cyclable.nodes, k=3)[[1]], "id"][[1]])
    nodes2 <- as.character(cyclable.nodes[st_nn(point2, cyclable.nodes, k=3)[[1]], "id"][[1]])
    
    # if nodes are different, proceed; if not (ie start and end are same) then return empty vectors
    if (!all(nodes1 == nodes2)) {
      #====== first direction ================================
      # if one of the 3 nearest has already been used as start/end point of another segment, then use it
      # (helps align segments end-to-end - this is what 'fromPoints' and 'toPoints' are for)
      if (length(intersect(nodes1, toPoints)) > 0) {
        fromNodes <- intersect(nodes1, toPoints)[1]
      } else {
        fromNodes <- nodes1
      }
      
      if (length(intersect(nodes2, fromPoints)) > 0) {
        toNodes <- intersect(nodes2, fromPoints)[1]
      } else {
        toNodes <- nodes2
      }
      
      distances1 <- distances(graph, fromNodes, toNodes, "out")  # from nodes 1 to nodes 2
      
      # find shortest path links for cycle path - first direction (but only if min dist not Inf)
      # 'function'shortest_paths' finds least cost path between nodes, and return edges ('epath')
      if (min(distances1) < Inf) {  # path exists ('Inf' means no path)
        # find row and column corresponding to min distances - row is [1] and column is [2]
        # if more than one - take the first
        shortestDist1 <- which(distances1 == min(distances1), arr.ind = TRUE)
        if (nrow(shortestDist1) > 1) {
          shortestDist1 <- shortestDist1[1,]
        }
        # find the vertex id's for the shortest distances
        shortestPair1 <- c(rownames(distances1)[shortestDist1[1]], 
                           colnames(distances1)[shortestDist1[2]])  
        
        if (shortestPair1[1] != shortestPair1[2]) { # selected nodes are not the same
          shortest1 <- shortest_paths(graph, 
                                      from = shortestPair1[1], 
                                      to = shortestPair1[2], 
                                      output = c("both"))
          
          # make vector of edges in the epath
          for (i in 1:length(shortest1$epath[[1]])) {
            # get the link_id for the link
            link <- directedlinks %>%
              filter(link_id == shortest1$epath[[1]][i][[1]]$link_id)  # this is required to return the link_id!!
            link_id <- link$link_id
            # add the link to the path.links vector
            path.links1 <- c(path.links1, link_id)
          }
          # remove duplicates (which arise from duplication for two-way links)
          path.links1 <- unique(path.links1)
          
          # add from and to points (so subsequent segments of the path can re-use them)
          fromPoints <- c(fromPoints, shortestPair1[1])
          toPoints <- c(toPoints, shortestPair1[2])
          
        } else {cat("No path found for cycle path #", path$scc_id, "segment", j, "direction 1 (selected nodes are same)\n")}
      } else {  # if min(distances1) is Inf - ie no path - then find partial paths from each end
        partial.path.link.outputs <- getPartialPathLinks(path,
                                                         graph,
                                                         directedlinks,
                                                         fromNodes,
                                                         toNodes,
                                                         path.links = path.links1,
                                                         sourcePoint = point1,
                                                         destPoint = point2,
                                                         fromPoints,
                                                         toPoints)
        path.links1 <- partial.path.link.outputs[[1]]
        fromPoints <- partial.path.link.outputs[[2]]
        toPoints <- partial.path.link.outputs[[3]]
        intersectionNodeA <- partial.path.link.outputs[[4]]
        intersectionNodeB <- partial.path.link.outputs[[5]]
        
        # create new links for the remaining gap
        new.links <- bridgeGaps(path, intersectionNodeA, intersectionNodeB, new.links) 
      }

      #====== second direction ================================
      # if one of the 3 nearest has already been used as start/end point of another segment, then use it
      # (helps align segments end-to-end - this is what 'fromPoints' and 'toPoints' are for)
      if (length(intersect(nodes2, toPoints)) > 1) {
        fromNodes <- intersect(nodes2, toPoints)[1]
      } else {
        fromNodes <- nodes2
      }
      
      if (length(intersect(nodes1, fromPoints)) > 1) {
        toNodes <- intersect(nodes1, fromPoints)[1]
      } else {
        toNodes <- nodes1
      }
      
      distances2 <- distances(graph, fromNodes, toNodes, "out")  # from nodes 2 to nodes 1
      
      # find shortest path links for cycle path - first direction (but only if min dist not Inf)
      # 'function'shortest_paths' finds least cost path between nodes, and return edges ('epath')
      if (min(distances2) < Inf) {  # path exists ('Inf' means no path)
        # find row and column corresponding to min distances - row is [1] and column is [2]
        # if more than one - take the first
        shortestDist2 <- which(distances2 == min(distances2), arr.ind = TRUE)
        if (nrow(shortestDist2) > 1) {
          shortestDist2 <- shortestDist2[1,]
        }
        # find the vertex id's for the shortest distances
        shortestPair2 <- c(rownames(distances2)[shortestDist2[1]], 
                           colnames(distances2)[shortestDist2[2]])  
        
        if (shortestPair2[1] != shortestPair2[2]) { # selected nodes are not the same
          shortest2 <- shortest_paths(graph, 
                                    from = shortestPair2[1], 
                                    to = shortestPair2[2], 
                                    output = c("both"))

          # make vector of edges in the epath
          for (i in 1:length(shortest2$epath[[1]])) {
            # get the link_id for the link
            link <- directedlinks %>%
              filter(link_id == shortest2$epath[[1]][i][[1]]$link_id)  # this is required to return the link_id!!
            link_id <- link$link_id
            # add the link to the path.links vector
            path.links2 <- c(path.links2, link_id)
          }
          # remove duplicates (which arise from duplication for two-way links)
          path.links2 <- unique(path.links2)
          
          # add from and to points (so subsequent segments of the path can re-use them)
          fromPoints <- c(fromPoints, shortestPair2[1])
          toPoints <- c(toPoints, shortestPair2[2])
          
        }  else {cat("No path found for cycle path #", path$scc_id, "segment", j, "direction 2 (selected nodes are same)\n")} # not expected to appear
      } else {  # if min(distances1) is Inf - ie no path - then find partial paths from each end
        partial.path.link.outputs <- getPartialPathLinks(path,
                                                         graph,
                                                         directedlinks,
                                                         fromNodes,
                                                         toNodes,
                                                         path.links = path.links2,
                                                         sourcePoint = point2,
                                                         destPoint = point1,
                                                         fromPoints,
                                                         toPoints)
        path.links2 <- partial.path.link.outputs[[1]]
        fromPoints <- partial.path.link.outputs[[2]]
        toPoints <- partial.path.link.outputs[[3]]
        intersectionNodeA <- partial.path.link.outputs[[4]]
        intersectionNodeB <- partial.path.link.outputs[[5]]
        
        # create new links for the remaining gap
        new.links <- bridgeGaps(path, intersectionNodeA, intersectionNodeB, new.links) 
      }
      
     } else {cat("No path for cycle path #", path$scc_id, "segment", j, "because start and end nodes are the same (too short)\n")}
  } else {cat("No local links found for cycle path #", path$scc_id, "segment", j, ", probably outside Melbourne\n")}

  # return the outputs as a list
  return(list(path.links1, path.links2, new.links, fromPoints, toPoints))
}
