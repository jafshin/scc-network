# Function - to find network links comprising shortest path for given cycle path

# getPathLinks <- function(path, bufferDistance, nextVertices) {
getPathLinks <- function(path, bufferDistance, startpoint, fromPoints, toPoints) {
  # getPathLinks <- function(path, bufferDistance) {
  # create empty vectors for outputs to be returned
  cat("Finding links for cycle path #", path$scc_id, "segment", j, "\n")
  path.links1 <- c()
  path.links2 <- c()
  
  # buffer cycle path to given distance (constrains to network links near cycle path route)
  path.buffered <- st_buffer(path, bufferDistance)
  
  # identify network links within the buffer
  local.links <- cyclable.links %>% 
    filter(st_intersects(GEOMETRY, path.buffered, sparse = FALSE))
  
  # if there are local links, proceed; if not (ie paths outside study region) then return empty vectors 
  if (nrow(local.links) > 0) {
    
    ### for simplified (where there is 'is_oneway' column) -
    # create graph to find shortest paths in network links for the given cycle path
    
    ### for simplified (where there is 'is_oneway' column) -
    # for two-way links, create an extra set of links with 'from' and 'to' reversed
    links2way <- local.links %>%
      st_drop_geometry() %>%
      filter(is_oneway == 0) %>%
      mutate(from = to_id, to = from_id) %>%  # reverse order
      dplyr::select(from, to, weight, link_id)  # 'weight' is created in #2 from 'length'
    
    # combine links with the extra set of reversed links for two-way links
    directedlinks <- local.links %>%
      st_drop_geometry() %>%
      dplyr::select(from = from_id, to = to_id, weight, link_id) %>%
      rbind(., links2way)
    
    ### for unsimplfied (where links are already directed)
    # directedlinks <- local.links %>%
    #   st_drop_geometry() %>%
    #   dplyr::select(from = from_id, to = to_id, weight, link_id)
    
    
    # create directed graph
    graph <- graph_from_data_frame(directedlinks, directed = T, vertices = cyclable.nodes)
    
    
    # get id's of nodes nearest to cycle path start and end points
    # point1 <- lwgeom::st_startpoint(path)
    # point2 <- lwgeom::st_endpoint(path)
    
    # alternative - get id's of nodes nearst to cycle path start and end points
    # ordered on basis of nearest to start of overall path
    ## assumes path segments are laid out end to end in increasing distance from start of
    ## overall path - wouldn't give good results for a U-shaped path
    pathstart <- lwgeom::st_startpoint(path)
    pathend <- lwgeom::st_endpoint(path)
    
    if (st_distance(startpoint, pathstart) < st_distance(startpoint, pathend)) {
      point1 <- pathstart
      point2 <- pathend
    } else {
      point1 <- pathend
      point2 <- pathstart
    }
    
    
    # node1 <- as.character(cyclable.nodes[st_nearest_feature(point1, cyclable.nodes), "id"][[1]])
    # node2 <- as.character(cyclable.nodes[st_nearest_feature(point2, cyclable.nodes), "id"][[1]])
    
    ## or alternative, to return nearest 3 - the igraph functions expect character input
    nodes1 <- as.character(cyclable.nodes[st_nn(point1, cyclable.nodes, k=3)[[1]], "id"][[1]])
    nodes2 <- as.character(cyclable.nodes[st_nn(point2, cyclable.nodes, k=3)[[1]], "id"][[1]])
    
    # if nodes are different, proceed; if not (ie start and end are same) then return empty vectors
    # if (node1 != node2) {
    if (nodes1 != nodes2) {
      # find best node pairs in each direction (shortest distance), trying 3 closest nodes at each end 
      # find shortest distance, in each direction 1 and 2, between nearest 3 nodes at start and end
      # distances1 <- distances(graph, nodes1, nodes2, "out")  # from nodes 1 to nodes 2
      # distances2 <- distances(graph, nodes2, nodes1, "out")  # from nodes 2 to nodes 1
      
      # alternative that uses already-used nodes
      # if any of nodes1 has already been used as a toPoint, then use it as the from node, and vice versa 
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
      
      

      
      
      # next vertices - the ones closest to the end points
      # paths are digitised so that 'startpoint' of path 1 = 'endpoint' of path 2
      # nextVertices <- c(shortestPair1[1], shortestPair2[2])
      
      
      
      # find shortest path links for cycle path - first direction (but only if min dist not Inf)
      # find shortest path between nodes, and return edges ('epath')
      if (min(distances1) < Inf) { # path exists
        # find row and column corresponding to min distances - row is [1] and column is [2]
        shortestDist1 <- which(distances1 == min(distances1), arr.ind = TRUE)
        # find the vertex id's for the shortest distances
        shortestPair1 <- c(rownames(distances1)[shortestDist1[1]], 
                           colnames(distances1)[shortestDist1[2]])  
        
        if (shortestPair1[1] != shortestPair1[2]) { # selected nodes are not the same
          shortest1 <- shortest_paths(graph, 
                                      from = shortestPair1[1], 
                                      to = shortestPair1[2], 
                                      output = c("both"))
          
          # ## alternative of shortest, 3x3
          # shortest1 <- shortest_paths(graph, 
          #                             from = node1, 
          #                             to = nodes2,  
          #                             output = c("both"))
          
          
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
          
          fromPoints <- c(fromPoints, shortestPair1[1])
          toPoints <- c(toPoints, shortestPair1[2])
          
        } else {cat("No path found for cycle path #", path$scc_id, "segment", j, "direction 1 (selected nodes are same)\n")}
      } else {cat("No path found for cycle path #", path$scc_id, "segment", j, "direction 1 (infinity exception)\n")}
      #========to be new function (maybe)==============================
      # ## to work on this - comment out from 'cat' above, and comment in } below
      # ## this is only going to bridge a single gap
      # 
      #   # for the from-nodes, find the one(s) from which you can get furtherst towards to-node
      #   # initialise candidates - will hold the from-nodes, and how close you can get from them
      #   candidate.nodes <- NULL
      #   for (i in 1:length(fromNodes)) {
      #     # find the vertices reachable from the from-node
      #     reachable <- subcomponent(graph, fromNodes[i], "out")
      #     reachable.nodes <- cyclable.nodes %>%
      #                                   filter(id %in% as.numeric(reachable$name))
      #     # convert those vertices into a convex hull, and into a linestring
      #     reachable.area <- st_union(reachable.nodes) %>%
      #       st_convex_hull(.) %>%
      #       st_cast(., "MULTILINESTRING")
      # 
      #     # find intersection of convex hull & path that is closest to the to-point (best)
      #     # and calculate distance, and find its 3 closest nodes in the reachable area
      #     path.intersecs <- st_intersection(path, reachable.area)
      #     best.path.intersec <- path.intersecs[st_nearest_feature(point2, path.intersecs)] %>%
      #       mutate(node = fromNodes[i],
      #              dist.to.dest = st_distance(., point2),
      #              toIntersecNodes = list(as.character(reachable.nodes[st_nn(., reachable.nodes, k=3)[[1]], "id"][[1]])))
      #     candidate.nodes <- rbind(candidate.nodes, best.path.intersec)
      #   }
      #   # keep only the pairs of the candidate(s) with the lowest distance to the to-node
      #   candidate.pairs <- candidate.nodes %>%
      #     filter(dist.to.dest == min(dist.to.dest)) %>%
      #     st_drop_geometry()
      # 
      #   # for each candidate, find which of the 'toIntersecNodes' produces the shortest distance, and find that distance
      #   for (i in 1:nrow(candidate.pairs)) {
      #     outNode <- candidate.pairs[i, "node"]
      #     toIntersecNodes <- candidate.pairs[i, "toIntersecNodes"]
      #     
      #     distancesOut <- distances(graph, outNode, toIntersecNodes[[1]], "out")  # from outNode to toIntersecNodes
      #     # find the 'toIntersecNode' with the lowest distance, and add to candidate pairs
      #     minDist <- min(distancesOut)
      #     
      #     candidate.pairs[i, "toIntersecNode"] <- colnames(distancesOut)[which(distancesOut == minDist, arr.ind = TRUE)[[2]]]
      #     candidate.pairs[i, "dist"] <- minDist
      #   }
      #   
      #   # from the candidates, find the from-node with the minimum distance
      #   best.pair <- candidate.pairs %>%
      #     filter(dist == min(dist))
      #   
      #   fromNode <- best.pair$node
      #   toIntersectNode <- best.pair$toIntersecNode
      #   
      #   ## next, find the shortest path between fromNode and toIntersectNode 
      #   
      #   
      #   
      #   
      #   ### then do it all again in the reverse direction
      #   
      #   
      #   ### then bridge the gap (probably first dividing into 500m steps)
      # 
      # 
      # }
      #=====================================
 
      # if any of nodes2 has already been used as a toPoint, then use it as the from node, and vice versa 
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
      
      # find shortest path links for cycle path - second direction (but only if min dist not Inf)
      # find shortest path between nodes, and return edges ('epath')
      if (min(distances2) < Inf) {  # path exists
        # find row and column corresponding to min distances - row is [1] and column is [2]
        shortestDist2 <- which(distances2 == min(distances2), arr.ind = TRUE)
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
          
          fromPoints <- c(fromPoints, shortestPair2[1])
          toPoints <- c(toPoints, shortestPair2[2])
          
        }  else {cat("No path found for cycle path #", path$scc_id, "segment", j, "direction 2 (selected nodes are same)\n")} # not expected to appear
      } else {cat("No path found for cycle path #", path$scc_id, "segment", j, "direction 2 (infinity exception)\n")}
       
      

      # find paths between nodes, and return edges ('epath')
      # shortest2 <- shortest_paths(graph, 
      #                             from = node2, 
      #                             to = node1, 
      #                             output = c("both"))
     } else {cat("No path for cycle path #", path$scc_id, "segment", j, "because start and end nodes are the same (too short)\n")}
  } else {cat("No local links found for cycle path #", path$scc_id, "segment", j, ", probably outside Melbourne\n")}

  # return the output as a list
  # return(list(path.links1, path.links2))
  # return(list(path.links1, path.links2, nextLinks))
  return(list(path.links1, path.links2, fromPoints, toPoints))
}
