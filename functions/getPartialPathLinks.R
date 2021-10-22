# Function - where can't find network links comprising shortest path for 
# whole of given given cycle path segment, finds links for part that can be found

# Returns 'intersectionNodeA' and 'intersectionNodeB', which are the point 
# at each end of the gap that remains to be bridged

# Note - only finds matching links at start and end of cycle path segment, leaving
# gap in the middle to be bridged with 'bridgeGaps' function; doesn't find
# any partly-routable areas in the middle that are disconnected from start and end

# [Note - the first- and second-direction parts of this code are very similar;
# could probably be re-written as a single function wiht better parameters]

getPartialPathLinks <- function(path, graph, directedlinks, fromNodes, toNodes,
                                path.links, sourcePoint, destPoint, fromPoints, toPoints) {

  # create empty vectors for outputs to be returned
  intersectionNodeA <- c()
  intersectionNodeB <- c()
  
  #====== from start of segment to the gap ================================
  # for the from-nodes, find the one(s) from which you can get furtherst towards to-node
  # initialise candidates - will hold the from-nodes, and how close you can get from them
  # candidates are the 3 nearest nodes to each of start and end points, 
  # because sometimes nearest is unsuitable (eg on one way road in wrong direction)
  candidate.nodes <- NULL
  for (i in 1:length(fromNodes)) {
    # find the vertices reachable from the from-node
    reachable <- subcomponent(graph, fromNodes[i], "out")
    reachable.nodes <- cyclable.nodes %>%
      filter(id %in% as.numeric(reachable$name))
    # if more than one, then convert those vertices into a convex hull, and into a linestring
    if (nrow(reachable.nodes) > 1) {
      reachable.area <- st_union(reachable.nodes) %>%
        st_convex_hull(.)
      
      # if path wholly within convex hull, then notional intersection point is destPoint
      # otherwise, find intersection of convex hull & path that is closest to the to-point (best)
      # and calculate distance, and find its 3 closest nodes in the reachable area
      if (st_within(path, reachable.area, sparse = FALSE)) {
        best.path.intersec <- destPoint %>%
          st_sf() %>%
          mutate(node = fromNodes[i],
                 dist.to.dest = 0,
                 toIntersecNodes = list(as.character(reachable.nodes[st_nn(., reachable.nodes, 
                                                                           k=min(nrow(reachable.nodes), 3))[[1]], "id"][[1]])))
        candidate.nodes <- rbind(candidate.nodes, best.path.intersec)
      } else {
        reachable.area <- reachable.area %>%
          st_cast(., "MULTILINESTRING")
        path.intersecs <- st_intersection(path, reachable.area)
        if (nrow(path.intersecs) > 0) {
          best.path.intersec <- path.intersecs[st_nearest_feature(destPoint, path.intersecs)] %>%
            mutate(node = fromNodes[i],
                   dist.to.dest = st_distance(., destPoint),
                   toIntersecNodes = list(as.character(reachable.nodes[st_nn(., reachable.nodes, 
                                                                             k=min(nrow(reachable.nodes), 3))[[1]], "id"][[1]])))
          candidate.nodes <- rbind(candidate.nodes, best.path.intersec)
        }
      }
    }
  }
  
  # keep only the pairs of the candidate(s) with the lowest distance to the to-node
  if (length(candidate.nodes) > 0) {  # only proceed if there are candidate nodes - otherwise, nearest node to from-point becomes intersectionNodeA
    candidate.pairs <- candidate.nodes %>%
      filter(dist.to.dest == min(dist.to.dest)) %>%
      st_drop_geometry()
    
    # for each candidate, find which of the 'toIntersecNodes' produces the shortest distance, and find that distance
    for (i in 1:nrow(candidate.pairs)) {
      outNode <- candidate.pairs[i, "node"]
      toIntersecNodes <- candidate.pairs[i, "toIntersecNodes"]
      
      distancesOut <- distances(graph, outNode, toIntersecNodes[[1]], "out")  # from outNode to toIntersecNodes
      # find the 'toIntersecNode' with the lowest distance, and add to candidate pairs
      minDist <- min(distancesOut)
      
      candidate.pairs[i, "toIntersecNode"] <- colnames(distancesOut)[which(distancesOut == minDist, arr.ind = TRUE)[[2]]]
      candidate.pairs[i, "dist"] <- minDist
    }
    
    # from the candidates, find the from-node with the minimum distance
    best.pair <- candidate.pairs %>%
      filter(dist == min(dist))
    
    fromNode <- best.pair$node[1] # if more than one 'best', then take the first
    toIntersectNode <- best.pair$toIntersecNode[1]  # if more than one 'best', then take the first
    
    ## find the shortest path between fromNode and toIntersectNode
    if (fromNode != toIntersectNode) { # selected nodes are not the same
      shortestFrom1 <- shortest_paths(graph, 
                                      from = fromNode, 
                                      to = toIntersectNode, 
                                      output = c("both"))
      
      # make vector of edges in the epath
      for (i in 1:length(shortestFrom1$epath[[1]])) {
        # get the link_id for the link
        link <- directedlinks %>%
          filter(link_id == shortestFrom1$epath[[1]][i][[1]]$link_id)  # this is required to return the link_id!!
        link_id <- link$link_id
        # add the link to the path.links vector
        path.links <- c(path.links, link_id)
      }
      # remove duplicates (which arise from duplication for two-way links)
      path.links <- unique(path.links)
      
      fromPoints <- c(fromPoints, fromNode)
      intersectionNodeA <- toIntersectNode
    } else {
      intersectionNodeA <- as.character(cyclable.nodes[st_nearest_feature(sourcePoint,cyclable.nodes), "id"][[1]])
    }
  } else {
    intersectionNodeA <- as.character(cyclable.nodes[st_nearest_feature(sourcePoint,cyclable.nodes), "id"][[1]])
  }
  #====== from the gap to end  of segment to the gap  ================================  
  # initialise candidates - will hold the to-nodes, and whats' the closest point to from-node
  # from which you can reach them
  # candidates are the 3 nearest nodes to each of start and end points, 
  # because sometimes nearest is unsuitable (eg on one way road in wrong direction)
  candidate.nodes <- NULL
  for (i in 1:length(toNodes)) {
    # find the vertices from which the to-node is reachable
    reachable <- subcomponent(graph, toNodes[i], "in")
    reachable.nodes <- cyclable.nodes %>%
      filter(id %in% as.numeric(reachable$name))
    # if more than one, then convert those vertices into a convex hull, and into a linestring
    if (nrow(reachable.nodes) > 1) {
      reachable.area <- st_union(reachable.nodes) %>%
        st_convex_hull(.)
      
      # if path wholly within convex hull, then notional intersection point is sourcePoint
      # otherwise, find intersection of convex hull & path that is closest to the from-point (best)
      # and calculate distance, and find its 3 closest nodes in the reachable area
      if (st_within(path, reachable.area, sparse = FALSE)) {
        best.path.intersec <- sourcePoint %>%
          st_sf() %>%
          mutate(node = toNodes[i],
                 dist.to.dest = 0,
                 fromIntersecNodes = list(as.character(reachable.nodes[st_nn(., reachable.nodes, 
                                                                             k=min(nrow(reachable.nodes), 3))[[1]], "id"][[1]])))
        candidate.nodes <- rbind(candidate.nodes, best.path.intersec)
      } else {
        reachable.area <- reachable.area %>%
          st_cast(., "MULTILINESTRING")
        path.intersecs <- st_intersection(path, reachable.area)
        if (nrow(path.intersecs) > 0) {
          best.path.intersec <- path.intersecs[st_nearest_feature(sourcePoint, path.intersecs)] %>%
            mutate(node = toNodes[i],
                   dist.to.dest = st_distance(., sourcePoint),
                   fromIntersecNodes = list(as.character(reachable.nodes[st_nn(., reachable.nodes, 
                                                                               k=min(nrow(reachable.nodes), 3))[[1]], "id"][[1]])))
          candidate.nodes <- rbind(candidate.nodes, best.path.intersec)
        }
      }
    }
  }
  
  # keep only the pairs of the candidate(s) with the lowest distance to the to-node
  if (length(candidate.nodes) > 0) {  # only proceed if there are candidate nodes - otherwise, nearest node to to-point becomes intersectionNodeB
    candidate.pairs <- candidate.nodes %>%
      filter(dist.to.dest == min(dist.to.dest)) %>%
      st_drop_geometry()
    
    # for each candidate, find which of the 'fromIntersecNodes' produces the shortest distance, and find that distance
    for (i in 1:nrow(candidate.pairs)) {
      inNode <- candidate.pairs[i, "node"]
      fromIntersecNodes <- candidate.pairs[i, "fromIntersecNodes"]
      
      distancesIn <- distances(graph, inNode, fromIntersecNodes[[1]], "in")  # to inNode from fromIntersecNodes
      # find the 'fromIntersecNode' with the lowest distance, and add to candidate pairs
      minDist <- min(distancesIn)
      
      candidate.pairs[i, "fromIntersecNode"] <- colnames(distancesIn)[which(distancesIn == minDist, arr.ind = TRUE)[[2]]]
      candidate.pairs[i, "dist"] <- minDist
    }
    
    # from the candidates, find the from-node with the minimum distance
    best.pair <- candidate.pairs %>%
      filter(dist == min(dist))
    
    toNode <- best.pair$node[1] # if more than one 'best', then take the first
    fromIntersectNode <- best.pair$fromIntersecNode[1]  # if more than one 'best', then take the first
    
    ## find the shortest path between toNode and fromIntersectNode
    if (toNode != fromIntersectNode) { # selected nodes are not the same
      shortestTo1 <- shortest_paths(graph, 
                                    from = fromIntersectNode, 
                                    to = toNode, 
                                    output = c("both"))
      
      # make vector of edges in the epath
      for (i in 1:length(shortestTo1$epath[[1]])) {
        # get the link_id for the link
        link <- directedlinks %>%
          filter(link_id == shortestTo1$epath[[1]][i][[1]]$link_id)  # this is required to return the link_id!!
        link_id <- link$link_id
        # add the link to the path.links vector
        path.links <- c(path.links, link_id)
      }
      # remove duplicates (which arise from duplication for two-way links)
      path.links <- unique(path.links)
      
      toPoints <- c(toPoints, toNode)
      intersectionNodeB <- fromIntersectNode
    } else {
      intersectionNodeB <- as.character(cyclable.nodes[st_nearest_feature(destPoint,cyclable.nodes), "id"][[1]])
    }
  } else {
    intersectionNodeB <- as.character(cyclable.nodes[st_nearest_feature(destPoint,cyclable.nodes), "id"][[1]])
  }
  
  # return the outputs as a list
  return(list(path.links, fromPoints, toPoints, intersectionNodeA, intersectionNodeB))
  
}

