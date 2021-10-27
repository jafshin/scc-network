# Function - bridge gaps between disconnected graphs at each end of cycle path segment,
# by creating links

# Note that each link is one-directional (but usually, the gap is the same in both directions,
# so a second link is created in the reverse direction, identical but for reversed 'to'and 'from' attributes)

bridgeGaps <- function(path, intersectionNodeA, intersectionNodeB, new.links) {
  # find nearest points along line https://gis.stackexchange.com/questions/288570/find-nearest-point-along-polyline-using-sf-package-in-r
  # convert path to points at 1m distance (density is in points per distance unit)
  pathPoints <- st_line_sample(path, density = 1) %>%
    st_sf() %>%
    st_cast("POINT") %>%
    mutate(group = 1)
  
  # find the nodes (as points) corresponding to the id's in intersectionNodeA and intersectionNodeB
  nodeA <- cyclable.nodes %>% 
    filter(id == intersectionNodeA)
  nodeB <- cyclable.nodes %>% 
    filter(id == intersectionNodeB)
  
  # find path points closest to intersection Nodes A and B
  pathPointA <- which.min(st_distance(nodeA, pathPoints))
  pathPointB <- which.min(st_distance(nodeB, pathPoints))
  
  # each point represents a 1m increment along the step, so 'break points' in the line
  # are pathPointA, then sequence in increments of 1000 (or other 'gapSegmentDistance') up to pathPointB
  # in each case, calculate steps from the minimum  value - where gaps in each direction
  # are aligned (which is the usual case), this will result in matching sections
  # eg one direction is at points 200, 1200, 2200, 2500; the reverse is 2500, 2200, 1200, 200
  breakpoints <- c()
  if (pathPointA < pathPointB) {
    breakpoints <- seq(pathPointA, pathPointB, by = gapSegmentDistance)
    if (!pathPointB %in% breakpoints) breakpoints <- c(breakpoints, pathPointB)
  } else {
    breakpoints <- seq(pathPointB, pathPointA, by = gapSegmentDistance)
    if (!pathPointA %in% breakpoints) breakpoints <- c(breakpoints, pathPointA)
    # and reverse order, so sequence begins at pathPointA
    breakpoints <- rev(breakpoints)
  }
  
  # find total length of path to be bridged (which is just the number of breakpoints)
  pathLength <- max(breakpoints) - 1
 
  # create links for section between each pair of breakpoints; first and last are Nodes A and B;
  # others are nearest nodes to the breakpoints (ie the pathpoints with the breakpoint numbers)
  if (length(breakpoints) > 1)  {
    # trial new links - subject to length adjustment
    trialNewLinks <- c()
    
    for (i in 1:(length(breakpoints)-1)) {
      # start node for section
      if (i == 1) {
        startNode <- nodeA
      } else {
        startNode <- cyclable.nodes[st_nearest_feature(pathPoints[breakpoints[i],], cyclable.nodes), ]
      }
      # end node for section
      if (i < length(breakpoints)-1) {
        endNode <- cyclable.nodes[st_nearest_feature(pathPoints[breakpoints[i+1],], cyclable.nodes), ]
      } else {
        endNode <- nodeB
      }
      
      # trial length of section - distance between start and end nodes (will be adjusted subsequently)
      trialSectionLength <- as.numeric(st_distance(startNode, endNode))
      
      # create the new link
      newlink <- rbind(startNode, endNode) %>%
        st_union(.) %>%
        st_cast("LINESTRING") %>%
        st_as_sf(crs = st_crs(links)) %>% # same crs as links
        rename(GEOMETRY = colnames(.)) %>%  # renames the sole column (inicially called x)
        # note - details below are for simplfied links file; column  names
        # may need to change to match structure of unsimplfied (or any other) links file
        mutate(from_id = startNode$id,
               to_id = endNode$id,
               fromx = st_coordinates(startNode)[1],
               fromy = st_coordinates(startNode)[2],
               tox = st_coordinates(endNode)[1],
               toy = st_coordinates(endNode)[2],
               length = trialSectionLength,
               freespeed = 30/3.6,
               permlanes = 1,
               capacity = 300,
               highway = "cycleway",
               is_oneway = 1,
               cycleway = "bikepath",
               is_cycle = 1,
               is_walk = 0,
               is_car = 0,
               modes = "bike",
               scc_id = path$scc_id,
               scc_type = path$type,
               scc_new = 1)  # '1' in this field indicates that the link is added by the bridge gap process
      
      trialNewLinks <- rbind(trialNewLinks, newlink)
    }
    
    # adjust link lengths so total equals 'pathLength', divided proportionately to link length
    trialSectionLengthsTotal <- sum(trialNewLinks$length)
    trialNewLinks <- trialNewLinks %>%
      mutate(length = length * (pathLength/trialSectionLengthsTotal)) %>%
      # omit zero-length links
      filter(length > 0)
    
    # add to new.links
    new.links <- rbind(new.links, trialNewLinks)
    
  }
  
  return(new.links)
  
}