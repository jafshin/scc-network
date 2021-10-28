# Function - if possile, to further links can be added at start and end (selected links will start and
# end at one of closest 3 - may be possible to add additional short links from closest point)

# Operates by checking whether start/end node is 2nd or 3rd closest, and if it is, sees whether
# (for 2nd) there is a link from 1st to 2nd, or
# (for 3rd) there are links from 1st to 2nd to 3rd, or from 2nd to 3rd
# and equivalent in reverse diraction

# Note that the script tests whether any of the 3 closest nodes are in 'fromPoints' or 'toPoints'
# However, where new links have been created, they will automaticlaly go from the closest node (and
# that closest node won't necessarily be in fromPoints or toPoints)


getStartEndLinks <- function(startpoint, endpoint, fromPoints, toPoints) {
  # create empty vector for outputs to be returned
  start.end.links <- c()  # to return short links at start and end
  
  # return 3 closest nodes to start and end points (returns in order, at least it seems so)
  startNodes <- cyclable.nodes[st_nn(startpoint, cyclable.nodes, k=3)[[1]], "id"][[1]]
  endNodes <- cyclable.nodes[st_nn(endpoint, cyclable.nodes, k=3)[[1]], "id"][[1]]
  
  #====== start of path ================================
  usable.start.links <- c()
  # check whether closest startNode is already in fromPoints (if it is, no further link needed at start)
  if (!startNodes[1] %in% fromPoints) { 
    # if not, check whether 2nd closest startNode is in fromPoints
    if (startNodes[2] %in% fromPoints) {
      #if it is, then test whether link from 1st to 2nd closest can be added
      # ('test' is whether from and to in correct order, or from and to in reverse order but link is two-way)
      usable.start.links <- cyclable.links %>%
        filter((from_id == startNodes[1] & to_id == startNodes[2]) |
                 (from_id == startNodes[2] & to_id == startNodes[1] & is_oneway == 0)) %>%
        .$link_id
    } else if (startNodes[3] %in% fromPoints) {  
      # if it is, then test sequentially whether links from 1st to 3rd, 1st to 2nd to 3rd, or 2nd to 3rd can be added
      # try 1st to 3rd
      usable.start.links <- cyclable.links %>% 
        filter((from_id == startNodes[1] & to_id == startNodes[3]) |
                 (from_id == startNodes[3] & to_id == startNodes[1] & is_oneway == 0)) %>%
        .$link_id
      # if not, try 1st to 2nd to 3rd - but if only 2nd to 3rd found, use that
      if (length(usable.start.links) == 0) {
        first.leg <- cyclable.links %>% 
          filter((from_id == startNodes[1] & to_id == startNodes[2]) |
                    (from_id == startNodes[2] & to_id == startNodes[1] & is_oneway == 0)) %>%
          .$link_id
        second.leg <- cyclable.links %>% 
          filter((from_id == startNodes[2] & to_id == startNodes[3]) |
                   (from_id == startNodes[3] & to_id == startNodes[2] & is_oneway == 0)) %>%
          .$link_id
        if(length(first.leg) > 0 & length(second.leg) > 0) {
          usable.start.links <- c(first.leg, second.leg)
        } else if (length(second.leg) > 0) {
          usable.start.links <- second.leg
        }
      }
    }
  }
    
    
  #====== end of path ================================
  usable.end.links <- c()
  # check whether closest endNode is already in toPoints (if it is, no further link needed at end)
  if (!endNodes[1] %in% toPoints) { 
     # if not, check whether 2nd closest endNode is in toPoints
    if (endNodes[2] %in% toPoints) {
      #if it is, then test whether link from 2nd to 1st closest can be added
      # ('test' is whether from and to in correct order, or from and to in reverse order but link is two-way)
      usable.end.links <- cyclable.links %>%
        filter((from_id == endNodes[2] & to_id == endNodes[1]) |
                 (from_id == endNodes[1] & to_id == endNodes[2] & is_oneway == 0)) %>%
        .$link_id
    } else if (endNodes[3] %in% toPoints) {  
      # if it is, then test sequentially whether links from 3rd to 1st, 3rd to 2nd to 1st, or 3rd to 2nd can be added
      # try 3rd to 1st
      usable.end.links <- cyclable.links %>% 
        filter((from_id == endNodes[3] & to_id == endNodes[1]) |
                 (from_id == endNodes[1] & to_id == endNodes[3] & is_oneway == 0)) %>%
        .$link_id
      # if not, try 3rd to 2nd to 1st - but if only 3rd to 2nd found, use that
      if (length(usable.end.links) == 0) {
        first.leg <- cyclable.links %>% 
          filter((from_id == endNodes[3] & to_id == endNodes[2]) |
                   (from_id == endNodes[2] & to_id == endNodes[3] & is_oneway == 0)) %>%
          .$link_id
        second.leg <- cyclable.links %>% 
          filter((from_id == endNodes[2] & to_id == endNodes[1]) |
                   (from_id == endNodes[1] & to_id == endNodes[2] & is_oneway == 0)) %>%
          .$link_id
        if (length(first.leg) > 0 & length(second.leg) > 0) {
          usable.end.links <- c(first.leg, second.leg)
        } else if (length(first.leg) > 0) {
          usable.end.links <- first.leg
        }
      }
    }

  }
  
  start.end.links <- c(usable.start.links, usable.end.links)
 
  return(start.end.links)
  
}
