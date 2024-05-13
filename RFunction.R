library('move2')
library('lubridate')
library('sf')
library('dplyr')
library('magrittr')
library('dbscan')

## The parameter "data" is reserved for the data object passed on from the previous app

# to display messages to the user in the log file of the App in MoveApps
# one can use the function from the logger.R file:
# logger.fatal(), logger.error(), logger.warn(), logger.info(), logger.debug(), logger.trace()

# Showcase injecting app setting (parameter `year`)
rFunction = function(data, 
                     eps1,
                     eps2, 
                     minPts,
                     startdate,
                     enddate
                     ) {
  
  
  
  
  # 1. Input Checks -------------------------------
  
  if (nrow(data) < minPts) {
    logger.warning("Input data contains fewer locations than minimum cluster size. No clusters can be generated")
  }
  
  if (!is.null(startdate)) {
    data %<>% filter(
      timestamp > as_datetime(startdate, tz = "UTC")
    )
  }
  if (!is.null(enddate)) {
    data %<>% filter(
      timestamp < as_datetime(enddate, tz = "UTC")
    )
  }
  
  
  # 2. Run ST-DBSCAN Clustering --------------------
  
  newclusts <- st_dbscan(
    data,
    eps1 = eps1,
    eps2 = eps2,
    minpts = minPts,
    epsg = st_crs(data)
  )
  
  
  return(newclusts)
}



#' Here, I'll implement a version of the ST-DBSCAN Algorithm (Birant & Kut, 2008) 
#' for use with spatiotemporal data. In particular, I'll modularise it for 
#' use with Move2 objects to allow re-application to other ecological data.

#' Inputs:
#' - data: a move2 object (see move2 documentation)
#' - eps1: the spatial connectivity radius (in metres)
#' - eps2: the temporal connectivity radius (in days)
#' - minpts: the minimum number of data points for a core point 
#' - epsg: the EPSG code for this data (to allow accurate calculations 
#' in metric units)
#' Outputs:
#' - data: the input data move2 object, with an additional column named 
#' "xy.clust" containing the cluster ID for each location

st_dbscan <- function(data, eps1, eps2, minpts, epsg = 29333) {
  
  # Verify that data is in move2 format - if not, return and restart
  if (!mt_is_move2(data)) {
    stop("Data is not in Move2 format. Please reformat and try again")
  }
  
  data %<>% mutate(
    index = 1:nrow(.),
    timestamp = mt_time(.)
  ) 
  
  # Extract coordinates 
  coords <- st_coordinates(data %>% st_transform(epsg))
  
  # 1. Run spatial DBSCAN on coordinate data
  spatdata <- data %>% mutate(
    spatclust = dbscan(coords, eps = eps1, minPts = minpts)$cluster
  ) 
  spatdata_filtered <- spatdata %>%
    # Remove non-clusters:
    filter(spatclust != 0)
  
  if (all(spatdata$spatclust == 0)) {
    warning("No clusters are identified within the dataset. Returning input")
    return(
      data %>%
        mutate(xy.clust = NA)
    )
  }
  
  spatdata_df <- spatdata_filtered %>%
    as.data.frame() %>%
    dplyr::select(c("index", "timestamp", "spatclust"))
  
  # Set up table for spatial-temporal cluster matching
  assigntable <- data.frame(
    index = 0,
    timeclust = 0
  )
  
  # Iterate over the spatial clusters and perform temporal DBSCAN on their contents:
  #' This sub-divides the spatial clusters into their temporally-differentiated
  #' constituent clusters.
  for (i in 1:max(spatdata$spatclust)) {
    
    # Retrieve the lowest possible cluster IDs to assign to new clusters:
    minclustID <- max(assigntable$timeclust, na.rm = T) 
    
    # Retrieve cluster data and perform internal clustering
    clustdat <- spatdata_df %>%
      filter(spatclust == i) %>%
      mutate(
        timeclust = dbscan(
          .$timestamp %>% as.numeric() %>% as.matrix(),
          eps = days(eps2) %>% as.numeric(),
          minPts = minpts
        )$cluster
        
      ) %>%
      filter(timeclust != 0) %>%
      mutate(timeclust = timeclust + minclustID) %>%
      # Adding the highest clustID brings the lowest possible cluster identifier up to the 
      # lowest available integer number. I.e. if we've already identified 5 
      # clusters, the next will be called 6.
      
      # Select only the cluster assignment data:
      dplyr::select(c("index", "timeclust")) 
    
    # Store outcome and move to next iteration
    assigntable %<>% rbind(clustdat)
    
  }
  
  if (nrow(assigntable) == 1) {
    warning("No temporal clusters have been generated. Returning")
    return(
      data %>%
        mutate(xy.clust = NA)
    )
  }
  
  # Add the temporal clusters to the data:
  allclusts <- spatdata %>% 
    left_join(assigntable %>% filter(index != 0), by = "index")
  
  # Generate new cluster IDs:
  newclusts <- allclusts %>% 
    as.data.frame() %>%
    filter(spatclust != 0 & !is.na(timeclust)) %>%
    group_by(spatclust, timeclust) %>%
    summarise(
      count = n(),
      .groups = "drop"
    ) %>%
    mutate(xy.clust = row_number()) %>%
    dplyr::select(-count)
  
  # Clean up and return:
  outdat <- allclusts %>% 
    left_join(newclusts, by = c("spatclust", "timeclust")) %>%
    dplyr::select(-c("spatclust", "timeclust", "index"))
  
  return(outdat)
  
}
