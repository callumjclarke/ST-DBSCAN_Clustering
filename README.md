# STDBSCAN Cluster Detection

MoveApps

Github repository: *github.com/callumjclarke/ST-DBSCAN_Clustering*

## Description

*This MoveApp is a limited adaptation of the ST-DBSCAN (Birant & Kut, 2007) spatiotemporal clustering process for event identification.*

Spatiotemporal Density-Based Spatial Clustering of Applications with Noise (ST-DBSCAN) is an extension of the DBSCAN clustering algorithm for spatiotemporal data. We implement a section of the work proposed by Birant & Kut (2007).

## Documentation

Two reachability distances are defined: a maximum *spatial* reachability distance (eps1) and a maximum *temporal* reachability distance (eps2). Generated clusters constitute points that are both spatially-reachable and temporally-reachable.

The data is first filtered to the user-selected window (defaulting to the full dataset). The clustering is then performed through a spatial DBSCAN process on the full dataset, followed by a secondary temporal DBSCAN on each individually-identified spatial cluster. The resulting clusters are, by definition, the union of spatial connectivity and temporal connectivity.

Each location's associated cluster ID is attached to the output data as an additional column called *xy.clust.*

### Input data

*move2* location object. This app performs best if the input data is in UTM format. See MoveApp *Standardise Formats and Calculate Basic Statistics* for conversion to UTM data.

### Output data

*Move2* location object with appended cluster ID column *xy.clust.*

### Artefacts

*None.*

### Settings

`eps1:` Maximum cluster spatial reachability. The maximum distance (in metres, if the data is in UTM format) between two locations defined as *spatially-reachable* from one another.

`eps2:` Maximum cluster temporal reachability. The maximum time difference (in days) between two locations defined as *temporally-reachable* from one another.

`minPts:` The minimum number of connected points to constitute a *Core Point* which can form a cluster.

`Startdate:` If applicable, the start date of the clustering period. Data will be filtered to this period. Defaults to the earliest timestamp within the dataset if not provided.

`Enddate:` If applicable, the end date of the clustering period. Data will be filtered to this period. Defaults to the final timestamp within the dataset if not provided.

### Most common errors

`Data not in UTM format:` If data provided is within WGS84, however `eps1` is still provided in metres, it is likely that all locations will form one cluster. Use a separate MoveApp to convert to UTM data before this MoveApp.

### Null or error handling

`Startdate` and `Enddate` default to the first and final timestamps within the dataset, respectively.
