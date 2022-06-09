# Forest Complexity Paper
Code used for LiDAR analysis in the paper: Unraveling Forest Complexity: Resource Use Efficiency, Disturbance, and the Structure-Function Relationship

Find the paper here: https://doi.org/10.1029/2021JG006748


Notes for using this code:

allmetrics.R is a library file containing all the functions we used to run metrics on LiDAR point clouds

allmetricsLOOP.R is a file for running all the metrics on a list of several LiDAR point clouds and saving the resulting data. Point clouds must be cleaned, normalized, and in .las format before running metrics. To use this file the filenames and locations will need to be editied.

metrics_settings.csv is mostly unused now. It just needs to be a .csv with a column of the different LiDAR data set names. The other columns don't matter.
