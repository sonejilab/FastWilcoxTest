# 0.2.0

Adding of a rgb2hex color function. This function works on integer matrixes from 0 to 255 values for r,g and b.
In addition there is a rgb2hexS function processing single colors.

# 0.1.17

The not logged FC has been added to make the FC easier to read by a human

# 0.1.15

Bugfix wilcox test to return empty matrix intead of failing.

# 0.1.14

Added a simple ZScoreAll function zscoring a sparse matrix returning a dense one.

# 0.1.13

CorMatrixIDS_N and CorMatrix_N now return the t value and the n cells expressing the gene correlated to for each cor result.


# 0.1.12

Added a function to estimate the entropy of a grouping vector based on euclidian sistance spheres around each cell. 
Based on specific drc 3D models to judge these.

# 0.1.11

Added a rolling sum function (rollSum(x,n)) and a linear correlation method for the rollSum results (CorNormalMatrix)

# 0.1.10

Added timeline functionallity to fit a timeline into 3D points.

# 0.1.9

Added a simple euclidian_distances function to allow timeline creation.

# 0.1.8

Added a ColNotZero function to quickly calculate the not zero entries in a sparse matrix. Column wise to e.g. quickly access the amount of expressed genes in a table.

# 0.1.7

Added the LinLang function to identify genes that show a linear raise over the detection limit in a given grouping.  

# 0.1.5

Added a NormalizeSamples function that simply scales the sparse matrix to a given value per sample.

# 0.1.4

Added a NormalizeCells function that mimics the seurat::RunUMISamplingPerCell function without up sampling.
In addition the function really enforces the usage of exactly nUMI reads per sample.

# 0.1.3

First almost public release.
