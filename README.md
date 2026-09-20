# qPCR-R
An R pipeline to automate data wrangling, calculation, and data plotting
Target genes: TNF-a and IL-6
Reference gene: GAPDH
Samples: Infected and Uninfected
Libraries: tidyverse, janitor, ggplot2

Example Dataset shows qPCR results of a class that is ungrouped. This script separates the rows into 6 separate groups accordingly for calculation.
Each group consists of triplicate data for target gene (TNF-a/IL-6) and a reference gene, GAPDH. The code finds the average of the 3 replicates of each gene of the group and the average Cq obtained from each group is taken as a replicate
- e.g. There are 3 groups with TNF-a and 3 other groups with IL-6, finding average Cq from each group and taking its average = there will be 3 replicates obtained.

### Calculation
* Filters GAPDH from the targets
* Merge variables by experiments, replicate, and sample
* Calculations are done using delta-delta-Ct method and the delta Ct, delta-delta-Ct, and Fold change results are represented in a matrix.

The relative fold change of infected and uninfected samples is also represented in a matrix and later in a bar graph by ggplot
