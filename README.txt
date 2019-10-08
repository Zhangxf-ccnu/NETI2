README file for R package supporting the paper "Joint reconstruction of multiple gene networks by simultaneously capturing inter-tumor and intra-tumor heterogeneity".


Contents of this archive
------------------------
This archive contains 
(1) pkg: subdirectory that contains the R package.
(2) NETI2-manual.pdf: reference manual.
(3) simulation: subdirectory that contains codes for carrying out simulation studies. Run the "demo.simulation.NETI2.R" file to perform simulation studies.
(4) TCGA.BRCA: subdirectory that contains codes for applying NETI2 to TCGA breast cancer datasets. Run the "TCGA.BRAC.Application.R" file to perform real data analysis.

The NETI2 package has the following R-package dependencies: MASS, Matrix,Matrix,QUIC,igraph, mvtnorm,foreach.
The dependents are automatically installed along with NETI2. You can use the following commands to install NETI2 from GitHub. 


# Step 1. Install the devtools package. Invoke R and then type
install.packages("devtools") 

# Step 2. Load the devtools package.
library("devtools") 

# Step 3. Install the NETI2 package from GitHub.
install_github("Zhangxf-ccnu/NETI2", subdir="pkg") 


Useage
Load the library NETI2 in R console, by running 
library(NETI2)

Taking the TCGA breast cancer datasets as an example, run the following code:
data("TCGA.BRCA")
TCGA.BRCA.NETI2 = NETI2(TCGA.BRCA$BRCA.data,TCGA.BRCA$BRCA.purity, lambda = 1.052, tau = 0.4,delta = 0.2).

For more examples about simulation studies, please refer to the "simulation" subdirectory. 
For more examples about real data application, please refer to the "TCGA.BRCA" subdirectory. 

Please do not hesitate to contact Dr. Xiao-Fei Zhang (zhangxf@mail.ccnu.edu.cn) or Miss Jia-Juan Tu (adaline_juan@mail.ccnu.edu.cn) to seek any clarifications regarding any  contents or operation of the archive.
