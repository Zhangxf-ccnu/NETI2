README file for R package supporting the paper "Joint reconstruction of multiple gene networks by simultaneously capturing inter-tumor and intra-tumor heterogeneity".


Contents of this archive
------------------------
This archive contains 
(1) pkg: subdirectory that contains the R package.
(2) NETI2-manual.pdf: reference manual.


The NETI2 package has the following R-package dependencies: MASS, Matrix,Matrix,QUIC,igraph, mvtnorm, foreach.
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
TCGA.BRCA.NETI2= NETI2(TCGA.BRCA$X,TCGA.BRCA$purity, lambda = 0.6, tau = 0.4,delta = 0.2)


Please do not hesitate to contact Dr. Xiao-Fei Zhang (zhangxf@mail.ccnu.edu.cn) or Miss Jia-Juan Tu (adaline_juan@mail.ccnu.edu.cn) to seek any clarifications regarding any  contents or operation of the archive.
