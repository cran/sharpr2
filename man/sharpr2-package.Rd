\name{sharpr2-package}
\alias{sharpr2-package}
\docType{package}
\title{
Estimating regularoty scores and identifying high resolution driver elements for ATAC-STARR data
}
\description{
The package develops an algorithm for identifying high-resolution driver elements for datasets from an ATAC-STARR library. 
}
\details{
\tabular{ll}{
Package: \tab sharpr2\cr
Type: \tab Package\cr
Version: \tab 1.1.1.0000 \cr
Date: \tab 2018-05-12\cr
License: \tab GPL \cr
}
}
\author{
Liang He

Maintainer: Liang He <lianghe@mit.edu>
}

\references{
High-resolution genome-wide functional dissection of transcriptional regulatory regions in human.
Xinchen Wang, Liang He, Sarah Goggin, Alham Saadat, Li Wang, Melina Claussnitzer, Manolis Kellis.
bioRxiv 193136; doi: https://doi.org/10.1101/193136
}
\keyword{ package }

\examples{
data(hidra_ex)
re <- sharpr2(hidra_ex[1:2000,], l_min = 150, l_max = 600, f_dna = 5, f_rna = 0, sig=FALSE)

}
