% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zfpkm.R
\name{zFPKM}
\alias{zFPKM}
\title{zFPKM Transformation}
\usage{
zFPKM(fpkmDF, assayName = "fpkm")
}
\arguments{
\item{fpkmDF}{A SummarizedExperiment or data frame containing raw FPKM (or TPM)
values. Each row corresponds to a gene/transcript and each column corresponds
to a sample.
NOTE: these are NOT log_2 transformed. Also, the rownames are gene/transcript
names and NOT included as a separate column}

\item{assayName}{When input is a SummarizedExperiment, names the specific
assay. Typically one of "fpkm" or "tpm" [default = "fpkm"]}
}
\value{
zFPKM data frame
}
\description{
Perform the zFPKM transform on RNA-seq FPKM data. This algorithm is
based on the publication by Hart et al., 2013 (Pubmed ID 24215113). Reference recommends
using zFPKM > -3 to select expressed genes.  Validated with encode open/closed promoter chromatin structure epigenetic
data on six of the ENCODE cell lines.  Works well for gene level data using FPKM or TPM. Does not appear to calibrate well
for transcript level data.
}
\examples{
library(dplyr)
gse94802 <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE94nnn/GSE94802/suppl/GSE94802_Minkina_etal_normalized_FPKM.csv.gz"
temp <- tempfile()
download.file(gse94802, temp)
fpkm <- read.csv(gzfile(temp), row.names=1)
MyFPKMdf <- select(fpkm, -MGI_Symbol)

zfpkm <- zFPKM(MyFPKMdf)

}
\references{
\url{http://www.ncbi.nlm.nih.gov/pubmed/24215113}
}
\author{
Ron Ammar, \email{ron.ammar@bms.com}
}
\keyword{zFPKM}
