#AUTHOR: Ron Ammar <ron.ammar@bms.com>, John R Thompson <john.thompson@bms.com>
#COMPANY: Bristol-Myers Squibb Co.
#DATE: 2015-07-08
#OBJECTIVE: Perform the zFPKM transform on RNA-seq FPKM data.
#SUMMARY:
#  Perform the zFPKM transform on RNA-seq FPKM data. This algorithm is based on
#  the publication by Hart et al., 2013 (Pubmed ID 24215113).
#CHANGELOG:
# JRT 8Sep2015:
#   Modified zFPKMTransformDF to optionally allow plotting and saving the
#   plot to a PNG file and set an output path.
#   Modified PlotGaussianFitDF to move legend to the top and optionally remove facet
#   titles on each individual plot (to allow more data to be plotted and still visually
#   interpretable).
#    Modfied to use double colon references to external libraries
# JRT 4Mar2016:
#   Added "floor" parameter to set lower limit on x axis log2FPKM value
# RA 10May2017:
#   Updated plotting options. Minor refactor. Using checkmate.
# RA 1Jun2017:
#   Separating zFPKM calc and plotting as per Bioconductor review suggestion.
# RA 10Jul2017:
#   Style changes for Bioconductor submission.
# JSL 30Jul2024
#   Added peak detection, useful for low input / (relatively dense) single cell data.

#' zFPKM Transformation
#'
#' Perform the zFPKM transform on RNA-seq FPKM data. This algorithm is
#' based on the publication by Hart et al., 2013 (Pubmed ID 24215113). Reference recommends
#' using zFPKM > -3 to select expressed genes.  Validated with encode open/closed promoter chromatin structure epigenetic
#' data on six of the ENCODE cell lines.  Works well for gene level data using FPKM or TPM. Does not appear to calibrate well
#' for transcript level data.
#'
#' @author Ron Ammar, \email{ron.ammar@bms.com}
#' @references \url{http://www.ncbi.nlm.nih.gov/pubmed/24215113}
#' @keywords zFPKM
#'
#' @param fpkmDF A SummarizedExperiment or data frame containing raw FPKM (or TPM)
#'  values. Each row corresponds to a gene/transcript and each column corresponds
#'  to a sample.
#'  NOTE: these are NOT log_2 transformed. Also, the rownames are gene/transcript
#'  names and NOT included as a separate column
#' @param assayName When input is a SummarizedExperiment, names the specific
#'  assay. Typically one of "fpkm" or "tpm" [default = "fpkm"]
#' @param peak_parameters A list containing two values for identifying peaks
#'  in the density plot: minpeakheight and minpeakdistance (from pracma::findpeaks).
#'  Default is \code{peak_parameters = list(minpeakheight = 0.02, minpeakdistance = 1)}
#'  \itemize{
#'  \item minpeakheight: the minimum (absolute) height a peak has to have to be recognized as such. Default is 0.02.
#'  \item minpeakdistance: the minimum distance (in indices) peaks have to have to be counted. Default is 1.
#'  }
#'
#' @return zFPKM data frame
#'
#' @examples
#' library(dplyr)
#' gse94802 <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE94nnn/GSE94802/suppl/GSE94802_Minkina_etal_normalized_FPKM.csv.gz"
#' temp <- tempfile()
#' download.file(gse94802, temp)
#' fpkm <- read.csv(gzfile(temp), row.names=1)
#' MyFPKMdf <- select(fpkm, -MGI_Symbol)
#'
#' zfpkm <- zFPKM(MyFPKMdf)
#'
#' # using different parameters to identify peaks
#' zfpkm <- zFPKM(MyFPKMdf, peak_parameters = list(0.01, 1))
#'
#' @import checkmate dplyr ggplot2 tidyr SummarizedExperiment pracma
#'
#' @export
zFPKM<- function(fpkmDF, assayName="fpkm", peak_parameters = list(0.02, 1)) {

  assert(checkDataFrame(fpkmDF), checkClass(fpkmDF, "SummarizedExperiment"),
         combine="or")

  message(paste("peak identification parameters are:",
        "minpeakheight =", peak_parameters[[1]], " and ",
        "minpeakdistance =", peak_parameters[[2]]))

  return(zFPKMTransform(fpkmDF, assayName, peak_parameters)[[2]])
}


removeNanInfRows <- function(fpkm) {
  # Remove FPKM rows containing all NaN values. These are most likely a result
  # of effective lengths = 0 when calculating FPKM.
  return(fpkm[which(!apply(fpkm, 1, function(r) all(is.nan(r) | is.infinite(r)))), ])
}


#' zFPKM Transformation
#'
#' Perform the zFPKM transform on RNA-seq FPKM data. This algorithm is
#' based on the publication by Hart et al., 2013 (Pubmed ID 24215113). Reference recommends
#' using zFPKM > -3 to select expressed genes.  Validated with encode open/closed promoter chromatin structure epigenetic
#' data on six of the ENCODE cell lines.  Works well for gene level data using FPKM or TPM. Does not appear to calibrate well
#' for transcript level data.
#'
#' @author Ron Ammar, \email{ron.ammar@bms.com}
#' @references \url{http://www.ncbi.nlm.nih.gov/pubmed/24215113}
#' @keywords zFPKM
#'
#' @param fpkmDF A SummarizedExperiment or data frame containing raw FPKM (or TPM)
#'  values. Each row corresponds to a gene/transcript and each column corresponds
#'  to a sample.
#'  NOTE: these are NOT log_2 transformed. Also, the rownames are gene/transcript
#'  names and NOT included as a separate column
#' @param assayName When input is a SummarizedExperiment, names the specific
#'  assay. Typically one of "fpkm" or "tpm" [default = "fpkm"]
#' @param FacetTitles use to label each facet with the sample name [default = FALSE]
#' @param PlotXfloor Lower limit for X axis (log2FPKM units) [default = -20] set to NULL to disable
#' @param peak_parameters A list containing two values for identifying peaks
#'  in the density plot: minpeakheight and minpeakdistance (from pracma::findpeaks).
#'  Default is \code{peak_parameters = list(minpeakheight = 0.02, minpeakdistance = 1)}
#'  \itemize{
#'  \item minpeakheight: the minimum (absolute) height a peak has to have to be recognized as such. Default is 0.02.
#'  \item minpeakdistance: the minimum distance (in indices) peaks have to have to be counted. Default is 1.
#'  }
#'
#' @return Displays plots of zFPKM distributions
#'
#' @examples
#' library(dplyr)
#' gse94802 <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE94nnn/GSE94802/suppl/GSE94802_Minkina_etal_normalized_FPKM.csv.gz"
#' temp <- tempfile()
#' download.file(gse94802, temp)
#' fpkm <- read.csv(gzfile(temp), row.names=1)
#' MyFPKMdf <- select(fpkm, -MGI_Symbol)
#'
#' zFPKMPlot(MyFPKMdf)
#'
#' # using different parameters to identify peaks
#' zFPKMPlot(MyFPKMdf, peak_parameters = list(0.01, 1))
#'
#' @import checkmate dplyr ggplot2 tidyr SummarizedExperiment pracma
#'
#' @export
zFPKMPlot <- function(fpkmDF, assayName="fpkm", FacetTitles=FALSE, PlotXfloor=-20, peak_parameters = list(0.02, 1)) {

  assert(checkDataFrame(fpkmDF), checkClass(fpkmDF, "SummarizedExperiment"),
         combine="or")

  p <- PlotGaussianFitDF(zFPKMTransform(fpkmDF, assayName, peak_parameters)[[1]], FacetTitles, PlotXfloor)
  return(p)
}


zFPKMTransform <- function(fpkmDF, assayName, peak_parameters = list(0.02, 1)) {
  # Helper function for zFPKM output and plotting. Do not call directly.
  #
  # Args:
  #   Defined by zFPKM() and zFPKMPlot()
  #   peak_parameters: a list containing two values
  #     1. minpeakheight (default = 0.02)
  #     2. minpeakdistance (default = 1)
  #
  # Returns:
  #   Internal objects used for zFPKM calculations and plotting.

  if (is(fpkmDF, "SummarizedExperiment")) {
    fpkmDF <- assay(fpkmDF, assayName)
  }

  fpkmDF <- removeNanInfRows(fpkmDF)

  zFPKMDF <- data.frame(row.names=row.names(fpkmDF))
  outputs <- list()
  for (c in colnames(fpkmDF)) {
    output <-
      zFPKMCalc(
        fpkmDF[, c],
        peak_parameters = peak_parameters)
    zFPKMDF[, c] <- output[["z"]]
    outputs[[c]] <- output
  }

  return(list(outputs, zFPKMDF))
}


zFPKMCalc <- function(fpkm, peak_parameters = list(0.02, 1)) {
  # Performs the zFPKM transform on RNA-seq FPKM data. This involves fitting a
  # Gaussian distribution based on the right side of the FPKM distribution, as
  # described by Hart et al., 2013 (Pubmed ID 24215113). The zFPKM transformed
  # FPKM values represent normalized FPKM data, and, according to Hart et al.,
  # a zFPKM value >= -3 can be considered expressed (has an active promoter).
  #
  # Args:
  #   fpkm: a vector of raw FPKM values. NOTE: these are NOT log_2 transformed.
  #   peak_parameters: a list containing two values
  #     1. minpeakheight (default = 0.02)
  #     2. minpeakdistance (default = 1)
  #
  # Returns:
  #   The zFPKM transformed vector of the input FPKM data

  if (!is.numeric(fpkm)) {
    stop("argument 'fpkm' must be numeric")
  }

  # log_2 transform the FPKM values
  fpkmLog2 <- log(fpkm, base=2)

  # Compute kernel density estimate
  d <- density(fpkmLog2)

  # Obtain the higher peak, in case the higher peak is lower than the peak
  # from the 'active' gene.

  peak_est <-
    pracma::findpeaks(
      d$y,
      minpeakheight = peak_parameters[[1]],
      minpeakdistance = peak_parameters[[2]])
  peak_pos <- d$x[peak_est[, 2]] # peak position
  max_x_index <- which.max(peak_pos) # peak with the highest x-value

  # Set the maximum point in the density as the mean for the fitted Gaussian
  mu <- peak_pos[max_x_index]

  # Estimate y value (i.e., maximum FPKM)
  maxFPKM <- d$y[peak_est[max_x_index, 2]]

  # Determine the standard deviation
  U <- mean(fpkmLog2[fpkmLog2 > mu])
  stdev <- (U - mu) * sqrt(pi / 2)

  # Compute zFPKM transform
  zFPKM <- (fpkmLog2 - mu) / stdev

  result <- ZFPKMResult(zFPKM, d, mu, stdev, maxFPKM)

  return(result)
}


ZFPKMResult <- function(zfpkmVector, density, mu, stdev, maxFPKM) {
  # S3 class to store zFPKM vector and related metrics to be used in plotting

  zfpkmRes <- list(
    z = zfpkmVector,
    d = density,
    m = mu,
    s = stdev,
    y = maxFPKM
  )

  class(zfpkmRes) <- append(class(zfpkmRes), "zFPKM")
  return(zfpkmRes)
}


PlotGaussianFitDF <- function(results, FacetTitles=FALSE, PlotXfloor) {
  # Plot a grid of the log_2(FPKM) density and the fitted Gaussian (scale both
  # to the same density scale so that the curves overlap).

  megaDF <- data.frame()

  for (name in names(results)) {
    result <- results[[name]]
    d <- result[["d"]]
    mu <- result[["m"]]
    stdev <- result[["s"]]
    maxFPKM <- result[["y"]]

    # Get max of each density and then compute the factor to multiply the
    # fitted Gaussian. NOTE: This is for plotting purposes only, not for zFPKM
    # transform.
    fitted <- dnorm(d[["x"]], mean=mu, sd=stdev)

    # maxFPKM <- max(d[["y"]])
    maxFitted <- max(fitted)

    scaleFitted <- fitted * (maxFPKM / maxFitted)

    df <- data.frame(sample_name=name, log2fpkm=d[["x"]], fpkm_density=d[["y"]],
                     fitted_density_scaled=scaleFitted)

    megaDF <- megaDF %>% dplyr::bind_rows(df)
  }

  megaDFG <- megaDF %>% tidyr::gather(source, density, -c(log2fpkm, sample_name))

  maxX = max(megaDFG[["log2fpkm"]])

  p <- ggplot2::ggplot(megaDFG, ggplot2::aes(x=log2fpkm, y=density, color=source)) +
    ggplot2::facet_wrap(~ sample_name) +
    ggplot2::geom_line(alpha=0.7) +
    ggplot2::theme_bw() +
    ggplot2::labs(x="log2(FPKM)", y="[scaled] density")  +
    ggplot2::theme(legend.position="top") +
    ggplot2::xlim(PlotXfloor, maxX)

  if (!FacetTitles) {  #remove the title on each facet
    p = p  + ggplot2::theme(strip.background = ggplot2::element_blank(),
                            strip.text = ggplot2::element_blank())
  }

  return(p)
}
