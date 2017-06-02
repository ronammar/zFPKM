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

### Libraries
# library(dplyr)
# library(ggplot2)
# library(tidyr)


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
#' @param fpkmDF A data frame of raw FPKM (or TPM) values. Each row corresponds
#'  to a gene/transcript and each column corresponds to a sample. NOTE: these
#'  are NOT log_2 transformed. Also, the rownames are gene/transcript names and
#'  NOT included as a separate column
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
#' @import checkmate dplyr ggplot2 tidyr
#'
#' @export
zFPKM<- function(fpkmDF) {

  assertDataFrame(fpkmDF)

  zFPKMDF <- data.frame(row.names=row.names(fpkmDF))
  outputs <- list()
  for (c in colnames(fpkmDF)) {
    output <- zFPKMCalc(fpkmDF[, c])
    zFPKMDF[, c] <- output$z
    outputs[[c]] <- output
  }

  return(zFPKMDF)
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
#' @param fpkmDF A data frame of raw FPKM (or TPM) values. Each row corresponds
#'  to a gene/transcript and each column corresponds to a sample. NOTE: these
#'  are NOT log_2 transformed. Also, the rownames are gene/transcript names and
#'  NOT included as a separate column
#' @param plotZFPKM Plot the scaled log_2(FPKM) and fitted Gaussian densities to screen [default = TRUE]
#' @param PlotFileName Plot the densities to specified file (.png) [default = NULL]
#' @param FacetTitles use to label each facet with the sample name [default = FALSE]
#' @param PlotXfloor Lower limit for X axis (log2FPKM units) [default = -20] set to NULL to disable
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
#' zFPKMPlot(MyFPKMdf)
#'
#' @import checkmate dplyr ggplot2 tidyr
#'
#' @export
zFPKMPlot <- function(fpkmDF, plotZFPKM=TRUE, PlotFileName=NULL,
                      FacetTitles=FALSE, PlotXfloor = -20) {

  assertDataFrame(fpkmDF)

  zFPKMDF <- data.frame(row.names=row.names(fpkmDF))
  outputs <- list()
  for (c in colnames(fpkmDF)) {
    output <- zFPKMCalc(fpkmDF[, c])
    zFPKMDF[, c] <- output$z
    outputs[[c]] <- output
  }

  if (plotZFPKM) {
    PlotGaussianFitDF(outputs, FacetTitles, PlotXfloor)
  }

  if (!is.null(PlotFileName)) {
    png(file=PlotFileName, width=8, height=8, units='in', res=300)
    PlotGaussianFitDF(outputs, FacetTitles, PlotXfloor)
    invisible(dev.off())
  }
}

zFPKMTransform <- function(fpkm) {
  # Performs the zFPKM transform on RNA-seq FPKM data.
  #
  # Args:
  #   fpkm: a vector of raw FPKM values. NOTE: these are NOT log_2 transformed.
  #
  # Returns:
  #   The zFPKM transformed data frame of the input FPKM data

  if (!is.numeric(fpkm)) {
    stop("argument 'fpkm' must be numeric")
  }

  output <- zFPKMCalc(fpkm)

  PlotGaussianFit(output)

  return(output$z)
}

zFPKMCalc <- function(fpkm) {
  # Performs the zFPKM transform on RNA-seq FPKM data. This involves fitting a
  # Gaussian distribution based on the right side of the FPKM distribution, as
  # described by Hart et al., 2013 (Pubmed ID 24215113). The zFPKM transformed
  # FPKM values represent normalized FPKM data, and, according to Hart et al.,
  # a zFPKM value >= -3 can be considered expressed (has an active promoter).
  #
  # Args:
  #   fpkm: a vector of raw FPKM values. NOTE: these are NOT log_2 transformed.
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

  # Set the maximum point in the density as the mean for the fitted Gaussian
  mu <- d$x[which.max(d$y)]

  # Determine the standard deviation
  U <- mean(fpkmLog2[fpkmLog2 > mu])
  stdev <- (U - mu) * sqrt(pi / 2)

  # Compute zFPKM transform
  zFPKM <- (fpkmLog2 - mu) / stdev

  result <- ZFPKMResult(zFPKM, d, mu, stdev)

  return(result)
}

ZFPKMResult <- function(zfpkmVector, density, mu, stdev) {
  # S3 class to store zFPKM vector and related metrics to be used in plotting

  this <- list(
    z = zfpkmVector,
    d = density,
    m = mu,
    s = stdev
  )

  class(this) <- append(class(this), "zFPKM")
  return(this)
}

PlotGaussianFit <- function(result) {
  # Plot the log_2(FPKM) density and the fitted Gaussian (scale both to the
  # same density scale so that the curves overlap).

  d <- result$d
  mu <- result$m
  stdev <- result$s

  # Get max of each density and then compute the factor to multiply the
  # fitted Gaussian. NOTE: This is for plotting purposes only, not for zFPKM
  # transform.
  fitted <- dnorm(d$x, mean=mu, sd=stdev)

  maxFPKM <- max(d$y)
  maxFitted <- max(fitted)

  scaleFitted <- fitted * (maxFPKM / maxFitted)

  df <- data.frame(log2fpkm=d$x, fpkm_density=d$y,
                   fitted_density_scaled=scaleFitted)
  ## No need to include the unscaled fitted Gaussian
  #df <- data.frame(log2fpkm=d$x, fpkm_density=d$y,
  #                  fitted_density_scaled=scaleFitted, fitted_density=fitted)
  dfg <- df %>% tidyr::gather(source, density, -log2fpkm)

  p <- ggplot2::ggplot(dfg, ggplot2::aes(x=log2fpkm, y=density, color=source)) +
    ggplot2::geom_line(size=1, alpha=0.7) +
    ggplot2::theme_bw() +
    ggplot2::labs(x="log2(FPKM)", y="[scaled] density",
                 title=paste("Check Gaussian fit (m=", round(mu, 2), " sd=",
                     round(stdev, 2), ")", sep=""))

  print(p)
}

PlotGaussianFitDF <- function(results, FacetTitles=FALSE, PlotXfloor) {
  # Plot a grid of the log_2(FPKM) density and the fitted Gaussian (scale both
  # to the same density scale so that the curves overlap).

  megaDF <- data.frame()

  for (name in names(results)) {
    result <- results[[name]]
    d <- result$d
    mu <- result$m
    stdev <- result$s

    # Get max of each density and then compute the factor to multiply the
    # fitted Gaussian. NOTE: This is for plotting purposes only, not for zFPKM
    # transform.
    fitted <- dnorm(d$x, mean=mu, sd=stdev)

    maxFPKM <- max(d$y)
    maxFitted <- max(fitted)

    scaleFitted <- fitted * (maxFPKM / maxFitted)

    df <- data.frame(sample_name=name, log2fpkm=d$x, fpkm_density=d$y,
                     fitted_density_scaled=scaleFitted)
#browser()
    megaDF <- megaDF %>% dplyr::bind_rows(df)
  }

  megaDFG <- megaDF %>% tidyr::gather(source, density, -c(log2fpkm, sample_name))

  maxX = max(megaDFG$log2fpkm)

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

  print(p)
}
