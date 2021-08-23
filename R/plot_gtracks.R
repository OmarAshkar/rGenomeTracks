#' This is a generic function used to plot `genome_track` objects.
#' @title Plotting genomic tracks
#' @param obj Object to be plotted.
#' @param obj genome_track object. Define all tracks to be plotted.
#' @param chr String or numeric value to indicate the chromosome desire.
#' @param start Numeric. Starting position of plotting on the defined
#'   chromosome.
#' @param end Numeric. Starting position of plotting on the defined chromosome.
#' @param dir String. Default is NULL. If defined, a string to directory and
#'   extension to which image is exported. Extension could be png, svg or pdf.
#' @param plot Boolean. Default if TRUE. If FALSE, plot will not be generated,
#'   only exported.
#' @param verbose If TRUE, print command that will be passed to pyGenomeTracks.
#' @param dpi Numeric. Default is 100
#' @param title String. Title of the generated plot. Default is NULL.
#' @param fontsize If set, global fontsize value overrides individual tracks.R .
#'   argument of all tracks passed.
#' @param width Numeric. The width of the plot. Default is 40
#' @param height Numeric. Height of the plot. Default is NULL to set is based on
#'   tracks height.
#' @param trackLabelFraction Numeric. Default is 0.05.
#' @param trackLabelHAlign String. Position of labels aligment. Options are
#'   "left", "right" or "center". Default is "left".
#' @param ... Extra arguments to be passed for generic plot().
#' @export
#' @docType methods
#' @rdname plot_gtracks
#' @return None
#' @import rGenomeTracksData
#' @note For this function to run, you need pyGenomeTracks
#' installed in R's loading enviroment. If not, please run install_pyGenomeTracks()
#' @examples
#' \dontrun{
#' # Get example data directories
#' # Download h5 example
#' ah <- AnnotationHub()
#' query(ah, "rGenomeTracksData")
#' h5_dir <- ah[["AH95901"]]
#' tads_dir <- system.file("extdata", "tad_classification.bed",
#'   package = "rGenomeTracks"
#' )
#' arcs_dir <- system.file("extdata", "links2.links", package = "rGenomeTracks")
#' bw_dir <- system.file("extdata", "bigwig2_X_2.5e6_3.5e6.bw", package = "rGenomeTracks")
#' #
#' # Create HiC track from HiC matrix
#' h5 <- track_hic_matrix(
#'   file = h5_dir, depth = 250000, min_value = 5, max_value = 200,
#'   transform = "log1p", show_masked_bins = FALSE
#' )
#'
#' # Create TADS track
#' tads <- track_domains(
#'   file = tads_dir, border_color = "black",
#'   color = "none", height = 5,
#'   line_width = 5,
#'   show_data_range = FALSE,
#'   overlay_previous = "share-y"
#' )
#'
#' # Create arcs track
#' arcs <- track_links(
#'   file = arcs_dir, links_type = "triangles", line_style = "dashed",
#'   overlay_previous = "share-y",
#'   color = "darkred",
#'   line_width = 3,
#'   show_data_range = FALSE
#' )
#'
#' # Create bigwig track
#' bw <- track_bigwig(
#'   file = bw_dir, color = "red",
#'   max_value = 50,
#'   min_value = 0,
#'   height = 4,
#'   overlay_previous = "yes",
#'   show_data_range = FALSE
#' )
#'
#' # Create one object from HiC, arcs and bigwid
#' tracks <- h5 + arcs + bw
#'
#' # Plot the tracks
#' plot_gtracks(tracks, chr = "X", start = 25 * 10^5, end = 31 * 10^5)
#' # Plot HiC, TADS and bigwig tracks
#' plot_gtracks(h5 + tads + bw, chr = "X", start = 25 * 10^5, end = 31 * 10^5)
#' }
#' @keywords install_pyGenomeTracks
#' @author Omar Elashkar
setGeneric("plot_gtracks", function(obj, chr, start, end,
                                    dir = NULL,
                                    plot = TRUE,
                                    verbose = FALSE,
                                    dpi = 100,
                                    title = NULL,
                                    fontsize = NULL,
                                    width = 40,
                                    height = NULL,
                                    trackLabelFraction = 0.05,
                                    trackLabelHAlign = "left", ...) {
  standardGeneric("plot_gtracks")
})


#' @rdname plot_gtracks
#' @aliases plot_gtracks, genome_track
#' @importFrom imager load.image
#' @importFrom reticulate import
#' @return None
#' @export
#' @author Omar Elashkar
setMethod("plot_gtracks", "genome_track", function(obj,
                                                   chr, start, end,
                                                   dir = NULL,
                                                   plot = TRUE,
                                                   verbose = FALSE,
                                                   dpi = 100,
                                                   title = NULL,
                                                   fontsize = NULL,
                                                   width = 40,
                                                   height = NULL,
                                                   trackLabelFraction = 0.05,
                                                   trackLabelHAlign = "left", ...) {
  obj <- obj@tracks
  start <- format(start, scientific = FALSE, trim = TRUE)
  end <- format(end, scientific = FALSE, trim = TRUE)
  if (verbose) {
    print(start)
    print(end)
  }
  # Take all list to create config file
  conf_file <- tempfile("tracks.ini")
  if (is.null(dir)) {
    dir <- tempfile(fileext = ".png")
  }

  cat("\n", file = conf_file)

  for (i in seq_along(obj)) {
    x <- obj[[i]]
    for (ii in seq_along(x)) {
      if (ii == 1) {
        cat(paste0("[", basename(x[[ii]]), "]", "\n"), file = conf_file, append = TRUE)
      }

      if (!x[ii] %in% c("", "spacer", "x-axis")) {
        s <- paste(names(x)[ii], "=", x[[ii]])
        cat(s, "\n", file = conf_file, append = TRUE)
      }
    }
  }

  # Run command
  cmd <- paste(
    "pyGenomeTracks",
    "--tracks", conf_file, "--region",
    paste0(chr, ":", start, "-", end),
    "--outFileName", normalizePath(dir, mustWork = FALSE),
    # "--trackLabelFraction 0.2 --width 38"
    "--dpi", dpi,
    if (!is.null(title)) {
      paste("--title", title)
    },
    if (!is.null(fontsize)) {
      paste("--fontSize", fontsize)
    },
    "--width", 40,
    if (!is.null(height)) {
      paste("--height", height)
    },
    "--trackLabelFraction", trackLabelFraction,
    "--trackLabelHAlign", trackLabelHAlign
  )
  if (verbose) print(cmd)

  py <- reticulate::import("os")
  py$system(cmd)
  if (plot) {
    img <- imager::load.image(dir)
    plot(img, axes = FALSE, ...)
  }
  NULL
})