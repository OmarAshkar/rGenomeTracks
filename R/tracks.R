#' @description Generate genome_track object from a bed file.
#' @details `track_bed()` supports all common bed files with minimal of 
#' 3 columns and maximum of 12 columns.
#' @title Generate bed track
#' @param file String. The location of the track file
#' @param title String. If specificed, the title of the track to be displayed.
#' @param height Numeric. The height of the plotted track in cm. Default is 2.
#' See notes.
#' @param overlay_previous String. Options are "no" (default) or "yes" or
#' "share-y".
#' @param fontsize Numeric value to font size of tracks's text.
#' @param orientation String. Set to "inverted" to make the track upside down.
#'   Default is NULL.
#' @param line_width Numeric. Default is 0.5.
#' @param color String. Hex color or string color. Default is "#1f78b4".
#' @param max_value Numeric. Default is NULL. The max value cut-off for the
#'   numeric column.
#' @param min_value Numeric. Default is NULL. The max value cut-off for the
#'   numeric column.
#' @param border_color String. default is "black"
#' @param prefered_name String. Denote which column to get elements names.
#'   Default is "transcript_name".
#' @param merge_transcripts Boolean. Default is FALSE.
#' @param labels Boolean. Default is FALSE.
#' @param style String. Options are "flybase" (default), or "UCSV" or
#'   "tassarrow".
#' @param display String. options are "stacked" (default) or "collapsed",
#'   "triangles" or "interleaved".
#' @param max_labels Numeric. Any integer about 1. Default is 60.
#' @param global_max_row Boolean. Default is FALSE.
#' @param gene_rows Numeric. Default is NULL.
#' @param arrow_interval Numeric. Should be above 1. Default is 2
#' @param arrowhead_included Boolean. Default is FALSE
#' @param color_utr String. Hex color or string. Default is "grey"
#' @param height_utr Numeric. Between 0 and 1. Default is 1.
#' @param arrow_length Numeric. Default is NULL.
#' @param all_labels_inside Boolean. Default is FALSE
#' @param labels_in_margin Boolean. Default is FALSE.
#' @return genome_track
#' @export
#' @examples
#' bed12_dir <- system.file("extdata", "dm3_genes.bed.gz",
#'   package = "rGenomeTracks"
#' )
#' bed4_dir <- system.file("extdata", "dm3_genes.bed4.gz",
#'   package = "rGenomeTracks"
#' )
#' bed6_dir <- system.file("extdata", "dm3_genes.bed6.gz",
#'   package = "rGenomeTracks"
#' )
#'
#' # Create bed track using bed4 file
#' bed4 <- track_bed(
#'   file = bed4_dir, height = 3, title = "bed4", color = "cyan", ,
#'   border_color = "#9ACD32", line_width = 1.5
#' )
#'
#' # Create bed track using bed6 file
#' bed6 <- track_bed(
#'   file = bed6_dir, height = 3, title = "bed4", fontsize = 8, color = "red",
#'   border_color = "yellow", arrowhead_included = TRUE
#' )
#'
#' # Create bed track using bed12 file
#' bed12 <- track_bed(
#'   file = bed12_dir, height = 3, title = "bed12", style = "UCSC",
#'   arrow_interval = 10, fontsize = 10
#' )
#'
#' # Create a spacer track
#' space <- track_spacer(height = 1)
#' \dontrun{
#' # Plotting the tracks
#' plot_gtracks(bed4 + space + bed6 + space + bed12 + space,
#'   chr = "X", start = 300 * 10^4, end = 330 * 10^4, verbose = TRUE
#' )
#' }
#' @note
#' `fontsize` argument can be overriden by the same
#' argument in `plot_gtracks()`
#' @importFrom methods getClass is new
#' @author Omar Elashkar
track_bed <- function(file,
                      title = NULL,
                      height = 2,
                      overlay_previous = "no",
                      fontsize = 12,
                      orientation = NULL,
                      line_width = 0.5,
                      color = "#1f78b4",
                      max_value = NULL,
                      min_value = NULL,
                      border_color = "black",
                      prefered_name = "transcript_name",
                      merge_transcripts = FALSE,
                      labels = TRUE,
                      style = "flybase",
                      display = "stacked",
                      max_labels = 60,
                      global_max_row = FALSE,
                      gene_rows = NULL,
                      arrow_interval = 2,
                      arrowhead_included = FALSE,
                      color_utr = 0,
                      height_utr = 1,
                      arrow_length = 0,
                      all_labels_inside = FALSE,
                      labels_in_margin = FALSE) {
  if (!file.exists(file)) stop("This file does not exit!")

  global_max_row <- ifelse(global_max_row == TRUE, "true", "false")
  labels <- ifelse(labels == TRUE, "true", "false")
  merge_transcripts <- ifelse(merge_transcripts == TRUE, "true", "false")
  arrowhead_included <- ifelse(arrowhead_included == TRUE, "true", "false")
  all_labels_inside <- ifelse(all_labels_inside == TRUE, "true", "false")
  labels_in_margin <- ifelse(labels_in_margin == TRUE, "true", "false")

  x <- list(
    file = normalizePath(file),
    height = height,
    overlay_previous = overlay_previous,
    fontsize = fontsize,
    line_width = line_width,
    color = color,
    border_color = border_color,
    prefered_name = prefered_name,
    merge_transcripts = merge_transcripts,
    labels = labels,
    style = style,
    display = display,
    max_labels = max_labels,
    global_max_row = global_max_row,
    arrow_interval = arrow_interval,
    arrowhead_included = arrowhead_included,
    color_utr = color_utr,
    height_utr = height_utr,
    arrow_length = arrow_length,
    all_labels_inside = all_labels_inside,
    labels_in_margin = labels_in_margin,
    file_type = "bed"
  )

  if (!is.null(title)) x["title"] <- title
  if (!is.null(orientation)) x["orientation"] <- orientation
  if (!is.null(max_value)) x["max_value"] <- max_value
  if (!is.null(min_value)) x["min_value"] <- max_value
  if (!is.null(gene_rows)) x["gene_rows"] <- gene_rows

  new("genome_track", tracks = list(x))
}

#' @description Generate genome_track object from bedgraph files.
#' @details summary_method parameter can be choosen to be by "mean", "average", "max",
#' "min", "stdev", "dev", "coverage", "cov" or "sum".
#' Tranform paramter options are "no" (default) or "log", "log1p", "-log", "log2" or "log10".
#' 'log1p': transformed_values = log(1 + initial_values)
#' 'log': transformed_values = log(log_pseudocount + initial_values)
#' 'log2': transformed_values = log2(log_pseudocount + initial_values)
#' 'log10': transformed_values = log10(log_pseudocount + initial_values)
#' '-log': transformed_values = log(log_pseudocount + initial_values)
#' To compute operations on the fly on the file or between 2 bedgraph files,
#' you can tweak operation parameter, it should contains file or
#' file and second_file. It is adviced to use `nans_to_zeros = TRUE` to
#' avoid unexpected results. Example value for operation are "0.89 * file",
#' "- file", "file - second_file", "log2((1 + file) / (1 + second_file))" and
#' "max(file, second_file)"
#'
#' to add the preferred line width or point size : type = "line:lw" where lw (linewidth) is numeric value.
#' Like `type = "line:0.5"` and `type = "points:0.5"`
#'
#' By default the bedgraph is plotted at the base pair resolution.
#' This can lead to very large pdf/svg files. If plotting large regions.
#' If you want to decrase the size of your file.
#' You can either rasterize the bedgraph profile by using: `rasterize = TRUE`
#'
#' @title Generate bedgraph track
#' @inheritParams track_bed
#' @param orientation String. Default is NULL. Other option is "inverted".
#' @param color String. Hex color or string color. Default is "#1f78b4".
#' @param alpha Numeric variable between 0 and 1 to indicate level of transparancy.
#' Default is 1.
#' @param use_middle Boolean. Default is FALSE.
#' @param show_data_range Boolean. Default is TRUE.
#' @param type String. Options are "fill" (default),"line", "points".
#' @param negative_color Hex color or string to indicate color of negative values. Default is NULL.
#' @param nans_to_zeros Boolean. To convert empty values to zeros,
#' set this to TRUE. Default is FALSE.
#' @param summary_method String. summary_method applied over bin range.
#' This parameter is set to NULL. See details for options.
#' @param number_of_bins Numeric value to indicate summary method used over the bin range. Default is 700
#' @param transform String to indicate type of transformation applied.
#' Default is "no".
#' @param log_pseudocount Numeric. Default is 0.
#' @param y_axis_values String with two options "transformed" (default) or "original".
#' @param second_file Path for another file to be included in operations.
#' This parameter is not set by default.
#' @param operation Default is set to "file". See details.
#' @param grid Boolean. Default is FALSE.
#' @param rasterize Boolean. Default is FALSE.
#' @return genome_track
#' @export
#' @examples
#' bg_dir <- system.file("extdata", "GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz",
#'   package = "rGenomeTracks"
#' )
#' bed_genes_dir <- system.file("extdata", "HoxD_cluster_regulatory_regions_mm10.bed",
#'   package = "rGenomeTracks"
#' )
#'
#' bg <- track_bedgraph(bg_dir, color = "green", height = 5, max_value = 10)
#' bg_middle <- track_bedgraph(bg_dir,
#'   use_middle = TRUE, color = "blue",
#'   height = 5, max_value = 10
#' )
#' bed_genes <- track_bed(bed_genes_dir,
#'   title = "Regulatory regions", ,
#'   color = "red", height = 3
#' )
#'
#' tracks <- track_x_axis(where = "top") + bg + bg_middle + bed_genes
#' \dontrun{
#' plot_gtracks(tracks,
#'   chr = 2, start = 738 * 10^5, end = 750 * 10^5,
#'   trackLabelFraction = 0.2
#' )
#' }
#' @note
#' `fontsize` parameter can be overriden by the same argument in `plot_gtracks()`
#' `height` parameter will be ignored if `overlay_previous` is set.
#' @importFrom methods getClass is new
#' @author Omar Elashkar
track_bedgraph <- function(file, title = NULL,
                           height = 2,
                           overlay_previous = "no",
                           orientation = NULL,
                           color = "#1f78b4",
                           alpha = 1,
                           max_value = NULL,
                           min_value = NULL,
                           use_middle = FALSE,
                           show_data_range = TRUE,
                           type = "fill",
                           negative_color = NULL,
                           nans_to_zeros = FALSE,
                           summary_method = NULL,
                           number_of_bins = 700,
                           transform = "no",
                           log_pseudocount = 0,
                           y_axis_values = "transformed",
                           second_file = NULL,
                           operation = "file",
                           grid = FALSE,
                           rasterize = FALSE) {
  if (!file.exists(file)) stop("This file does not exit!")

  use_middle <- ifelse(use_middle == TRUE, "true", "false")
  show_data_range <- ifelse(show_data_range == TRUE, "true", "false")
  nans_to_zeros <- ifelse(nans_to_zeros == TRUE, "true", "false")
  grid <- ifelse(grid == TRUE, "true", "false")
  rasterize <- ifelse(rasterize == TRUE, "true", "false")

  x <- list(
    file = normalizePath(file),
    height = height,
    overlay_previous = overlay_previous,
    color = color,
    alpha = alpha,
    use_middle = use_middle,
    show_data_range = show_data_range,
    type = type,
    nans_to_zeros = nans_to_zeros,
    number_of_bins = number_of_bins,
    transform = transform,
    log_pseudocount = log_pseudocount,
    y_axis_values = y_axis_values,
    operation = operation,
    grid = grid,
    rasterize = rasterize,
    file_type = "bedgraph"
  )

  if (!is.null(title)) x["title"] <- title
  if (!is.null(orientation)) x["orientation"] <- orientation
  if (!is.null(max_value)) x["max_value"] <- max_value
  if (!is.null(min_value)) x["min_value"] <- min_value
  if (!is.null(second_file)) x["second_file"] <- second_file
  if (!is.null(summary_method)) x["summary_method"] <- summary_method
  if (!is.null(negative_color)) x["negative_color"] <- negative_color

  new("genome_track", tracks = list(x))
}




#' A track for file like bedgraph but with more than 4 columns, like the insulation score from hicPlotTADs
#'
#' The different options for color maps can be found here: https://matplotlib.org/users/colormaps.html.
#' @title Generate bedgraph matrix track
#' @param file String. The location of the track file
#' @param title String. If specificed, the title of the track to be displayed.
#' @param height Numeric. The height of the plotted track in cm.
# Default is 2. See notes.
#' @param overlay_previous String. Options are "no" (default) or "yes" or
# "share-y".
#' @param orientation String. Set to "inverted" to make the track upside down. Default is NULL.
#' @param max_value  Numeric. Default is NULL.
#' The max value cut-off for the numeric column.
#' @param min_value  Numeric. Default is NULL.
#' The min value cut-off for the numeric column.
#' @param show_data_range Boolean. Default is FALSE.
#' @param type "matrix" (default) or "lines".
#' @param rasterize Boolean. Default is TRUE
#' @param pos_score_in_bin String value to indicate the position of score with respect to bin start and end.
#' Possible values are either "center" (default) or "block".
#' @param plot_horizontal_lines Boolean. Can be used only if type parameter is set to "lines".
# Default is FALSE.
#' @param colormap String with matplotlib-compatible colormap. Default is set to "viridis".
#' @return genome_track
#' @export
#' @note
#' `fontsize` argument can be overriden by the same argument in `plot_gtracks()`
#' @importFrom methods getClass is new
#' @author Omar Elashkar
track_bedgraph_matrix <- function(file,
                                  title = NULL,
                                  height = 2,
                                  overlay_previous = "no",
                                  orientation = NULL,
                                  max_value = NULL,
                                  min_value = NULL,
                                  show_data_range = FALSE,
                                  type = "matrix",
                                  rasterize = TRUE,
                                  pos_score_in_bin = "center",
                                  plot_horizontal_lines = FALSE,
                                  colormap = "virdis") {
  if (!file.exists(file)) stop("This file does not exit!")

  show_data_range <- ifelse(show_data_range == TRUE, "true", "false")
  rasterize <- ifelse(rasterize == TRUE, "true", "false")
  plot_horizontal_lines <- ifelse(plot_horizontal_lines == TRUE, "true", "false")


  x <- list(
    file = normalizePath(file),
    height = height,
    overlay_previous = overlay_previous,
    show_data_range = show_data_range,
    type = type,
    rasterize = rasterize,
    pos_score_in_bin = pos_score_in_bin,
    plot_horizontal_lines = plot_horizontal_lines,
    colormap = colormap
  )

  if (!is.null(title)) x["title"] <- title
  if (!is.null(orientation)) x["orientation"] <- orientation
  if (!is.null(max_value)) x["max_value"] <- max_value
  if (!is.null(min_value)) x["min_value"] <- min_value
  new("genome_track", tracks = list(x))
}



#' @description Create genome_track object from bigwig file.
#' @details summary_method parameter can be choosen to be by "mean", "average", "max",
#' "min", "stdev", "dev", "coverage", "cov" or "sum".
#' Tranform paramter options are "no" (default) or "log", "log1p", "-log", "log2" or "log10".
#' 'log1p': transformed_values = log(1 + initial_values)
#' 'log': transformed_values = log(log_pseudocount + initial_values)
#' 'log2': transformed_values = log2(log_pseudocount + initial_values)
#' 'log10': transformed_values = log10(log_pseudocount + initial_values)
#' '-log': transformed_values = log(log_pseudocount + initial_values)
#' To compute operations on the fly on the file or between 2 bedgraph files,
#' you can tweak operation parameter, it should contains file or
#' file and second_file. It is adviced to use `nans_to_zeros = TRUE` to
#' avoid unexpected results. Example value for operation are "0.89 * file",
#' "- file", "file - second_file", "log2((1 + file) / (1 + second_file))" and
#' "max(file, second_file)"
#'
#' # to add the preferred line width or point size : type = "line:lw" where lw (linewidth) is numeric value.
#' Like `type = "line:0.5"` and `type = "points:0.5"`
#' @title Generate bigwig track
#' @param file String. The location of the track file
#' @param title String. If specificed, the title of the track to be displayed.
#' @param height Numeric. The height of the plotted track in cm.
# Default is 2.
#' @param overlay_previous String. Options are "no" (default) or "yes" or
# "share-y".
#' @param orientation String. Set to "inverted" to make the track upside down. Default is NULL.
#' @param color String. Hex color or string color. Default is "#1f78b4".
#' @param alpha Numeric variable between 0 and 1 to indicate level of transparancy.
#' Default is 1.
#' @param max_value Numeric. Default is NULL.
#' The max value cut-off for the numeric column.
#' @param min_value Numeric. Default is NULL.
#' The max value cut-off for the numeric column.
#' @param show_data_range Boolean. Default is TRUE.
#' @param type String. Options are "fill" (default),"line", "points".
#' @param negative_color Hex color or string to indicate color of negative values. Default is NULL.
#' @param nans_to_zeros Boolean. To convert empty values to zeros,
#' set this to TRUE. Default is FALSE.
#' @param summary_method String. summary_method applied over bin range.
#' This parameter is set to NULL. See details for options.
#' @param number_of_bins Numeric value to indicate summary method used over the bin range. Default is 700
#' @param transform String to indicate type of transformation applied.
#' Default is "no".
#' @param log_pseudocount Numeric. Default is 0.
#' @param y_axis_values String with two options "transformed" (default) or "original".
#' @param second_file Path for another file to be included in operations.
#' This parameter is not set by default.
#' @param operation Default is set to "file". See details.
#' @param grid Boolean. Default is FALSE.
#' @return None
#' @export
#' @examples
#' bw_dir <- system.file("extdata", "bigwig2_X_2.5e6_3.5e6.bw",
#'   package = "rGenomeTracks"
#' )
#' mean_bw <- track_bigwig(
#'   file = bw_dir, color = "gray",
#'   type = "point:1", summary_method = "mean", number_of_bins = 300, max_value = 200, min_value = -5
#' )
#' min_bw <- track_bigwig(
#'   file = bw_dir, color = "blue", type = "line:1", summary_method = "min", number_of_bins = 300,
#'   overlay_previous = "share-y", show_data_range = FALSE,
#'   max_value = 200, min_value = -5
#' )
#' max_bw <- track_bigwig(
#'   file = bw_dir, color = "red", type = "line:1", summary_method = "max", number_of_bins = 300,
#'   overlay_previous = "share-y", show_data_range = FALSE,
#'   max_value = 200, min_value = -5
#' )
#' hlines <- track_hlines(
#'   y_values = "10, 150",
#'   overlay_previous = "share-y",
#'   color = "blue", line_style = "dotted"
#' )
#' \dontrun{
#' plot_gtracks(mean_bw + min_bw + max_bw + hlines, chr = "X", start = 27 * 10^5, end = 31 * 10^5)
#' }
#' @author Omar Elashkar
track_bigwig <- function(file,
                         title = NULL,
                         height = 2,
                         overlay_previous = "no",
                         orientation = NULL,
                         color = "#1f78b4",
                         alpha = 1,
                         max_value = NULL,
                         min_value = NULL,
                         show_data_range = TRUE,
                         type = "fill",
                         negative_color = NULL,
                         nans_to_zeros = FALSE,
                         summary_method = "mean",
                         number_of_bins = 700,
                         transform = "no",
                         log_pseudocount = 0,
                         y_axis_values = "transformed",
                         second_file = NULL,
                         operation = "file",
                         grid = FALSE) {
  if (!file.exists(file)) stop("This file does not exit!")

  show_data_range <- ifelse(show_data_range == TRUE, "true", "false")
  nans_to_zeros <- ifelse(nans_to_zeros == TRUE, "true", "false")
  grid <- ifelse(grid == TRUE, "true", "false")

  x <- list(
    file = normalizePath(file),
    height = height,
    overlay_previous = overlay_previous,
    color = color,
    alpha = alpha,
    show_data_range = show_data_range,
    type = type,
    nans_to_zeros = nans_to_zeros,
    summary_method = summary_method,
    number_of_bins = number_of_bins,
    transform = transform,
    log_pseudocount = log_pseudocount,
    y_axis_values = y_axis_values,
    operation = operation,
    grid = grid,
    file_type = "bigwig"
  )


  if (!is.null(title)) x["title"] <- title
  if (!is.null(orientation)) x["orientation"] <- orientation
  if (!is.null(max_value)) x["max_value"] <- max_value
  if (!is.null(min_value)) x["min_value"] <- min_value
  if (!is.null(negative_color)) x["negative_color"] <- negative_color
  if (!is.null(second_file)) x["second_file"] <- second_file
  new("genome_track", tracks = list(x))
}


#' @description Domain files are bed files represents TADS in the case of HiC analysis.
#' @details To remove the border, set 'border_color' parameter to "none".
#' @title Generate domains track
#' @inheritParams track_bed
#' @inheritParams track_bedgraph
#' @return genome_track
#' @export
#' @examples
#' tads_dir <- system.file("extdata", "tad_classification.bed",
#'   package = "rGenomeTracks"
#' )
#' tads <- track_domains(
#'   file = tads_dir, border_color = "black",
#'   color = "#11FF34", height = 5
#' )
#' tads_i <- track_domains(
#'   file = tads_dir, border_color = "red",
#'   color = "#cccccc", height = 3, orientation = "inverted"
#' )
#' tracks <- track_x_axis(where = "top") +
#'   tads + tads_i
#'
#' \dontrun{
#' plot_gtracks(tracks, chr = "X", start = 30 * 10^5, end = 35 * 10^5)
#' }
#' @importFrom methods getClass is new
#' @author Omar Elashkar
track_domains <- function(file,
                          title = NULL,
                          height = 2,
                          overlay_previous = "no",
                          orientation = NULL,
                          line_width = 0.5,
                          color = "#1f78b4",
                          max_value = NULL,
                          show_data_range = TRUE,
                          min_value = NULL,
                          border_color = "black",
                          prefered_name = "transcript_name",
                          merge_transcripts = FALSE) {
  if (!file.exists(file)) stop("This file does not exit!")
  show_data_range <- ifelse(show_data_range == TRUE, "true", "false")
  merge_transcripts <- ifelse(merge_transcripts == TRUE, "true", "false")

  x <- list(
    file = normalizePath(file),
    height = height,
    overlay_previous = overlay_previous,
    line_width = line_width,
    show_data_range = show_data_range,
    color = color,
    border_color = border_color,
    prefered_name = prefered_name,
    merge_transcripts = merge_transcripts,
    file_type = "domains"
  )

  if (!is.null(title)) x["title"] <- title
  if (!is.null(orientation)) x["orientation"] <- orientation
  if (!is.null(max_value)) x["max_value"] <- max_value
  if (!is.null(min_value)) x["min_value"] <- min_value

  new("genome_track", tracks = list(x))
}


#' @description A convience function to generate epilogo json configuration file to be passed for `epi_logos()`
#' @details The only argument passed to this function is data.frame or data.frame similar object.
#' It should have 3 column: First is the state number of epilogos.
#' The second is the label of the state. Finally, the desired colored of such state.
#' Check the example provided for the structure of this data.frame.
#' @title Generate epilogo json configuration file
#' @param cat_df Dataframe with 3 columns of categories, names and colors
#' @return Directory
#' @export
#' @inherit track_epilogos examples
#' @export
#' @keywords track_epilogos
#' @author Omar Elashkar
epilogos_json <- function(cat_df) {
  if (ncol(cat_df) != 3) stop("Object passed has less or more than 3 columns.
Please check documentation.")

  jfile <- tempfile(fileext = ".json", pattern = "epilogo_config_")
  cat('{"categories":{', file = jfile)
  for (i in seq_len(nrow(cat_df))) {
    if (i != nrow(cat_df)) {
      cat(paste0(
        '"', cat_df[i, 1],
        '"', ":[", '"', cat_df[i, 2],
        '"', ",", '"', cat_df[i, 3], '"', "],"
      ), file = jfile, append = TRUE)
    } else {
      cat(paste0(
        '"', cat_df[i, 1],
        '"', ":[", '"', cat_df[i, 2],
        '"', ",", '"', cat_df[i, 3], '"', "]"
      ), file = jfile, append = TRUE)
    }
  }
  cat("}}", file = jfile, append = TRUE)
  jfile
}

#' @description Generate epilogos genome_track from qcat file.
#' @details Epilogos is used widely to represent multiple "states" across genome,
#' like ChromHMM states. More details \href{https://epilogos.altius.org/}{here}
#' `qcat` file is needed which can be generated using \href{https://github.com/Altius/epilogos}{epilogos}
#' `track_epiolog` can optionally take categories_file parameter which specify the color scheme for
#' the states present in `qcat` file. Check the example section for demonestration.
#' @title Generate epilogos track
#' @inheritParams track_bed
#' @param categories_file Optionally pass a string of JSON custom colors configuration file directory. Default is NULL.
#' @return None
#' @export
#' @examples
#' epilog_dir <- system.file("extdata", "epilog.qcat.bgz", package = "rGenomeTracks")
#' epi_cat <- data.frame(
#'   category = 1:15,
#'   label = c(
#'     "Active TSS",
#'     "Flanking Active TSS",
#'     "Transcr at gene 5 and 3",
#'     "Strong transcription",
#'     "Weak transcription",
#'     "Genic enhancers",
#'     "Enhancers",
#'     "ZNF genes & repeats",
#'     "Heterochromatin",
#'     "Bivalent/Poised TSS",
#'     "Flanking Bivalent TSS/Enh",
#'     "Bivalent Enhancer",
#'     "Repressed PolyComb",
#'     "Weak Repressed PolyComb",
#'     "Quiescent/Low"
#'   ),
#'   color = c(
#'     "#ff0000", "#ff4500", "#32cd32", "#008000",
#'     "#006400", "#c2e105", "#ffff00", "#66cdaa",
#'     "#8a91d0", "#cd5c5c", "#e9967a", "#bdb76b",
#'     "#808080", "#c0c0c0", "#ffffff"
#'   )
#' )
#' epilog <- track_epilogos(file = epilog_dir, categories_file = epilogos_json(epi_cat))
#' \dontrun{
#' plot_gtracks(epilog, chr = "X", start = 3100000, 3150000)
#' }
#' @note
#' `fontsize` argument can be overriden by the same argument in `plot_gtracks()`
#' @keywords epilogos_json
#' @importFrom methods getClass is new
#' @author Omar Elashkar
track_epilogos <- function(file,
                           title = NULL,
                           height = 2,
                           overlay_previous = "no",
                           categories_file = NULL,
                           orientation = NULL) {
  x <- list(
    file = normalizePath(file),
    height = height,
    overlay_previous = overlay_previous
  )

  if (!is.null(title)) x["title"] <- title
  if (!is.null(categories_file)) x["categories_file"] <- categories_file
  if (!is.null(orientation)) x["orientation"] <- orientation

  new("genome_track", tracks = list(x))
}


#' @description Create genome_track object for gtf annotation files.
#' @details gtf files, unlike bed file, can provide richer annotation regarding levels of
#' annotation where genomic features can be grouped based on the composing entity.
#' @title Generate gtf track
#' @inheritParams track_bed
#' @inheritParams track_bedgraph
#' @return genome_track
#' @export
#' @examples
#' gtf_dir <- system.file("extdata", "dm3_subset_BDGP5.78.gtf.gz",
#'   package = "rGenomeTracks"
#' )
#' gtf <- track_gtf(
#'   file = gtf_dir, height = 10,
#'   prefered_name = "gene_name", merge_transcripts = TRUE, fontsize = 12
#' )
#' \dontrun{
#' plot_gtracks(gtf + track_spacer() +
#'   track_x_axis(), chr = "X", start = 30 * 10^5, end = 33 * 10^5)
#' }
#' @note `fontsize` argument can be overriden by the same argument in `plot_gtracks()`
#' @importFrom methods getClass is new
#' @author Omar Elashkar
track_gtf <- function(file,
                      title = NULL,
                      height = 2,
                      overlay_previous = "no",
                      fontsize = 12,
                      orientation = NULL,
                      line_width = 0.5,
                      color = "#1f78b4",
                      border_color = "black",
                      prefered_name = "transcript_name",
                      merge_transcripts = FALSE,
                      labels = FALSE,
                      display = "stacked",
                      max_labels = 60,
                      global_max_row = FALSE,
                      gene_rows = NULL,
                      arrow_interval = 2,
                      arrowhead_included = FALSE,
                      color_utr = "grey",
                      height_utr = 1,
                      arrow_length = NULL,
                      all_labels_inside = FALSE,
                      labels_in_margin = FALSE) {
  merge_transcripts <- ifelse(merge_transcripts == TRUE, "true", "false")
  labels <- ifelse(labels == TRUE, "true", "false")
  global_max_row <- ifelse(global_max_row == TRUE, "true", "false")
  arrowhead_included <- ifelse(arrowhead_included == TRUE, "true", "false")
  all_labels_inside <- ifelse(all_labels_inside == TRUE, "true", "false")
  labels_in_margin <- ifelse(labels_in_margin == TRUE, "true", "false")

  x <- list(
    file = normalizePath(file),
    height = height,
    overlay_previous = overlay_previous,
    fontsize = fontsize,
    line_width = line_width,
    color = color,
    border_color = border_color,
    prefered_name = prefered_name,
    merge_transcripts = merge_transcripts,
    labels = labels,
    display = display,
    max_labels = max_labels,
    global_max_row = global_max_row,
    arrow_interval = arrow_interval,
    arrowhead_included = arrowhead_included,
    color_utr = color_utr,
    height_utr = height_utr,
    all_labels_inside = all_labels_inside,
    labels_in_margin = labels_in_margin
  )

  if (!is.null(title)) x["title"] <- title
  if (!is.null(orientation)) x["orientation"] <- orientation
  if (!is.null(gene_rows)) x["gene_rows"] <- gene_rows
  if (!is.null(arrow_length)) x["arrow_length"] <- arrow_length

  new("genome_track", tracks = list(x))
}


#' @description Create a genome_track for matrix files. Currently, only cool format and h5 format.
#' @details
#' This function expect cool or h5 format. Format converter like
#' \href{https://hicexplorer.readthedocs.io/en/latest/content/tools/hicConvertFormat.html#hicconvertformat}{hicConvertFormat}
#' can help converting to supported formats.
#' depth is the maximum distance that should be plotted.
#' If it is more than 125% of the plotted region, it will be adjsted to this maximum value.
#' colormap argument should be compatible with \href{https://matplotlib.org/2.0.2/users/colormaps.html}{matplotlib}.
#' show_masked_bins plots bins not used during the corrections as white lines.
#' Setting this argument to FALSE (default) extends neighboring bins to
#' obtain an aesthetically pleasant output.
#' scale argument scales the matrix by specific factor.
#' This is useful if plotting multiple hic-matrices to be on the same scale.
#' @title Generate HiC track
#' @inheritParams track_bed
#' @inheritParams track_bedgraph
#' @inheritParams track_bedgraph_matrix
#' @param depth Numeric value above 1 to indicate the maximum distance that should be plotted. Default is 100000.
#' @param show_masked_bins Boolean. If TRUE, showing masked bins as white lines. Default is FALSE.
# those bins that were not used during the correction
#' @param scale_factor Numeric factor by which matrix is to be scaled.
#' @return genom_track
#' @export
#' @inherit plot_gtracks examples
#' @importFrom methods getClass is new
#' @author Omar Elashkar
track_hic_matrix <- function(file,
                             title = NULL,
                             height = NULL,
                             overlay_previous = "no",
                             orientation = NULL,
                             max_value = NULL,
                             min_value = NULL,
                             transform = "no",
                             rasterize = TRUE,
                             colormap = "RdYlBu_r",
                             depth = 100000,
                             show_masked_bins = FALSE,
                             scale_factor = 1) {
  rasterize <- ifelse(rasterize == TRUE, "true", "false")
  show_masked_bins <- ifelse(show_masked_bins == TRUE, "true", "false")

  x <- list(
    file = normalizePath(file),
    overlay_previous = overlay_previous,
    transform = transform,
    rasterize = rasterize,
    colormap = colormap,
    depth = depth,
    show_masked_bins = show_masked_bins,
    scale_factor = scale_factor,
    file_type = "hic_matrix"
  )

  if (!is.null(title)) x["title"] <- title
  if (!is.null(height)) x["height"] <- height
  if (!is.null(orientation)) x["orientation"] <- orientation
  if (!is.null(max_value)) x["max_value"] <- max_value
  if (!is.null(min_value)) x["min_value"] <- min_value
  new("genome_track", tracks = list(x))
}


#' @description `track_hlines()` creates a genome_track with horizonal lines that
#' can be overlayed on the previous track or, by default, track the lines in separate track.
#' @details y_values argument specify locations on the genome where where
#' horizontal lines should be plotted separated by comma, like "50, 90"
#' @title Generate a track with horizontal lines
#' @inheritParams track_bedgraph
#' @param y_values String for y-values where horizontal lines should be plotted separated by comma.
#' @param line_width Numeric value for line width.
#' @param line_style String with options of either "solid", "dashed", "dotted", and "dashdot".
#' @return genome_track
#' @export
#' @inherit track_bigwig examples
#' @importFrom methods getClass is new
#' @author Omar Elashkar
track_hlines <- function(y_values,
                         title = NULL,
                         height = 0.5,
                         overlay_previous = NULL,
                         orientation = NULL,
                         line_width = 0.5,
                         line_style = "solid",
                         color = "black",
                         alpha = 1,
                         max_value = NULL,
                         min_value = NULL,
                         show_data_range = TRUE) {
  show_data_range <- ifelse(show_data_range == TRUE, "true", "false")

  x <- list(
    y_values = y_values,
    height = height,
    line_width = line_width,
    line_style = line_style,
    color = color,
    alpha = alpha,
    show_data_range = show_data_range,
    file_type = "hlines"
  )

  if (!is.null(title)) x["title"] <- title
  if (!is.null(overlay_previous)) x["overlay_previous"] <- overlay_previous
  if (!is.null(orientation)) x["orientation"] <- orientation
  if (!is.null(max_value)) x["max_value"] <- max_value
  if (!is.null(min_value)) x["min_value"] <- min_value

  new("genome_track", tracks = list(x))
}


#' @description Generate links track from arc file.
#' @details
#' Level of compactness relative to arcs' length can be manipulated using the argument compact_arcs_level where:
#' \itemize{
#' \item{compact_arcs_level = 0, The default where the height is proportional to distance}
#' \item{compact_arcs_level = 1, the height is proportional to the square root of the distance}
#' \item{compact_arcs_level = 2, the height is the same for all distances}
#' }
#' ylim argument sets the cutoff for arcs' height.
#' This could be handy if you have small arc overridden by larger arc.
#' @title Generate links track
#' @inheritParams track_bedgraph
#' @inheritParams track_hlines
#' @param links_type String value with options "arcs" (default) or "triangles" or "loops".
#' @param ylim Numeric value above 0 to set arcs' height cutoff. Default is NULL
#' @param compact_arcs_level Numeric value of either 0, 1 or 2
#' to indicate level of arcs' compactness by distance it travels.
#' @return genome_track
#' @export
#' @note ylim argument is incompatible with compact_arcs_level = 2
#' @examples
#' tads_dir <- system.file("extdata", "tad_classification.bed",
#'   package = "rGenomeTracks"
#' )
#' genes_dir <- system.file("extdata", "dm3_genes.bed.gz",
#'   package = "rGenomeTracks"
#' )
#' links_dir <- system.file("extdata", "test.arcs",
#'   package = "rGenomeTracks"
#' )
#' tads <- track_domains(tads_dir, color = "#cccccc", border_color = "red")
#' links_overlay <- track_links(links_dir,
#'   color = "red",
#'   line_width = 3, links_type = "loop",
#'   overlay_previous = "share-y"
#' )
#' links <- track_links(links_dir,
#'   color = "blue",
#'   line_width = 3, height = 3
#' )
#' genes <- track_bed(genes_dir,
#'   height = 7, style = "flybase",
#'   fontsize = 10
#' )
#' \dontrun{
#' plot_gtracks(tads + links_overlay + links + genes, chr = "X", start = 30 * 10^5, end = 35 * 10^5)
#' }
#' @importFrom methods getClass is new
#' @author Omar Elashkar
track_links <- function(file,
                        title = NULL,
                        height = 2,
                        overlay_previous = "no",
                        orientation = NULL,
                        links_type = "arcs",
                        line_width = NULL,
                        line_style = "solid",
                        color = "blue",
                        alpha = 0.8,
                        max_value = NULL,
                        min_value = NULL,
                        ylim = NULL,
                        show_data_range = FALSE,
                        compact_arcs_level = 0,
                        use_middle = FALSE) {
  show_data_range <- ifelse(show_data_range == TRUE, "true", "false")
  use_middle <- ifelse(use_middle == TRUE, "true", "false")
  x <- list(
    file = normalizePath(file),
    height = height,
    overlay_previous = overlay_previous,
    links_type = links_type,
    line_style = line_style,
    color = color,
    alpha = alpha,
    compact_arcs_level = compact_arcs_level,
    use_middle = use_middle
  )

  if (!is.null(title)) x["title"] <- title
  if (!is.null(orientation)) x["orientation"] <- orientation
  if (!is.null(line_width)) x["line_width"] <- line_width
  if (!is.null(max_value)) x["max_value"] <- max_value
  if (!is.null(min_value)) x["min_value"] <- min_value
  if (!is.null(ylim)) x["ylim"] <- ylim

  new("genome_track", tracks = list(x))
}


#' @description Create genome_track object from narrow peak bed format.
#' @details
#' narrowPeak file is bed file (4+3), where the 5th column is peak name, 6th column in p-value and 7th column in q-value. You might increase height it increased font size.
#' narrowPeak format is very common with analysis pipelines involving MACS2.
#' narrowPeak format provides the information of the peak summit.
#' use_summit argument is used to deterimine if this information should be used.
#' By default this information is used (use_summit = TRUE) although some peaks may look crooked.
#' type argument specify if the plot will be:
#' \itemize{
#' \item{"box" which will plot a rectangle of the peak width}
#' \item{or "peak" which will plot the shape of the peak, whose height is the
#' narrowPeak file signal value (usually peak coverage) }
#' }
#' @title Generate narrow peaks track
#' @inheritParams track_bedgraph
#' @inheritParams track_links
#' @param use_summit  Boolean. If TRUE, peak summit data will be plotted.
#' @param width_adjust Numeric value above 0 to adjust peaks' width. Default is 1.5.
#' @param type  String with options either "peak" or "box".
#' @param show_labels Boolean. If TRUE, display labels on plotting which include peak tag, p-val and q-val.
#' @return genome_track
#' @export
#' @examples
#' np_bed_dir <- system.file("extdata", "test2.narrowPeak", package = "rGenomeTracks")
#'
#' tracks <-
#'   track_scalebar() +
#'   track_narrow_peak(np_bed_dir,
#'     title = "peak type with summit",
#'     height = 3,
#'     type = "peak",
#'     color = "green"
#'   ) +
#'
#'   track_spacer(height = 2) +
#'   track_narrow_peak(np_bed_dir,
#'     title = "peak type without summit",
#'     height = 3,
#'     type = "peak",
#'     color = "green",
#'     use_summit = FALSE
#'   ) +
#'   track_spacer(height = 2) +
#'   track_narrow_peak(np_bed_dir,
#'     title = "Box type with summit",
#'     height = 3,
#'     type = "box",
#'     color = "blue"
#'   ) +
#'   track_spacer(height = 2) +
#'   track_narrow_peak(np_bed_dir,
#'     title = "Box type without summit",
#'     height = 3,
#'     type = "box",
#'     color = "blue",
#'     use_summit = FALSE
#'   ) +
#'   track_x_axis()
#' \dontrun{
#' plot_gtracks(tracks, chr = "X", start = 276 * 10^4, end = 280 * 10^4, trackLabelFraction = 0.2)
#' }
#' @importFrom methods getClass is new
#' @author Omar Elashkar
track_narrow_peak <- function(file,
                              title = NULL,
                              height = 3,
                              overlay_previous = "no",
                              orientation = NULL,
                              line_width = 1,
                              color = "#FF000080",
                              max_value = NULL,
                              show_data_range = TRUE,
                              show_labels = TRUE,
                              use_summit = TRUE,
                              width_adjust = 1.5,
                              type = "peak") {
  show_data_range <- ifelse(show_data_range == TRUE, "true", "false")
  show_labels <- ifelse(show_labels == TRUE, "true", "false")
  use_summit <- ifelse(use_summit == TRUE, "true", "false")

  x <- list(
    file = normalizePath(file),
    height = height,
    overlay_previous = overlay_previous,
    line_width = line_width,
    color = color,
    show_data_range = show_data_range,
    show_labels = show_labels,
    use_summit = use_summit,
    width_adjust = width_adjust,
    type = type,
    file_type = "narrow_peak"
  )

  if (!is.null(title)) x["title"] <- title
  if (!is.null(orientation)) x["orientation"] <- orientation
  if (!is.null(max_value)) x["max_value"] <- max_value

  new("genome_track", tracks = list(x))
}


#' @description scalebar track is a track with a stretch that highlights
#' specific distance on the genomic coordiantes
#' @title Generate scalebar track
#' @inheritParams track_bed
#' @inheritParams track_bedgraph
#' @param where "left" (default), "right", "top" or "bottom".
#' @param line_width 0.5 (default) or any float above 0.
#' @param x_center Numeric value above 0. Default is NULL.
#' @param size Numeric value above 0. Default is NULL.
#' @param scalebar_start_position Numeric value above 0. Default is NULL.
#' @param scalebar_end_position Numeric value above 0. Default is NULL.
#' @return genome_track
#' @export
#' @examples
#' np_bed_dir <- system.file("extdata", "test2.narrowPeak", package = "rGenomeTracks")
#'
#' tracks <-
#'   track_scalebar(
#'     scalebar_start_position = 2785 * 10^3,
#'     scalebar_end_position = 2799 * 10^3
#'   ) +
#'   track_narrow_peak(np_bed_dir,
#'     title = "peak type with summit",
#'     height = 3,
#'     type = "peak",
#'     color = "green"
#'   ) + track_x_axis()
#' \dontrun{
#' plot_gtracks(tracks, chr = "X", start = 276 * 10^4, end = 280 * 10^4, trackLabelFraction = 0.2)
#' }
#' @note
#' `fontsize` argument can be overriden by the same argument in `plot_gtracks()`
#' @importFrom methods getClass is new
#' @author Omar Elashkar
track_scalebar <- function(title = NULL,
                           height = 2,
                           overlay_previous = "no",
                           where = "left",
                           fontsize = 12,
                           line_width = 0.5,
                           color = "black",
                           alpha = 1,
                           x_center = NULL,
                           size = NULL,
                           scalebar_start_position = NULL,
                           scalebar_end_position = NULL) {
  x <- list(
    file_type = "scalebar",
    height = height,
    overlay_previous = overlay_previous,
    where = where,
    fontsize = fontsize,
    line_width = line_width,
    color = color,
    alpha = alpha
  )

  if (!is.null(title)) x["title"] <- title
  if (!is.null(x_center)) x["x_center"] <- x_center
  if (!is.null(size)) x["size"] <- size
  if (!is.null(scalebar_start_position)) x["scalebar_start_position"] <- scalebar_start_position
  if (!is.null(scalebar_end_position)) x["scalebar_end_position"] <- scalebar_end_position
  new("genome_track", tracks = list(x))
}


#' @description Create spacing track with custom height.
#' @title Generate spacing track
#' @inheritParams track_bed
#' @return None
#' @export
#' @inherit track_bed examples
#' @importFrom methods getClass is new
#' @author Omar Elashkar
track_spacer <- function(title = NULL,
                         height = 2,
                         overlay_previous = "no") {
  x <- list(
    file = "spacer",
    title = title,
    height = height,
    overlay_previous = overlay_previous
  )

  if (!is.null(title)) x["title"] <- title

  new("genome_track", tracks = list(x))
}


#' @description This track will specifiy the options for x-axis for location,
#' height, font size and wheather to overlay previous track.
#' @title Specify x_axis option for genome_track.
#' @inheritParams track_bed
#' @param where String. Either "bottom" (default) or "top"
#' @return genome_track
#' @export
#' @inherit track_domains examples
#' @note
#' `fontsize` argument can be overriden by the same argument in `plot_gtracks()`
#' @importFrom methods getClass is new
#' @author Omar Elashkar
track_x_axis <- function(title = NULL,
                         height = 2,
                         overlay_previous = "no",
                         where = "bottom",
                         fontsize = 15) {
  x <- list(
    file = "x-axis",
    height = height,
    overlay_previous = overlay_previous,
    where = where,
    fontsize = fontsize
  )

  if (!is.null(title)) x["title"] <- title

  new("genome_track", tracks = list(x))
}
