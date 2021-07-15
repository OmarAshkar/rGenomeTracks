# No constructor, only generated via track_
setClass(
  "genome_track",
  slots = list(tracks = "list")
)

#' This method adds two "genome_track" objects together.
#' @title Adding genome_track Objects
#' @param e1 genome_track object.
#' @param e2 genome_track object.
#' @return genome_track object
#' @export
#' @inherit track_links examples
#' @author Omar Elashkar
setMethod(
  "+", c("genome_track", "genome_track"),
  function(e1, e2) {
    if (is(getClass(e2), "genome_track")) {
      e1@tracks[length(e1@tracks) + 1] <- e2@tracks
      new("genome_track", tracks = e1@tracks)
    }
  }
)
