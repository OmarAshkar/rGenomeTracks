
test_that("track_bed", {
  bed12_dir <- system.file("extdata",
    "dm3_genes.bed.gz",
    package = "rGenomeTracks"
  )
  test <- track_bed(bed12_dir)
  expect_s4_class(test, "genome_track")
})

test_that("track_bedgraph", {
  bg_dir <- system.file("extdata", "GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz",
    package = "rGenomeTracks"
  )
  test <- track_bedgraph(bg_dir)
  expect_s4_class(test, "genome_track")
})


test_that("track_gtf", {
  gtf_dir <- system.file("extdata", "dm3_subset_BDGP5.78.gtf.gz",
    package = "rGenomeTracks"
  )
  test <- track_gtf(file = gtf_dir)
  expect_s4_class(test, "genome_track")
})

test_that("track_links", {
  links_dir <- system.file("extdata",
    "test.arcs",
    package = "rGenomeTracks"
  )
  test <- track_links(links_dir)
  expect_s4_class(test, "genome_track")
})
test_that("track_hlines", {
  test <- track_hlines(y_values = "12 15")
  expect_s4_class(test, "genome_track")
})
test_that("track_spacer", {
  test <- track_spacer()
  expect_s4_class(test, "genome_track")
})
test_that("track_bigwig", {
  bw_dir <- system.file("extdata", "bigwig2_X_2.5e6_3.5e6.bw", package = "rGenomeTracks")
  test <- track_bigwig(file = bw_dir)
  expect_s4_class(test, "genome_track")
})
test_that("track_epilogos", {
  epilog_dir <- system.file("extdata", "epilog.qcat.bgz", package = "rGenomeTracks")
  test <- track_epilogos(file = epilog_dir)
  expect_s4_class(test, "genome_track")
})

test_that("track_narrow_peak", {
  np_bed_dir <- system.file("extdata", "test2.narrowPeak", package = "rGenomeTracks")
  test <- track_narrow_peak(np_bed_dir)
  expect_s4_class(test, "genome_track")
})

test_that("track_domains", {
  tads_dir <- system.file("extdata", "tad_classification.bed", package = "rGenomeTracks")
  test <- track_domains(
    file = tads_dir
  )

  expect_s4_class(test, "genome_track")
})

test_that("track_hic_matrix", {
  h5_dir <- system.file("extdata", "Li_et_al_2015.h5", package = "rGenomeTracks")
  test <- track_hic_matrix(h5_dir)
  expect_s4_class(test, "genome_track")
})

test_that("track_links", {
  links_dir <- system.file("extdata", "test.arcs",
    package = "rGenomeTracks"
  )
  test <- track_links(links_dir)
  expect_s4_class(test, "genome_track")
})

test_that("track_scalebar", {
  test <- track_scalebar()
  expect_s4_class(test, "genome_track")
})
