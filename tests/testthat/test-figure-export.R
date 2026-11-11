source(file.path(proj_root, "R", "figure_export.R"))

test_that("the CRC-32 matches the standard check value", {
  expect_equal(png_crc32(charToRaw("123456789")), as.raw(c(0xcb, 0xf4, 0x39, 0x26)))
  expect_equal(png_crc32(raw(0)), as.raw(c(0, 0, 0, 0)))
})

test_that("a PNG's resolution is stamped, replaced and read back, and the file is still a valid image", {
  f <- tempfile(fileext = ".png"); grDevices::png(f, width = 300, height = 200, res = 72); plot(1:3); grDevices::dev.off()
  png_set_dpi(f, 300)
  i <- png_info(f)
  expect_equal(unname(i[["dpi"]]), 300, tolerance = 0.01); expect_equal(unname(i[c("width_px", "height_px")]), c(300, 200))
  png_set_dpi(f, 600)
  expect_equal(unname(png_info(f)[["dpi"]]), 600, tolerance = 0.01)
  expect_equal(sum(vapply(png_chunks(readBin(f, "raw", file.size(f))), function(c) c$type == "pHYs", TRUE)), 1L)
  skip_if_not_installed("rsvg")
  expect_equal(dim(rsvg::rsvg_raw(charToRaw("<svg xmlns='http://www.w3.org/2000/svg' width='10' height='10'><rect width='10' height='10'/></svg>"), width = 20))[2:3], c(20, 20))
})

test_that("an SVG is rendered at the stated width and resolution", {
  skip_if_not_installed("rsvg")
  f <- tempfile(fileext = ".png")
  svg_to_png("<svg xmlns='http://www.w3.org/2000/svg' width='100' height='50'><rect width='100' height='50' fill='#1B365D'/></svg>", f, width_in = 2, dpi = 300)
  i <- png_info(f)
  expect_equal(unname(i[c("width_px", "height_px", "dpi")]), c(600, 300, 300), tolerance = 0.01)
})

test_that("the PRISMA side labels are inserted into each blue box of the SVG", {
  box <- function(nm) sprintf("<g><a xlink:title=\"Blue %s box\">\n<path fill=\"#bcd2ee\" d=\"M19.3,-666C19.3,-666 9.7,-666 9.7,-666 0,-656 0,-550 0,-550 0,-540 19.3,-540 19.3,-540\"/>\n</a></g>", nm)
  svg <- paste(box("identification"), box("screening"), box("included"), sep = "\n")
  out <- prisma_side_labels(svg)
  expect_equal(lengths(regmatches(out, gregexpr("<text ", out))), 3L)
  expect_true(all(vapply(c("Identification", "Screening", "Included"), grepl, TRUE, x = out)))
  expect_match(out, "translate\\(13\\.65,-603\\.00\\) rotate\\(-90\\)")
  expect_error(prisma_side_labels("<svg></svg>"), "no identification box")
})
