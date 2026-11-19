# Helpers for exporting figures at a stated resolution (scripts 44 and 53). Base R plus, for SVG, the rsvg package.
#
# R's png device and rsvg write the pixels but not always the resolution, so png_set_dpi() writes the pHYs chunk (pixels per metre) that tells a viewer or a
# journal's checker the intended print size. png_info() reads it back.

# CRC-32 of a raw vector (the PNG chunk checksum), in 16-bit halves so that it stays inside R's 32-bit signed integers
png_crc32 <- function(bytes) {
  tab_hi <- tab_lo <- integer(256)
  for (n in 0:255) {
    hi <- 0L; lo <- n
    for (k in 1:8) {
      carry <- bitwAnd(lo, 1L)
      lo <- bitwOr(bitwShiftR(lo, 1L), bitwShiftL(bitwAnd(hi, 1L), 15L)); hi <- bitwShiftR(hi, 1L)
      if (carry == 1L) { hi <- bitwXor(hi, 0xEDB8L); lo <- bitwXor(lo, 0x8320L) }
    }
    tab_hi[n + 1L] <- hi; tab_lo[n + 1L] <- lo
  }
  hi <- 0xFFFFL; lo <- 0xFFFFL
  for (b in as.integer(bytes)) {
    i <- bitwAnd(bitwXor(lo, b), 255L) + 1L
    nlo <- bitwXor(bitwOr(bitwShiftR(lo, 8L), bitwShiftL(bitwAnd(hi, 255L), 8L)), tab_lo[i])
    hi <- bitwXor(bitwShiftR(hi, 8L), tab_hi[i]); lo <- nlo
  }
  hi <- bitwXor(hi, 0xFFFFL); lo <- bitwXor(lo, 0xFFFFL)
  as.raw(c(bitwShiftR(hi, 8L), bitwAnd(hi, 255L), bitwShiftR(lo, 8L), bitwAnd(lo, 255L)))
}

be32 <- function(x) { x <- as.numeric(x); as.raw(c(x %/% 16777216, (x %/% 65536) %% 256, (x %/% 256) %% 256, x %% 256)) }

# Read the chunks of a PNG as a list of (type, start, length)
png_chunks <- function(b) {
  stopifnot(identical(b[1:8], as.raw(c(0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A))))
  pos <- 9L; out <- list()
  while (pos < length(b)) {
    len <- sum(as.numeric(b[pos:(pos + 3L)]) * c(16777216, 65536, 256, 1))
    out[[length(out) + 1L]] <- list(type = rawToChar(b[(pos + 4L):(pos + 7L)]), start = pos, size = len + 12)
    pos <- pos + as.integer(len) + 12L
  }
  out
}

# Set (or replace) the resolution of a PNG file, in dots per inch
png_set_dpi <- function(file, dpi = 300) {
  b <- readBin(file, "raw", file.size(file))
  ch <- png_chunks(b)
  keep <- !vapply(ch, function(c) c$type == "pHYs", logical(1))
  ppm <- round(dpi / 0.0254)
  data <- c(be32(ppm), be32(ppm), as.raw(1L))
  phys <- c(be32(9), charToRaw("pHYs"), data, png_crc32(c(charToRaw("pHYs"), data)))
  pieces <- lapply(ch[keep], function(c) b[c$start:(c$start + c$size - 1)])
  writeBin(c(b[1:8], pieces[[1]], phys, unlist(pieces[-1])), file)
  invisible(file)
}

# Pixel size and resolution (dpi, from the pHYs chunk; NA if absent) of a PNG file
png_info <- function(file) {
  b <- readBin(file, "raw", file.size(file))
  ch <- png_chunks(b)
  px <- function(x) sum(as.numeric(x) * c(16777216, 65536, 256, 1))
  p <- Filter(function(c) c$type == "pHYs", ch)
  dpi <- if (length(p)) px(b[(p[[1]]$start + 8):(p[[1]]$start + 11)]) * 0.0254 else NA_real_
  c(width_px = px(b[17:20]), height_px = px(b[21:24]), dpi = dpi)
}

# Render an SVG string to a PNG `width_in` inches wide at `dpi`, and stamp the resolution
svg_to_png <- function(svg, file, width_in, dpi = 300) {
  rsvg::rsvg_png(charToRaw(svg), file, width = round(width_in * dpi))
  png_set_dpi(file, dpi)
}

# The PRISMA2020 package draws the rotated side labels ("Identification", "Screening", "Included") with JavaScript in the HTML widget, so they are missing from its SVG.
# This writes them into the SVG, centred in the three blue boxes, so that the SVG renders the same as the widget.
prisma_side_labels <- function(svg, labels = c(identification = "Identification", screening = "Screening", included = "Included")) {
  for (nm in names(labels)) {
    m <- regmatches(svg, regexec(sprintf('<a xlink:title="Blue %s box">[[:space:]]*<path [^>]*d="([^"]*)"[^>]*/>', nm), svg))[[1]]
    if (length(m) != 2L) stop("no ", nm, " box in the SVG", call. = FALSE)
    nums <- as.numeric(strsplit(gsub("[A-Za-z]", " ", m[2]), "[ ,]+")[[1]]); nums <- nums[!is.na(nums)]
    x <- nums[c(TRUE, FALSE)]; y <- nums[c(FALSE, TRUE)]
    txt <- sprintf('<text transform="translate(%.2f,%.2f) rotate(-90)" text-anchor="middle" font-family="Helvetica,sans-Serif" font-size="12.00" fill="#000000">%s</text>',
                   mean(range(x)) + 4, mean(range(y)), labels[[nm]])
    svg <- sub(m[1], paste0(m[1], "\n", txt), svg, fixed = TRUE)
  }
  svg
}

# The PRISMA2020 package stops at "New studies included in review". This adds a box below it for the extraction stage, joined by an arrow, in the package's
# own box style, and lengthens the page and the blue "Included" band to hold it. Call it before prisma_side_labels() so that the label is centred in the longer band.
# `lines` are the box's text lines; `extra` is the added height in points.
prisma_extraction_box <- function(svg, lines, extra = 120) {
  m <- regmatches(svg, regexec('<a xlink:title="New studies included in review">[[:space:]]*<polygon [^>]*points="([^"]*)"', svg))[[1]]
  if (length(m) != 2L) stop("no 'New studies included' box in the SVG", call. = FALSE)
  pts <- matrix(as.numeric(strsplit(m[2], "[ ,]+")[[1]]), ncol = 2, byrow = TRUE)
  x0 <- min(pts[, 1]); x1 <- max(pts[, 1]); yb <- max(pts[, 2]); xm <- (x0 + x1) / 2
  top <- 12; lh <- 14.4; h <- length(lines) * lh + 14
  txt <- sprintf('<text text-anchor="middle" x="%.2f" y="%.2f" font-family="Helvetica,sans-Serif" font-size="12.00" fill="#000000">%s</text>',
                 xm, top + 7 + lh * (seq_along(lines) - 0.3), lines)
  node <- c('<g id="extraction" class="node"><a xlink:title="Data extraction">',
            sprintf('<polygon fill="none" stroke="#36648b" points="%.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f"/>',
                    x1, top, x0, top, x0, top + h, x1, top + h, x1, top), txt, '</a></g>',
            sprintf('<g class="edge"><path fill="none" stroke="#36648b" d="M%.1f,%.1f L%.1f,%.1f"/><polygon fill="#36648b" stroke="#36648b" points="%.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f"/></g>',
                    xm, yb, xm, top - 4, xm - 3.5, top - 11, xm + 3.5, top - 11, xm, top - 4, xm - 3.5, top - 11))
  if (top + h > extra - 6) stop("`extra` is too small for the box", call. = FALSE)
  # Lengthen the page, its white background and the blue Included band
  h0 <- as.numeric(regmatches(svg, regexec('viewBox="0.00 0.00 ([0-9.]+) ([0-9.]+)"', svg))[[1]][3])
  svg <- sub('viewBox="0.00 0.00 ([0-9.]+) [0-9.]+"', sprintf('viewBox="0.00 0.00 \\1 %.1f"', h0 + extra), svg)
  svg <- sub('height="[0-9]+pt"', sprintf('height="%.0fpt"', h0 + extra), svg)
  bg <- regmatches(svg, regexec('<polygon fill="#ffffff" stroke="transparent" points="([^"]*)"', svg))[[1]]
  svg <- sub(bg[1], sub(bg[2], gsub("(,)4( |$)", paste0("\\1", 4 + extra, "\\2"), bg[2]), bg[1], fixed = TRUE), svg, fixed = TRUE)
  band <- regmatches(svg, regexec('(<a xlink:title="Blue included box">[[:space:]]*<path [^>]*d=")([^"]*)(")', svg))[[1]]
  tok <- regmatches(band[3], gregexpr("-?[0-9.]+,-?[0-9.]+", band[3]))[[1]]
  xy <- do.call(rbind, strsplit(tok, ",")); y <- as.numeric(xy[, 2]); y <- ifelse(y > -90, y + extra, y)
  regmatches(band[3], gregexpr("-?[0-9.]+,-?[0-9.]+", band[3])) <- list(paste0(xy[, 1], ",", y))
  svg <- sub(band[1], paste0(band[2], band[3], band[4]), svg, fixed = TRUE)
  sub("</g>\n</svg>", paste0(paste(node, collapse = "\n"), "\n</g>\n</svg>"), svg, fixed = TRUE)
}
