# Build submission/reference.docx, the Word style template used to render the submission documents (manuscript and supplement).
# Run from the project root. Needs Quarto (for Pandoc's default reference document) and the system `zip`.
#
# Starts from Pandoc's default reference document and sets what review manuscripts at JID/OUP expect: US Letter, 1-inch margins, 12-point Times New Roman,
# double-spaced body text, continuous line numbers, page numbers centred in the footer, black headings. Table cells (Pandoc's "Compact" style) stay single-spaced.

quarto <- Sys.which("quarto")
stopifnot(nzchar(quarto), nzchar(Sys.which("zip")))
tmp <- tempfile("refdocx"); dir.create(tmp)
src <- file.path(tmp, "default.docx")
system2(quarto, c("pandoc", "-o", shQuote(src), "--print-default-data-file", "reference.docx"))
unz <- file.path(tmp, "x"); utils::unzip(src, exdir = unz)

rd <- function(f) paste(readLines(file.path(unz, f), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
wr <- function(x, f) writeLines(x, file.path(unz, f), useBytes = TRUE)
sub1 <- function(pattern, replacement, x) { if (!grepl(pattern, x, perl = TRUE)) stop("pattern not found: ", pattern, call. = FALSE); sub(pattern, replacement, x, perl = TRUE) }

# ---- styles ----
st <- rd("word/styles.xml")
st <- sub1("<w:rFonts [^>]*/>\\s*<w:sz w:val=\"24\" />",
           "<w:rFonts w:ascii=\"Times New Roman\" w:hAnsi=\"Times New Roman\" w:eastAsia=\"Times New Roman\" w:cs=\"Times New Roman\" />\n        <w:sz w:val=\"24\" />", st)
st <- sub1("<w:pPrDefault>\\s*<w:pPr>\\s*<w:spacing w:after=\"200\" />",
           "<w:pPrDefault>\n      <w:pPr>\n        <w:spacing w:after=\"0\" w:line=\"240\" w:lineRule=\"auto\" />", st)
# Single spacing by default (table cells drawn by gt have no paragraph style); body text, first paragraphs and the reference list double-spaced
# Replace (or add) the paragraph properties of one style, working inside that style's own block only
style_ppr <- function(x, id, ppr) {
  pat <- sprintf("(?s)<w:style [^>]*w:styleId=\"%s\"[^>]*>.*?</w:style>", id)
  m <- regmatches(x, regexpr(pat, x, perl = TRUE))
  if (length(m) != 1L) stop("no style ", id, call. = FALSE)
  block <- if (grepl("<w:pPr>", m, fixed = TRUE)) sub("(?s)<w:pPr>.*?</w:pPr>", ppr, m, perl = TRUE)
           else if (grepl("<w:rPr>", m, fixed = TRUE)) sub("<w:rPr>", paste0(ppr, "<w:rPr>"), m, fixed = TRUE)
           else sub("</w:style>", paste0(ppr, "</w:style>"), m, fixed = TRUE)
  sub(m, block, x, fixed = TRUE)
}
st <- style_ppr(st, "BodyText", "<w:pPr><w:spacing w:before=\"0\" w:after=\"120\" w:line=\"480\" w:lineRule=\"auto\" /></w:pPr>")
st <- style_ppr(st, "FirstParagraph", "<w:pPr><w:spacing w:before=\"0\" w:after=\"120\" w:line=\"480\" w:lineRule=\"auto\" /></w:pPr>")
st <- style_ppr(st, "Compact", "<w:pPr><w:spacing w:before=\"0\" w:after=\"0\" w:line=\"240\" w:lineRule=\"auto\" /></w:pPr>")
st <- style_ppr(st, "Bibliography", "<w:pPr><w:spacing w:before=\"0\" w:after=\"120\" w:line=\"480\" w:lineRule=\"auto\" /></w:pPr>")
# Headings, title and captions in black (Pandoc's default is theme blue)
st <- gsub("<w:color [^>]*w:themeColor=\"accent1\"[^>]*/>", "<w:color w:val=\"000000\" />", st, perl = TRUE)
st <- gsub("<w:rFonts [^>]*(majorHAnsi|majorBidi)[^>]*/>",
           "<w:rFonts w:ascii=\"Times New Roman\" w:hAnsi=\"Times New Roman\" w:eastAsia=\"Times New Roman\" w:cs=\"Times New Roman\" />", st, perl = TRUE)
# Title and subtitle: 16 and 13 point bold, left-aligned (Pandoc's default title is 28 point and centred)
title_rpr <- function(x, id, half_points) {
  pat <- sprintf("(?s)<w:style [^>]*w:styleId=\"%s\"[^>]*>.*?</w:style>", id)
  m <- regmatches(x, regexpr(pat, x, perl = TRUE))
  block <- gsub("<w:sz w:val=\"[0-9]+\" />", sprintf("<w:sz w:val=\"%d\" />", half_points), m, perl = TRUE)
  block <- gsub("<w:szCs w:val=\"[0-9]+\" />", sprintf("<w:szCs w:val=\"%d\" />", half_points), block, perl = TRUE)
  block <- gsub("<w:jc w:val=\"center\" />", "<w:jc w:val=\"left\" />", block, fixed = TRUE)
  if (!grepl("<w:b />", block, fixed = TRUE)) block <- sub("<w:rPr>", "<w:rPr><w:b />", block, fixed = TRUE)
  sub(m, block, x, fixed = TRUE)
}
st <- title_rpr(st, "Title", 32)
st <- title_rpr(st, "Subtitle", 26)
wr(st, "word/styles.xml")

# ---- footer with the page number ----
wr(paste0('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n',
  '<w:ftr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">',
  '<w:p><w:pPr><w:jc w:val="center" /><w:spacing w:after="0" w:line="240" w:lineRule="auto" /></w:pPr>',
  '<w:r><w:fldChar w:fldCharType="begin" /></w:r><w:r><w:instrText xml:space="preserve"> PAGE </w:instrText></w:r>',
  '<w:r><w:fldChar w:fldCharType="separate" /></w:r><w:r><w:t>1</w:t></w:r><w:r><w:fldChar w:fldCharType="end" /></w:r></w:p></w:ftr>'),
  "word/footer1.xml")
rels <- rd("word/_rels/document.xml.rels")
rels <- sub1("</Relationships>", paste0('<Relationship Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/footer" Id="rIdFtr1" Target="footer1.xml" /></Relationships>'), rels)
wr(rels, "word/_rels/document.xml.rels")
ct <- rd("[Content_Types].xml")
ct <- sub1("</Types>", '<Override PartName="/word/footer1.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.footer+xml" /></Types>', ct)
wr(ct, "[Content_Types].xml")

# ---- section: page size, margins, footer, line numbers ----
doc <- rd("word/document.xml")
if (!grepl("xmlns:r=", doc, fixed = TRUE)) doc <- sub1("<w:document ", "<w:document xmlns:r=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships\" ", doc)
doc <- sub1("(?s)<w:sectPr>.*?</w:sectPr>", paste0(
  "<w:sectPr><w:footerReference w:type=\"default\" r:id=\"rIdFtr1\" /><w:footnotePr><w:numRestart w:val=\"eachSect\" /></w:footnotePr>",
  "<w:pgSz w:w=\"12240\" w:h=\"15840\" /><w:pgMar w:top=\"1440\" w:right=\"1440\" w:bottom=\"1440\" w:left=\"1440\" w:header=\"720\" w:footer=\"720\" w:gutter=\"0\" />",
  "<w:lnNumType w:countBy=\"1\" w:restart=\"continuous\" /></w:sectPr>"), doc)
wr(doc, "word/document.xml")

# Zip the template folder only. Inside a function, so that on.exit() restores the working directory after the zip (at the top level of a sourced
# script it would fire at once, and the zip would run in the project root).
zip_dir <- function(dir, out) {
  owd <- setwd(dir); on.exit(setwd(owd))
  files <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
  stopifnot(all(c("[Content_Types].xml", "word/document.xml", "word/styles.xml") %in% files), length(files) < 40)
  utils::zip(out, files = files, flags = "-X -q")
}
out <- file.path(normalizePath("submission"), "reference.docx")
if (file.exists(out)) file.remove(out)
zip_dir(unz, out)
stopifnot(file.size(out) < 1e6)
message("wrote ", out, " (", round(file.size(out) / 1024), " KB)")
