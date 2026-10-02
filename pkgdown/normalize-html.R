# Work around two pkgdown 2.2.1 HTML-generation defects. These exact replacements
# are harmless if a later pkgdown version has already corrected the templates.
site_url <- yaml::read_yaml("pkgdown/_pkgdown.yml")$url
for (path in list.files("docs", pattern = "\\.html$", recursive = TRUE,
                       full.names = TRUE)) {
  original <- readLines(path, warn = FALSE, encoding = "UTF-8")
  corrected <- gsub('type="\u201dimage/svg+xml\u201d"', 'type="image/svg+xml"',
                    original, fixed = TRUE)
  if (basename(path) == "404.html") {
    # Absolute-link conversion invents a src for inline scripts (missing = NA).
    # Only remove that exact artificial attribute from scripts with a body.
    html <- xml2::read_html(paste(corrected, collapse = "\n"))
    scripts <- xml2::xml_find_all(html, "//script[@src]")
    fake_src <- paste0(site_url, "NA")
    inline <- xml2::xml_attr(scripts, "src") == fake_src &
      nzchar(trimws(xml2::xml_text(scripts)))
    if (any(inline)) {
      corrected <- gsub(paste0('<script src="', fake_src, '">'),
                        '<script>', corrected, fixed = TRUE)
    }
  }
  if (!identical(original, corrected)) {
    writeLines(enc2utf8(corrected), path, useBytes = TRUE)
  }
}
