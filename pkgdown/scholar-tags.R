# Google Scholar (Highwire Press) tags for the built article pages.
# Sourced by make_site.R and by .github/workflows/pkgdown.yaml.

#' Write Google Scholar (Highwire Press) meta tags into built article pages
#'
#' @param docs The built site folder.
#' @param site The site's address, for the PDF links.
#' @return The articles tagged, invisibly. A page that already carries
#'   `citation_title` is left as it is, so running twice changes nothing.
add_scholar_tags <- function(docs, site) {
  sources <- c(list.files("vignettes", pattern = "\\.(Rmd|qmd)$", full.names = TRUE),
               list.files(file.path("vignettes", "articles"),
                          pattern = "\\.(Rmd|qmd)$", full.names = TRUE))
  slugs <- tools::file_path_sans_ext(basename(sources))
  has_pdf <- file.exists(file.path("pkgdown", "assets", "pdf", paste0(slugs, ".pdf")))
  pages <- file.path(docs, "articles", paste0(slugs, ".html"))
  keep <- has_pdf & file.exists(pages)
  tag_one <- function(source, slug, page) {
    front <- rmarkdown::yaml_front_matter(source)
    title <- paste(c(front$title, front$subtitle), collapse = ": ")
    authors <- vapply(front$author, function(a) {
      name <- if (is.list(a)) a$name else a
      parts <- strsplit(name, " ", fixed = TRUE)[[1L]]
      # "Mohammed Saqr" -> "Saqr, Mohammed", the form Scholar reads
      paste0(parts[length(parts)], ", ", paste(parts[-length(parts)], collapse = " "))
    }, character(1L))
    escape <- function(x) {
      x <- gsub("&", "&amp;", x, fixed = TRUE)
      x <- gsub("\"", "&quot;", x, fixed = TRUE)
      gsub("<", "&lt;", x, fixed = TRUE)
    }
    meta <- function(name, content) {
      sprintf('<meta name="%s" content="%s">', name, escape(content))
    }
    tags <- c(meta("citation_title", title),
              vapply(authors, function(a) meta("citation_author", a), character(1L)),
              meta("citation_publication_date", "2026"),
              meta("citation_publisher", "Nestimate"),
              meta("citation_pdf_url", sprintf("%s/pdf/%s.pdf", site, slug)),
              meta("citation_abstract_html_url",
                   sprintf("%s/articles/%s.html", site, slug)))
    html <- readLines(page, warn = FALSE, encoding = "UTF-8")
    if (any(grepl('name="citation_title"', html, fixed = TRUE))) return(FALSE)
    at <- grep("</head>", html, fixed = TRUE)[1L]
    stopifnot("the page has a </head>" = !is.na(at))
    html <- append(html, tags, after = at - 1L)
    writeLines(html, page, useBytes = TRUE)
    TRUE
  }
  done <- unlist(Map(tag_one, sources[keep], slugs[keep], pages[keep]))
  invisible(slugs[keep][done])
}
