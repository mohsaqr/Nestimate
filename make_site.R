#!/usr/bin/env Rscript
# Build the pkgdown site, then enforce .Rbuildignore on the published pages.
#
# pkgdown renders every root-level *.md as a site page (see
# pkgdown:::package_mds) and does NOT consult .Rbuildignore, so internal
# working files (CLAUDE.md, HANDOFF.md, LEARNINGS.md, CHANGES.md, FEATURES.md,
# TODO.md, AGENT-NOTE-*.md, Nestimate-paper.md, ...) leak into docs/ and onto
# the public site. This wrapper applies the rule pkgdown skips: if R CMD build
# ignores a root *.md, its rendered page is removed from docs/ and the search
# index is rebuilt so it carries no stale entries.
#
# Run from the package root:  Rscript make_site.R

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

local({
  # Use the locally extracted Quarto if present (needed for the .qmd articles).
  quarto_bin <- path.expand("~/.local/quarto/bin")
  if (dir.exists(quarto_bin)) {
    Sys.setenv(PATH = paste(quarto_bin, Sys.getenv("PATH"),
                            sep = .Platform$path.sep))
    Sys.setenv(QUARTO_PATH = file.path(quarto_bin, "quarto"))
  }

  pkgdown::build_site(".", preview = FALSE, install = FALSE)

  # Root *.md files that R CMD build would ignore (the same patterns CRAN uses).
  ignore <- if (file.exists(".Rbuildignore")) {
    Filter(nzchar, readLines(".Rbuildignore", warn = FALSE))
  } else {
    character(0)
  }
  root_md <- list.files(".", pattern = "\\.md$", ignore.case = TRUE)
  is_ignored <- function(f) {
    any(vapply(ignore, function(p) grepl(p, f, perl = TRUE), logical(1)))
  }
  leaked_md <- Filter(is_ignored, root_md)

  pages <- file.path("docs", sub("\\.md$", ".html", leaked_md))
  removed <- pages[file.exists(pages)]
  unlink(removed)

  if (length(removed)) {
    message("Removed leaked internal pages: ",
            paste(basename(removed), collapse = ", "))
    pkgdown::build_search(".")  # re-index so search.json drops the removed pages
  } else {
    message("No leaked internal pages found in docs/.")
  }

  # Google Scholar indexes an article from Highwire Press meta tags in its
  # head; pkgdown cannot give each article its own, so they are written into
  # the built page of every article that has a PDF in pkgdown/assets/pdf/.
  tagged <- add_scholar_tags("docs", site = "https://pak.dynasite.org/Nestimate")
  message("Google Scholar tags added to ", length(tagged), " article(s).")
})
