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
})
