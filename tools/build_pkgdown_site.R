build_pkgdown_site <- function() {
  cnd_build_error <- NULL

  # Patch README.md: Quarto renders fig.alt as data-fig-alt, but pkgdown's
  # accessibility checker requires the standard HTML alt attribute.
  # Convert every occurrence before pkgdown reads the file.
  # Read and write as raw bytes to avoid iconv truncation on Windows:
  # R text connections convert through the native locale (e.g. windows-1252),
  # and characters like U+2139 (ℹ) that have no native equivalent cause
  # writeLines() to stop mid-file. The substitution targets pure ASCII
  # attribute names, so byte-level matching is correct.
  raw_readme <- readBin(
    con = "README.md",
    what = "raw",
    n = file.info("README.md")$size
  )
  str_readme <- rawToChar(raw_readme)
  str_readme <- gsub(
    pattern = " data-fig-alt=\"",
    replacement = " alt=\"",
    x = str_readme,
    fixed = TRUE,
    useBytes = TRUE
  )
  writeBin(
    object = charToRaw(str_readme),
    con = "README.md"
  )

  tryCatch(
    expr = {
      pkgdown::build_site_github_pages(
        pkg = ".",
        new_process = TRUE
      )
    },
    error = function(cnd_error) {
      cnd_build_error <<- cnd_error
    }
  )

  vec_keep_github_pages <-
    c("404", "CODE_OF_CONDUCT", "CONTRIBUTING", "SUPPORT")

  vec_github_md <-
    list.files(
      path = ".github",
      pattern = "\\.md$",
      full.names = FALSE
    )

  vec_drop_pages <-
    setdiff(
      tools::file_path_sans_ext(vec_github_md),
      vec_keep_github_pages
    )

  if (length(vec_drop_pages) > 0L) {
    unlink(
      file.path("docs", paste0(vec_drop_pages, ".html")),
      force = TRUE
    )
  }

  vec_markdown_mirrors <-
    list.files(
      path = "docs",
      pattern = "\\.md$",
      recursive = TRUE,
      full.names = TRUE
    )

  if (length(vec_markdown_mirrors) > 0L) {
    unlink(vec_markdown_mirrors, force = TRUE)
  }

  unlink(file.path("docs", "llms.txt"), force = TRUE)

  if (is.null(cnd_build_error)) {
    getFromNamespace("build_sitemap", ns = "pkgdown")()
    pkgdown::build_search()
  }

  if (!is.null(cnd_build_error)) {
    stop(cnd_build_error)
  }

  invisible(NULL)
}

if (sys.nframe() == 0L) {
  build_pkgdown_site()
}
