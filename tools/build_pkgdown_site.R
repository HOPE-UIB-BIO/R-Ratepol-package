build_pkgdown_site <- function() {
  cnd_build_error <- NULL

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
