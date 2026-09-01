# Verify that traverse builds and loads from an isolated installation.

find_project_root <- function(start = getwd()) {
  path <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION")) &&
        dir.exists(file.path(path, "R"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Run this script inside the traverse project.", call. = FALSE)
    }
    path <- parent
  }
}

run_command <- function(command, args) {
  output <- system2(command, args = args, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && !identical(as.integer(status), 0L)) {
    stop(
      paste(c(sprintf("Command failed with status %s:", status), output), collapse = "\n"),
      call. = FALSE
    )
  }
  invisible(output)
}

root <- find_project_root()
build_dir <- tempfile("traverse-build-")
install_dir <- tempfile("traverse-install-")
source_dir <- tempfile("traverse-source-")
dir.create(build_dir)
dir.create(install_dir)
dir.create(source_dir)

# Build from a clean source tree so local Git metadata and development outputs
# cannot affect the package source or make the result depend on the checkout.
excluded_names <- c(
  ".git", ".Rproj.user", ".quarto", "validation-output", ".Rhistory",
  ".RData", ".RDataTmp", ".Ruserdata", ".Rcheck", "..Rcheck",
  "traverse.Rcheck"
)
source_items <- list.files(root, all.files = TRUE, no.. = TRUE, full.names = TRUE)
source_items <- source_items[!basename(source_items) %in% excluded_names]
source_items <- source_items[!grepl("[.]tar[.]gz$|[.]Rout$", basename(source_items))]
copied <- vapply(
  source_items,
  function(path) file.copy(path, source_dir, recursive = TRUE, copy.date = TRUE),
  logical(1L)
)
if (!all(copied)) {
  stop("Could not create the clean source tree for installation verification.", call. = FALSE)
}

r_exec <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "R.exe" else "R")
old_dir <- setwd(build_dir)
on.exit(setwd(old_dir), add = TRUE)

run_command(r_exec, c("CMD", "build", "--no-build-vignettes", shQuote(source_dir)))
tarballs <- list.files(build_dir, pattern = "^traverse_.*[.]tar[.]gz$", full.names = TRUE)
if (length(tarballs) != 1L) {
  stop("Expected exactly one traverse source tarball after R CMD build.", call. = FALSE)
}

run_command(
  r_exec,
  c("CMD", "INSTALL", "-l", shQuote(install_dir), shQuote(tarballs[[1L]]))
)

old_libpaths <- .libPaths()
on.exit(.libPaths(old_libpaths), add = TRUE)
.libPaths(c(install_dir, old_libpaths))
library(traverse)

installed_path <- normalizePath(find.package("traverse"), mustWork = TRUE)
expected_prefix <- normalizePath(install_dir, mustWork = TRUE)
if (!startsWith(installed_path, expected_prefix)) {
  stop("traverse was not loaded from the isolated installation.", call. = FALSE)
}

if (!inherits(citation("traverse"), "citation")) {
  stop("Installed package citation metadata could not be loaded.", call. = FALSE)
}

expected_exports <- c(
  "traverse", "traverse_domain", "traverse_mask", "traverse_nodes",
  "traverse_pairs", "traverse_passage", "traverse_surface"
)
missing_exports <- setdiff(expected_exports, getNamespaceExports("traverse"))
if (length(missing_exports) > 0L) {
  stop(
    sprintf("Installed package is missing exports: %s", paste(missing_exports, collapse = ", ")),
    call. = FALSE
  )
}

cat("Clean installation verified from", installed_path, "\n")
