# Repository publication guard; runs with base R before installing dependencies.
# This checks tracked file paths, formats, size and binary content.
# Scientific provenance and redistribution rights still need human review.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1L]] else ".", winslash = "/", mustWork = TRUE)
git <- function(arguments) {
  result <- system2("git", c("-C", shQuote(root), arguments), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(result, "status"))) stop(paste(result, collapse = "\n"), call. = FALSE)
  result
}
top <- normalizePath(git(c("rev-parse", "--show-toplevel")), winslash = "/", mustWork = TRUE)
if (!identical(root, top)) stop("Run this guard at the package repository root.", call. = FALSE)
files <- git(c("-c", "core.quotepath=false", "ls-files"))
if (!length(files)) stop("No tracked files to review. Stage the public sources first.", call. = FALSE)
if (any(grepl("[\r\n]", files) | startsWith(files, '"')))
  stop("Unsupported tracked file names; review them explicitly.", call. = FALSE)

problems <- character()
flag <- function(bad, reason) {
  if (any(bad)) problems <<- c(problems, paste0(reason, ": ", paste(files[bad], collapse = ", ")))
}
lower <- tolower(files)
flag(grepl("(^|/)(private|benchmarks?|outputs?|cache|caches|data|local-library|artifacts)(/|$)", lower),
     "Private data, outputs or generated directories")
flag(grepl("(^|/)(snps\\.txt|solutions|pev_pec|hinv\\.rds|sigma_focals_complete\\.rds|gima22i_ren\\.txt)$", lower),
     "Private evaluation inputs or outputs")
flag(grepl("(^|/)(\\.renviron|\\.env([.]|$)|id_rsa|id_ed25519)(/|$)", lower),
     "Local configuration or key files")

allowed_exact <- c("DESCRIPTION", "NAMESPACE", "LICENSE", "inst/CITATION",
                   ".gitignore", ".gitattributes", ".Rbuildignore")
allowed_type <- grepl("\\.(R|Rd|md|cff|ya?ml)$", files)
flag(!(files %in% allowed_exact | allowed_type), "Unreviewed file format")

paths <- file.path(root, files)
info <- file.info(paths)
flag(is.na(info$size) | info$isdir, "Missing file or directory")
flag(!is.na(info$size) & info$size > 1024^2, "File exceeds the 1 MiB publication review limit")
if (!length(problems)) {
  binary <- vapply(paths, function(path) {
    con <- file(path, "rb")
    on.exit(close(con))
    any(readBin(con, "raw", n = file.info(path)$size) == as.raw(0))
  }, logical(1))
  flag(binary, "Binary content")
}
if (length(problems)) stop(paste(problems, collapse = "\n"), call. = FALSE)
cat("PUBLIC FILE CHECK PASSED:", length(files), "tracked text files;",
    sum(info$size), "bytes.\n")
