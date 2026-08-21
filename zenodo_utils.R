if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop(
    "Package 'jsonlite' is required. Install it with install.packages('jsonlite').",
    call. = FALSE
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

normalize_zenodo_identifier <- function(x) {
  trimws(x)
}

extract_zenodo_query <- function(x) {
  x <- normalize_zenodo_identifier(x)
  if (!grepl("?", x, fixed = TRUE)) {
    return(NA_character_)
  }

  query <- sub("^[^?]*[?]", "", x)
  if (!nzchar(query)) NA_character_ else query
}

append_zenodo_query <- function(url, query) {
  if (is.na(query) || !nzchar(query)) {
    return(url)
  }

  separator <- if (grepl("?", url, fixed = TRUE)) "&" else "?"
  paste0(url, separator, query)
}

extract_zenodo_doi <- function(x) {
  x <- normalize_zenodo_identifier(x)
  x <- sub("^https?://(dx\\.)?doi\\.org/", "", x, ignore.case = TRUE)
  if (grepl("^10\\.5281/zenodo\\.[[:digit:]]+$", x, ignore.case = TRUE)) {
    return(x)
  }
  NA_character_
}

extract_zenodo_record_id <- function(x) {
  x <- normalize_zenodo_identifier(x)

  if (grepl("^[[:digit:]]+$", x)) {
    return(x)
  }

  record_match <- regexec("/records/([[:digit:]]+)", x, ignore.case = TRUE)
  record_parts <- regmatches(x, record_match)[[1]]
  if (length(record_parts) >= 2) {
    return(record_parts[2])
  }

  doi_match <- regexec("zenodo\\.([[:digit:]]+)$", x, ignore.case = TRUE)
  doi_parts <- regmatches(x, doi_match)[[1]]
  if (length(doi_parts) >= 2) {
    return(doi_parts[2])
  }

  NA_character_
}

read_zenodo_json <- function(url) {
  jsonlite::fromJSON(url, simplifyVector = FALSE)
}

fetch_zenodo_record_by_id <- function(record_id, identifier = record_id) {
  query <- extract_zenodo_query(identifier)
  endpoint <- if (is.na(query)) {
    sprintf("https://zenodo.org/api/records/%s", record_id)
  } else {
    sprintf("https://zenodo.org/api/records/%s/draft", record_id)
  }
  url <- append_zenodo_query(endpoint, query)
  tryCatch(read_zenodo_json(url), error = function(e) NULL)
}

fetch_zenodo_record_by_doi <- function(doi) {
  query <- URLencode(sprintf('doi:"%s"', doi), reserved = TRUE)
  url <- sprintf("https://zenodo.org/api/records?q=%s&size=1", query)
  result <- tryCatch(read_zenodo_json(url), error = function(e) NULL)

  hits <- result$hits$hits %||% list()
  if (length(hits) == 0) {
    return(NULL)
  }

  hits[[1]]
}

resolve_zenodo_record <- function(identifier) {
  doi <- extract_zenodo_doi(identifier)
  if (!is.na(doi)) {
    record <- fetch_zenodo_record_by_doi(doi)
    if (!is.null(record)) {
      return(record)
    }
  }

  record_id <- extract_zenodo_record_id(identifier)
  if (!is.na(record_id)) {
    record <- fetch_zenodo_record_by_id(record_id, identifier)
    if (!is.null(record)) {
      return(record)
    }
  }

  stop(
    "Could not resolve the Zenodo identifier. Use a Zenodo record ID, version-specific DOI, or record URL.",
    call. = FALSE
  )
}

download_zenodo_inputs <- function(
  identifier,
  files = NULL,
  destination_dir = "input_folder",
  overwrite = FALSE
) {
  record <- resolve_zenodo_record(identifier)
  access_query <- extract_zenodo_query(identifier)
  record_files <- record$files %||% list()

  if (length(record_files) == 0) {
    stop("The Zenodo record does not contain downloadable files.", call. = FALSE)
  }

  available_files <- vapply(
    record_files,
    function(file_rec) file_rec$key %||% NA_character_,
    character(1)
  )

  if (!is.null(files)) {
    missing_files <- setdiff(files, available_files)
    if (length(missing_files) > 0) {
      stop(
        sprintf(
          "These requested files are not present in the Zenodo record: %s",
          paste(missing_files, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    keep_idx <- which(available_files %in% files)
    record_files <- record_files[keep_idx]
    available_files <- available_files[keep_idx]
  }

  dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)

  record_id <- as.character(record$id %||% extract_zenodo_record_id(identifier))
  record_doi <- record$doi %||% record$metadata$doi %||% NA_character_
  concept_doi <- record$conceptdoi %||% record$metadata$conceptdoi %||% NA_character_
  record_title <- record$title %||% record$metadata$title %||% "Untitled Zenodo record"

  message("Resolved Zenodo record: ", record_title)
  message("Record ID: ", record_id)

  if (!is.na(record_doi)) {
    message("Version DOI: ", record_doi)
  }

  if (!is.na(concept_doi) && !identical(concept_doi, record_doi)) {
    message("Concept DOI: ", concept_doi)
  }

  downloaded <- setNames(character(length(record_files)), available_files)

  for (i in seq_along(record_files)) {
    file_rec <- record_files[[i]]
    file_key <- available_files[[i]]
    target_path <- file.path(destination_dir, file_key)
    download_url <- file_rec$links$content %||%
      file_rec$links$download %||%
      file_rec$links$self

    if (is.null(download_url)) {
      stop("Missing download URL for file: ", file_key, call. = FALSE)
    }

    download_url <- append_zenodo_query(download_url, access_query)

    dir.create(dirname(target_path), recursive = TRUE, showWarnings = FALSE)

    if (file.exists(target_path) && !overwrite) {
      message("Skipping existing file: ", target_path)
      downloaded[[file_key]] <- normalizePath(target_path, winslash = "/", mustWork = FALSE)
      next
    }

    message("Downloading ", file_key, " ...")
    # Keep preview tokens out of console logs; progress is reported above.
    utils::download.file(download_url, destfile = target_path, mode = "wb", quiet = TRUE)
    downloaded[[file_key]] <- normalizePath(target_path, winslash = "/", mustWork = FALSE)
  }

  message("")
  message(
    "Downloaded ",
    length(downloaded),
    " file(s) into: ",
    normalizePath(destination_dir, winslash = "/", mustWork = FALSE)
  )

  invisible(downloaded)
}
