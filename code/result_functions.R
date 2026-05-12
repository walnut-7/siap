diff_label <- function(flag) {
  if (flag) "diff" else "nodiff"
}

result_path_simu <- function(stem, pdt = NULL, gp_diff_flag = NULL, marss_diff_flag = NULL, method = NULL, ext = ".rds") {
  if (!is.null(gp_diff_flag) & !is.null(marss_diff_flag)) {
    result_suffix <- paste0(
      "_gp_", diff_label(gp_diff_flag),
      "_marss_", diff_label(marss_diff_flag)
    )
  } else {
    result_suffix <- ""
  }
  
  if (is.null(pdt)) {
    return(
      file.path(
        "./output/simulation",
        paste0(
          if (is.null(method)) stem else paste0(method, "_", stem),
          result_suffix,
          ext
        )
      )
    )
  }
  
  file.path(
    "./output/simulation",
    paste0(
      if (is.null(method)) stem else paste0(method, "_", stem),
      result_suffix,
      "_", 
      pdt * 10,
      ext
    )
  )
}
