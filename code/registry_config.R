registry_locations <- list(
  registries = list(
    bt_simulation = list(
      work.dir = getwd(),
      file.dir = "bt_simulation"
    ),
    bt_simu_gp_uq = list(
      work.dir = getwd(),
      file.dir = "bt_simu_gp_uq"
    ),
    bt_simu_marss = list(
      work.dir = getwd(),
      file.dir = "bt_simu_marss"
    ),
    bt_realdata = list(
      work.dir = getwd(),
      file.dir = "bt_realdata"
    )
  )
)

get_registry_location <- function(name) {
  registries <- registry_locations$registries
  if (!name %in% names(registries)) {
    stop(sprintf("Unknown registry '%s'.", name))
  }

  location <- registries[[name]]
  required_fields <- c("work.dir", "file.dir")
  missing_fields <- setdiff(required_fields, names(location))
  if (length(missing_fields) > 0) {
    stop(sprintf(
      "Registry '%s' is missing required field(s): %s",
      name,
      paste(missing_fields, collapse = ", ")
    ))
  }

  location
}

get_registry_file_dir <- function(name) {
  get_registry_location(name)$file.dir
}
