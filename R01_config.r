read_config <- function(path = here::here("config", "config.yaml")) {
  cfg <- yaml::read_yaml(path)
  cfg$project$input_dir <- here::here(cfg$project$input_dir)
  cfg$project$out_dir   <- here::here(cfg$project$out_dir)
  dir.create(cfg$project$out_dir, showWarnings = FALSE, recursive = TRUE)
  cfg
}