save_plot <- function(cfg, filename, plot = ggplot2::last_plot(), width = 8, height = 6, dpi = 300) {
  if (isTRUE(cfg$qc$save_plots)) {
    ggplot2::ggsave(filename = file.path(cfg$project$out_dir, filename),
                    plot = plot, width = width, height = height, dpi = dpi)
  }
}

write_csv <- function(cfg, df, filename) {
  utils::write.csv(df, file.path(cfg$project$out_dir, filename), row.names = FALSE)
}

stop_if_missing <- function(paths, label = "file") {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop(sprintf("Missing %s(s):\n%s", label, paste(missing, collapse = "\n")))
  }
}