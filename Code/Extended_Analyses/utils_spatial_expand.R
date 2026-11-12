.st_expand_adjacent_1 <- function(df_location,
                                  x = "col",
                                  y = "row",
                                  mode = 1) {
  df <- df_location %>%
    dplyr::relocate(setNames(c(x, y), c("x", "y")))
  df_adjacent <- dplyr::bind_rows(
    df %>% dplyr::mutate(x = x + 1, y = y + 0),
    df %>% dplyr::mutate(x = x - 1, y = y + 0),
    df %>% dplyr::mutate(x = x + 0, y = y + 1),
    df %>% dplyr::mutate(x = x + 0, y = y - 1)
  )
  if (mode == 1) {
    df_adjacent <- df_adjacent
  } else if (mode == 2) {
    df_adjacent <- dplyr::bind_rows(
      df_adjacent,
      df %>% dplyr::mutate(x = x + 1, y = y + 1),
      df %>% dplyr::mutate(x = x + 1, y = y - 1),
      df %>% dplyr::mutate(x = x - 1, y = y + 1),
      df %>% dplyr::mutate(x = x - 1, y = y - 1)
    )
  } else {
    stop("mode must be 1 or 2!", call. = FALSE)
  }
  df_adjacent <- df_adjacent %>%
    dplyr::relocate(setNames(c("x", "y"), c(x, y)))
  return(df_adjacent)
}

.st_expand_1 <- function(df_location,
                         x = "col",
                         y = "row",
                         expand_n = "expand_n",
                         mode = 1) {
  if (!expand_n %in% colnames(df_location)) {df_location[[expand_n]] <- 0}
  df <- df_location %>%
    dplyr::relocate(setNames(c(x, y, expand_n), c("x", "y", "expand_n")))
  df_expand <- .st_expand_adjacent_1(
    df_location = df %>% dplyr::filter(expand_n == max(expand_n)),
    x = "x",
    y = "y",
    mode = mode
  ) %>%
    dplyr::anti_join(df, by = c("x", "y")) %>%
    dplyr::mutate(expand_n = expand_n + 1)
  df_expand <- df_expand %>%
    dplyr::relocate(setNames(c("x", "y", "expand_n"), c(x, y, expand_n)))
  return(df_expand)
}

.st_expand_add_1 <- function(df_location,
                             x = "col",
                             y = "row",
                             expand_n = "expand_n",
                             group.by = NULL,
                             mode = 1) {
  if (is.null(group.by)) {
    group.by <- "group.by"
    df_location$group.by <- 1
  }
  if (!expand_n %in% colnames(df_location)) {df_location[[expand_n]] <- 0}
  df <- df_location %>%
    dplyr::relocate(setNames(c(x, y, expand_n,group.by),
                             c("x", "y", "expand_n", "group.by")))
  df_expand_add <- dplyr::bind_rows(
    df,
    df %>%
      dplyr::group_split(group.by) %>%
      purrr::map(
        ~ .st_expand_1(df_location = .x,
                       x = "x",
                       y = "y",
                       expand_n = "expand_n",
                       mode = mode),
        .progress = TRUE
      ) %>%
      dplyr::bind_rows()
  )
  df_expand_add <- df_expand_add %>%
    dplyr::relocate(setNames(c("x", "y", "expand_n", "group.by"),
                             c(x, y, expand_n, group.by)))
  return(df_expand_add)
}

st_expand_n <- function(df_location,
                        n,
                        x = "col",
                        y = "row",
                        expand_n = "expand_n",
                        group.by = NULL,
                        mode = 1) {
  if (is.null(group.by)) {
    group.by <- "group.by"
    df_location$group.by <- 1
  }
  if (!expand_n %in% colnames(df_location)) {df_location[[expand_n]] <- 0}
  df <- df_location %>%
    dplyr::select(setNames(c(x, y, expand_n,group.by),
                           c("x", "y", "expand_n", "group.by"))) %>%
    dplyr::distinct()
  df_center <- df %>%
    dplyr::group_by(group.by) %>%
    dplyr::summarise(xc = mean(range(x)), yc = mean(range(y)))
  cli::cli_progress_bar("Expanding Layers ...",
                        total = n, clear = FALSE)
  for (i in seq_len(n)) {
    df <- .st_expand_add_1(df_location = df,
                           x = "x",
                           y = "y",
                           expand_n = "expand_n",
                           group.by = "group.by",
                           mode = mode) %>%
      dplyr::distinct()
    cli::cli_progress_update()
  }
  cli::cli_progress_done()
  df_expand <- df %>%
    dplyr::group_by(x, y) %>%
    dplyr::filter(expand_n == min(expand_n)) %>%
    dplyr::mutate(n = dplyr::n_distinct(group.by)) %>%
    dplyr::ungroup()
  if (any(df_expand$n > 1)) {
    df_expand_single <- df_expand %>%
      dplyr::filter(n == 1)
    df_expand_multiple <- df_expand %>%
      dplyr::filter(n > 1) %>%
      dplyr::left_join(df_center, by = "group.by") %>%
      dplyr::mutate(dist = (x - xc)^2 + (y - yc) ^ 2) %>%
      dplyr::group_by(x, y) %>%
      dplyr::filter(dist == min(dist)) %>%
      dplyr::filter(dplyr::n_distinct(group.by) == 1)
    df_expand <- dplyr::bind_rows(df_expand_single, df_expand_multiple)
  }
  df_expand <- df_expand %>%
    dplyr::distinct(x, y, expand_n, group.by) %>%
    dplyr::rename(setNames(c("x", "y"), c(x, y)))
  return(df_expand)
}
