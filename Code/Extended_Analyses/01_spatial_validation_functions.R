read_truth <- function(path) {
  truth <- read.delim(
    path, check.names = FALSE, stringsAsFactors = FALSE
  )
  if (!all(c("Hclone", "Lclone") %in% names(truth))) {
    stop("Unexpected ground-truth columns in ", path, call. = FALSE)
  }
  truth
}

make_adjacency <- function(coords, radius) {
  if (!nrow(coords)) return(matrix(FALSE, 0L, 0L))
  pmax(
    abs(outer(coords$row, coords$row, "-")),
    abs(outer(coords$col, coords$col, "-"))
  ) <= radius
}

count_components <- function(selected, adjacency) {
  selected <- unique(as.integer(selected))
  if (!length(selected)) return(0L)
  if (length(selected) == 1L) return(1L)
  graph <- adjacency[selected, selected, drop = FALSE]
  visited <- rep(FALSE, length(selected))
  count <- 0L
  for (start in seq_along(selected)) {
    if (visited[[start]]) next
    count <- count + 1L
    queue <- start
    visited[[start]] <- TRUE
    while (length(queue)) {
      current <- queue[[1L]]
      queue <- queue[-1L]
      neighbours <- which(graph[current, ] & !visited)
      if (length(neighbours)) {
        visited[neighbours] <- TRUE
        queue <- c(queue, neighbours)
      }
    }
  }
  count
}

repeat_null <- function(k_h, k_l, n, adjacency, observed, n_perm) {
  overlap <- rhyper(n_perm, m = k_h, n = n - k_h, k = k_l)
  simulated <- matrix(
    0L, nrow = n_perm, ncol = length(adjacency),
    dimnames = list(NULL, names(adjacency))
  )
  one <- which(overlap == 1L)
  if (length(one)) simulated[one, ] <- 1L
  many <- which(overlap >= 2L)
  if (length(many)) {
    for (j in many) {
      selected <- sample.int(k_h, overlap[[j]], replace = FALSE)
      for (radius_name in names(adjacency)) {
        simulated[j, radius_name] <- count_components(
          selected, adjacency[[radius_name]]
        )
      }
    }
  }
  list(
    p = vapply(names(adjacency), function(radius_name) {
      (1 + sum(simulated[, radius_name] >= observed[[radius_name]])) /
        (n_perm + 1)
    }, numeric(1)),
    expected = colMeans(simulated)
  )
}

read_sample <- function(sample_name) {
  paths <- sample_paths[[sample_name]]
  spatial <- fread(paths[["info"]], select = c("type", "locus", "clone_id", "binid"))
  spatial <- spatial[grepl("Spatial", type)]
  spatial[, clone_name := paste(locus, as.character(clone_id), sep = "_")]
  spatial <- unique(spatial[, .(clone_name, binid = as.character(binid))])

  metadata <- fread(paths[["metadata"]], select = c("BinID", "row", "col", "tissue"))
  metadata[, binid := sub(paste0("^", sample_name, "_"), "", as.character(BinID))]
  metadata <- unique(metadata[, .(binid, row, col, tissue)])
  if (anyDuplicated(metadata$binid) || anyDuplicated(metadata[, .(row, col)])) {
    stop("Duplicated bin coordinates in ", sample_name, call. = FALSE)
  }

  final <- fread(paths[["final"]])
  final <- final[!is.na(Lclone_t), .(H = Hclone, L = Lclone)]
  test <- fread(paths[["test"]], select = c("Hclone", "Lclone"))
  test <- unique(test[, .(H = Hclone, L = Lclone)])
  if (anyDuplicated(final) || anyDuplicated(test) ||
      nrow(final) != nrow(test) ||
      !setequal(paste(final$H, final$L), paste(test$H, test$L))) {
    stop("final_e100 and test_e100 disagree for ", sample_name, call. = FALSE)
  }

  setorder(test, H, L)

  if (!all(unique(spatial$binid) %in% metadata$binid)) {
    stop("Spatial bins missing from metadata for ", sample_name, call. = FALSE)
  }
  list(spatial = spatial, metadata = metadata, pairs = test)
}

clip_probability <- function(p, epsilon = 1e-8) {
  pmin(pmax(p, epsilon), 1 - epsilon)
}

equal_sample_weights <- function(sample_name) {
  counts <- table(sample_name)
  weights <- 1 / as.numeric(counts[as.character(sample_name)])
  weights / mean(weights)
}

make_grouped_folds <- function(sample_name, k, seed) {
  group_sizes <- sort(table(sample_name), decreasing = TRUE)
  groups <- names(group_sizes)
  k <- min(as.integer(k), length(groups))
  set.seed(seed)
  groups <- groups[order(-as.numeric(group_sizes), runif(length(groups)))]
  group_sizes <- as.numeric(table(sample_name)[groups])
  fold_load <- rep(0, k)
  group_fold <- setNames(integer(length(groups)), groups)
  for (i in seq_along(groups)) {
    available <- which(fold_load == min(fold_load))
    selected <- sample(available, 1L)
    group_fold[[groups[[i]]]] <- selected
    fold_load[[selected]] <- fold_load[[selected]] + group_sizes[[i]]
  }
  unname(group_fold[as.character(sample_name)])
}

standardize_train_test <- function(x_train, x_test = NULL) {
  center <- colMeans(x_train)
  scale <- apply(x_train, 2, sd)
  scale[!is.finite(scale) | scale < 1e-8] <- 1
  list(
    train = sweep(sweep(x_train, 2, center, "-"), 2, scale, "/"),
    test = if (is.null(x_test)) NULL else
      sweep(sweep(x_test, 2, center, "-"), 2, scale, "/"),
    center = center, scale = scale
  )
}

fit_ridge <- function(x, y, weights, lambda, start = NULL) {
  design <- cbind("(Intercept)" = 1, x)
  if (is.null(start)) start <- rep(0, ncol(design))
  weight_sum <- sum(weights)
  objective <- function(beta) {
    eta <- drop(design %*% beta)
    softplus <- pmax(eta, 0) + log1p(exp(-abs(eta)))
    sum(weights * (softplus - y * eta)) / weight_sum +
      0.5 * lambda * sum(beta[-1L]^2)
  }
  gradient <- function(beta) {
    probability <- plogis(drop(design %*% beta))
    drop(crossprod(design, weights * (probability - y))) / weight_sum +
      c(0, lambda * beta[-1L])
  }
  fit <- optim(
    start, objective, gradient, method = "BFGS",
    control = list(maxit = 1000L, reltol = 1e-9)
  )
  list(coefficients = fit$par, convergence = fit$convergence)
}

predict_ridge <- function(fit, x) {
  if (!nrow(x)) return(numeric())
  plogis(drop(cbind("(Intercept)" = 1, x) %*% fit$coefficients))
}

fit_ridge_path <- function(x, y, weights, lambdas) {
  fits <- vector("list", length(lambdas))
  start <- NULL
  for (j in order(lambdas, decreasing = TRUE)) {
    fits[[j]] <- fit_ridge(x, y, weights, lambdas[[j]], start)
    start <- fits[[j]]$coefficients
  }
  fits
}

select_lambda <- function(x, y, sample_name, lambdas, k, seed) {
  fold_id <- make_grouped_folds(sample_name, k, seed)
  oof <- matrix(NA_real_, nrow(x), length(lambdas))
  convergence <- matrix(NA_integer_, max(fold_id), length(lambdas))
  for (fold in seq_len(max(fold_id))) {
    train <- fold_id != fold
    validation <- !train
    scaled <- standardize_train_test(
      x[train, , drop = FALSE], x[validation, , drop = FALSE]
    )
    path <- fit_ridge_path(
      scaled$train, y[train], equal_sample_weights(sample_name[train]), lambdas
    )
    for (j in seq_along(lambdas)) {
      oof[validation, j] <- predict_ridge(path[[j]], scaled$test)
      convergence[fold, j] <- path[[j]]$convergence
    }
  }
  loss <- vapply(seq_along(lambdas), function(j) {
    row_loss <- -(
      y * log(clip_probability(oof[, j])) +
        (1 - y) * log(1 - clip_probability(oof[, j]))
    )
    mean(vapply(split(row_loss, sample_name), mean, numeric(1)))
  }, numeric(1))
  selected <- which.min(loss)
  list(
    lambda = lambdas[[selected]], loss = loss[[selected]],
    converged = all(convergence[, selected] == 0L)
  )
}

theme_report <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      plot.title = element_text(face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

save_figure <- function(plot, path, width, height) {
  ggsave(paste0(path, ".pdf"), plot, width = width, height = height,
         device = grDevices::cairo_pdf, bg = "white")
  ggsave(paste0(path, ".png"), plot, width = width, height = height,
         dpi = 600, bg = "white", limitsize = FALSE)
}

format_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

make_schematic <- function(label, h_cells, l_cells) {
  grid <- expand.grid(x = 1:10, y = 1:8, KEEP.OUT.ATTRS = FALSE)
  key <- paste(grid$x, grid$y, sep = "_")
  h_key <- paste(h_cells$x, h_cells$y, sep = "_")
  l_key <- paste(l_cells$x, l_cells$y, sep = "_")
  grid$state <- case_when(
    key %in% intersect(h_key, l_key) ~ "H–L overlap",
    key %in% h_key ~ "H only",
    key %in% l_key ~ "L only",
    TRUE ~ "No signal"
  )
  grid$scenario <- label
  grid
}

make_confusion_data <- function(counts) {
  tibble::tibble(
    evidence = factor(c("Support", "Support", "No support", "No support"),
                      levels = c("No support", "Support")),
    correctness = factor(rep(c("Correct", "Incorrect"), 2),
                         levels = c("Correct", "Incorrect")),
    n = counts
  ) %>%
    mutate(
      percent = 100 * n / sum(n),
      label = sprintf("%d\n(%.1f%%)", n, percent),
      label_color = if_else(n >= 320, "white", "black")
    )
}

make_confusion_plot <- function(counts, title, subtitle, y_title) {
  ggplot(make_confusion_data(counts), aes(correctness, evidence, fill = n)) +
    geom_tile(color = "white", linewidth = 1.1) +
    geom_text(aes(label = label, color = label_color),
              size = 4.2 * font_scale, lineheight = 0.95,
              fontface = "bold", show.legend = FALSE) +
    scale_fill_gradientn(
      colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
      limits = c(0, 563), oob = scales::squish, guide = "none"
    ) +
    scale_color_identity() +
    coord_cartesian(expand = FALSE) +
    labs(title = title, subtitle = subtitle,
         x = "SIPEAC prediction vs ground truth", y = y_title) +
    theme_report(9) +
    theme(
      axis.line = element_blank(), axis.ticks = element_blank(),
      plot.title = element_text(size = 11 * font_scale),
      plot.subtitle = element_text(size = 8.2 * font_scale, color = "#3F3F3F"),
      axis.title = element_text(size = 8.5 * font_scale),
      axis.text = element_text(size = 8.3 * font_scale),
      plot.background = element_rect(fill = "#FCFCFC", color = "#D9D9D9", linewidth = 0.45),
      plot.margin = margin(4, 4, 4, 4)
    )
}

component_membership <- function(coords, radius = 2L) {
  if (!nrow(coords)) return(integer())
  adjacency <- make_adjacency(coords, radius)
  component <- rep(NA_integer_, nrow(coords))
  count <- 0L
  for (start in seq_len(nrow(coords))) {
    if (!is.na(component[[start]])) next
    count <- count + 1L
    component[[start]] <- count
    queue <- start
    while (length(queue)) {
      current <- queue[[1L]]
      queue <- queue[-1L]
      neighbours <- which(adjacency[current, ] & is.na(component))
      if (length(neighbours)) {
        component[neighbours] <- count
        queue <- c(queue, neighbours)
      }
    }
  }
  component
}

square_extent <- function(extent) {
  side <- max(extent[["xmax"]] - extent[["xmin"]],
              extent[["ymax"]] - extent[["ymin"]])
  xmid <- mean(extent[c("xmin", "xmax")])
  ymid <- mean(extent[c("ymin", "ymax")])
  c(xmin = xmid - side / 2, xmax = xmid + side / 2,
    ymin = ymid - side / 2, ymax = ymid + side / 2)
}

spatial_pair <- function(example) {
  dat <- read_sample(example$sample_name)
  h_bins <- unique(dat$spatial[clone_name == example$H, binid])
  l_bins <- unique(dat$spatial[clone_name == example$L, binid])
  overlap_bins <- intersect(h_bins, l_bins)
  signal <- copy(dat$metadata[binid %in% union(h_bins, l_bins)])
  signal[, state := case_when(
    binid %in% overlap_bins ~ "H–L overlap",
    binid %in% h_bins ~ "H only",
    TRUE ~ "L only"
  )]
  signal[, state := factor(state, levels = c("H only", "L only", "H–L overlap"))]
  overlap <- copy(dat$metadata[binid %in% overlap_bins])
  overlap[, event := component_membership(as.data.frame(overlap[, .(row, col)]), 2L)]
  event_centres <- overlap[, .(col = mean(col), row = mean(row), n_overlap_bins = .N),
                           by = event][order(event)]
  event_centres[, event_label := paste0("E", event)]
  if (nrow(event_centres) != example$R_r2) {
    stop("Spatial event check failed for ", example$sample_name, call. = FALSE)
  }

  min_row <- min(dat$metadata$row)
  max_row <- max(dat$metadata$row)
  min_col <- min(dat$metadata$col)
  max_col <- max(dat$metadata$col)
  raster <- matrix("#FFFFFF00", nrow = max_row - min_row + 1L,
                   ncol = max_col - min_col + 1L)
  raster[cbind(dat$metadata$row - min_row + 1L,
               dat$metadata$col - min_col + 1L)] <- "#E8E8E8"
  crop <- square_extent(c(
    xmin = min(overlap$col) - 14, xmax = max(overlap$col) + 14,
    ymin = min(overlap$row) - 14, ymax = max(overlap$row) + 14
  ))
  raster_extent <- c(
    xmin = min_col - 0.5, xmax = max_col + 0.5,
    ymin = -max_row - 0.5, ymax = -min_row + 0.5
  )
  list(
    signal = signal, event_centres = event_centres, raster = raster,
    raster_extent = raster_extent, plot_extent = square_extent(raster_extent), crop = crop
  )
}

spatial_base <- function(spatial, crop = FALSE) {
  plot <- ggplot() +
    annotation_raster(
      spatial$raster,
      xmin = spatial$raster_extent[["xmin"]], xmax = spatial$raster_extent[["xmax"]],
      ymin = spatial$raster_extent[["ymin"]], ymax = spatial$raster_extent[["ymax"]]
    ) +
    geom_point(
      data = spatial$signal, aes(col, -row, color = state),
      size = if (crop) 2.0 * font_scale else 1.15 * font_scale, alpha = 0.95
    ) +
    scale_color_manual(values = signal_colors, breaks = names(signal_colors), drop = FALSE) +
    theme_void(base_size = 8) +
    theme(legend.position = "none", plot.margin = margin(1, 1, 1, 1))

  if (crop) {
    plot +
      geom_point(
        data = spatial$event_centres, aes(col, -row), shape = 21,
        size = 5 * font_scale, stroke = 0.85, color = "#B7791F", fill = NA
      ) +
      geom_text(
        data = spatial$event_centres, aes(col, -row, label = event_label),
        color = "#7A4E00", fontface = "bold", size = 3.1 * font_scale, nudge_y = 5
      ) +
      coord_equal(xlim = spatial$crop[c("xmin", "xmax")],
                  ylim = -rev(spatial$crop[c("ymin", "ymax")]), expand = FALSE) +
      theme(panel.border = element_rect(color = "#B7791F", fill = NA, linewidth = 0.55))
  } else {
    plot +
      geom_rect(
        aes(xmin = spatial$crop[["xmin"]], xmax = spatial$crop[["xmax"]],
            ymin = -spatial$crop[["ymax"]], ymax = -spatial$crop[["ymin"]]),
        inherit.aes = FALSE, fill = NA, color = "#B7791F",
        linewidth = 0.55, linetype = "dashed"
      ) +
      coord_equal(xlim = spatial$plot_extent[c("xmin", "xmax")],
                  ylim = spatial$plot_extent[c("ymin", "ymax")], expand = FALSE)
  }
}

make_spatial_panel <- function(example, panel_letter) {
  spatial <- spatial_pair(example)
  full <- spatial_base(spatial, FALSE) +
    labs(title = "Whole tissue") +
    theme(plot.title = element_text(size = 7.5 * font_scale, face = "bold", hjust = 0.5))
  local <- spatial_base(spatial, TRUE) +
    labs(title = "Same region enlarged") +
    theme(plot.title = element_text(size = 7.5 * font_scale, face = "bold", hjust = 0.5))
  bridge <- ggplot() +
    geom_segment(
      aes(x = 0.08, y = 0.5, xend = 0.92, yend = 0.5),
      linewidth = 0.75, color = "#B7791F",
      arrow = grid::arrow(length = grid::unit(5, "pt"), type = "closed")
    ) +
    annotate("text", x = 0.5, y = 0.64, label = "same\nregion",
             size = 2.4 * font_scale, color = "#7A4E00", lineheight = 0.9) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void()
  maps <- full + bridge + local + plot_layout(widths = c(1, 0.13, 1))
  legend_data <- data.frame(
    x = c(1, 2.15, 3.25), y = 1,
    state = factor(names(signal_colors), levels = names(signal_colors)),
    label = names(signal_colors)
  )
  legend <- ggplot(legend_data, aes(x, y, color = state)) +
    geom_point(size = 2.5) +
    geom_text(aes(label = label), nudge_x = 0.13, hjust = 0,
              color = "#2F2F2F", size = 2.8 * font_scale) +
    scale_color_manual(values = signal_colors, guide = "none") +
    coord_cartesian(xlim = c(0.8, 4.2), ylim = c(0.8, 1.2)) +
    theme_void()

  probability <- if (is.na(example$calibration_probability)) {
    "Tier 2 not evaluated"
  } else {
    sprintf("outer-LOSO P(correct) = %.3f", example$calibration_probability)
  }
  title <- switch(
    as.character(example$example_group),
    "Tier 1 support / Correct" = paste0(panel_letter, "  Tier 1 support"),
    "Tier 2 rescue / Correct" = paste0(panel_letter, "  Tier 2 rescue"),
    "Final no support / Correct" = paste0(panel_letter, "  No support")
  )
  caption <- switch(
    as.character(example$example_group),
    "Tier 1 support / Correct" = "Four overlap bins form three separated events",
    "Tier 2 rescue / Correct" = "One event fails Tier 1; the spatial profile is rescued",
    "Final no support / Correct" =
      "Correct SIPEAC pair, but spatial evidence remains insufficient"
  )
  (maps / legend + plot_layout(heights = c(1, 0.085))) +
    plot_annotation(
      title = title,
      subtitle = sprintf(
        "%s | O = %d; R2 = %d; q = %s; %s\n%s × %s | %s",
        example$sample_name, example$O_ij, example$R_r2,
        format_p(example$q_repeat_r2), probability,
        example$H, example$L,
        ifelse(example$final_support, "Final support", "Final no support")
      ),
      caption = caption,
      theme = theme(
        plot.title = element_text(face = "bold", size = 10.5 * font_scale),
        plot.subtitle = element_text(size = 6.8 * font_scale, color = "#3F3F3F"),
        plot.caption = element_text(size = 7.5 * font_scale, hjust = 0,
                                    color = "#3F3F3F"),
        plot.background = element_rect(fill = "#FCFCFC", color = "#D9D9D9", linewidth = 0.45),
        plot.margin = margin(4, 4, 3, 4)
      )
    )
}
