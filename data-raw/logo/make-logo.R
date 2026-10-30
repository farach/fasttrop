# Make the fasttrop hex logo and the favicons of the website.
#
# Run this script in the package folder:
#   Rscript data-raw/logo/make-logo.R
#
# The script writes these files:
#   man/figures/logo.png      the hex logo for the README and the website
#   pkgdown/favicon/*         the icons for the browser tabs
#
# The design is a small panel in a hexagon. The rows are units and the columns
# are periods. The coral block shows the treated cells. An untreated cell has
# a stronger color when TROP gives it more weight: its unit is similar to the
# treated units, or its period is near the treated periods. The three short
# lines show speed.
#
# The script uses grid, ragg, systemfonts, svglite, and magick. The package
# does not import these packages. The font is Segoe UI Semibold, a Windows
# font. On other computers, the script uses the default sans font.

library(grid)

pal <- list(
  rim = "#F4EBDD",     # warm cream
  fill = "#0B4F5C",    # deep teal
  cell = "#F4EBDD",
  treated = "#F07A5A", # coral
  text = "#FBF6EC",
  streak = "#7FD3D6"   # light teal
)

font_family <- tryCatch({
  systemfonts::register_variant("fasttrop-logo", family = "Segoe UI", weight = "semibold")
  "fasttrop-logo"
}, error = function(e) "sans")

# Corners of a hexagon with a point at the top.
hex_xy <- function(r) {
  a <- (90 + 60 * 0:5) * pi / 180
  list(x = r * cos(a), y = r * sin(a))
}

# Rows are units (top to bottom) and columns are periods (left to right). The
# last two units are treated in the last three periods.
panel_cells <- function(n_row = 6, n_col = 8, treated_rows = 5:6, treated_cols = 6:8) {
  cells <- expand.grid(row = seq_len(n_row), col = seq_len(n_col))
  cells$treated <- cells$row %in% treated_rows & cells$col %in% treated_cols
  unit_weight <- c(0.45, 0.95, 0.60, 0.85, 0.70, 0.90)[cells$row]
  unit_weight[cells$row %in% treated_rows] <- 1
  time_weight <- exp(-0.35 * abs(cells$col - mean(treated_cols)))
  cells$alpha <- ifelse(cells$treated, 1, 0.16 + 0.84 * unit_weight * time_weight)
  cells
}

rounded_square <- function(cx, cy, size, radius, n = 8) {
  h <- size / 2 - radius
  corners <- rbind(c(h, h), c(-h, h), c(-h, -h), c(h, -h))
  pts <- lapply(1:4, function(k) {
    a <- seq(90 * (k - 1), 90 * k, length.out = n) * pi / 180
    cbind(cx + corners[k, 1] + radius * cos(a), cy + corners[k, 2] + radius * sin(a))
  })
  do.call(rbind, pts)
}

draw_panel <- function(cx, cy, size, gap, cells) {
  n_row <- max(cells$row)
  n_col <- max(cells$col)
  width <- n_col * size + (n_col - 1) * gap
  height <- n_row * size + (n_row - 1) * gap
  x <- cx - width / 2 + (cells$col - 1) * (size + gap) + size / 2
  y <- cy + height / 2 - (cells$row - 1) * (size + gap) - size / 2
  for (k in seq_len(nrow(cells))) {
    colour <- if (cells$treated[k]) pal$treated else pal$cell
    pts <- rounded_square(x[k], y[k], size, 0.22 * size)
    grid.polygon(pts[, 1], pts[, 2], default.units = "native",
                 gp = gpar(fill = grDevices::adjustcolor(colour, cells$alpha[k]), col = NA))
  }
  invisible(list(left = cx - width / 2, top = cy + height / 2, height = height))
}

draw_streaks <- function(right, cy, length, spacing, lwd) {
  for (k in 1:3) {
    len <- length * c(0.75, 1, 0.55)[k]
    y <- cy + (2 - k) * spacing
    grid.lines(c(right - len, right), c(y, y), default.units = "native",
               gp = gpar(col = pal$streak, lwd = lwd, lineend = "round", alpha = 0.9))
  }
}

# The full logo. Native units are pixels, with the origin at the center.
draw_logo <- function(width, height, res) {
  r <- height / 2 - 1
  pushViewport(viewport(xscale = c(-width / 2, width / 2), yscale = c(-height / 2, height / 2)))
  outer <- hex_xy(r)
  inner <- hex_xy(r - 0.06 * r)
  grid.polygon(outer$x, outer$y, default.units = "native", gp = gpar(fill = pal$rim, col = NA))
  grid.polygon(inner$x, inner$y, default.units = "native", gp = gpar(fill = pal$fill, col = NA))
  s <- r / 600
  px_to_pt <- 72 / res
  box <- draw_panel(cx = 40 * s, cy = 105 * s, size = 58 * s, gap = 13 * s, cells = panel_cells())
  draw_streaks(right = box$left - 28 * s, cy = 105 * s, length = 120 * s, spacing = 52 * s,
               lwd = 13 * s * px_to_pt / 0.75)
  grid.text("fasttrop", x = 0, y = -285 * s, default.units = "native",
            gp = gpar(col = pal$text, fontfamily = font_family, fontsize = 165 * s * px_to_pt))
  popViewport()
}

# The favicon: the hexagon and a smaller panel, without text.
draw_icon <- function(width, height) {
  r <- min(width / sqrt(3), height / 2) * 0.98
  pushViewport(viewport(xscale = c(-width / 2, width / 2), yscale = c(-height / 2, height / 2)))
  outer <- hex_xy(r)
  inner <- hex_xy(r - 0.08 * r)
  grid.polygon(outer$x, outer$y, default.units = "native", gp = gpar(fill = pal$rim, col = NA))
  grid.polygon(inner$x, inner$y, default.units = "native", gp = gpar(fill = pal$fill, col = NA))
  cells <- panel_cells(n_row = 4, n_col = 4, treated_rows = 3:4, treated_cols = 3:4)
  draw_panel(cx = 0, cy = 0, size = 0.22 * r, gap = 0.06 * r, cells = cells)
  popViewport()
}

png_file <- function(path, width, height, background, draw) {
  ragg::agg_png(path, width = width, height = height, background = background, res = 144)
  grid.newpage()
  draw()
  invisible(grDevices::dev.off())
  magick::image_read(path)
}

dir.create(file.path("man", "figures"), showWarnings = FALSE, recursive = TRUE)
fav <- file.path("pkgdown", "favicon")
dir.create(fav, showWarnings = FALSE, recursive = TRUE)

# Logo: draw at two times the size, then scale to 518 x 600 pixels.
logo <- png_file(tempfile(fileext = ".png"), 1036, 1200, "transparent",
                 function() draw_logo(1036, 1200, res = 144))
magick::image_write(magick::image_scale(logo, "518x600"), file.path("man", "figures", "logo.png"),
                    format = "png")

# Icon master at 512 x 512 pixels with a transparent background.
icon <- png_file(tempfile(fileext = ".png"), 512, 512, "transparent",
                 function() draw_icon(512, 512))

# Icons for the home screen of phones need a background and a margin.
icon_on_cream <- function(px, scale) {
  png_file(tempfile(fileext = ".png"), px, px, pal$rim, function() {
    pushViewport(viewport(width = scale, height = scale))
    draw_icon(px * scale, px * scale)
    popViewport()
  })
}

magick::image_write(magick::image_scale(icon, "96x96"), file.path(fav, "favicon-96x96.png"))
magick::image_write(icon_on_cream(180, 0.84), file.path(fav, "apple-touch-icon.png"))
magick::image_write(icon_on_cream(192, 0.72), file.path(fav, "web-app-manifest-192x192.png"))
magick::image_write(icon_on_cream(512, 0.72), file.path(fav, "web-app-manifest-512x512.png"))
ico <- magick::image_join(lapply(c("16x16", "32x32", "48x48"),
                                 function(g) magick::image_scale(icon, g)))
magick::image_write(ico, file.path(fav, "favicon.ico"), format = "ico")

svglite::svglite(file.path(fav, "favicon.svg"), width = 4, height = 4, bg = "transparent")
grid.newpage()
draw_icon(288, 288)
invisible(grDevices::dev.off())

writeLines(c(
  "{",
  "  \"name\": \"fasttrop\",",
  "  \"short_name\": \"fasttrop\",",
  "  \"icons\": [",
  "    {\"src\": \"web-app-manifest-192x192.png\", \"sizes\": \"192x192\", \"type\": \"image/png\", \"purpose\": \"maskable\"},",
  "    {\"src\": \"web-app-manifest-512x512.png\", \"sizes\": \"512x512\", \"type\": \"image/png\", \"purpose\": \"maskable\"}",
  "  ],",
  sprintf("  \"theme_color\": \"%s\",", pal$fill),
  sprintf("  \"background_color\": \"%s\",", pal$rim),
  "  \"display\": \"standalone\"",
  "}"
), file.path(fav, "site.webmanifest"))
