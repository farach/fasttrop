# Shared settings for the charts in the website articles of fasttrop. Each
# article runs this file in a hidden chunk. The website does not show the
# chart code. The code is in the source files of the articles.

library(ggplot2)

# The charts use the body font of the website. If the download of the font
# is not possible, the charts use the default sans font.
chart_font <- tryCatch({
  systemfonts::require_font("Atkinson Hyperlegible Next", verbose = FALSE)
  "Atkinson Hyperlegible Next"
}, error = function(e) "sans")

# Colors of the website. They are safe for the common types of color
# blindness when the charts use them as pairs: teal with coral, or dark with
# light.
ft_col <- c(
  ink = "#1E2A32", muted = "#55636A", grid = "#E3DCCF",
  teal = "#0B6E79", teal_light = "#8FCFD4", teal_dark = "#0B4F5C",
  coral = "#E4572E", coral_light = "#F4A58A", sand = "#D8CDBB", grey = "#9AA6AB"
)

theme_ft <- function(base_size = 15) {
  theme_minimal(base_size = base_size, base_family = chart_font) +
    theme(
      text = element_text(colour = ft_col[["ink"]]),
      axis.text = element_text(colour = ft_col[["muted"]]),
      axis.title = element_text(colour = ft_col[["muted"]], size = rel(0.9)),
      panel.grid.major = element_line(colour = ft_col[["grid"]], linewidth = 0.4),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_blank(),
      strip.text = element_text(colour = ft_col[["ink"]], hjust = 0, size = rel(1)),
      plot.margin = margin(8, 14, 8, 8),
      plot.background = element_blank(),
      panel.background = element_blank()
    )
}

theme_set(theme_ft())
update_geom_defaults("text", list(family = chart_font, colour = ft_col[["ink"]]))
update_geom_defaults("label", list(family = chart_font, colour = ft_col[["ink"]]))

# Weights of the local method for one treated cell (i, t), from Eq. 3 of the
# paper. The time weight is exp(-lambda_time * |t - s|). The unit weight is
# exp(-lambda_unit * d), where d is the root mean square difference between
# the outcomes of unit j and unit i. The difference uses the periods in which
# both units are untreated, but not period t. Treated cells get no weight.
local_weights <- function(Y, W, i, t, lambda_time, lambda_unit) {
  untreated <- W == 0 & is.finite(Y)
  theta <- exp(-lambda_time * abs(t - seq_len(ncol(Y))))
  omega <- vapply(seq_len(nrow(Y)), function(j) {
    if (j == i || lambda_unit == 0) return(1)
    s <- setdiff(which(untreated[i, ] & untreated[j, ]), t)
    if (length(s) == 0L) return(0)
    exp(-lambda_unit * sqrt(mean((Y[i, s] - Y[j, s])^2)))
  }, double(1))
  w <- outer(omega, theta) * untreated
  list(weights = w, unit = omega, time = theta)
}
