# =============================================================================
# Create hex sticker for cfperformance
# =============================================================================

# Install hexSticker if needed
# install.packages("hexSticker")
# install.packages("showtext")

library(hexSticker)
library(ggplot2)
library(showtext)

# Enable showtext for better fonts
font_add_google("Fira Sans", "firasans")
showtext_auto()

# -----------------------------------------------------------------------------
# Design 1: ROC Curve with Counterfactual Split
# Shows an ROC curve splitting into two paths (treatment vs control)
# -----------------------------------------------------------------------------

create_roc_design <- function() {
  # Create data for the curves
  x <- seq(0, 1, length.out = 100)
  
  # Main ROC curve (observed)
  roc_main <- data.frame(
    x = x,
    y = x^0.35,
    group = "main"
  )
  
  # Counterfactual ROC curve (what could have been)
  roc_cf <- data.frame(
    x = x,
    y = x^0.5,
    group = "cf"
  )
  
  p <- ggplot() +
    # Reference diagonal
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                 linetype = "dashed", color = "gray60", linewidth = 0.5) +
    # Counterfactual ROC (dashed, representing what we estimate)
    geom_line(data = roc_cf, aes(x = x, y = y), 
              color = "#7ECCE8", linewidth = 1.2, linetype = "dashed") +
    # Observed ROC
    geom_line(data = roc_main, aes(x = x, y = y), 
              color = "#FFFFFF", linewidth = 1.5) +
    # Add small arrows to suggest causality
    annotate("segment", x = 0.5, y = 0.5^0.35, xend = 0.5, yend = 0.5^0.5,
             arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
             color = "#F4D35E", linewidth = 0.8) +
    coord_fixed(ratio = 1) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(legend.position = "none")
  
  return(p)
}

# -----------------------------------------------------------------------------
# Design 2: Stylized "cf" with Performance Curve
# -----------------------------------------------------------------------------

create_cf_design <- function() {
  # Create an elegant curve representing performance
  x <- seq(0, 1, length.out = 100)
  
  p <- ggplot() +
    # Calibration-style curve
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                 linetype = "dashed", color = "gray50", linewidth = 0.4) +
    # Performance curve
    stat_function(fun = function(x) pnorm(qnorm(x) + 0.5), 
                  color = "#FFFFFF", linewidth = 1.3, n = 200) +
    # Second curve for counterfactual
    stat_function(fun = function(x) pnorm(qnorm(x) + 1.0), 
                  color = "#7ECCE8", linewidth = 1.0, linetype = "dashed", n = 200) +
    coord_fixed(ratio = 1) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void()
  
  return(p)
}

# -----------------------------------------------------------------------------
# Design 3: Abstract Causal Paths
# Two diverging paths representing treatment/control outcomes
# -----------------------------------------------------------------------------

create_causal_design <- function() {
  # Create diverging paths
  t <- seq(0, 1, length.out = 50)
  
  # Path 1 (treatment outcome trajectory)
  path1 <- data.frame(
    x = t,
    y = 0.5 + 0.3 * sin(t * pi) + 0.15 * t,
    group = "treated"
  )
  
  # Path 2 (control outcome trajectory)
  path2 <- data.frame(
    x = t,
    y = 0.5 + 0.15 * sin(t * pi) - 0.1 * t,
    group = "control"
  )
  
  # Starting point
  start_point <- data.frame(x = 0, y = 0.5)
  
  p <- ggplot() +
    geom_line(data = path1, aes(x = x, y = y), 
              color = "#FFFFFF", linewidth = 1.5) +
    geom_line(data = path2, aes(x = x, y = y), 
              color = "#7ECCE8", linewidth = 1.2, linetype = "longdash") +
    geom_point(data = start_point, aes(x = x, y = y), 
               color = "#F4D35E", size = 3) +
    xlim(-0.1, 1.1) + ylim(0, 1) +
    theme_void()
  
  return(p)
}

# -----------------------------------------------------------------------------
# Create the hex stickers
# -----------------------------------------------------------------------------

# Color palette (medical/clinical blues)
bg_color <- "#1A365D"      # Deep navy
border_color <- "#2E86AB"  # Bright blue
accent_color <- "#F4D35E"  # Gold accent

# Design 1: ROC curve (recommended)
p1 <- create_roc_design()

sticker(p1, 
        package = "cfperformance", 
        p_size = 16, 
        p_y = 1.45,
        p_color = "#FFFFFF",
        p_family = "firasans",
        s_x = 1, 
        s_y = 0.75, 
        s_width = 1.4, 
        s_height = 1.1,
        h_fill = bg_color, 
        h_color = border_color,
        h_size = 1.5,
        url = "github.com/boyercb/cfperformance",
        u_size = 3.5,
        u_color = "#7ECCE8",
        filename = "inst/hex/cfperformance_hex_roc.png",
        dpi = 300)

# Design 2: Calibration/performance curve
p2 <- create_cf_design()

sticker(p2, 
        package = "cfperformance", 
        p_size = 16, 
        p_y = 1.45,
        p_color = "#FFFFFF",
        p_family = "firasans",
        s_x = 1, 
        s_y = 0.75, 
        s_width = 1.4, 
        s_height = 1.1,
        h_fill = bg_color, 
        h_color = border_color,
        h_size = 1.5,
        url = "github.com/boyercb/cfperformance",
        u_size = 3.5,
        u_color = "#7ECCE8",
        filename = "inst/hex/cfperformance_hex_calib.png",
        dpi = 300)

# Design 3: Causal paths
p3 <- create_causal_design()

sticker(p3, 
        package = "cfperformance", 
        p_size = 16, 
        p_y = 1.45,
        p_color = "#FFFFFF",
        p_family = "firasans",
        s_x = 1, 
        s_y = 0.75, 
        s_width = 1.4, 
        s_height = 1.1,
        h_fill = bg_color, 
        h_color = border_color,
        h_size = 1.5,
        url = "github.com/boyercb/cfperformance",
        u_size = 3.5,
        u_color = "#7ECCE8",
        filename = "inst/hex/cfperformance_hex_causal.png",
        dpi = 300)

# -----------------------------------------------------------------------------
# Print summary
# -----------------------------------------------------------------------------

cat("\n")
cat("==============================================\n")
cat("Hex stickers created in inst/hex/\n")
cat("==============================================\n")
cat("\n")
cat("Files:\n")
cat("  - cfperformance_hex_roc.png   (ROC curve design)\n")
cat("  - cfperformance_hex_calib.png (Calibration curve design)\n")
cat("  - cfperformance_hex_causal.png (Causal paths design)\n")
cat("\n")
cat("To use in README, add to man/figures/ and reference:\n")
cat('  <img src="man/figures/cfperformance_hex.png" align="right" height="139" />\n')
cat("\n")
