library(readxl)
library(dplyr)
library(ggplot2)

# Read data ----
# Sheet layout (col A = labels, cols B-I = 8 cell lines):
#   row 1   = cell line names
#   row 2   = INKA rank
#   rows 3-5 = three viability replicates at 1.1 uM
#   row 7   = pre-computed mean (with a blank row 6 above it)
raw <- read_excel("redo 8C.xlsx", sheet = 1, col_names = FALSE,
                  col_types = "text")

df <- tibble(
  cell_line = as.character(raw[1, -1]),
  inka_rank = as.numeric(raw[2, -1]),
  viability = as.numeric(raw[7, -1])
)

# Spearman rank correlation (INKA rank is ordinal; exact = FALSE because of ties)
ct      <- cor.test(df$inka_rank, df$viability, method = "spearman", exact = FALSE)
rho_lab <- sprintf("rho == %.2f", ct$estimate)
p_lab   <- sprintf("italic(p) == %.4f", ct$p.value)

# Plot ----
# 058B and OE33 have nearly identical viability at rank 1 (~59.7 and ~59.1),
# so a tiny x-jitter keeps both dots visible without distorting the ranks.
x_ann <- max(df$inka_rank) - 0.1
y_lo  <- min(df$viability)

p <- ggplot(df, aes(inka_rank, viability)) +
  geom_smooth(method = "lm", formula = y ~ x,
              color = "#3B7AC8", fill = "grey80", linewidth = 0.7) +
  geom_point(position = position_jitter(width = 0.05, height = 0, seed = 1),
             size = 2) +
  annotate("text", x = x_ann, y = y_lo + 10, label = rho_lab,
           parse = TRUE, hjust = 1, vjust = 0, size = 4.5) +
  annotate("text", x = x_ann, y = y_lo + 2,  label = p_lab,
           parse = TRUE, hjust = 1, vjust = 0, size = 4.5) +
  scale_x_continuous(breaks = 1:4) +
  labs(x = "INKA rank", y = "Viability") +
  theme_classic(base_size = 14) +
  theme(aspect.ratio = 1)

print(p)

ggsave("INKA_vs_viability_spearman.pdf", p, width = 3.2, height = 3.2)
ggsave("INKA_vs_viability_spearman.png", p, width = 3.2, height = 3.2, dpi = 300)