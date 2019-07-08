library(tidyverse)
library(plotly)

# Plot for errors per model
p <- read_csv("../data/processed/x_normalized_errors.csv") %>%
  select(-X1) %>%
  filter(author != "Ishii" & author != "Ec_core" & author != "Chassagnole") %>%
  mutate(author = factor(author)) %>%
  mutate(author = fct_relevel(author, levels=c("Khodayari", "Kurata", "Millard", "iML1515", "ECC2","Exp_ECC2", "Exp_iML1515", "McCloskey"))) %>%
  ggplot(aes(x=author, y=normalized_error, color=author, label=sample_id)) 

errors_p <- p + geom_jitter(width=0.3) + ggtitle("Normalized error knockout simulations")
# It's possible to save it as html using export tool in graphics pane
errors_plotly = ggplotly(errors_p, tooltip = c("x", "label", "y"))

# Plot for ANOVA-like analysis
p <- read_csv("../data/processed/subsystems.csv") %>%
  select(-X1) %>%
  mutate(subsystem = fct_relevel(subsystem, levels = c("Uptake", "EMP", "EDA", "PPP", "Pyr", "TCA"))) %>%
  filter(author != "Ishii" & author != "Ec_core") %>%
  ggplot(aes(x = subsystem, y = relative_error, color = subsystem))

anova_p <- p + geom_jitter(height = 0, width = 0.1)