# install.packages('tidyverse')
# install.packages('ggplot2')
# install,packages('readr')

library(tidyverse)
library(ggplot2)
library(readr)

# loading in files
score_calls <- as.data.frame(read_tsv('./results/scoring/scores_callsweep.tsv'))

# first -> making confusion matrices
cm_maker <- function(score_calls, threshold, tool_name, cmu_value) {
  data <- score_calls %>% filter(tool == tool_name & truth_min_umis == 1 & call_min_umis == cmu_value) # accessing row
  # extracting the TP, FP, FN, and TN from score calls
  tp <- data$TP  # True Positive
  fp <- data$FP  # False Positive
  fn <- data$FN  # False Negative
  tn <- data$TN  # True Negative

  cm_df <- data.frame(
    Actual = factor(c("Positive", "Negative", "Positive", "Negative"), 
                    levels = c("Negative", "Positive")),
    Predicted = factor(c("Positive", "Positive", "Negative", "Negative"), 
                       levels = c("Positive", "Negative")),
    Count = c(tp, fp, fn, tn),
    Label = c("TP", "FP", "FN", "TN")
  )
  print(cm_df)
  plot <- ggplot(cm_df, aes(x = Predicted, y = Actual, fill = Count)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = paste0(Label, "\n", Count)), size = 6, fontface = "bold") +
    scale_fill_gradient(low = "#e8f1fa", high = "#1f78b4") +
    coord_equal() +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank()
      #legend.position = "none"
    ) +
    labs(
      title = paste0("Confusion Matrix: ",tool_name, '-', data$label),
      x = "Predicted Viral Read Presence per Cell",
      y = "Actual Viral Read Presence per Cell"
    ) + 
  theme(
    axis.title = element_text(size = 16),   
    plot.title = element_text(size = 18),   
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),  
    axis.title.y = element_text(margin = margin(r = 15)), # bringing the axes slightly farther away
    axis.title.x = element_text(margin = margin(t = 15))
  )
  return(plot)
}

# making cm for each tool
spq_cm <- cm_maker(score_calls, threshold = 1, tool = 'ScPathoQuant', cmu_value = 1)
dvsc_cm <- cm_maker(score_calls, threshold = 1, tool = 'DVsc', cmu_value = 1)
vt_cm <- cm_maker(score_calls, threshold = 1, tool = 'ViralTrack', cmu_value=1)
pt_cm <- cm_maker(score_calls, threshold = 1, tool = 'PathogenTrack',cmu_value = 1)
