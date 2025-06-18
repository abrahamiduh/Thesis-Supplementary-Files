install.packages("promises")

library('tidyverse')
library('dtwclust')
library(ggcorrplot)
install.packages("dtwclust")  # If not already installed
library(dtwclust)
library(emmeans)
library(ggpubr)

pheno_data <- read_rds('pheno_data.rds') %>% 
  as_tibble() %>% 
  mutate(timestamp = as_datetime(timestamp)) %>% 
  mutate(variety = as.character(variety)) %>% 
  mutate(treatment = as.character(treatment)) %>% 
  mutate(timestamp_cluster = as.character(timestamp_cluster))

pheno_data %>% 
  group_by(timestamp_cluster) %>% 
  mutate(timestamp = median(timestamp)) %>% 
  ungroup() %>% 
  group_by(variety, treatment, timestamp_cluster) %>% 
  summarise(across(where(is.numeric) | where(is.POSIXct), 
                   ~median(.x, na.rm = TRUE))) %>% 
  ungroup() -> pheno_data_summarized

pheno_data_summarized %>% 
  group_by(timestamp_cluster) %>% 
  count() # remove cluster 2 as incomplete


pheno_data_summarized <- pheno_data_summarized %>% 
  filter(timestamp_cluster != 2)

cluster_to_timestamp <- pheno_data_summarized %>% 
  select(timestamp_cluster, timestamp) %>% 
  unique()

pheno_data_summarized %>%
  ggplot(aes(x=timestamp, y=NDVI.Average)) + 
  geom_point() +
  geom_smooth() +
  facet_grid(variety~treatment, labeller = label_both)









test_var <- 'Digital.Biomass.mm.' #changable


time_series_df <- pheno_data_summarized %>% 
  select(variety, treatment, timestamp_cluster, all_of(test_var)) %>% 
  pivot_wider(names_from = timestamp_cluster, values_from = all_of(test_var))

# can be improved by better algorithm
clusters <- time_series_df %>% select(where(is.numeric)) %>% 
  tsclust(k = 2L,
          distance = "L2", centroid = "pam",
          seed = 3247, trace = TRUE,
          control = partitional_control(nrep = 1L))


printable_df <- time_series_df %>% 
  mutate(shape_cluster = as.character(clusters@cluster)) %>% 
  pivot_longer(where(is.numeric), 
               names_to = 'timestamp_cluster', values_to = test_var) %>% 
  left_join(cluster_to_timestamp)

printable_df %>% 
  ggplot(aes(x=timestamp, y=!!sym(test_var), colour = shape_cluster)) + 
  geom_point() +
  geom_smooth() +
  facet_grid(variety~treatment, labeller = label_both)

cluster_shape_data <- printable_df %>% 
  select(variety, treatment, shape_cluster) %>% 
  unique() %>% 
  mutate(variable = test_var)

corr_table <- cluster_shape_data %>% 
  left_join(select(pheno_data_summarized, 
                   where(is.character) | contains('CM00')) %>% 
              select(-timestamp_cluster) %>% unique)


model.matrix(~0+., data=corr_table %>% select(-variable, -variety)) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag=FALSE, type="lower", lab=TRUE, lab_size=2)


#picked timepoint 5 as one with highest effect on shapes

modeling_table <- printable_df %>% filter(timestamp_cluster == "5") %>% 
  left_join(select(pheno_data_summarized, 
                   where(is.character) | contains('CM00')) %>% 
              select(-timestamp_cluster) %>% unique)
  

modeling_table %>% select(Digital.Biomass.mm., 
                          treatment,
                          CM000847.3_511875_A_G,
                          CM000852.4_48659698_G_A) %>% 
  aov(Digital.Biomass.mm. ~ treatment + CM000847.3_511875_A_G + CM000852.4_48659698_G_A, .) %>% 
  summary()
#########################################################################################

library(rstatix)
library(emmeans)
library(tidyverse)
library(ggcorrplot)
library(dtwclust)
library(emmeans)
library(ggpubr)
set.seed(123)
data <- read.csv("new_data1.csv", 
                 header = TRUE,    # First row as column names
                 stringsAsFactors = FALSE,  # Don't convert strings to factors
                 na.strings = c("NA", "", " "))
data$CM000846.3_39596360_G_A <- factor(data$CM000846.3_39596360_G_A, levels = c(0, 1, 2))
new_df <- data %>% 
  select(treatment, CM000846.3_39596360_G_A, Digital.Biomass.mm.)
new_df
new_df_c <- new_df %>% 
  group_by(treatment) %>%
  emmeans_test(Digital.Biomass.mm. ~ CM000846.3_39596360_G_A, p.adjust.method = "bonferroni") 
new_df_c

res.aov <- new_df %>% anova_test(Digital.Biomass.mm. ~ treatment * CM000846.3_39596360_G_A)

# Visualization: box plots with p-values
new_df_c <- new_df_c %>% add_xy_position(x = "treatment")
bxp +
  stat_pvalue_manual(new_df_c) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(new_df_c)
  )
head(new_df_c)
#################################################################################
set.seed(123)
data("jobsatisfaction", package = "datarium")
jobsatisfaction %>% sample_n_by(gender, education_level, size = 1)


pwc <- jobsatisfaction %>% 
  group_by(gender) %>%
  emmeans_test(score ~ education_level, p.adjust.method = "bonferroni") 
pwc
res.aov <- jobsatisfaction %>% anova_test(score ~ gender * education_level)
pwc <- pwc %>% add_xy_position(x = "gender")
bxp +
  stat_pvalue_manual(pwc) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

##################################################################################
# Alternative using ggpubr's built-in functions
ggboxplot(new_df, x = "treatment", y = "Digital.Biomass.mm.",
          color = "CM000847.3_511875_A_G") +
  stat_pvalue_manual(new_df_c, y.position = max(new_df$Digital.Biomass.mm.) * 1.1) +
  labs(subtitle = get_test_label(res.aov, detailed = TRUE))

bxp <- ggboxplot(
  data = new_df,
  x = "treatment",
  y = "Digital.Biomass.mm.",
  color = "CM000847.3_511875_A_G"  # Add color if you want grouped boxplots
) +
  theme_classic()
final_plot <- bxp +
  stat_pvalue_manual(
    new_df_c,
    label = "p.adj.signif",  # Show significance stars
    tip.length = 0.01         # Adjust bracket tip length if needed
  ) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(new_df_c)
  )

# Display the plot
print(final_plot)




###############################################
###############################################
modeling_table %>%
  select(Digital.Biomass.mm., 
         treatment,
         CM000847.3_511875_A_G,
         CM000852.4_48659698_G_A) %>%
  aov(Digital.Biomass.mm. ~ treatment * CM000847.3_511875_A_G +
        treatment * CM000852.4_48659698_G_A +
        CM000847.3_511875_A_G * CM000852.4_48659698_G_A, .) %>%
  summary()
########################################
#######################################
# Load necessary package
library(ggplot2)

# Run the ANOVA and print the summary
anova_result <- modeling_table %>%
  select(Digital.Biomass.mm., 
         treatment,
         CM000847.3_511875_A_G,
         CM000852.4_48659698_G_A) %>%
  aov(Digital.Biomass.mm. ~ treatment * CM000847.3_511875_A_G +
        treatment * CM000852.4_48659698_G_A +
        CM000847.3_511875_A_G * CM000852.4_48659698_G_A, .)

summary(anova_result)

# Create a new grouping variable to visualize combinations
modeling_table$group <- with(modeling_table, 
                             paste(treatment, CM000847.3_511875_A_G, CM000852.4_48659698_G_A, sep = "_"))

# Plot: Violin + Boxplot of Digital Biomass by combined group
ggplot(modeling_table, aes(x = treatment)
###########################################
###########################################

corr_table %>% 
  group_by(treatment) %>% 
  summarize(cor = cor(CM000847.3_511875_A_G, Digital.Biomass.mm.))
