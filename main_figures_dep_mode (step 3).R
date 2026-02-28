# Set working directory
setwd ("C:/Users/Belle/Documents/Belle - Nanosight/dados_raw")
nanosight_intersect_dep_mode <- read.csv("nanosight_intersect_dep_mode.csv")

######## FIGURE 1 ########

### Number of individuals with each disorder at each time point
## dcany (any disorder - DSM IV)
##Tabela só com w1 e só w2
nanosight_w1_dep_mode <- subset(nanosight_intersect_dep_mode, wave == "t1")
nanosight_w2_dep_mode <- subset(nanosight_intersect_dep_mode, wave == "t2")
nanosight_intersect_pares_dep_mode <- semi_join(nanosight_w1_dep_mode, nanosight_w2_dep_mode, by = ("subjectid"))
#w1
dcany_w1_dep_mode <- subset(nanosight_w1_dep_mode, grepl("T", dcany))
num_dcany_w1_dep_mode <- nrow(dcany_w1_dep_mode)
num_dcany_w1_dep_mode
#w2
dcany_w2_dep_mode <- subset(nanosight_w2_dep_mode, grepl("T", dcany))
num_dcany_w2_dep_mode <- nrow(dcany_w2_dep_mode)
num_dcany_w2_dep_mode
## dcmadep (major depression)
#w1
dcmadep_w1_dep_mode <- subset(nanosight_w1_dep_mode, grepl("T", dcmadep))
num_dcmadep_w1_dep_mode <- nrow(dcmadep_w1_dep_mode)
num_dcmadep_w1_dep_mode
#w2
dcmadep_w2_dep_mode <- subset(nanosight_w2_dep_mode, grepl("T", dcmadep))
num_dcmadep_w2_dep_mode <- nrow(dcmadep_w2_dep_mode)
num_dcmadep_w2_dep_mode
## dcanyanx (any anxiety)
#w1
dcanyanx_w1_dep_mode <- subset(nanosight_w1_dep_mode, grepl("T", dcanyanx))
num_dcanyanx_w1_dep_mode <- nrow(dcanyanx_w1_dep_mode)
num_dcanyanx_w1_dep_mode
#w2
dcanyanx_w2_dep_mode <- subset(nanosight_w2_dep_mode, grepl("T", dcanyanx))
num_dcanyanx_w2_dep_mode <- nrow(dcanyanx_w2_dep_mode)
num_dcanyanx_w2_dep_mode
## dcgena (generalised anxiety)
#w1
dcgena_w1_dep_mode <- subset(nanosight_w1_dep_mode, grepl("T", dcgena))
num_dcgena_w1_dep_mode <- nrow(dcgena_w1_dep_mode)
num_dcgena_w1_dep_mode
#w2
dcgena_w2_dep_mode <- subset(nanosight_w2_dep_mode, grepl("T", dcgena))
num_dcgena_w2_dep_mode <- nrow(dcgena_w2_dep_mode)
num_dcgena_w2_dep_mode
## dcanyhk (ADHD)
#w1
dcanyhk_w1_dep_mode <- subset(nanosight_w1_dep_mode, grepl("T", dcanyhk))
num_dcanyhk_w1_dep_mode <- nrow(dcanyhk_w1_dep_mode)
num_dcanyhk_w1_dep_mode
#w2
dcanyhk_w2_dep_mode <- subset(nanosight_w2_dep_mode, grepl("T", dcanyhk))
num_dcanyhk_w2_dep_mode <- nrow(dcanyhk_w2_dep_mode)
num_dcanyhk_w2_dep_mode
## dcpsych (psychosis)
#w1
dcpsych_w1_dep_mode <- subset(nanosight_w1_dep_mode, grepl("T", dcpsych))
num_dcpsych_w1_dep_mode <- nrow(dcpsych_w1_dep_mode)
num_dcpsych_w1_dep_mode
#w2
dcpsych_w2_dep_mode <- subset(nanosight_w2_dep_mode, grepl("T", dcpsych))
num_dcpsych_w2_dep_mode <- nrow(dcpsych_w2_dep_mode)
num_dcpsych_w2_dep_mode
## dcmania (mania/bip)
#w1
dcmania_w1_dep_mode <- subset(nanosight_w1_dep_mode, grepl("T", dcmania))
num_dcmania_w1_dep_mode <- nrow(dcmania_w1_dep_mode)
num_dcmania_w1_dep_mode
#w2
dcmania_w2_dep_mode <- subset(nanosight_w2_dep_mode, grepl("T", dcmania))
num_dcmania_w2_dep_mode <- nrow(dcmania_w2_dep_mode)
num_dcmania_w2_dep_mode
## dcptsd 
#w1
dcptsd_w1_dep_mode <- subset(nanosight_w1_dep_mode, grepl("T", dcptsd))
num_dcptsd_w1_dep_mode <- nrow(dcptsd_w1_dep_mode)
num_dcptsd_w1_dep_mode
#w2
dcptsd_w2_dep_mode <- subset(nanosight_w2_dep_mode, grepl("T", dcptsd))
num_dcptsd_w2_dep_mode <- nrow(dcptsd_w2_dep_mode)
num_dcptsd_w2_dep_mode

# Packages
library(ggplot2)
library(dplyr)
library(patchwork)

# Creating dataframe 
df_transtornos <- data.frame(
  Disorder = rep(c("Any disorder", "Major depression", "Anxiety disorder", 
                   "Generalized anxiety", "Hyperactivity and Deficit Attention Disorder", "Psychosis", "Mania/Bipolar disorder", "Post Traumatic Stress Disorder"), each = 2),
  Timepoint = rep(c("t1", "t2"), times = 8),
  N = c(num_dcany_w1_dep_mode, num_dcany_w2_dep_mode,
        num_dcmadep_w1_dep_mode, num_dcmadep_w2_dep_mode,
        num_dcanyanx_w1_dep_mode, num_dcanyanx_w2_dep_mode,
        num_dcgena_w1_dep_mode, num_dcgena_w2_dep_mode,
        num_dcanyhk_w1_dep_mode, num_dcanyhk_w2_dep_mode,
        num_dcpsych_w1_dep_mode, num_dcpsych_w2_dep_mode,
        num_dcmania_w1_dep_mode, num_dcmania_w2_dep_mode,
        num_dcptsd_w1_dep_mode, num_dcptsd_w2_dep_mode)
)

# Organiza fatores para ordem correta
df_transtornos$Disorder <- factor(df_transtornos$Disorder,
                                  levels = c("Any disorder", "Major depression", "Anxiety disorder", 
                                             "Generalized anxiety", "Hyperactivity and Deficit Attention Disorder", "Psychosis", "Mania/Bipolar disorder", "Post Traumatic Stress Disorder"))
# Bar Chart
ggplot(df_transtornos, aes(x = Disorder, y = N, fill = Timepoint)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
  geom_text(aes(label = N), position = position_dodge(width = 0.8), vjust = -0.5, size = 4.5) +
  scale_fill_manual(values = c("t1" = "#FFB6C1", "t2" = "#81D8D0"), labels = c("t 1", "t 2")) +
  labs(title = "Number of individuals with each disorder by time point (depression analysis)",
       x = "Mental condition",
       y = "Number of individuals",
       fill = "Timepoint") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 13),
    axis.text.y = element_text(size = 13),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  )
# Save
ggsave("contagem_transtornos_por_timepoint_dep_mode.png", width = 10, height = 8, dpi = 300)


######## FIGURE 2 ########

### Descriptive graphs for size, concentration and percentage
# Graph A - Size
plot_a <- ggplot(nanosight_intersect_dep_mode, aes(x = grupo_analise_dep, y = tamanho_mode_average, fill = wave)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(aes(color = wave), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 2.2, alpha = 0.8) +
  labs(title = "A) EV's Average of Mode's Size (depression)", x = "Trajectory", y = "Size (nm)", fill = "Time point", color = "Time point") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40")
# Graph B - Concentration
plot_b <- ggplot(nanosight_intersect_dep_mode, aes(x = grupo_analise_dep, y = concentracao_real, fill = wave)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(aes(color = wave), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 2.2, alpha = 0.8) +
  labs(title = "B) EV's Concentration", x = "Trajectory", y = "Concentration (particles/mL)", fill = "Time point", color = "Time point") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40")
# Graph C - Small EV %
plot_c <- ggplot(nanosight_intersect_dep_mode, aes(x = grupo_analise_dep, y = EV_pequenas_porcentagem, fill = wave)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(aes(color = wave), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 2.2, alpha = 0.8) +
  labs(title = "C) Small EV's Percentage", x = "Trajectory", y = "Percentage (%)", fill = "Time point", color = "Time point") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40")
# Combine graphs
painel_final_dep_mode <- (plot_a / plot_b / plot_c) + plot_layout(guides = "collect")
# Save
ggsave("painel_EV_boxplots_dep_mode.png", painel_final_dep_mode, width = 11, height = 14, dpi = 300)


######## FIGURE 3 ########

# Packages
library(ggplot2)
library(patchwork)

##Sex
# Graphs by sex
g1 <- ggplot(nanosight_intersect_dep_mode, aes(x = sex, y = tamanho_mode_average, fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(aes(color = sex), width = 0.15, size = 2, alpha = 0.7) +
  labs(title = "A) Size vs Sex", x = "Sex", y = "Size (nm)") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
g2 <- ggplot(nanosight_intersect_dep_mode, aes(x = sex, y = concentracao_real, fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(aes(color = sex), width = 0.15, size = 2, alpha = 0.7) +
  labs(title = "B) Concentration vs Sex", x = "Sex", y = "Concentration (particles/mL)") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
g3 <- ggplot(nanosight_intersect_dep_mode, aes(x = sex, y = EV_pequenas_porcentagem, fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(aes(color = sex), width = 0.15, size = 2, alpha = 0.7) +
  labs(title = "C) Small EVs % vs Sex", x = "Sex", y = "Percentage (%)") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
# Graphs by site
g4 <- ggplot(nanosight_intersect_dep_mode, aes(x = site, y = tamanho_mode_average, fill = site)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(aes(color = site), width = 0.15, size = 2, alpha = 0.7) +
  labs(title = "D) Size vs Site", x = "Site", y = "Size (nm)") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
g5 <- ggplot(nanosight_intersect_dep_mode, aes(x = site, y = concentracao_real, fill = site)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(aes(color = site), width = 0.15, size = 2, alpha = 0.7) +
  labs(title = "E) Concentration vs Site", x = "Site", y = "Concentration (particles/mL)") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
g6 <- ggplot(nanosight_intersect_dep_mode, aes(x = site, y = EV_pequenas_porcentagem, fill = site)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(aes(color = site), width = 0.15, size = 2, alpha = 0.7) +
  labs(title = "F) Small EVs % vs Site", x = "Site", y = "Percentage (%)") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
# Graphs by age
g7 <- ggplot(nanosight_intersect_dep_mode, aes(x = bage, y = tamanho_mode_average)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  labs(title = "G) Size vs Age", x = "Age", y = "Size (nm)") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
g8 <- ggplot(nanosight_intersect_dep_mode, aes(x = bage, y = concentracao_real)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  labs(title = "H) Concentration vs Age", x = "Age", y = "Concentration (particles/mL)") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
g9 <- ggplot(nanosight_intersect_dep_mode, aes(x = bage, y = EV_pequenas_porcentagem)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  labs(title = "I) Small EVs % vs Age", x = "Age", y = "Percentage (%)") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
# Compor o painel
painel_completo_dep_mode <- (g1 | g2 | g3) / (g4 | g5 | g6) / (g7 | g8 | g9)
# Visualizar
painel_completo_dep_mode
# Salvar o painel
ggsave("painel_vesiculas_completo_dep_mode.png", painel_completo_dep_mode, width = 16, height = 14, dpi = 300)


