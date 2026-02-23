#SET WORKING DIRECTORY
# Ajuste o caminho abaixo conforme necessário
project_directory <- "C:/Users/Belle/Documents/Belle - Nanosight/dados_raw"
tryCatch({ setwd(project_directory); cat("Working directory successfully set to:\n", getwd(), "\n\n")}, 
         error = function(e) { stop("ERROR: The specified directory was not found. Please check the path.") })

nanosight_intersect_gena <- read.csv("nanosight_intersect_gena.csv")

if (!file.exists("nanosight_intersect_gena.csv")) {
  stop("Arquivo nanosight_intersect_gena.csv não encontrado.")
}

# ENVIRONMENT AND PACKAGES
cat("--- SECTION 0: SETTING UP THE ENVIRONMENT ---\n")
pacotes_necessarios <- c("lme4", "lmerTest", "glmmTMB", "car", "emmeans", 
                         "performance", "DHARMa", "dplyr") # Adicionei dplyr por conveniência
for (pacote in pacotes_necessarios) {
  if (!require(pacote, character.only = TRUE)) {
    install.packages(pacote, dependencies = TRUE); library(pacote, character.only = TRUE)
  }
}

# Criar novo diretório para os resultados customizados
output_dir <- "analysis_results_ansiedade_controleoriginal"
if (!dir.exists(output_dir)) { dir.create(output_dir) }
cat("Packages loaded and results directory is ready.\n\n")

# LOG
log_file_name <- paste0(output_dir, "/log_analise_custom_", Sys.Date(), ".txt")
sink(log_file_name, append = FALSE, split = TRUE) 

cat("### INÍCIO DO LOG DE ANÁLISE (GRUPOS CUSTOMIZADOS) ###\n")
cat("Data e Hora:", as.character(Sys.time()), "\n")
cat("Diretório de Trabalho:", getwd(), "\n\n")

tryCatch({
  
  # LOADING AND PREPARING DATA
  cat("--- SECTION 1: LOADING AND PREPARING DATA ---\n")
  nanosight_intersect_gena <- read.csv("nanosight_intersect_gena.csv")
  cat("File loaded successfully.\n")
  
  # Definir fatores e níveis (Control como referência)
  nanosight_intersect_gena$grupo_analise_gena <- factor(nanosight_intersect_gena$grupo_analise_gena, 
                                                                  levels = c("Control", "Incident", "Remitted", "Persistent"))
  nanosight_intersect_gena$wave <- as.factor(nanosight_intersect_gena$wave)
  nanosight_intersect_gena$subjectid <- as.factor(nanosight_intersect_gena$subjectid)
  
  # Transformação da porcentagem (mesma lógica anterior)
  cat("Verifying and transforming the percentage variable...\n")
  if (any(nanosight_intersect_gena$EV_pequenas_porcentagem == 0, na.rm = TRUE) || any(nanosight_intersect_gena$EV_pequenas_porcentagem == 100, na.rm = TRUE)) {
    cat("Values of 0 or 100 detected. Applying transformation.\n")
    n <- nrow(nanosight_intersect_gena)
    nanosight_intersect_gena$percentage_prop <- (nanosight_intersect_gena$EV_pequenas_porcentagem * (n - 1) + 0.5) / n
  } else {
    nanosight_intersect_gena$percentage_prop <- nanosight_intersect_gena$EV_pequenas_porcentagem / 100
  }
  cat("Data preparation complete.\n\n")
  
  # CREATE DHARMa diagonosis table 
  criar_tabela_diagnostico_dharma_qq <- function(simulationOutput) {
    alfa <- 0.05
    teste_uniformidade <- testUniformity(simulationOutput, plot = FALSE); teste_dispersao <- testDispersion(simulationOutput, plot = FALSE); teste_outliers <- testOutliers(simulationOutput, plot = FALSE)
    tabela <- data.frame(Diagnostic_Test = c("Overall Uniformity (KS)", "Dispersion", "Outliers"), p_value = c(round(teste_uniformidade$p.value, 3), round(teste_dispersao$p.value, 3), round(teste_outliers$p.value, 3)))
    tabela$Significance <- ifelse(tabela$p_value < alfa, "Significant Violation", "OK")
    return(tabela)
  }
  
  resultados_finais_aic_bic <- data.frame()
  alfa <- 0.05 
  
  # MODEL 1: EV MEAN SIZE
  cat("\n\n=========================================================\n"); cat("  ANALYSIS 1: EV MEAN SIZE (tamanho_mean_average)\n"); cat("=========================================================\n\n")
  cat("--- 1.1. Fitting Gaussian Model with Dispersion Formula ---\n")
  equacao_size <- "glmmTMB(tamanho_mean_average ~ wave * grupo_analise_gena + (1 | subjectid), family = gaussian(), dispformula = grupo_analise_gena)"
  cat("Model Equation: ", equacao_size, "\n")
  modelo_size <- tryCatch({
    glmmTMB(tamanho_mean_average ~ wave * grupo_analise_gena + (1 | subjectid), 
            data = nanosight_intersect_gena, 
            family = gaussian())
  }, error = function(e) { stop("Erro ao ajustar modelo de tamanho: ", conditionMessage(e)) })
  resultados_finais_aic_bic <- rbind(resultados_finais_aic_bic, data.frame(Variable = "Size", Model = "Gaussian (glmmTMB + disp)", AIC = AIC(modelo_size), BIC = BIC(modelo_size)))
  cat("\n--- 1.2. DHARMa Residual Diagnostics ---\n")
  residuos_size <- simulateResiduals(fittedModel = modelo_size)
  tabela_diag_size <- criar_tabela_diagnostico_dharma_qq(residuos_size)
  png(paste0(output_dir, "/plot_diagnostics_size_GAUSSIAN_DISP.png"), width = 3000, height = 1000, res = 300)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 1, 0), cex = 1)
  plotQQunif(residuos_size, main = "A) QQ Plot - Size (Gaussian+Disp)", testUniformity = FALSE, testDispersion = FALSE, testOutliers = FALSE)
  plotResiduals(residuos_size, main = "B) Residuals vs. Predicted - Size")
  dev.off()
  print(tabela_diag_size); write.csv(tabela_diag_size, paste0(output_dir, "/table_diagnostics_dharma_size.csv"), row.names = FALSE)
  cat("\n--- 1.3. ANOVA Table (Type III) for Size ---\n")
  # glmmTMB usa Anova do pacote car (Wald chisquare), não anova F-test do lmerTest
  anova_size <- Anova(modelo_size, type = "III")
  df_anova_size <- as.data.frame(anova_size)
  df_anova_size$Significance <- ifelse(df_anova_size$`Pr(>Chisq)` < alfa, "Significant", "Not Significant")
  print(df_anova_size); write.csv(df_anova_size, paste0(output_dir, "/table_anova_size.csv"))
  cat("\n--- 1.4. Conditional Pairwise Comparisons (Size) ---\n")
  p_interacao_size <- df_anova_size["wave:grupo_analise_gena", "Pr(>Chisq)"]
  if (!is.na(p_interacao_size) && p_interacao_size < 0.05) {
    cat("Interaction is significant. Performing contrasts for BOTH interaction perspectives.\n")
    emm_size_p1 <- emmeans(modelo_size, specs = pairwise ~ grupo_analise_gena | wave, adjust = "tukey")
    write.csv(as.data.frame(emm_size_p1$contrasts), paste0(output_dir, "/table_contrasts_perspective1_size.csv"), row.names = FALSE)
    emm_size_p2 <- emmeans(modelo_size, specs = pairwise ~ wave | grupo_analise_gena, adjust = "tukey")
    write.csv(as.data.frame(emm_size_p2$contrasts), paste0(output_dir, "/table_contrasts_perspective2_size.csv"), row.names = FALSE)
  } else {
    cat("Interaction is not significant. Proceeding to analyze main effects (Type II).\n")
    anova_size_type2 <- Anova(modelo_size, type = "II")
    print(anova_size_type2)
    
    # Check Wave
    if (anova_size_type2["wave", "Pr(>Chisq)"] < 0.05) {
      emm_size_wave <- emmeans(modelo_size, specs = pairwise ~ wave, adjust = "tukey")
      write.csv(as.data.frame(emm_size_wave$contrasts), paste0(output_dir, "/table_contrasts_wave_size.csv"), row.names = FALSE)
    }
    # Check Group
    if (anova_size_type2["grupo_analise_gena", "Pr(>Chisq)"] < 0.05) {
      emm_size_traj <- emmeans(modelo_size, specs = pairwise ~ grupo_analise_gena, adjust = "tukey")
      write.csv(as.data.frame(emm_size_traj$contrasts), paste0(output_dir, "/table_contrasts_trajectory_size.csv"), row.names = FALSE)
    }
  }
  
  cat("\n--- 1.5. VIF (Size) ---\n")
  # VIF model without interaction
  modelo_size_vif <- glmmTMB(tamanho_mean_average ~ wave + grupo_analise_gena + (1 | subjectid), 
                             data = nanosight_intersect_gena, 
                             family = gaussian(), 
                             dispformula = ~ grupo_analise_gena)
  
  vif_df_size <- as.data.frame(check_collinearity(modelo_size_vif))
  print(vif_df_size); write.csv(vif_df_size, paste0(output_dir, "/table_vif_size.csv"), row.names = FALSE)
  
  icc_size <- icc(modelo_size)
  cat("\nICC:\n"); print(icc_size)
  
  # MODEL 2: EV's concentration
  cat("\n\n====================================================================\n"); cat("  ANALYSIS 2: EV CONCENTRATION (concentracao_real)\n"); cat("====================================================================\n\n")
  cat("--- 2.1. Model Simplification Process for Concentration ---\n")
  modelo_conc <- NULL;
  equacao_conc_1 <- "glmmTMB(concentracao_real ~ wave + grupo_analise_gena + (1 | subjectid), data = nanosight_intersect_gena, family = nbinom2(link = log)"
  cat("\nAttempt 1: Full Negative Binomial GLMM...\n", "equacao: ", equacao_conc_1, "\n", "\n")
  modelo_conc_tentativa1 <- tryCatch({glmmTMB(concentracao_real ~ wave + grupo_analise_gena + (1 | subjectid), data = nanosight_intersect_gena, family = nbinom2(link = "log"))}, warning = function(w) {cat("CONVERGENCE WARNING in Attempt 1:\n", conditionMessage(w), "\n"); return(NULL)}, error = function(e) {cat("CONVERGENCE ERROR in Attempt 1:\n", conditionMessage(e), "\n"); return(NULL)})
  if (is.null(modelo_conc_tentativa1)) {equacao_conc_2 <- "glmmTMB(concentracao_real ~ wave + grupo_analise_gena + (1 | subjectid), data = nanosight_intersect_gena, family = nbinom1(link = log)"; cat("\nAttempt 2: Negative Binomial GLMM (main effects only)...\n", "equacao: ", equacao_conc_2, "\n", "\n"); modelo_conc_tentativa2 <- tryCatch({glmmTMB(concentracao_real ~ wave + grupo_analise_gena + (1 | subjectid), data = nanosight_intersect_gena, family = nbinom1(link = "log"))}, warning = function(w) {cat("CONVERGENCE WARNING in Attempt 2:\n", conditionMessage(w), "\n"); return(NULL)}, error = function(e) {cat("CONVERGENCE ERROR in Attempt 2:\n", conditionMessage(e), "\n"); return(NULL)})}
  if(is.null(modelo_conc_tentativa1) && is.null(modelo_conc_tentativa2)) {equacao_conc_3 <- "glmmTMB(concentracao_real ~ wave + grupo_analise_gena, data = nanosight_intersect_gena, family = nbinom1(link = log)"; cat("\nBoth mixed models failed. Selecting final model: GLM with main effects only.\n", "equacao: ", equacao_conc_3,"\n", "\n"); modelo_conc <- glmmTMB(concentracao_real ~ wave + grupo_analise_gena, data = nanosight_intersect_gena, family = nbinom1(link = "log"))}
  resultados_finais_aic_bic <- rbind(resultados_finais_aic_bic, data.frame(Variable = "Concentration", Model = "Negative Binomial (GLM)", AIC = AIC(modelo_conc), BIC = BIC(modelo_conc)))
  cat("\n--- 2.2. ANOVA Table (Type II) for Concentration ---\n")
  anova_conc <- Anova(modelo_conc, type = "II"); df_anova_conc <- as.data.frame(anova_conc); df_anova_conc$Significance <- ifelse(df_anova_conc$`Pr(>Chisq)` < alfa, "Significant", "Not Significant"); print(df_anova_conc); write.csv(df_anova_conc, paste0("analysis_results_ansiedade_controleoriginal/table_anova_concentration_",  ".csv"))
  cat("\n--- 2.3. Pairwise Comparisons for Main Effects (Concentration) ---\n")
  if (df_anova_conc["wave", "Pr(>Chisq)"] < 0.05) {cat("\n Main effect of 'wave' is significant. Performing contrasts.\n\n"); emm_conc_wave <- emmeans(modelo_conc, specs = pairwise ~ wave, adjust = "tukey"); df_emm_conc_wave <- as.data.frame(emm_conc_wave$contrasts); df_emm_conc_wave$Significance <- ifelse(df_emm_conc_wave$p.value < alfa, "Significant", "Not Significant"); print(df_emm_conc_wave); write.csv(df_emm_conc_wave, paste0("analysis_results_ansiedade_controleoriginal/",  "/table_contrasts_wave_concentration_",  ".csv"), row.names = FALSE)}
  if (df_anova_conc["grupo_analise_gena", "Pr(>Chisq)"] < 0.05) {cat("\n Main effect of 'grupo_analise_gena' is significant. Performing contrasts.\n\n"); emm_conc_traj <- emmeans(modelo_conc, specs = pairwise ~ grupo_analise_gena, adjust = "tukey"); df_emm_conc_traj <- as.data.frame(emm_conc_traj$contrasts); df_emm_conc_traj$Significance <- ifelse(df_emm_conc_traj$p.value < alfa, "Significant", "Not Significant"); print(df_emm_conc_traj); write.csv(df_emm_conc_traj, paste0("analysis_results_ansiedade_controleoriginal/",  "/table_contrasts_trajectory_concentration_",  ".csv"), row.names = FALSE)}
  cat("\n--- 2.4. Model Diagnostics (Concentration) ---\n")
  vif_df_conc <- as.data.frame(check_collinearity(modelo_conc)); vif_df_conc$Interpretation <- ifelse(vif_df_conc$VIF >= 10, "High", ifelse(vif_df_conc$VIF >= 5, "Moderate", "Low"));
  cat("VIF:\n"); print(vif_df_conc); write.csv(vif_df_conc, paste0("analysis_results_ansiedade_controleoriginal/",  "table_vif_concentration_",  ".csv"), row.names = FALSE); cat("\nICC:\nNot applicable (GLM without random effects).\n")
  cat("\n--- 2.5. DHARMa Residual Diagnostics (Concentration) ---\n")
  residuos_conc <- simulateResiduals(fittedModel = modelo_conc); tabela_diag_conc <- criar_tabela_diagnostico_dharma_qq(residuos_conc)
  png(paste0("analysis_results_ansiedade_controleoriginal/",  "plot_diagnostics_concentration_",  ".png"), width = 3000, height = 1000, res = 300); par(mfrow = c(1, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 1, 0), cex = 1); plotQQunif(residuos_conc, main = "A) QQ Plot - Concentration", testUniformity = FALSE, testDispersion = FALSE, testOutliers = FALSE); plotResiduals(residuos_conc, main = "B) Residuals vs. Predicted - Concentration"); dev.off()
  print(tabela_diag_conc); write.csv(tabela_diag_conc, paste0("analysis_results_ansiedade_controleoriginal/",  "table_diagnostics_dharma_concentration_",  ".csv"), row.names = FALSE)
  
  # MODEL 3: EVs percentages
  cat("\n\n======================================================================\n"); cat("  ANALYSIS 3: PERCENTAGE OF SMALL EVs (p90_porcentagem)\n"); cat("======================================================================\n\n")
  equacao_perc <- "glmmTMB(percentage_prop ~ wave * grupo_analise_gena + (1 | subjectid), data = nanosight_intersect_gena, family = beta_family(link = logit), dispformula = grupo_analise_gena"
  cat("--- 3.1. Fitting the Beta GLMM for Percentage ---\n", "equacao: ", equacao_perc, "\n")
  modelo_perc <- tryCatch({glmmTMB(percentage_prop ~ wave * grupo_analise_gena + (1 | subjectid), data = nanosight_intersect_gena, family = beta_family(link = "logit"))}) 
  resultados_finais_aic_bic <- rbind(resultados_finais_aic_bic, data.frame(Variable = "Percentage", Model = "Beta (GLMM)", AIC = AIC(modelo_perc), BIC = BIC(modelo_perc)))
  cat("\n--- 3.2. ANOVA Table (Type III) for Percentage ---\n")
  anova_perc_type3 <- Anova(modelo_perc, type = "III"); df_anova_perc3 <- as.data.frame(anova_perc_type3); df_anova_perc3$Significance <- ifelse(df_anova_perc3$`Pr(>Chisq)` < alfa, "Significant", "Not Significant"); print(df_anova_perc3); write.csv(df_anova_perc3, paste0("analysis_results_ansiedade_controleoriginal/",  "table_anova_percentage_",  ".csv"))
  cat("\n--- 3.3. Conditional Pairwise Comparisons (Percentage) ---\n")
  p_interacao_perc <- df_anova_perc3["wave:grupo_analise_gena", "Pr(>Chisq)"]
  if (p_interacao_perc < 0.05) {
    cat("Interaction is significant. Analyzing both perspectives.\n")
    emm_size_p1 <- emmeans(modelo_perc, specs = pairwise ~ grupo_analise_gena | wave, adjust = "tukey"); df_emm_conc_p1 <- as.data.frame(emm_size_p1$contrasts); df_emm_conc_p1$Significance <- ifelse(df_emm_conc_p1$p.value < alfa, "Significant", "Not Significant"); print(df_emm_conc_p1); write.csv(df_emm_conc_p1, paste0("analysis_results_ansiedade_controleoriginal/",  "table_contrasts_wave_perc_",  ".csv"), row.names = FALSE)
    emm_size_p2 <- emmeans(modelo_perc, specs = pairwise ~ wave | grupo_analise_gena, adjust = "tukey"); df_emm_conc_p2 <- as.data.frame(emm_size_p2$contrasts); df_emm_conc_p2$Significance <- ifelse(df_emm_conc_p2$p.value < alfa, "Significant", "Not Significant"); print(df_emm_conc_p2); write.csv(df_emm_conc_p2, paste0("analysis_results_ansiedade_controleoriginal/",  "table_contrasts_grupo_analise_gena_perc_",  ".csv"), row.names = FALSE)
  } else {
    cat("Interaction is not significant. Analyzing main effects with a Type II ANOVA.\n")
    anova_perc_type2 <- Anova(modelo_perc, type = "II"); df_anova_perc2 <- as.data.frame(anova_perc_type2); df_anova_perc2$Significance <- ifelse(df_anova_perc2$`Pr(>Chisq)` < alfa, "Significant", "Not Significant"); print(df_anova_perc2)
    if (df_anova_perc2["wave", "Pr(>Chisq)"] < 0.05) {cat("\n Main effect of 'wave' is significant. Performing contrasts.\n\n"); emm_perc_wave <- emmeans(modelo_perc, specs = pairwise ~ wave, adjust = "tukey"); df_emm_perc_wave <- as.data.frame(emm_perc_wave$contrasts); df_emm_perc_wave$Significance <- ifelse(df_emm_perc_wave$p.value < alfa, "Significant", "Not Significant"); print(df_emm_perc_wave); write.csv(df_emm_perc_wave, paste0("analysis_results_ansiedade_controleoriginal/",  "/table_contrasts_wave_percentage_",  ".csv"), row.names = FALSE) }
    if (df_anova_perc2["grupo_analise_gena", "Pr(>Chisq)"] < 0.05) {cat("\n Main effect of 'grupo_analise_gena' is significant. Performing contrasts.\n\n");emm_perc_traj <- emmeans(modelo_perc, specs = pairwise ~ grupo_analise_gena, adjust = "tukey"); df_emm_perc_traj <- as.data.frame(emm_perc_traj$contrasts); df_emm_perc_traj$Significance <- ifelse(df_emm_perc_traj$p.value < alfa, "Significant", "Not Significant"); print(df_emm_perc_traj); write.csv(df_emm_perc_traj, paste0("analysis_results_ansiedade_controleoriginal/",  "/table_contrasts_trajectory_percentage_",  ".csv"), row.names = FALSE) }
  }
  cat("\n--- 3.4. Model Diagnostics (Percentage) ---\n")
  cat("Assessing VIF on a model without the interaction term to prevent artificial inflation:\n")
  modelo_perc_vif <- glmmTMB(percentage_prop ~ wave + grupo_analise_gena + (1 | subjectid), data = nanosight_intersect_gena, family = beta_family(link = "logit"), dispformula = ~ grupo_analise_gena)
  vif_df_perc <- as.data.frame(check_collinearity(modelo_perc_vif)); vif_df_perc$Interpretation <- ifelse(vif_df_perc$VIF >= 10, "High", ifelse(vif_df_perc$VIF >= 5, "Moderate", "Low")); 
  icc_perc <- icc(modelo_perc)
  cat("VIF:\n"); print(vif_df_perc); write.csv(vif_df_perc, paste0("analysis_results_ansiedade_controleoriginal/",  "table_vif_percentage_",  ".csv"), row.names = FALSE)
  cat("\nICC:\n"); print(icc_perc)
  cat("\n--- 3.5. DHARMa Residual Diagnostics (Percentage) ---\n")
  residuos_perc <- simulateResiduals(fittedModel = modelo_perc); tabela_diag_perc <- criar_tabela_diagnostico_dharma_qq(residuos_perc)
  png(paste0("analysis_results_ansiedade_controleoriginal/",  "plot_diagnostics_dharma_percentage_",  ".png"), width = 3000, height = 1000, res = 300); par(mfrow = c(1, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 1, 0), cex = 1); plotQQunif(residuos_perc, main = paste0("A) QQ Plot - Percentage "), testUniformity = FALSE, testDispersion = FALSE, testOutliers = FALSE); plotResiduals(residuos_perc, main = paste0("B) Residuals vs. Predicted - Percentage")); dev.off()
  print(tabela_diag_perc); write.csv(tabela_diag_perc, paste0("analysis_results_ansiedade_controleoriginal/",  "table_diagnostics_dharma_percentage_",  ".csv"), row.names = FALSE)
  # AIC/BIC
  cat("\n\n=========================================================\n"); cat("  SUMMARY TABLE: AIC/BIC OF FINAL MODELS\n"); cat("=========================================================\n\n")
  print(resultados_finais_aic_bic)
  write.csv(resultados_finais_aic_bic, paste0("analysis_results_ansiedade_controleoriginal/",  "table_summary_aic_bic_final_",  ".csv"), row.names = FALSE)
  
}, finally = {
  cat("\n\n### FIM DO LOG DE ANÁLISE ###\n")
  sink() 
})

###### --- BOXPLOT ---

library(ggplot2)

# Definir ordem dos fatores
nanosight_intersect_gena$grupo_analise_gena <- factor(nanosight_intersect_gena$grupo_analise_gena, levels = c("Control", "Incident", "Remitted"))

# GRÁFICO 1: TAMANHO MÉDIO (Size)
p1 <- ggplot(nanosight_intersect_gena, aes(x = grupo_analise_gena, y = tamanho_mean_average, fill = wave)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) + # Boxplot transparente
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.6) + # Pontos individuais
  theme_minimal() +
  labs(title = "Size Dispersion by Group and Wave",
       subtitle = "Generalized Anxiety Disorder (GAD)",
       x = "Group",
       y = "Mean Size (nm)") +
  scale_fill_brewer(palette = "Set2")

print(p1)
ggsave("boxplot_size_wave_gena_sem_persistent_remitted.png", plot = p1, width = 8, height = 6)

# GRÁFICO 2: CONCENTRAÇÃO
p2 <- ggplot(nanosight_intersect_gena, aes(x = grupo_analise_gena, y = concentracao_real, fill = wave)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.6) +
  theme_minimal() +
  labs(title = "Concentration Dispersion by Group and Wave",
       x = "Group",
       y = "Concentration (Particles/mL)") +
  scale_y_log10() + # Escala Log para ver melhor
  scale_fill_brewer(palette = "Set2")

print(p2)
ggsave("boxplot_conc_wave_gena_sem_persistent_remitted.png", plot = p2, width = 8, height = 6)

# GRÁFICO 3: PORCENTAGEM EVs PEQUENAS
p2 <- ggplot(nanosight_intersect_gena, aes(x = grupo_analise_gena, y = EV_pequenas_porcentagem, fill = wave)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.6) +
  theme_minimal() +
  labs(title = "Percentage Dispersion by Group and Wave",
       x = "Group",
       y = "Small EVs' Percentage (%)") +
  scale_fill_brewer(palette = "Set2")

print(p2)
ggsave("boxplot_perc_wave_sem_persistent_remitted.png", plot = p2, width = 8, height = 6)
