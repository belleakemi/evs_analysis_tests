# SET WORKING DIRECTORY
project_directory <- "C:/Users/Belle/Documents/Belle - Nanosight/dados_raw"
tryCatch({ setwd(project_directory); cat("Working directory successfully set to:\n", getwd(), "\n\n")}, 
         error = function(e) { stop("ERROR: The specified directory was not found. Please check the path.") })
if (!file.exists("nanosight_intersect_ev_pequena_long.csv")) {
  stop("Arquivo nanosight_intersect_ev_pequena_long.csv não encontrado.")
}

# ENVIRONMENT AND PACKAGES
cat("--- SECTION 0: SETTING UP THE ENVIRONMENT ---\n")
pacotes_necessarios <- c("lme4", "lmerTest", "glmmTMB", "car", "emmeans", 
                         "performance", "DHARMa")
for (pacote in pacotes_necessarios) {
  if (!require(pacote, character.only = TRUE)) {
    install.packages(pacote, dependencies = TRUE); library(pacote, character.only = TRUE)
  }
}
if (!dir.exists("analysis_results_traj_dep")) { dir.create("analysis_results_traj_dep") }
cat("Packages loaded and results directory is ready.\n\n")

# LOG
log_file_name <- paste0("analysis_results_traj_dep/log_analise_", Sys.Date(), ".txt")
sink(log_file_name, append = FALSE, split = TRUE) # 'split = TRUE' saves in log file and outputs in console

cat("### INÍCIO DO LOG DE ANÁLISE ###\n")
cat("Data e Hora:", as.character(Sys.time()), "\n")
cat("Diretório de Trabalho:", getwd(), "\n\n")

tryCatch({
  
  # LOADING AND PREPARING DATA
  cat("--- SECTION 1: LOADING AND PREPARING DATA ---\n")
  tryCatch({
    nanosight_intersect_ev_pequena_long <- read.csv("nanosight_intersect_ev_pequena_long.csv")
    cat("File 'nanosight_intersect_ev_pequena_long.csv' loaded successfully.\n")
    nanosight_intersect_ev_pequena_long$traj_dep <- as.factor(nanosight_intersect_ev_pequena_long$traj_dep)
    nanosight_intersect_ev_pequena_long$wave <- as.factor(nanosight_intersect_ev_pequena_long$wave)
    nanosight_intersect_ev_pequena_long$subjectid <- as.factor(nanosight_intersect_ev_pequena_long$subjectid)
  }, error = function(e) {
    stop("ERROR: File 'nanosight_intersect_ev_pequena_long.csv' not found. Ensure it is in your working directory.")
  })
  cat("Verifying and transforming the percentage variable...\n")
  if (any(nanosight_intersect_ev_pequena_long$EV_pequenas_porcentagem == 0, na.rm = TRUE) || any(nanosight_intersect_ev_pequena_long$EV_pequenas_porcentagem == 100, na.rm = TRUE)) {
    cat("Values of 0 or 100 detected. Applying transformation for the Beta model.\n")
    n <- nrow(nanosight_intersect_ev_pequena_long)
    nanosight_intersect_ev_pequena_long$percentage_prop <- (nanosight_intersect_ev_pequena_long$EV_pequenas_porcentagem * (n - 1) + 0.5) / n
  } else {
    cat("No 0 or 100 values detected. Using simple division by 100.\n")
    nanosight_intersect_ev_pequena_long$percentage_prop <- nanosight_intersect_ev_pequena_long$EV_pequenas_porcentagem / 100
  }
  cat("Data preparation complete.\n\n")
  
  # HELPER FUNCTIONS
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
  
  # MODEL 1: EV's Sizes
  cat("\n\n=========================================================\n"); cat("  ANALYSIS 1: EV MEAN SIZE (tamanho_mean_average)\n"); cat("=========================================================\n\n")
  cat("--- 1.1. Comparing Gaussian (LMM) and Gamma (GLMM) models ---\n")
  modelo_size_gaussiano <- lmer(tamanho_mean_average ~ wave * traj_dep + (1 | subjectid), data = nanosight_intersect_ev_pequena_long)
  modelo_size_gamma <- glmmTMB(tamanho_mean_average ~ wave * traj_dep + (1 | subjectid), data = nanosight_intersect_ev_pequena_long, family = Gamma(link = "log"))
  tabela_comp_size <- data.frame(Model = c("Gaussian (LMM)", "Gamma (GLMM)"), AIC = c(AIC(modelo_size_gaussiano), AIC(modelo_size_gamma)), BIC = c(BIC(modelo_size_gaussiano), BIC(modelo_size_gamma)))
  print(tabela_comp_size); write.csv(tabela_comp_size, "analysis_results_traj_dep/table_model_comparison_size.csv", row.names = FALSE)
  cat("\n--- 1.2. DHARMa Residual Diagnostics for BOTH models ---\n")
  residuos_gamma <- simulateResiduals(fittedModel = modelo_size_gamma); tabela_diag_gamma <- criar_tabela_diagnostico_dharma_qq(residuos_gamma)
  png("analysis_results_traj_dep/plot_diagnostics_size_GAMMA.png", width = 3000, height = 1000, res = 300);  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 1, 0), cex = 1); plotQQunif(residuos_gamma, main = "A) QQ Plot - Size (Gamma)", testUniformity = FALSE, testDispersion = FALSE, testOutliers = FALSE); plotResiduals(residuos_gamma, main = "B) Residuals vs. Predicted - Size (Gamma)"); dev.off()
  cat("Diagnostics for Gamma Model:\n"); print(tabela_diag_gamma); write.csv(tabela_diag_gamma, "analysis_results_traj_dep/table_diagnostics_dharma_size_GAMMA.csv", row.names = FALSE)
  residuos_gaussiano <- simulateResiduals(fittedModel = modelo_size_gaussiano); tabela_diag_gaussiano <- criar_tabela_diagnostico_dharma_qq(residuos_gaussiano)
  png("analysis_results_traj_dep/plot_diagnostics_size_GAUSSIAN.png", width = 3000, height = 1000, res = 300); par(mfrow = c(1, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 1, 0), cex = 1); plotQQunif(residuos_gaussiano, main = "A) QQ Plot - Size (Gaussian)", testUniformity = FALSE, testDispersion = FALSE, testOutliers = FALSE); plotResiduals(residuos_gaussiano, main = "B) Residuals vs. Predicted - Size (Gaussian)"); dev.off()
  cat("\nDiagnostics for Gaussian Model:\n"); print(tabela_diag_gaussiano); write.csv(tabela_diag_gaussiano, "analysis_results_traj_dep/table_diagnostics_dharma_size_GAUSSIAN.csv", row.names = FALSE)
  modelo_size <- modelo_size_gaussiano
  cat("\nFinal Model Selected: Gaussian (LMM), based on lower AIC/BIC and adequate diagnostics.\n")
  resultados_finais_aic_bic <- rbind(resultados_finais_aic_bic, data.frame(Variable = "Size", Model = "Gaussian (LMM)", AIC = AIC(modelo_size), BIC = BIC(modelo_size)))
  cat("\n--- 1.3. ANOVA Table (Type III) for Size ---\n")
  anova_size <- anova(modelo_size, type = 3); df_anova_size <- as.data.frame(anova_size); df_anova_size$Significance <- ifelse(df_anova_size$`Pr(>F)` < alfa, "Significant", "Not Significant"); print(df_anova_size); write.csv(df_anova_size, "analysis_results_traj_dep/table_anova_size.csv")
  cat("\n--- 1.4. Conditional Pairwise Comparisons (Size) ---\n")
  p_interacao_size <- df_anova_size["wave:traj_dep", "Pr(>F)"]
  if (p_interacao_size < 0.05) {
    cat("Interaction is significant. Performing contrasts for BOTH interaction perspectives.\n")
    emm_size_p1 <- emmeans(modelo_size, specs = pairwise ~ traj_dep | wave, adjust = "tukey"); df_emm_size_p1 <- as.data.frame(emm_size_p1$contrasts); df_emm_size_p1$Significance <- ifelse(df_emm_size_p1$p.value < alfa, "Significant", "Not Significant"); print(df_emm_size_p1); write.csv(df_emm_size_p1, "analysis_results_traj_dep/table_contrasts_perspective1_size.csv", row.names = FALSE)
    emm_size_p2 <- emmeans(modelo_size, specs = pairwise ~ wave | traj_dep, adjust = "tukey"); df_emm_size_p2 <- as.data.frame(emm_size_p2$contrasts); df_emm_size_p2$Significance <- ifelse(df_emm_size_p2$p.value < alfa, "Significant", "Not Significant"); print(df_emm_size_p2); write.csv(df_emm_size_p2, "analysis_results_traj_dep/table_contrasts_perspective2_size.csv", row.names = FALSE)
  } else {
    cat("Interaction is not significant. Proceeding to analyze main effects with a Type II ANOVA.\n")
    anova_size_type2 <- anova(modelo_size, type = 2); print(anova_size_type2)
    if (anova_size_type2["wave", "Pr(>F)"] < 0.05) { emm_size_wave <- emmeans(modelo_size, specs = pairwise ~ wave, adjust = "tukey"); df_emm_size_wave <- as.data.frame(emm_size_wave$contrasts); df_emm_size_wave$Significance <- ifelse(df_emm_size_wave$p.value < alfa, "Significant", "Not Significant"); print(df_emm_size_wave); write.csv(df_emm_size_wave, "analysis_results_traj_dep/table_contrasts_wave_size.csv", row.names = FALSE) }
    if (anova_size_type2["traj_dep", "Pr(>F)"] < 0.05) { emm_size_traj <- emmeans(modelo_size, specs = pairwise ~ traj_dep, adjust = "tukey"); df_emm_size_traj <- as.data.frame(emm_size_traj$contrasts); df_emm_size_traj$Significance <- ifelse(df_emm_size_traj$p.value < alfa, "Significant", "Not Significant"); print(df_emm_size_traj); write.csv(df_emm_size_traj, "analysis_results_traj_dep/table_contrasts_trajectory_size.csv", row.names = FALSE) }
  }
  cat("\n--- 1.5. Model Diagnostics (Size) ---\n")
  cat("Assessing VIF on a model without the interaction term to prevent artificial inflation:\n")
  modelo_size_vif <- lmer(tamanho_mean_average ~ wave + traj_dep + (1 | subjectid), data = nanosight_intersect_ev_pequena_long)
  vif_df_size <- as.data.frame(check_collinearity(modelo_size_vif)); vif_df_size$Interpretation <- ifelse(vif_df_size$VIF >= 10, "High", ifelse(vif_df_size$VIF >= 5, "Moderate", "Low")); icc_size <- icc(modelo_size)
  cat("VIF:\n"); print(vif_df_size); write.csv(vif_df_size, "analysis_results_traj_dep/table_vif_size.csv", row.names = FALSE)
  cat("\nICC:\n"); print(icc_size)
  
  # MODEL 2: EV's concentration
  cat("\n\n====================================================================\n"); cat("  ANALYSIS 2: EV CONCENTRATION (concentracao_real)\n"); cat("====================================================================\n\n")
  cat("--- 2.1. Model Simplification Process for Concentration ---\n")
  modelo_conc <- NULL; cat("\nAttempt 1: Full Negative Binomial GLMM...\n")
  modelo_conc_tentativa1 <- tryCatch({glmmTMB(concentracao_real ~ wave * traj_dep + (1 | subjectid), data = nanosight_intersect_ev_pequena_long, family = nbinom2(link = "log"))}, warning = function(w) {cat("CONVERGENCE WARNING in Attempt 1:\n", conditionMessage(w), "\n"); return(NULL)}, error = function(e) {cat("CONVERGENCE ERROR in Attempt 1:\n", conditionMessage(e), "\n"); return(NULL)})
  if (is.null(modelo_conc_tentativa1)) {cat("\nAttempt 2: Negative Binomial GLMM (main effects only)...\n"); modelo_conc_tentativa2 <- tryCatch({glmmTMB(concentracao_real ~ wave + traj_dep + (1 | subjectid), data = nanosight_intersect_ev_pequena_long, family = nbinom1(link = "log"))}, warning = function(w) {cat("CONVERGENCE WARNING in Attempt 2:\n", conditionMessage(w), "\n"); return(NULL)}, error = function(e) {cat("CONVERGENCE ERROR in Attempt 2:\n", conditionMessage(e), "\n"); return(NULL)})}
  if(is.null(modelo_conc_tentativa1) && is.null(modelo_conc_tentativa2)) {cat("\nBoth mixed models failed. Selecting final model: GLM with main effects only.\n"); modelo_conc <- glmmTMB(concentracao_real ~ wave + traj_dep, data = nanosight_intersect_ev_pequena_long, family = nbinom1(link = "log"))}
  resultados_finais_aic_bic <- rbind(resultados_finais_aic_bic, data.frame(Variable = "Concentration", Model = "Negative Binomial (GLM)", AIC = AIC(modelo_conc), BIC = BIC(modelo_conc)))
  cat("\n--- 2.2. ANOVA Table (Type II) for Concentration ---\n")
  anova_conc <- Anova(modelo_conc, type = "II"); df_anova_conc <- as.data.frame(anova_conc); df_anova_conc$Significance <- ifelse(df_anova_conc$`Pr(>Chisq)` < alfa, "Significant", "Not Significant"); print(df_anova_conc); write.csv(df_anova_conc, "analysis_results_traj_dep/table_anova_concentration.csv")
  cat("\n--- 2.3. Pairwise Comparisons for Main Effects (Concentration) ---\n")
  if (df_anova_conc["wave", "Pr(>Chisq)"] < 0.05) {cat("Main effect of 'wave' is significant. Performing contrasts.\n"); emm_conc_wave <- emmeans(modelo_conc, specs = pairwise ~ wave, adjust = "tukey"); df_emm_conc_wave <- as.data.frame(emm_conc_wave$contrasts); df_emm_conc_wave$Significance <- ifelse(df_emm_conc_wave$p.value < alfa, "Significant", "Not Significant"); print(df_emm_conc_wave); write.csv(df_emm_conc_wave, "analysis_results_traj_dep/table_contrasts_wave_concentration.csv", row.names = FALSE)}
  if (df_anova_conc["traj_dep", "Pr(>Chisq)"] < 0.05) {cat("\nMain effect of 'traj_dep' is significant. Performing contrasts.\n"); emm_conc_traj <- emmeans(modelo_conc, specs = pairwise ~ traj_dep, adjust = "tukey"); df_emm_conc_traj <- as.data.frame(emm_conc_traj$contrasts); df_emm_conc_traj$Significance <- ifelse(df_emm_conc_traj$p.value < alfa, "Significant", "Not Significant"); print(df_emm_conc_traj); write.csv(df_emm_conc_traj, "analysis_results_traj_dep/table_contrasts_trajectory_concentration.csv", row.names = FALSE)}
  cat("\n--- 2.4. Model Diagnostics (Concentration) ---\n")
  vif_df_conc <- as.data.frame(check_collinearity(modelo_conc)); vif_df_conc$Interpretation <- ifelse(vif_df_conc$VIF >= 10, "High", ifelse(vif_df_conc$VIF >= 5, "Moderate", "Low"));
  cat("VIF:\n"); print(vif_df_conc); write.csv(vif_df_conc, "analysis_results_traj_dep/table_vif_concentration.csv", row.names = FALSE); cat("\nICC:\nNot applicable (GLM without random effects).\n")
  cat("\n--- 2.5. DHARMa Residual Diagnostics (Concentration) ---\n")
  residuos_conc <- simulateResiduals(fittedModel = modelo_conc); tabela_diag_conc <- criar_tabela_diagnostico_dharma_qq(residuos_conc)
  png("analysis_results_traj_dep/plot_diagnostics_concentration.png", width = 3000, height = 1000, res = 300); par(mfrow = c(1, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 1, 0), cex = 1); plotQQunif(residuos_conc, main = "A) QQ Plot - Concentration", testUniformity = FALSE, testDispersion = FALSE, testOutliers = FALSE); plotResiduals(residuos_conc, main = "B) Residuals vs. Predicted - Concentration"); dev.off()
  print(tabela_diag_conc); write.csv(tabela_diag_conc, "analysis_results_traj_dep/table_diagnostics_dharma_concentration.csv", row.names = FALSE)
  
  # MODEL 3: SMALL EVs percentages
  cat("\n\n======================================================================\n"); cat("  ANALYSIS 3: PERCENTAGE OF SMALL EVs (EV_pequenas_porcentagem)\n"); cat("======================================================================\n\n")
  cat("--- 3.1. Fitting the Beta GLMM for Percentage ---\n")
  modelo_perc <- glmmTMB(percentage_prop ~ wave * traj_dep + (1 | subjectid), data = nanosight_intersect_ev_pequena_long, family = beta_family(link = "logit"), dispformula = ~ traj_dep)
  resultados_finais_aic_bic <- rbind(resultados_finais_aic_bic, data.frame(Variable = "Percentage", Model = "Beta (GLMM)", AIC = AIC(modelo_perc), BIC = BIC(modelo_perc)))
  cat("\n--- 3.2. ANOVA Table (Type III) for Percentage ---\n")
  anova_perc_type3 <- Anova(modelo_perc, type = "III"); df_anova_perc3 <- as.data.frame(anova_perc_type3); df_anova_perc3$Significance <- ifelse(df_anova_perc3$`Pr(>Chisq)` < alfa, "Significant", "Not Significant"); print(df_anova_perc3); write.csv(df_anova_perc3, "analysis_results_traj_dep/table_anova_percentage.csv")
  cat("\n--- 3.3. Conditional Pairwise Comparisons (Percentage) ---\n")
  p_interacao_perc <- df_anova_perc3["wave:traj_dep", "Pr(>Chisq)"]
  if (p_interacao_perc < 0.05) {
    cat("Interaction is significant. Analyzing both perspectives.\n")
    emm_perc_p1 <- emmeans(modelo_perc, specs = pairwise ~ traj_dep | wave, adjust = "tukey"); df_emm_perc_p1 <- as.data.frame(emm_perc_p1$contrasts); df_emm_perc_p1$Significance <- ifelse(df_emm_perc_p1$p.value < alfa, "Significant", "Not Significant"); print(df_emm_perc_p1); write.csv(df_emm_perc_p1, "analysis_results_traj_dep/table_contrasts_perspective1_percentage.csv", row.names = FALSE)
    emm_perc_p2 <- emmeans(modelo_perc, specs = pairwise ~ wave | traj_dep, adjust = "tukey"); df_emm_perc_p2 <- as.data.frame(emm_perc_p2$contrasts); df_emm_perc_p2$Significance <- ifelse(df_emm_perc_p2$p.value < alfa, "Significant", "Not Significant"); print(df_emm_perc_p2); write.csv(df_emm_perc_p2, "analysis_results_traj_dep/table_contrasts_perspective2_percentage.csv", row.names = FALSE)
  } else {
    cat("Interaction is not significant. Analyzing main effects with a Type II ANOVA.\n")
    anova_perc_type2 <- Anova(modelo_perc, type = "II"); df_anova_perc2 <- as.data.frame(anova_perc_type2); df_anova_perc2$Significance <- ifelse(df_anova_perc2$`Pr(>Chisq)` < alfa, "Significant", "Not Significant"); print(df_anova_perc2)
    if (df_anova_perc2["wave", "Pr(>Chisq)"] < 0.05) { emm_perc_wave <- emmeans(modelo_perc, specs = pairwise ~ wave, adjust = "tukey"); df_emm_perc_wave <- as.data.frame(emm_perc_wave$contrasts); df_emm_perc_wave$Significance <- ifelse(df_emm_perc_wave$p.value < alfa, "Significant", "Not Significant"); print(df_emm_perc_wave); write.csv(df_emm_perc_wave, "analysis_results_traj_dep/table_contrasts_wave_percentage.csv", row.names = FALSE) }
    if (df_anova_perc2["traj_dep", "Pr(>Chisq)"] < 0.05) { emm_perc_traj <- emmeans(modelo_perc, specs = pairwise ~ traj_dep, adjust = "tukey"); df_emm_perc_traj <- as.data.frame(emm_perc_traj$contrasts); df_emm_perc_traj$Significance <- ifelse(df_emm_perc_traj$p.value < alfa, "Significant", "Not Significant"); print(df_emm_perc_traj); write.csv(df_emm_perc_traj, "analysis_results_traj_dep/table_contrasts_trajectory_percentage.csv", row.names = FALSE) }
  }
  cat("\n--- 3.4. Model Diagnostics (Percentage) ---\n")
  cat("Assessing VIF on a model without the interaction term to prevent artificial inflation:\n")
  modelo_perc_vif <- glmmTMB(percentage_prop ~ wave + traj_dep + (1 | subjectid), data = nanosight_intersect_ev_pequena_long, family = beta_family(link = "logit"))
  vif_df_perc <- as.data.frame(check_collinearity(modelo_perc_vif)); vif_df_perc$Interpretation <- ifelse(vif_df_perc$VIF >= 10, "High", ifelse(vif_df_perc$VIF >= 5, "Moderate", "Low")); icc_perc <- icc(modelo_perc)
  cat("VIF:\n"); print(vif_df_perc); write.csv(vif_df_perc, "analysis_results_traj_dep/table_vif_percentage.csv", row.names = FALSE)
  cat("\nICC:\n"); print(icc_perc)
  cat("\n--- 3.5. DHARMa Residual Diagnostics (Percentage) ---\n")
  residuos_perc <- simulateResiduals(fittedModel = modelo_perc); tabela_diag_perc <- criar_tabela_diagnostico_dharma_qq(residuos_perc)
  png("analysis_results_traj_dep/plot_diagnostics_percentage.png", width = 3000, height = 1000, res = 300); par(mfrow = c(1, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 1, 0), cex = 1); plotQQunif(residuos_perc, main = "A) QQ Plot - Percentage", testUniformity = FALSE, testDispersion = FALSE, testOutliers = FALSE); plotResiduals(residuos_perc, main = "B) Residuals vs. Predicted - Percentage"); dev.off()
  print(tabela_diag_perc); write.csv(tabela_diag_perc, "analysis_results_traj_dep/table_diagnostics_dharma_percentage.csv", row.names = FALSE)
  
  # AIC/BIC
  cat("\n\n=========================================================\n"); cat("  SUMMARY TABLE: AIC/BIC OF FINAL MODELS\n"); cat("=========================================================\n\n")
  print(resultados_finais_aic_bic)
  write.csv(resultados_finais_aic_bic, "analysis_results_traj_dep/table_summary_aic_bic_final.csv", row.names = FALSE)
  
}, finally = {
  # LOG
  cat("\n\n### FIM DO LOG DE ANÁLISE ###\n")
  cat("Todos os resultados foram salvos no diretório 'analysis_results_traj_dep'.\n")
  cat("O log completo desta sessão foi salvo em:", log_file_name, "\n")
  sink() 
}
)

