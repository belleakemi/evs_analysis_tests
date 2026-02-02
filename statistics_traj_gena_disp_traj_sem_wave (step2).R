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
if (!dir.exists("analysis_results_traj_gena_disp_traj_sem_wave")) { dir.create("analysis_results_traj_gena_disp_traj_sem_wave") }
cat("Packages loaded and results directory is ready.\n\n")

# LOG
log_file_name <- paste0("analysis_results_traj_gena_disp_traj_sem_wave/log_analise_", Sys.Date(), ".txt")
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
    nanosight_intersect_ev_pequena_long$traj_gena <- as.factor(nanosight_intersect_ev_pequena_long$traj_gena)
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
  
  # MODEL 2: EV's concentration
  cat("\n\n====================================================================\n"); cat("  ANALYSIS 2: EV CONCENTRATION (concentracao_real)\n"); cat("====================================================================\n\n")
  cat("--- 2.1. Model Simplification Process for Concentration ---\n")
  modelo_conc <- NULL; cat("\nAttempt 1: Full Negative Binomial GLMM...\n")
  modelo_conc_tentativa1 <- tryCatch({glmmTMB(concentracao_real ~ traj_gena + (1 | subjectid), data = nanosight_intersect_ev_pequena_long, family = nbinom2(link = "log"), dispformula = "traj_gena")}, warning = function(w) {cat("CONVERGENCE WARNING in Attempt 1:\n", conditionMessage(w), "\n"); return(NULL)}, error = function(e) {cat("CONVERGENCE ERROR in Attempt 1:\n", conditionMessage(e), "\n"); return(NULL)})
  if (is.null(modelo_conc_tentativa1)) {cat("\nAttempt 2: Negative Binomial GLMM (main effects only)...\n"); modelo_conc_tentativa2 <- tryCatch({glmmTMB(concentracao_real ~ traj_gena + (1 | subjectid), data = nanosight_intersect_ev_pequena_long, family = nbinom1(link = "log"))}, warning = function(w) {cat("CONVERGENCE WARNING in Attempt 2:\n", conditionMessage(w), "\n"); return(NULL)}, error = function(e) {cat("CONVERGENCE ERROR in Attempt 2:\n", conditionMessage(e), "\n"); return(NULL)})}
  if(is.null(modelo_conc_tentativa1) && is.null(modelo_conc_tentativa2)) {cat("\nBoth mixed models failed. Selecting final model: GLM with main effects only.\n"); modelo_conc <- glmmTMB(concentracao_real ~ traj_gena, data = nanosight_intersect_ev_pequena_long, family = nbinom1(link = "log"))}
  resultados_finais_aic_bic <- rbind(resultados_finais_aic_bic, data.frame(Variable = "Concentration", Model = "Negative Binomial (GLM)", AIC = AIC(modelo_conc), BIC = BIC(modelo_conc)))
  cat("\n--- 2.2. ANOVA Table (Type II) for Concentration ---\n")
  anova_conc <- Anova(modelo_conc, type = "II"); df_anova_conc <- as.data.frame(anova_conc); df_anova_conc$Significance <- ifelse(df_anova_conc$`Pr(>Chisq)` < alfa, "Significant", "Not Significant"); print(df_anova_conc); write.csv(df_anova_conc, "analysis_results_traj_gena_disp_traj_sem_wave/table_anova_concentration.csv")
  cat("\n--- 2.3. Pairwise Comparisons for Main Effects (Concentration) ---\n")
  if (df_anova_conc["traj_gena", "Pr(>Chisq)"] < 0.05) {cat("\nMain effect of 'traj_gena' is significant. Performing contrasts.\n"); emm_conc_traj <- emmeans(modelo_conc, specs = pairwise ~ traj_gena, adjust = "tukey"); df_emm_conc_traj <- as.data.frame(emm_conc_traj$contrasts); df_emm_conc_traj$Significance <- ifelse(df_emm_conc_traj$p.value < alfa, "Significant", "Not Significant"); print(df_emm_conc_traj); write.csv(df_emm_conc_traj, "analysis_results_traj_gena_disp_traj_sem_wave/table_contrasts_trajectory_concentration.csv", row.names = FALSE)}
  cat("\n--- 2.4. Model Diagnostics (Concentration) ---\n")
  vif_df_conc <- as.data.frame(check_collinearity(modelo_conc)); vif_df_conc$Interpretation <- ifelse(vif_df_conc$VIF >= 10, "High", ifelse(vif_df_conc$VIF >= 5, "Moderate", "Low"));
  cat("VIF:\n"); print(vif_df_conc); write.csv(vif_df_conc, "analysis_results_traj_gena_disp_traj_sem_wave/table_vif_concentration.csv", row.names = FALSE); cat("\nICC:\nNot applicable (GLM without random effects).\n")
  cat("\n--- 2.5. DHARMa Residual Diagnostics (Concentration) ---\n")
  residuos_conc <- simulateResiduals(fittedModel = modelo_conc); tabela_diag_conc <- criar_tabela_diagnostico_dharma_qq(residuos_conc)
  png("analysis_results_traj_gena_disp_traj_sem_wave/plot_diagnostics_concentration.png", width = 3000, height = 1000, res = 300); par(mfrow = c(1, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 1, 0), cex = 1); plotQQunif(residuos_conc, main = "A) QQ Plot - Concentration", testUniformity = FALSE, testDispersion = FALSE, testOutliers = FALSE); plotResiduals(residuos_conc, main = "B) Residuals vs. Predicted - Concentration"); dev.off()
  print(tabela_diag_conc); write.csv(tabela_diag_conc, "analysis_results_traj_gena_disp_traj_sem_wave/table_diagnostics_dharma_concentration.csv", row.names = FALSE)
  
  # AIC/BIC
  cat("\n\n=========================================================\n"); cat("  SUMMARY TABLE: AIC/BIC OF FINAL MODELS\n"); cat("=========================================================\n\n")
  print(resultados_finais_aic_bic)
  write.csv(resultados_finais_aic_bic, "analysis_results_traj_gena_disp_traj_sem_wave/table_summary_aic_bic_final.csv", row.names = FALSE)
  
}, finally = {
  # LOG
  cat("\n\n### FIM DO LOG DE ANÁLISE ###\n")
  cat("Todos os resultados foram salvos no diretório 'analysis_results_traj_gena_disp_traj_sem_wave'.\n")
  cat("O log completo desta sessão foi salvo em:", log_file_name, "\n")
  sink() 
}
)

