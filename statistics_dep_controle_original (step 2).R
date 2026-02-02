# SET WORKING DIRECTORY
# Ajuste o caminho abaixo conforme necessário
project_directory <- "C:/Users/Belle/Documents/Belle - Nanosight/dados_raw"
tryCatch({ setwd(project_directory); cat("Working directory successfully set to:\n", getwd(), "\n\n")}, 
         error = function(e) { stop("ERROR: The specified directory was not found. Please check the path.") })

if (!file.exists("nanosight_intersect_ev_pequena_long.csv")) {
  stop("Arquivo nanosight_intersect_ev_pequena_long.csv não encontrado.")
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
output_dir <- "analysis_results_depressao_controleoriginal"
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
  nanosight_intersect_ev_pequena_long <- read.csv("nanosight_intersect_ev_pequena_long.csv")
  cat("File loaded successfully.\n")
  
  # --- CUSTOMIZAÇÃO DOS GRUPOS ---
  cat("Criando a variável de grupo customizada (grupo_analise_dep)...\n")
  
  # 1. Criar coluna vazia
  nanosight_intersect_ev_pequena_long$grupo_analise_dep <- NA
  
  # 2. Atribuir grupos baseados em 'traj_dep' (Assumindo que traj_dcmadep se refere a traj_dep)
  # Normalizando para maiúscula inicial para ficar bonito no gráfico
  nanosight_intersect_ev_pequena_long$grupo_analise_dep[nanosight_intersect_ev_pequena_long$traj_dep == "incidente"] <- "Incidente"
  nanosight_intersect_ev_pequena_long$grupo_analise_dep[nanosight_intersect_ev_pequena_long$traj_dep == "remitente"] <- "Remitente"
  nanosight_intersect_ev_pequena_long$grupo_analise_dep[nanosight_intersect_ev_pequena_long$traj_dep == "persistente"] <- "Persistente"
  
  # 3. Atribuir grupo Control baseado na coluna 'Trajetoria' (SOBRESCREVENDO qualquer valor anterior se for o caso)
  nanosight_intersect_ev_pequena_long$grupo_analise_dep[nanosight_intersect_ev_pequena_long$Trajetoria == "Control"] <- "Control"
  
  # 4. Filtrar NAs (remover quem não se encaixa em nenhum dos 4 grupos)
  n_antes <- nrow(nanosight_intersect_ev_pequena_long)
  nanosight_intersect_ev_pequena_long <- nanosight_intersect_ev_pequena_long[!is.na(nanosight_intersect_ev_pequena_long$grupo_analise_dep), ]
  n_depois <- nrow(nanosight_intersect_ev_pequena_long)
  
  cat(paste("Filtragem realizada. Linhas antes:", n_antes, "| Linhas depois:", n_depois, "\n"))
  cat("Contagem por grupo:\n")
  print(table(nanosight_intersect_ev_pequena_long$grupo_analise_dep))
  
  # 5. Definir fatores e níveis (Control como referência)
  nanosight_intersect_ev_pequena_long$grupo_analise_dep <- factor(nanosight_intersect_ev_pequena_long$grupo_analise_dep, 
                                                              levels = c("Control", "Incidente", "Remitente", "Persistente"))
  nanosight_intersect_ev_pequena_long$wave <- as.factor(nanosight_intersect_ev_pequena_long$wave)
  nanosight_intersect_ev_pequena_long$subjectid <- as.factor(nanosight_intersect_ev_pequena_long$subjectid)
  
  # Transformação da porcentagem (mesma lógica anterior)
  cat("Verifying and transforming the percentage variable...\n")
  if (any(nanosight_intersect_ev_pequena_long$EV_pequenas_porcentagem == 0, na.rm = TRUE) || any(nanosight_intersect_ev_pequena_long$EV_pequenas_porcentagem == 100, na.rm = TRUE)) {
    cat("Values of 0 or 100 detected. Applying transformation.\n")
    n <- nrow(nanosight_intersect_ev_pequena_long)
    nanosight_intersect_ev_pequena_long$percentage_prop <- (nanosight_intersect_ev_pequena_long$EV_pequenas_porcentagem * (n - 1) + 0.5) / n
  } else {
    nanosight_intersect_ev_pequena_long$percentage_prop <- nanosight_intersect_ev_pequena_long$EV_pequenas_porcentagem / 100
  }
  cat("Data preparation complete.\n\n")
  
  # HELPER FUNCTIONS
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
  
cat("\n\n=========================================================\n"); cat("  ANALYSIS 1: EV MEAN SIZE (tamanho_mean_average)\n"); cat("=========================================================\n\n")
  cat("--- 1.1. Comparing Gaussian (LMM) and Gamma (GLMM) models ---\n")
  equacao_gaussiano <- "lmer(tamanho_mean_average ~ wave * grupo_analise_dep + (1 | subjectid), data = nanosight_intersect_ev_pequena_long)"
  cat("Gaussiano: ", equacao_gaussiano,"\n"  )
  modelo_size_gaussiano <- tryCatch({lmer(tamanho_mean_average ~ wave * grupo_analise_dep + (1 | subjectid), data = nanosight_intersect_ev_pequena_long)})
  equacao_gamma <- "glmmTMB(tamanho_mean_average ~ wave * grupo_analise_dep + (1 | subjectid), data = nanosight_intersect_ev_pequena_long, family = Gamma(link = log)"
  cat("Gamma: ", equacao_gamma, "\n")
  modelo_size_gamma <- tryCatch({glmmTMB(tamanho_mean_average ~ wave * grupo_analise_dep + (1 | subjectid), data = nanosight_intersect_ev_pequena_long, family = Gamma(link = "log"))})
  tabela_comp_size <- data.frame(Model = c("Gaussian (LMM)", "Gamma (GLMM)"), AIC = c(AIC(modelo_size_gaussiano), AIC(modelo_size_gamma)), BIC = c(BIC(modelo_size_gaussiano), BIC(modelo_size_gamma)))
  print(tabela_comp_size); write.csv(tabela_comp_size, paste0("analysis_results/",  "table_model_comparison_size_",  ".csv"), row.names = FALSE)
  cat("\n--- 1.2. DHARMa Residual Diagnostics ---\n")
  residuos_gamma <- simulateResiduals(fittedModel = modelo_size_gamma); tabela_diag_gamma <- criar_tabela_diagnostico_dharma_qq(residuos_gamma)
  png(paste0("analysis_results/",  "/plot_diagnostics_size_GAMMA_",  ".png"), width = 3000, height = 1000, res = 300);  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 1, 0), cex = 1); plotQQunif(residuos_gamma, main = "A) QQ Plot - Size (Gamma)", testUniformity = FALSE, testDispersion = FALSE, testOutliers = FALSE); plotResiduals(residuos_gamma, main = "B) Residuals vs. Predicted - Size (Gamma)"); dev.off()
  cat("Diagnostics for Gamma Model:\n"); print(tabela_diag_gamma); write.csv(tabela_diag_gamma, paste0("analysis_results/",  "table_diagnostics_dharma_size_GAMMA_",  ".csv"), row.names = FALSE)
  residuos_gaussiano <- simulateResiduals(fittedModel = modelo_size_gaussiano); tabela_diag_gaussiano <- criar_tabela_diagnostico_dharma_qq(residuos_gaussiano)
  png(paste0("analysis_results/",  "/plot_diagnostics_size_GAUSSIAN_",  ".png"), width = 3000, height = 1000, res = 300); par(mfrow = c(1, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 1, 0), cex = 1); plotQQunif(residuos_gaussiano, main = "A) QQ Plot - Size (Gaussian)", testUniformity = FALSE, testDispersion = FALSE, testOutliers = FALSE); plotResiduals(residuos_gaussiano, main = "B) Residuals vs. Predicted - Size (Gaussian)"); dev.off()
  cat("\nDiagnostics for Gaussian Model:\n"); print(tabela_diag_gaussiano); write.csv(tabela_diag_gaussiano, paste0("analysis_results/",  "table_diagnostics_dharma_size_GAUSSIAN_",  ".csv"), row.names = FALSE)
  cat("\nFinal Model Selected: Gaussian (LMM), based on lower AIC/BIC and adequate diagnostics.\n")
  modelo_size <- modelo_size_gaussiano
  resultados_finais_aic_bic <- rbind(resultados_finais_aic_bic, data.frame(Variable = "Size", Model = "Gaussian (LMM)", AIC = AIC(modelo_size), BIC = BIC(modelo_size)))
  cat("\n--- 1.3. ANOVA Table (Type III) for Size ---\n")
  anova_size <- anova(modelo_size, type = 3); df_anova_size <- as.data.frame(anova_size); df_anova_size$Significance <- ifelse(df_anova_size$`Pr(>F)` < alfa, "Significant", "Not Significant"); print(df_anova_size); write.csv(df_anova_size, paste0("analysis_results/",  "/table_anova_size_",  ".csv"))
  cat("\n--- 1.4. Conditional Pairwise Comparisons (Size) ---\n")
  p_interacao_size <- df_anova_size["wave:grupo_analise_dep", "Pr(>F)"]
  if (p_interacao_size < 0.05) {
    cat("Interaction is significant. Performing contrasts for BOTH interaction perspectives.\n")
    emm_size_p1 <- emmeans(modelo_size, specs = pairwise ~ grupo_analise_dep | wave, adjust = "tukey"); df_emm_size_p1 <- as.data.frame(emm_size_p1$contrasts); df_emm_size_p1$Significance <- ifelse(df_emm_size_p1$p.value < alfa, "Significant", "Not Significant"); print(df_emm_size_p1); write.csv(df_emm_size_p1, paste0("analysis_results/",  "/table_contrasts_wave_size_",  ".csv"), row.names = FALSE)
    emm_size_p2 <- emmeans(modelo_size, specs = pairwise ~ wave | grupo_analise_dep, adjust = "tukey"); df_emm_size_p2 <- as.data.frame(emm_size_p2$contrasts); df_emm_size_p2$Significance <- ifelse(df_emm_size_p2$p.value < alfa, "Significant", "Not Significant"); print(df_emm_size_p2); write.csv(df_emm_size_p2, paste0("analysis_results/",  "/table_contrasts_grupo_analise_dep_size_",  ".csv"), row.names = FALSE)
  } else {
    cat("Interaction is not significant. Proceeding to analyze main effects with a Type II ANOVA.\n")
    anova_size_type2 <- anova(modelo_size, type = 2); print(anova_size_type2)
    if (anova_size_type2["wave", "Pr(>F)"] < 0.05) {cat("\n Main effect of 'wave' is significant. Performing contrasts.\n\n"); emm_size_wave <- emmeans(modelo_size, specs = pairwise ~ wave, adjust = "tukey"); df_emm_size_wave <- as.data.frame(emm_size_wave$contrasts); df_emm_size_wave$Significance <- ifelse(df_emm_size_wave$p.value < alfa, "Significant", "Not Significant"); print(df_emm_size_wave); write.csv(df_emm_size_wave, paste0("analysis_results/",  "/table_contrasts_wave_size_",  ".csv"), row.names = FALSE) }
    if (anova_size_type2["grupo_analise_dep", "Pr(>F)"] < 0.05) {cat("\n Main effect of 'grupo_analise_dep' is significant. Performing contrasts.\n\n"); emm_size_traj <- emmeans(modelo_size, specs = pairwise ~ grupo_analise_dep, adjust = "tukey"); df_emm_size_traj <- as.data.frame(emm_size_traj$contrasts); df_emm_size_traj$Significance <- ifelse(df_emm_size_traj$p.value < alfa, "Significant", "Not Significant"); print(df_emm_size_traj); write.csv(df_emm_size_traj, paste0("analysis_results/",  "/table_contrasts_trajectory_size", ".csv"), row.names = FALSE) }
  }
  cat("\n--- 1.5. Model Diagnostics (Size) ---\n")
  cat("Assessing VIF on a model without the interaction term to prevent artificial inflation:\n")
  modelo_size_vif <- lmer(tamanho_mean_average ~ wave + grupo_analise_dep + (1 | subjectid), data = nanosight_intersect_ev_pequena_long)
  vif_df_size <- as.data.frame(check_collinearity(modelo_size_vif)); vif_df_size$Interpretation <- ifelse(vif_df_size$VIF >= 10, "High", ifelse(vif_df_size$VIF >= 5, "Moderate", "Low")); icc_size <- icc(modelo_size)
  cat("VIF:\n"); print(vif_df_size); write.csv(vif_df_size, paste0("analysis_results/",  "table_vif_size_",  ".csv"), row.names = FALSE)
  cat("\nICC:\n"); print(icc_size)
  
 # MODEL 2: EV's concentration
  cat("\n\n====================================================================\n"); cat("  ANALYSIS 2: EV CONCENTRATION (concentracao_real)\n"); cat("====================================================================\n\n")
  cat("--- 2.1. Model Simplification Process for Concentration ---\n")
  modelo_conc <- NULL;
  equacao_conc_1 <- "glmmTMB(concentracao_real ~ wave * grupo_analise_dep + (1 | subjectid), data = nanosight_intersect_ev_pequena_long, family = nbinom2(link = log)"
  cat("\nAttempt 1: Full Negative Binomial GLMM...\n", "equacao: ", equacao_conc_1, "\n", "\n")
  modelo_conc_tentativa1 <- tryCatch({glmmTMB(concentracao_real ~ wave * grupo_analise_dep + (1 | subjectid), data = nanosight_intersect_ev_pequena_long, family = nbinom2(link = "log"))}, warning = function(w) {cat("CONVERGENCE WARNING in Attempt 1:\n", conditionMessage(w), "\n"); return(NULL)}, error = function(e) {cat("CONVERGENCE ERROR in Attempt 1:\n", conditionMessage(e), "\n"); return(NULL)})
  if (is.null(modelo_conc_tentativa1)) {equacao_conc_2 <- "glmmTMB(concentracao_real ~ wave + grupo_analise_dep + (1 | subjectid), data = nanosight_intersect_ev_pequena_long, family = nbinom1(link = log)"; cat("\nAttempt 2: Negative Binomial GLMM (main effects only)...\n", "equacao: ", equacao_conc_2, "\n", "\n"); modelo_conc_tentativa2 <- tryCatch({glmmTMB(concentracao_real ~ wave + grupo_analise_dep + (1 | subjectid), data = nanosight_intersect_ev_pequena_long, family = nbinom1(link = "log"))}, warning = function(w) {cat("CONVERGENCE WARNING in Attempt 2:\n", conditionMessage(w), "\n"); return(NULL)}, error = function(e) {cat("CONVERGENCE ERROR in Attempt 2:\n", conditionMessage(e), "\n"); return(NULL)})}
  if(is.null(modelo_conc_tentativa1) && is.null(modelo_conc_tentativa2)) {equacao_conc_3 <- "glmmTMB(concentracao_real ~ wave + grupo_analise_dep, data = nanosight_intersect_ev_pequena_long, family = nbinom1(link = log)"; cat("\nBoth mixed models failed. Selecting final model: GLM with main effects only.\n", "equacao: ", equacao_conc_3,"\n", "\n"); modelo_conc <- glmmTMB(concentracao_real ~ wave + grupo_analise_dep, data = nanosight_intersect_ev_pequena_long, family = nbinom1(link = "log"), dispformula = ~ grupo_analise_dep)}
  resultados_finais_aic_bic <- rbind(resultados_finais_aic_bic, data.frame(Variable = "Concentration", Model = "Negative Binomial (GLM)", AIC = AIC(modelo_conc), BIC = BIC(modelo_conc)))
  cat("\n--- 2.2. ANOVA Table (Type II) for Concentration ---\n")
  anova_conc <- Anova(modelo_conc, type = "II"); df_anova_conc <- as.data.frame(anova_conc); df_anova_conc$Significance <- ifelse(df_anova_conc$`Pr(>Chisq)` < alfa, "Significant", "Not Significant"); print(df_anova_conc); write.csv(df_anova_conc, paste0("analysis_results/table_anova_concentration_",  ".csv"))
  cat("\n--- 2.3. Pairwise Comparisons for Main Effects (Concentration) ---\n")
  if (df_anova_conc["wave", "Pr(>Chisq)"] < 0.05) {cat("\n Main effect of 'wave' is significant. Performing contrasts.\n\n"); emm_conc_wave <- emmeans(modelo_conc, specs = pairwise ~ wave, adjust = "tukey"); df_emm_conc_wave <- as.data.frame(emm_conc_wave$contrasts); df_emm_conc_wave$Significance <- ifelse(df_emm_conc_wave$p.value < alfa, "Significant", "Not Significant"); print(df_emm_conc_wave); write.csv(df_emm_conc_wave, paste0("analysis_results/",  "/table_contrasts_wave_concentration_",  ".csv"), row.names = FALSE)}
  if (df_anova_conc["grupo_analise_dep", "Pr(>Chisq)"] < 0.05) {cat("\n Main effect of 'grupo_analise_dep' is significant. Performing contrasts.\n\n"); emm_conc_traj <- emmeans(modelo_conc, specs = pairwise ~ grupo_analise_dep, adjust = "tukey"); df_emm_conc_traj <- as.data.frame(emm_conc_traj$contrasts); df_emm_conc_traj$Significance <- ifelse(df_emm_conc_traj$p.value < alfa, "Significant", "Not Significant"); print(df_emm_conc_traj); write.csv(df_emm_conc_traj, paste0("analysis_results/",  "/table_contrasts_trajectory_concentration_",  ".csv"), row.names = FALSE)}
  cat("\n--- 2.4. Model Diagnostics (Concentration) ---\n")
  vif_df_conc <- as.data.frame(check_collinearity(modelo_conc)); vif_df_conc$Interpretation <- ifelse(vif_df_conc$VIF >= 10, "High", ifelse(vif_df_conc$VIF >= 5, "Moderate", "Low"));
  cat("VIF:\n"); print(vif_df_conc); write.csv(vif_df_conc, paste0("analysis_results/",  "table_vif_concentration_",  ".csv"), row.names = FALSE); cat("\nICC:\nNot applicable (GLM without random effects).\n")
  cat("\n--- 2.5. DHARMa Residual Diagnostics (Concentration) ---\n")
  residuos_conc <- simulateResiduals(fittedModel = modelo_conc); tabela_diag_conc <- criar_tabela_diagnostico_dharma_qq(residuos_conc)
  png(paste0("analysis_results/",  "/plot_diagnostics_concentration_",  ".png"), width = 3000, height = 1000, res = 300); par(mfrow = c(1, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 1, 0), cex = 1); plotQQunif(residuos_conc, main = "A) QQ Plot - Concentration", testUniformity = FALSE, testDispersion = FALSE, testOutliers = FALSE); plotResiduals(residuos_conc, main = "B) Residuals vs. Predicted - Concentration"); dev.off()
  print(tabela_diag_conc); write.csv(tabela_diag_conc, paste0("analysis_results/",  "table_diagnostics_dharma_concentration_",  ".csv"), row.names = FALSE)
  # MODEL 3: EVs percentages
  cat("\n\n======================================================================\n"); cat("  ANALYSIS 3: PERCENTAGE OF SMALL EVs (p90_porcentagem)\n"); cat("======================================================================\n\n")
  equacao_perc <- "glmmTMB(percentage_prop ~ wave * grupo_analise_dep + (1 | subjectid), data = nanosight_intersect_ev_pequena_long, family = beta_family(link = logit)"
  cat("--- 3.1. Fitting the Beta GLMM for Percentage ---\n", "equacao: ", equacao_perc, "\n")
  modelo_perc <- tryCatch({glmmTMB(percentage_prop ~ wave * grupo_analise_dep + (1 | subjectid), data = nanosight_intersect_ev_pequena_long, family = beta_family(link = "logit"))}) 
  resultados_finais_aic_bic <- rbind(resultados_finais_aic_bic, data.frame(Variable = "Percentage", Model = "Beta (GLMM)", AIC = AIC(modelo_perc), BIC = BIC(modelo_perc)))
  cat("\n--- 3.2. ANOVA Table (Type III) for Percentage ---\n")
  anova_perc_type3 <- Anova(modelo_perc, type = "III"); df_anova_perc3 <- as.data.frame(anova_perc_type3); df_anova_perc3$Significance <- ifelse(df_anova_perc3$`Pr(>Chisq)` < alfa, "Significant", "Not Significant"); print(df_anova_perc3); write.csv(df_anova_perc3, paste0("analysis_results/",  "table_anova_percentage_",  ".csv"))
  cat("\n--- 3.3. Conditional Pairwise Comparisons (Percentage) ---\n")
  p_interacao_perc <- df_anova_perc3["wave:grupo_analise_dep", "Pr(>Chisq)"]
   if (p_interacao_perc < 0.05) {
      cat("Interaction is significant. Analyzing both perspectives.\n")
     emm_size_p1 <- emmeans(modelo_perc, specs = pairwise ~ grupo_analise_dep | wave, adjust = "tukey"); df_emm_conc_p1 <- as.data.frame(emm_size_p1$contrasts); df_emm_conc_p1$Significance <- ifelse(df_emm_conc_p1$p.value < alfa, "Significant", "Not Significant"); print(df_emm_conc_p1); write.csv(df_emm_conc_p1, paste0("analysis_results/",  "table_contrasts_wave_perc_",  ".csv"), row.names = FALSE)
     emm_size_p2 <- emmeans(modelo_perc, specs = pairwise ~ wave | grupo_analise_dep, adjust = "tukey"); df_emm_conc_p2 <- as.data.frame(emm_size_p2$contrasts); df_emm_conc_p2$Significance <- ifelse(df_emm_conc_p2$p.value < alfa, "Significant", "Not Significant"); print(df_emm_conc_p2); write.csv(df_emm_conc_p2, paste0("analysis_results/",  "table_contrasts_grupo_analise_dep_perc_",  ".csv"), row.names = FALSE)
    } else {
      cat("Interaction is not significant. Analyzing main effects with a Type II ANOVA.\n")
      anova_perc_type2 <- Anova(modelo_perc, type = "II"); df_anova_perc2 <- as.data.frame(anova_perc_type2); df_anova_perc2$Significance <- ifelse(df_anova_perc2$`Pr(>Chisq)` < alfa, "Significant", "Not Significant"); print(df_anova_perc2)
      if (df_anova_perc2["wave", "Pr(>Chisq)"] < 0.05) {cat("\n Main effect of 'wave' is significant. Performing contrasts.\n\n"); emm_perc_wave <- emmeans(modelo_perc, specs = pairwise ~ wave, adjust = "tukey"); df_emm_perc_wave <- as.data.frame(emm_perc_wave$contrasts); df_emm_perc_wave$Significance <- ifelse(df_emm_perc_wave$p.value < alfa, "Significant", "Not Significant"); print(df_emm_perc_wave); write.csv(df_emm_perc_wave, paste0("analysis_results/",  "/table_contrasts_wave_percentage_",  ".csv"), row.names = FALSE) }
      if (df_anova_perc2["grupo_analise_dep", "Pr(>Chisq)"] < 0.05) {cat("\n Main effect of 'grupo_analise_dep' is significant. Performing contrasts.\n\n");emm_perc_traj <- emmeans(modelo_perc, specs = pairwise ~ grupo_analise_dep, adjust = "tukey"); df_emm_perc_traj <- as.data.frame(emm_perc_traj$contrasts); df_emm_perc_traj$Significance <- ifelse(df_emm_perc_traj$p.value < alfa, "Significant", "Not Significant"); print(df_emm_perc_traj); write.csv(df_emm_perc_traj, paste0("analysis_results/",  "/table_contrasts_trajectory_percentage_",  ".csv"), row.names = FALSE) }
   }
    cat("\n--- 3.4. Model Diagnostics (Percentage) ---\n")
    cat("Assessing VIF on a model without the interaction term to prevent artificial inflation:\n")
    modelo_perc_vif <- glmmTMB(percentage_prop ~ wave + grupo_analise_dep + (1 | subjectid), data = nanosight_intersect_ev_pequena_long, family = beta_family(link = "logit"), dispformula = ~ grupo_analise_dep)
    vif_df_perc <- as.data.frame(check_collinearity(modelo_perc_vif)); vif_df_perc$Interpretation <- ifelse(vif_df_perc$VIF >= 10, "High", ifelse(vif_df_perc$VIF >= 5, "Moderate", "Low")); 
    icc_perc <- icc(modelo_perc)
    cat("VIF:\n"); print(vif_df_perc); write.csv(vif_df_perc, paste0("analysis_results/",  "table_vif_percentage_",  ".csv"), row.names = FALSE)
    cat("\nICC:\n"); print(icc_perc)
    cat("\n--- 3.5. DHARMa Residual Diagnostics (Percentage) ---\n")
    residuos_perc <- simulateResiduals(fittedModel = modelo_perc); tabela_diag_perc <- criar_tabela_diagnostico_dharma_qq(residuos_perc)
    png(paste0("analysis_results/",  "/plot_diagnostics_dharma_percentage_",  ".png"), width = 3000, height = 1000, res = 300); par(mfrow = c(1, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 1, 0), cex = 1); plotQQunif(residuos_perc, main = paste0("A) QQ Plot - Percentage ", tag), testUniformity = FALSE, testDispersion = FALSE, testOutliers = FALSE); plotResiduals(residuos_perc, main = paste0("B) Residuals vs. Predicted - Percentage", tag)); dev.off()
    print(tabela_diag_perc); write.csv(tabela_diag_perc, paste0("analysis_results/",  "table_diagnostics_dharma_percentage_",  ".csv"), row.names = FALSE)
  # AIC/BIC
  cat("\n\n=========================================================\n"); cat("  SUMMARY TABLE: AIC/BIC OF FINAL MODELS\n"); cat("=========================================================\n\n")
  print(resultados_finais_aic_bic)
  write.csv(resultados_finais_aic_bic, paste0("analysis_results/",  "table_summary_aic_bic_final_",  ".csv"), row.names = FALSE)
  
}, finally = {
  cat("\n\n### FIM DO LOG DE ANÁLISE ###\n")
  sink() 
})


