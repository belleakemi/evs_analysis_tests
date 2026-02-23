# SET WORKING DIRECTORY
project_directory <- "C:/Users/Belle/Documents/Belle - Nanosight/dados_raw"
tryCatch({ setwd(project_directory); cat("Working directory successfully set to:\n", getwd(), "\n\n")},
         error = function(e) { stop("ERROR: The specified directory was not found. Please check the path.") })
options(warn = 1) # Força avisos a serem impressos no exato momento em que ocorre

# ENVIRONMENT AND PACKAGES
cat("--- SECTION 0: SETTING UP THE ENVIRONMENT ---\n")
pacotes_necessarios <- c("lme4", "lmerTest", "glmmTMB", "car", "emmeans",
                         "performance", "DHARMa")
for (pacote in pacotes_necessarios) {
  if (!require(pacote, character.only = TRUE)) {
    install.packages(pacote, dependencies = TRUE); library(pacote, character.only = TRUE)
  }
}
if (!dir.exists("analysis_results_traj_anyhk")) { dir.create("analysis_results_traj_anyhk") }
cat("Packages loaded and results directory is ready.\n\n")

config_analises <- list(
  list(arquivo = "nanosight_intersect_ev_pequena_long.csv", coluna_perc = "EV_pequenas_porcentagem", tag = "EV_PEQUENA")
  # list(arquivo = "nanosight_intersect_p90.csv",        coluna_perc = "p90_porcentagem",        tag = "P90"),
  # list(arquivo = "nanosight_intersect_p95.csv",        coluna_perc = "p95_porcentagem",        tag = "P95")
)

for (config in config_analises) {
  
  tag <- config$tag
  arquivo_csv <- config$arquivo
  coluna_alvo <- config$coluna_perc
  if (!dir.exists(paste0("analysis_results_traj_anyhk/", tag))) { dir.create(paste0("analysis_results_traj_anyhk/", tag)) }
  
  # LOG
  log_file_name <- paste0("analysis_results_traj_anyhk/", tag, "/log_analise_", tag, "_", Sys.Date(), ".txt")
  ## ABRINDO CONEXÃO DE LOG
  con <- file(log_file_name, open = "wt")
  sink(con, split = TRUE)            # Saída padrão (cat, print) vai para Console + Log
  sink(con, type = "message")        # Mensagens de Erro/Aviso vão APENAS para o Log (limpa o console)
  
  cat("### INÍCIO DO LOG DE ANÁLISE - GRUPO:", tag, " ###\n")
  cat("Data e Hora:", as.character(Sys.time()), "\n")
  cat("Diretório de Trabalho:", getwd(), "\n\n")
  cat("Arquivo de Origem:", arquivo_csv, "\n")
  cat("Coluna de Porcentagem Analisada:", coluna_alvo, "\n\n")
  
  tryCatch({
    
    # LOADING AND PREPARING DATA
    cat("--- SECTION 1: LOADING AND PREPARING DATA ---\n")
    tryCatch({
      nanosight_intersect <- read.csv(arquivo_csv)
      cat("File ", arquivo_csv, " loaded successfully.\n")
      nanosight_intersect$traj_anyhk <- as.factor(nanosight_intersect$traj_anyhk)
      nanosight_intersect$wave <- as.factor(nanosight_intersect$wave)
      nanosight_intersect$subjectid <- as.factor(nanosight_intersect$subjectid)
    }, error = function(e) {
      stop("ERROR: Initial input file not found. Ensure it is in your working directory.")
    })
    cat("Verifying and transforming the percentage variable...\n")
    nanosight_intersect$coluna_porcentagem_em_uso <- nanosight_intersect[[coluna_alvo]]
    proporcao_pura <- nanosight_intersect$coluna_porcentagem_em_uso / 100
    if (any(nanosight_intersect$coluna_porcentagem_em_uso == 0, na.rm = TRUE) || any(nanosight_intersect$coluna_porcentagem_em_uso == 100, na.rm = TRUE)) {
      cat("Values of 0 or 100 detected. Applying transformation for the Beta model: (valor/100)*((n-1) + 0.5)/n.\n")
      n <- nrow(nanosight_intersect)
      nanosight_intersect$percentage_prop <- (proporcao_pura * (n - 1) + 0.5) / n
    } else {
      cat("No 0 or 100 values detected. Using simple division by 100.\n")
      nanosight_intersect$percentage_prop <- proporcao_pura
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
    equacao_gaussiano <- "lmer(tamanho_mean_average ~ wave * traj_anyhk + (1 | subjectid), data = nanosight_intersect)"
    cat("Gaussiano: ", equacao_gaussiano,"\n"  )
    modelo_size_gaussiano <- tryCatch({lmer(tamanho_mean_average ~ wave * traj_anyhk + (1 | subjectid), data = nanosight_intersect)})
    equacao_gamma <- "glmmTMB(tamanho_mean_average ~ wave * traj_anyhk + (1 | subjectid), data = nanosight_intersect, family = Gamma(link = log)"
    cat("Gamma: ", equacao_gamma, "\n")
    modelo_size_gamma <- tryCatch({glmmTMB(tamanho_mean_average ~ wave * traj_anyhk + (1 | subjectid), data = nanosight_intersect, family = Gamma(link = "log"))})
    tabela_comp_size <- data.frame(Model = c("Gaussian (LMM)", "Gamma (GLMM)"), AIC = c(AIC(modelo_size_gaussiano), AIC(modelo_size_gamma)), BIC = c(BIC(modelo_size_gaussiano), BIC(modelo_size_gamma)))
    print(tabela_comp_size); write.csv(tabela_comp_size, paste0("analysis_results_traj_anyhk/", tag, "/table_model_comparison_size_", tag, ".csv"), row.names = FALSE)
    cat("\n--- 1.2. DHARMa Residual Diagnostics ---\n")
    residuos_gamma <- simulateResiduals(fittedModel = modelo_size_gamma); tabela_diag_gamma <- criar_tabela_diagnostico_dharma_qq(residuos_gamma)
    png(paste0("analysis_results_traj_anyhk/", tag, "/plot_diagnostics_size_GAMMA_", tag, ".png"), width = 3000, height = 1000, res = 300);  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 1, 0), cex = 1); plotQQunif(residuos_gamma, main = "A) QQ Plot - Size (Gamma)", testUniformity = FALSE, testDispersion = FALSE, testOutliers = FALSE); plotResiduals(residuos_gamma, main = "B) Residuals vs. Predicted - Size (Gamma)"); dev.off()
    cat("Diagnostics for Gamma Model:\n"); print(tabela_diag_gamma); write.csv(tabela_diag_gamma, paste0("analysis_results_traj_anyhk/", tag, "/table_diagnostics_dharma_size_GAMMA_", tag, ".csv"), row.names = FALSE)
    residuos_gaussiano <- simulateResiduals(fittedModel = modelo_size_gaussiano); tabela_diag_gaussiano <- criar_tabela_diagnostico_dharma_qq(residuos_gaussiano)
    png(paste0("analysis_results_traj_anyhk/", tag, "/plot_diagnostics_size_GAUSSIAN_", tag, ".png"), width = 3000, height = 1000, res = 300); par(mfrow = c(1, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 1, 0), cex = 1); plotQQunif(residuos_gaussiano, main = "A) QQ Plot - Size (Gaussian)", testUniformity = FALSE, testDispersion = FALSE, testOutliers = FALSE); plotResiduals(residuos_gaussiano, main = "B) Residuals vs. Predicted - Size (Gaussian)"); dev.off()
    cat("\nDiagnostics for Gaussian Model:\n"); print(tabela_diag_gaussiano); write.csv(tabela_diag_gaussiano, paste0("analysis_results_traj_anyhk/", tag, "/table_diagnostics_dharma_size_GAUSSIAN_", tag, ".csv"), row.names = FALSE)
    cat("\nFinal Model Selected: Gaussian (LMM), based on lower AIC/BIC and adequate diagnostics.\n")
    modelo_size <- modelo_size_gaussiano
    resultados_finais_aic_bic <- rbind(resultados_finais_aic_bic, data.frame(Variable = "Size", Model = "Gaussian (LMM)", AIC = AIC(modelo_size), BIC = BIC(modelo_size)))
    cat("\n--- 1.3. ANOVA Table (Type III) for Size ---\n")
    anova_size <- anova(modelo_size, type = 3); df_anova_size <- as.data.frame(anova_size); df_anova_size$Significance <- ifelse(df_anova_size$`Pr(>F)` < alfa, "Significant", "Not Significant"); print(df_anova_size); write.csv(df_anova_size, paste0("analysis_results_traj_anyhk/", tag, "/table_anova_size_", tag, ".csv"))
    cat("\n--- 1.4. Conditional Pairwise Comparisons (Size) ---\n")
    p_interacao_size <- df_anova_size["wave:traj_anyhk", "Pr(>F)"]
    if (p_interacao_size < 0.05) {
      cat("Interaction is significant. Performing contrasts for BOTH interaction perspectives.\n")
      emm_size_p1 <- emmeans(modelo_size, specs = pairwise ~ traj_anyhk | wave, adjust = "tukey"); df_emm_size_p1 <- as.data.frame(emm_size_p1$contrasts); df_emm_size_p1$Significance <- ifelse(df_emm_size_p1$p.value < alfa, "Significant", "Not Significant"); print(df_emm_size_p1); write.csv(df_emm_size_p1, paste0("analysis_results_traj_anyhk/", tag, "/table_contrasts_wave_size_", tag, ".csv"), row.names = FALSE)
      emm_size_p2 <- emmeans(modelo_size, specs = pairwise ~ wave | traj_anyhk, adjust = "tukey"); df_emm_size_p2 <- as.data.frame(emm_size_p2$contrasts); df_emm_size_p2$Significance <- ifelse(df_emm_size_p2$p.value < alfa, "Significant", "Not Significant"); print(df_emm_size_p2); write.csv(df_emm_size_p2, paste0("analysis_results_traj_anyhk/", tag, "/table_contrasts_traj_anyhk_size_", tag, ".csv"), row.names = FALSE)
    } else {
      cat("Interaction is not significant. Proceeding to analyze main effects with a Type II ANOVA.\n")
      anova_size_type2 <- anova(modelo_size, type = 2); print(anova_size_type2)
      if (anova_size_type2["wave", "Pr(>F)"] < 0.05) {cat("\n Main effect of 'wave' is significant. Performing contrasts.\n\n"); emm_size_wave <- emmeans(modelo_size, specs = pairwise ~ wave, adjust = "tukey"); df_emm_size_wave <- as.data.frame(emm_size_wave$contrasts); df_emm_size_wave$Significance <- ifelse(df_emm_size_wave$p.value < alfa, "Significant", "Not Significant"); print(df_emm_size_wave); write.csv(df_emm_size_wave, paste0("analysis_results_traj_anyhk/", tag, "/table_contrasts_wave_size_", tag, ".csv"), row.names = FALSE) }
      if (anova_size_type2["traj_anyhk", "Pr(>F)"] < 0.05) {cat("\n Main effect of 'traj_anyhk' is significant. Performing contrasts.\n\n"); emm_size_traj <- emmeans(modelo_size, specs = pairwise ~ traj_anyhk, adjust = "tukey"); df_emm_size_traj <- as.data.frame(emm_size_traj$contrasts); df_emm_size_traj$Significance <- ifelse(df_emm_size_traj$p.value < alfa, "Significant", "Not Significant"); print(df_emm_size_traj); write.csv(df_emm_size_traj, paste0("analysis_results_traj_anyhk/", tag, "/table_contrasts_trajectory_size", tag,".csv"), row.names = FALSE) }
    }
    cat("\n--- 1.5. Model Diagnostics (Size) ---\n")
    cat("Assessing VIF on a model without the interaction term to prevent artificial inflation:\n")
    modelo_size_vif <- lmer(tamanho_mean_average ~ wave + traj_anyhk + (1 | subjectid), data = nanosight_intersect)
    vif_df_size <- as.data.frame(check_collinearity(modelo_size_vif)); vif_df_size$Interpretation <- ifelse(vif_df_size$VIF >= 10, "High", ifelse(vif_df_size$VIF >= 5, "Moderate", "Low")); icc_size <- icc(modelo_size)
    cat("VIF:\n"); print(vif_df_size); write.csv(vif_df_size, paste0("analysis_results_traj_anyhk/", tag, "/table_vif_size_", tag, ".csv"), row.names = FALSE)
    cat("\nICC:\n"); print(icc_size)
    
    # MODEL 2: EV's concentration
    cat("\n\n====================================================================\n"); cat("  ANALYSIS 2: EV CONCENTRATION (concentracao_real)\n"); cat("====================================================================\n\n")
    cat("--- 2.1. Model Simplification Process for Concentration ---\n")
    modelo_conc <- NULL;
    equacao_conc_1 <- "glmmTMB(concentracao_real ~ wave * traj_anyhk + (1 | subjectid), data = nanosight_intersect, family = nbinom2(link = log)"
    cat("\nAttempt 1: Full Negative Binomial GLMM...\n", "equacao: ", equacao_conc_1, "\n", "\n")
    modelo_conc_tentativa1 <- tryCatch({glmmTMB(concentracao_real ~ wave * traj_anyhk + (1 | subjectid), data = nanosight_intersect, family = nbinom2(link = "log"))}, warning = function(w) {cat("CONVERGENCE WARNING in Attempt 1:\n", conditionMessage(w), "\n"); return(NULL)}, error = function(e) {cat("CONVERGENCE ERROR in Attempt 1:\n", conditionMessage(e), "\n"); return(NULL)})
    if (is.null(modelo_conc_tentativa1)) {equacao_conc_2 <- "glmmTMB(concentracao_real ~ wave + traj_anyhk + (1 | subjectid), data = nanosight_intersect, family = nbinom1(link = log)"; cat("\nAttempt 2: Negative Binomial GLMM (main effects only)...\n", "equacao: ", equacao_conc_2, "\n", "\n"); modelo_conc_tentativa2 <- tryCatch({glmmTMB(concentracao_real ~ wave + traj_anyhk + (1 | subjectid), data = nanosight_intersect, family = nbinom1(link = "log"))}, warning = function(w) {cat("CONVERGENCE WARNING in Attempt 2:\n", conditionMessage(w), "\n"); return(NULL)}, error = function(e) {cat("CONVERGENCE ERROR in Attempt 2:\n", conditionMessage(e), "\n"); return(NULL)})}
    if(is.null(modelo_conc_tentativa1) && is.null(modelo_conc_tentativa2)) {equacao_conc_3 <- "glmmTMB(concentracao_real ~ wave + traj_anyhk, data = nanosight_intersect, family = nbinom1(link = log)"; cat("\nBoth mixed models failed. Selecting final model: GLM with main effects only.\n", "equacao: ", equacao_conc_3,"\n", "\n"); modelo_conc <- glmmTMB(concentracao_real ~ wave + traj_anyhk, data = nanosight_intersect, family = nbinom1(link = "log"))}
    resultados_finais_aic_bic <- rbind(resultados_finais_aic_bic, data.frame(Variable = "Concentration", Model = "Negative Binomial (GLM)", AIC = AIC(modelo_conc), BIC = BIC(modelo_conc)))
    cat("\n--- 2.2. ANOVA Table (Type II) for Concentration ---\n")
    anova_conc <- Anova(modelo_conc, type = "II"); df_anova_conc <- as.data.frame(anova_conc); df_anova_conc$Significance <- ifelse(df_anova_conc$`Pr(>Chisq)` < alfa, "Significant", "Not Significant"); print(df_anova_conc); write.csv(df_anova_conc, paste0("analysis_results_traj_anyhk/table_anova_concentration_", tag, ".csv"))
    cat("\n--- 2.3. Pairwise Comparisons for Main Effects (Concentration) ---\n")
    if (df_anova_conc["wave", "Pr(>Chisq)"] < 0.05) {cat("\n Main effect of 'wave' is significant. Performing contrasts.\n\n"); emm_conc_wave <- emmeans(modelo_conc, specs = pairwise ~ wave, adjust = "tukey"); df_emm_conc_wave <- as.data.frame(emm_conc_wave$contrasts); df_emm_conc_wave$Significance <- ifelse(df_emm_conc_wave$p.value < alfa, "Significant", "Not Significant"); print(df_emm_conc_wave); write.csv(df_emm_conc_wave, paste0("analysis_results_traj_anyhk/", tag, "/table_contrasts_wave_concentration_", tag, ".csv"), row.names = FALSE)}
    if (df_anova_conc["traj_anyhk", "Pr(>Chisq)"] < 0.05) {cat("\n Main effect of 'traj_anyhk' is significant. Performing contrasts.\n\n"); emm_conc_traj <- emmeans(modelo_conc, specs = pairwise ~ traj_anyhk, adjust = "tukey"); df_emm_conc_traj <- as.data.frame(emm_conc_traj$contrasts); df_emm_conc_traj$Significance <- ifelse(df_emm_conc_traj$p.value < alfa, "Significant", "Not Significant"); print(df_emm_conc_traj); write.csv(df_emm_conc_traj, paste0("analysis_results_traj_anyhk/", tag, "/table_contrasts_trajectory_concentration_", tag, ".csv"), row.names = FALSE)}
    cat("\n--- 2.4. Model Diagnostics (Concentration) ---\n")
    vif_df_conc <- as.data.frame(check_collinearity(modelo_conc)); vif_df_conc$Interpretation <- ifelse(vif_df_conc$VIF >= 10, "High", ifelse(vif_df_conc$VIF >= 5, "Moderate", "Low"));
    cat("VIF:\n"); print(vif_df_conc); write.csv(vif_df_conc, paste0("analysis_results_traj_anyhk/", tag, "/table_vif_concentration_", tag, ".csv"), row.names = FALSE); cat("\nICC:\nNot applicable (GLM without random effects).\n")
    cat("\n--- 2.5. DHARMa Residual Diagnostics (Concentration) ---\n")
    residuos_conc <- simulateResiduals(fittedModel = modelo_conc); tabela_diag_conc <- criar_tabela_diagnostico_dharma_qq(residuos_conc)
    png(paste0("analysis_results_traj_anyhk/", tag, "/plot_diagnostics_concentration_", tag, ".png"), width = 3000, height = 1000, res = 300); par(mfrow = c(1, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 1, 0), cex = 1); plotQQunif(residuos_conc, main = "A) QQ Plot - Concentration", testUniformity = FALSE, testDispersion = FALSE, testOutliers = FALSE); plotResiduals(residuos_conc, main = "B) Residuals vs. Predicted - Concentration"); dev.off()
    print(tabela_diag_conc); write.csv(tabela_diag_conc, paste0("analysis_results_traj_anyhk/", tag, "/table_diagnostics_dharma_concentration_", tag, ".csv"), row.names = FALSE)
    
    # MODEL 3: EVs percentages
    cat("\n\n======================================================================\n"); cat("  ANALYSIS 3: PERCENTAGE OF SMALL EVs (p90_porcentagem)\n"); cat("======================================================================\n\n")
    equacao_perc <- "glmmTMB(percentage_prop ~ wave + traj_anyhk + (1 | subjectid), data = nanosight_intersect, family = beta_family(link = logit)"
    cat("--- 3.1. Fitting the Beta GLMM for Percentage ---\n", "equacao: ", equacao_perc, "\n")
    modelo_perc <- tryCatch({glmmTMB(percentage_prop ~ wave + traj_anyhk + (1 | subjectid), data = nanosight_intersect, family = beta_family(link = "logit"))})
    resultados_finais_aic_bic <- rbind(resultados_finais_aic_bic, data.frame(Variable = "Percentage", Model = "Beta (GLMM)", AIC = AIC(modelo_perc), BIC = BIC(modelo_perc)))
    cat("\n--- 3.2. ANOVA Table (Type III) for Percentage ---\n")
    anova_perc_type3 <- Anova(modelo_perc, type = "III"); df_anova_perc3 <- as.data.frame(anova_perc_type3); df_anova_perc3$Significance <- ifelse(df_anova_perc3$`Pr(>Chisq)` < alfa, "Significant", "Not Significant"); print(df_anova_perc3); write.csv(df_anova_perc3, paste0("analysis_results_traj_anyhk/", tag, "/table_anova_percentage_", tag, ".csv"))
    cat("\n--- 3.3. Conditional Pairwise Comparisons (Percentage) ---\n")
    #p_interacao_perc <- df_anova_perc3["wave:traj_anyhk", "Pr(>Chisq)"]
    #if (p_interacao_perc < 0.05) {
    #  cat("Interaction is significant. Analyzing both perspectives.\n")
    #emm_size_p1 <- emmeans(modelo_size, specs = pairwise ~ traj_anyhk | wave, adjust = "tukey"); df_emm_size_p1 <- as.data.frame(emm_size_p1$contrasts); df_emm_size_p1$Significance <- ifelse(df_emm_size_p1$p.value < alfa, "Significant", "Not Significant"); print(df_emm_size_p1); write.csv(df_emm_size_p1, paste0("analysis_results_traj_anyhk/", tag, "/table_contrasts_wave_perc_", tag, ".csv"), row.names = FALSE)
    #emm_size_p2 <- emmeans(modelo_size, specs = pairwise ~ wave | traj_anyhk, adjust = "tukey"); df_emm_size_p2 <- as.data.frame(emm_size_p2$contrasts); df_emm_size_p2$Significance <- ifelse(df_emm_size_p2$p.value < alfa, "Significant", "Not Significant"); print(df_emm_size_p2); write.csv(df_emm_size_p2, paste0("analysis_results_traj_anyhk/", tag, "/table_contrasts_traj_anyhk_perc_", tag, ".csv"), row.names = FALSE)
    #} else {
    #cat("Interaction is not significant. Analyzing main effects with a Type II ANOVA.\n")
    #anova_perc_type2 <- Anova(modelo_perc, type = "II"); df_anova_perc2 <- as.data.frame(anova_perc_type2); df_anova_perc2$Significance <- ifelse(df_anova_perc2$`Pr(>Chisq)` < alfa, "Significant", "Not Significant"); print(df_anova_perc2)
    #if (df_anova_perc2["wave", "Pr(>Chisq)"] < 0.05) {cat("\n Main effect of 'wave' is significant. Performing contrasts.\n\n"); emm_perc_wave <- emmeans(modelo_perc, specs = pairwise ~ wave, adjust = "tukey"); df_emm_perc_wave <- as.data.frame(emm_perc_wave$contrasts); df_emm_perc_wave$Significance <- ifelse(df_emm_perc_wave$p.value < alfa, "Significant", "Not Significant"); print(df_emm_perc_wave); write.csv(df_emm_perc_wave, paste0("analysis_results_traj_anyhk/", tag, "/table_contrasts_wave_percentage_", tag, ".csv"), row.names = FALSE) }
    #if (df_anova_perc2["traj_anyhk", "Pr(>Chisq)"] < 0.05) {cat("\n Main effect of 'traj_anyhk' is significant. Performing contrasts.\n\n");emm_perc_traj <- emmeans(modelo_perc, specs = pairwise ~ traj_anyhk, adjust = "tukey"); df_emm_perc_traj <- as.data.frame(emm_perc_traj$contrasts); df_emm_perc_traj$Significance <- ifelse(df_emm_perc_traj$p.value < alfa, "Significant", "Not Significant"); print(df_emm_perc_traj); write.csv(df_emm_perc_traj, paste0("analysis_results_traj_anyhk/", tag, "/table_contrasts_trajectory_percentage_", tag, ".csv"), row.names = FALSE) }
    #}
    cat("\n--- 3.4. Model Diagnostics (Percentage) ---\n")
    cat("Assessing VIF on a model without the interaction term to prevent artificial inflation:\n")
    modelo_perc_vif <- glmmTMB(percentage_prop ~ wave + traj_anyhk + (1 | subjectid), data = nanosight_intersect, family = beta_family(link = "logit"), dispformula = ~ traj_anyhk)
    vif_df_perc <- as.data.frame(check_collinearity(modelo_perc_vif)); vif_df_perc$Interpretation <- ifelse(vif_df_perc$VIF >= 10, "High", ifelse(vif_df_perc$VIF >= 5, "Moderate", "Low"));
    icc_perc <- icc(modelo_perc)
    cat("VIF:\n"); print(vif_df_perc); write.csv(vif_df_perc, paste0("analysis_results_traj_anyhk/", tag, "/table_vif_percentage_", tag, ".csv"), row.names = FALSE)
    cat("\nICC:\n"); print(icc_perc)
    cat("\n--- 3.5. DHARMa Residual Diagnostics (Percentage) ---\n")
    residuos_perc <- simulateResiduals(fittedModel = modelo_perc); tabela_diag_perc <- criar_tabela_diagnostico_dharma_qq(residuos_perc)
    png(paste0("analysis_results_traj_anyhk/", tag, "/plot_diagnostics_dharma_percentage_", tag, ".png"), width = 3000, height = 1000, res = 300); par(mfrow = c(1, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 1, 0), cex = 1); plotQQunif(residuos_perc, main = paste0("A) QQ Plot - Percentage ", tag), testUniformity = FALSE, testDispersion = FALSE, testOutliers = FALSE); plotResiduals(residuos_perc, main = paste0("B) Residuals vs. Predicted - Percentage", tag)); dev.off()
    print(tabela_diag_perc); write.csv(tabela_diag_perc, paste0("analysis_results_traj_anyhk/", tag, "/table_diagnostics_dharma_percentage_", tag, ".csv"), row.names = FALSE)
    # AIC/BIC
    cat("\n\n=========================================================\n"); cat("  SUMMARY TABLE: AIC/BIC OF FINAL MODELS\n"); cat("=========================================================\n\n")
    print(resultados_finais_aic_bic)
    write.csv(resultados_finais_aic_bic, paste0("analysis_results_traj_anyhk/", tag, "/table_summary_aic_bic_final_", tag, ".csv"), row.names = FALSE)
    
  }, finally = {
    # LOG
    cat("\n\n### FIM DO LOG DE ANÁLISE ###\n")
    cat("Todos os resultados foram salvos no diretório 'analysis_results_traj_anyhk'.\n")
    cat("O log completo desta sessão foi salvo em:", log_file_name, "\n")
    sink(type = "message")
    while(sink.number() > 0) sink()
    if(exists("con")) close(con)
  }
  )}

