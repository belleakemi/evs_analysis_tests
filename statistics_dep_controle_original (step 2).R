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
  criar_tabela_diagnostico_dharma_qq <- function(simulationOutput) {
    alfa <- 0.05
    teste_uniformidade <- testUniformity(simulationOutput, plot = FALSE)
    teste_dispersao <- testDispersion(simulationOutput, plot = FALSE)
    teste_outliers <- testOutliers(simulationOutput, plot = FALSE)
    tabela <- data.frame(Diagnostic_Test = c("Overall Uniformity (KS)", "Dispersion", "Outliers"), 
                         p_value = c(round(teste_uniformidade$p.value, 3), round(teste_dispersao$p.value, 3), round(teste_outliers$p.value, 3)))
    tabela$Significance <- ifelse(tabela$p_value < alfa, "Significant Violation", "OK")
    return(tabela)
  }
  
  resultados_finais_aic_bic <- data.frame()
  alfa <- 0.05 
  
  # MODEL 1: EV's Sizes (USANDO grupo_analise_dep)
  cat("\n\n=========================================================\n"); cat("  ANALYSIS 1: EV MEAN SIZE (tamanho_mean_average)\n"); cat("=========================================================\n\n")
  cat("--- 1.1. Comparing Gaussian (LMM) and Gamma (GLMM) models ---\n")
  modelo_size_gaussiano <- lmer(tamanho_mean_average ~ wave * grupo_analise_dep + (1 | subjectid), data = nanosight_intersect_ev_pequena_long)
  modelo_size_gamma <- glmmTMB(tamanho_mean_average ~ wave * grupo_analise_dep + (1 | subjectid), data = nanosight_intersect_ev_pequena_long, family = Gamma(link = "log"))
  
  tabela_comp_size <- data.frame(Model = c("Gaussian (LMM)", "Gamma (GLMM)"), AIC = c(AIC(modelo_size_gaussiano), AIC(modelo_size_gamma)), BIC = c(BIC(modelo_size_gaussiano), BIC(modelo_size_gamma)))
  print(tabela_comp_size); write.csv(tabela_comp_size, paste0(output_dir, "/table_model_comparison_size.csv"), row.names = FALSE)
  
  cat("\n--- 1.2. DHARMa Residual Diagnostics for BOTH models ---\n")
  residuos_gamma <- simulateResiduals(fittedModel = modelo_size_gamma)
  tabela_diag_gamma <- criar_tabela_diagnostico_dharma_qq(residuos_gamma)
  png(paste0(output_dir, "/plot_diagnostics_size_GAMMA.png"), width = 3000, height = 1000, res = 300)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 1, 0), cex = 1)
  plotQQunif(residuos_gamma, main = "A) QQ Plot - Size (Gamma)", testUniformity = FALSE, testDispersion = FALSE, testOutliers = FALSE)
  plotResiduals(residuos_gamma, main = "B) Residuals vs. Predicted - Size (Gamma)")
  dev.off()
  
  residuos_gaussiano <- simulateResiduals(fittedModel = modelo_size_gaussiano)
  tabela_diag_gaussiano <- criar_tabela_diagnostico_dharma_qq(residuos_gaussiano)
  png(paste0(output_dir, "/plot_diagnostics_size_GAUSSIAN.png"), width = 3000, height = 1000, res = 300)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 1, 0), cex = 1)
  plotQQunif(residuos_gaussiano, main = "A) QQ Plot - Size (Gaussian)", testUniformity = FALSE, testDispersion = FALSE, testOutliers = FALSE)
  plotResiduals(residuos_gaussiano, main = "B) Residuals vs. Predicted - Size (Gaussian)")
  dev.off()
  
  modelo_size <- modelo_size_gaussiano
  cat("\nFinal Model Selected: Gaussian (LMM).\n")
  resultados_finais_aic_bic <- rbind(resultados_finais_aic_bic, data.frame(Variable = "Size", Model = "Gaussian (LMM)", AIC = AIC(modelo_size), BIC = BIC(modelo_size)))
  
  cat("\n--- 1.3. ANOVA Table (Type III) for Size ---\n")
  anova_size <- anova(modelo_size, type = 3)
  df_anova_size <- as.data.frame(anova_size)
  df_anova_size$Significance <- ifelse(df_anova_size$`Pr(>F)` < alfa, "Significant", "Not Significant")
  print(df_anova_size); write.csv(df_anova_size, paste0(output_dir, "/table_anova_size.csv"))
  
  cat("\n--- 1.4. Conditional Pairwise Comparisons (Size) ---\n")
  # Verificar se existe interação significativa. Se não, analisar efeitos principais.
  # Nota: lmerTest anova retorna Pr(>F).
  p_interacao_size <- tryCatch(df_anova_size["wave:grupo_analise_dep", "Pr(>F)"], error = function(e) 1)
  if(is.na(p_interacao_size)) p_interacao_size <- 1
  
  if (p_interacao_size < 0.05) {
    cat("Interaction is significant. Performing contrasts.\n")
    emm_size_p1 <- emmeans(modelo_size, specs = pairwise ~ grupo_analise_dep | wave, adjust = "tukey")
    write.csv(as.data.frame(emm_size_p1$contrasts), paste0(output_dir, "/table_contrasts_perspective1_size.csv"), row.names = FALSE)
    emm_size_p2 <- emmeans(modelo_size, specs = pairwise ~ wave | grupo_analise_dep, adjust = "tukey")
    write.csv(as.data.frame(emm_size_p2$contrasts), paste0(output_dir, "/table_contrasts_perspective2_size.csv"), row.names = FALSE)
  } else {
    cat("Interaction NOT significant. Analyzing Main Effects.\n")
    emm_size_group <- emmeans(modelo_size, specs = pairwise ~ grupo_analise_dep, adjust = "tukey")
    write.csv(as.data.frame(emm_size_group$contrasts), paste0(output_dir, "/table_contrasts_group_size.csv"), row.names = FALSE)
    print(emm_size_group$contrasts)
  }
  
  cat("\n--- 1.5. VIF (Size) ---\n")
  modelo_size_vif <- lmer(tamanho_mean_average ~ wave + grupo_analise_dep + (1 | subjectid), data = nanosight_intersect_ev_pequena_long)
  vif_df_size <- as.data.frame(check_collinearity(modelo_size_vif))
  print(vif_df_size); write.csv(vif_df_size, paste0(output_dir, "/table_vif_size.csv"), row.names = FALSE)
  
  # MODEL 2: EV's concentration
  cat("\n\n====================================================================\n"); cat("  ANALYSIS 2: EV CONCENTRATION (concentracao_real)\n"); cat("====================================================================\n\n")
  modelo_conc <- NULL
  cat("\nAttempt 1: Full Negative Binomial GLMM...\n")
  modelo_conc_tentativa1 <- tryCatch({glmmTMB(concentracao_real ~ wave * grupo_analise_dep + (1 | subjectid), data = nanosight_intersect_ev_pequena_long, family = nbinom2(link = "log"))}, warning = function(w) NULL, error = function(e) NULL)
  
  if (is.null(modelo_conc_tentativa1)) {
    cat("\nAttempt 2: Main effects only...\n")
    modelo_conc <- glmmTMB(concentracao_real ~ wave + grupo_analise_dep + (1 | subjectid), data = nanosight_intersect_ev_pequena_long, family = nbinom1(link = "log"))
  } else {
    modelo_conc <- modelo_conc_tentativa1
  }
  
  resultados_finais_aic_bic <- rbind(resultados_finais_aic_bic, data.frame(Variable = "Concentration", Model = "Negative Binomial", AIC = AIC(modelo_conc), BIC = BIC(modelo_conc)))
  
  cat("\n--- 2.2. ANOVA Table (Type II) for Concentration ---\n")
  anova_conc <- Anova(modelo_conc, type = "II")
  print(anova_conc); write.csv(as.data.frame(anova_conc), paste0(output_dir, "/table_anova_concentration.csv"))
  
  # Contrastes para Concentração (Exemplo Main Effects se interação não for foco do modelo simplificado)
  emm_conc_group <- emmeans(modelo_conc, specs = pairwise ~ grupo_analise_dep, adjust = "tukey")
  write.csv(as.data.frame(emm_conc_group$contrasts), paste0(output_dir, "/table_contrasts_group_concentration.csv"), row.names = FALSE)
  
  cat("\n--- 2.5. DHARMa Residual Diagnostics (Concentration) ---\n")
  residuos_conc <- simulateResiduals(fittedModel = modelo_conc)
  png(paste0(output_dir, "/plot_diagnostics_concentration.png"), width = 3000, height = 1000, res = 300)
  par(mfrow = c(1, 2)); plotQQunif(residuos_conc); plotResiduals(residuos_conc); dev.off()
  
  # MODEL 3: SMALL EVs percentages
  cat("\n\n======================================================================\n"); cat("  ANALYSIS 3: PERCENTAGE OF SMALL EVs (EV_pequenas_porcentagem)\n"); cat("======================================================================\n\n")
  cat("--- 3.1. Fitting the Beta GLMM for Percentage ---\n")
  # NOTA: Ajuste dispformula para usar grupo_analise_dep
  modelo_perc <- glmmTMB(percentage_prop ~ wave * grupo_analise_dep + (1 | subjectid), 
                         data = nanosight_intersect_ev_pequena_long, 
                         family = beta_family(link = "logit"), 
                         dispformula = ~ grupo_analise_dep)
  
  resultados_finais_aic_bic <- rbind(resultados_finais_aic_bic, data.frame(Variable = "Percentage", Model = "Beta (GLMM)", AIC = AIC(modelo_perc), BIC = BIC(modelo_perc)))
  
  cat("\n--- 3.2. ANOVA Table (Type III) for Percentage ---\n")
  anova_perc_type3 <- Anova(modelo_perc, type = "III")
  print(anova_perc_type3); write.csv(as.data.frame(anova_perc_type3), paste0(output_dir, "/table_anova_percentage.csv"))
  
  cat("\n--- 3.3. Contrasts Percentage ---\n")
  # Verificando interação
  p_interacao_perc <- as.data.frame(anova_perc_type3)["wave:grupo_analise_dep", "Pr(>Chisq)"]
  if (!is.na(p_interacao_perc) && p_interacao_perc < 0.05) {
    emm_perc_p1 <- emmeans(modelo_perc, specs = pairwise ~ grupo_analise_dep | wave, adjust = "tukey")
    write.csv(as.data.frame(emm_perc_p1$contrasts), paste0(output_dir, "/table_contrasts_perspective1_percentage.csv"), row.names = FALSE)
  } else {
    emm_perc_main <- emmeans(modelo_perc, specs = pairwise ~ grupo_analise_dep, adjust = "tukey")
    write.csv(as.data.frame(emm_perc_main$contrasts), paste0(output_dir, "/table_contrasts_group_percentage.csv"), row.names = FALSE)
  }
  
  cat("\n--- 3.5. DHARMa Residual Diagnostics (Percentage) ---\n")
  residuos_perc <- simulateResiduals(fittedModel = modelo_perc)
  png(paste0(output_dir, "/plot_diagnostics_percentage.png"), width = 3000, height = 1000, res = 300)
  par(mfrow = c(1, 2)); plotQQunif(residuos_perc); plotResiduals(residuos_perc); dev.off()
  
  # SUMMARY
  print(resultados_finais_aic_bic)
  write.csv(resultados_finais_aic_bic, paste0(output_dir, "/table_summary_aic_bic_final.csv"), row.names = FALSE)
  
}, finally = {
  cat("\n\n### FIM DO LOG DE ANÁLISE ###\n")
  sink() 
})

