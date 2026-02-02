# =========================================================
# MODEL 1: EV MEAN SIZE (COM CORREÇÃO DE DISPERSÃO)
# =========================================================
cat("\n\n=========================================================\n"); cat("  ANALYSIS 1: EV MEAN SIZE (tamanho_mean_average)\n"); cat("=========================================================\n\n")

cat("--- 1.1. Fitting Gaussian Model with Dispersion Formula ---\n")

# Tenta rodar o modelo com glmmTMB + dispformula
modelo_size <- tryCatch({
  glmmTMB(tamanho_mean_average ~ wave * grupo_analise_gena + (1 | subjectid), 
          data = nanosight_intersect_ev_pequena_long, 
          family = gaussian(), 
          dispformula = ~ grupo_analise_gena)
}, error = function(e) {
  cat("Erro no modelo misto de tamanho. Tentando GLM simples (sem random effects)...\n")
  glmmTMB(tamanho_mean_average ~ wave * grupo_analise_gena, 
          data = nanosight_intersect_ev_pequena_long, 
          family = gaussian(), 
          dispformula = ~ grupo_analise_gena)
})

if(!is.null(modelo_size)){
  resultados_finais_aic_bic <- rbind(resultados_finais_aic_bic, data.frame(Variable = "Size", Model = "Gaussian (glmmTMB)", AIC = AIC(modelo_size), BIC = BIC(modelo_size)))
  
  cat("\n--- 1.2. DHARMa Residual Diagnostics ---\n")
  residuos_size <- simulateResiduals(fittedModel = modelo_size)
  tabela_diag_size <- criar_tabela_diagnostico_dharma_qq(residuos_size)
  
  png(paste0(output_dir, "/plot_diagnostics_size_GAUSSIAN_DISP.png"), width = 3000, height = 1000, res = 300)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 1, 0), cex = 1)
  plotQQunif(residuos_size, main = "A) QQ Plot - Size", testUniformity = FALSE, testDispersion = FALSE, testOutliers = FALSE)
  plotResiduals(residuos_size, main = "B) Residuals vs. Predicted - Size")
  dev.off()
  print(tabela_diag_size); write.csv(tabela_diag_size, paste0(output_dir, "/table_diagnostics_dharma_size.csv"), row.names = FALSE)
  
  cat("\n--- 1.3. ANOVA Table (Type III) for Size ---\n")
  anova_size <- Anova(modelo_size, type = "III")
  df_anova_size <- as.data.frame(anova_size)
  df_anova_size$Significance <- ifelse(df_anova_size$`Pr(>Chisq)` < alfa, "Significant", "Not Significant")
  print(df_anova_size); write.csv(df_anova_size, paste0(output_dir, "/table_anova_size.csv"))
  
  cat("\n--- 1.4. Conditional Pairwise Comparisons (Size) ---\n")
  p_interacao_size <- df_anova_size["wave:grupo_analise_gena", "Pr(>Chisq)"]
  if (!is.na(p_interacao_size) && p_interacao_size < 0.05) {
    cat("Interaction is significant.\n")
    emm_size_p1 <- emmeans(modelo_size, specs = pairwise ~ grupo_analise_gena | wave, adjust = "tukey")
    write.csv(as.data.frame(emm_size_p1$contrasts), paste0(output_dir, "/table_contrasts_perspective1_size.csv"), row.names = FALSE)
    emm_size_p2 <- emmeans(modelo_size, specs = pairwise ~ wave | grupo_analise_gena, adjust = "tukey")
    write.csv(as.data.frame(emm_size_p2$contrasts), paste0(output_dir, "/table_contrasts_perspective2_size.csv"), row.names = FALSE)
  } else {
    cat("Interaction NOT significant. Analyzing Main Effects (Type II).\n")
    anova_size_type2 <- Anova(modelo_size, type = "II")
    print(anova_size_type2)
    
    if (anova_size_type2["wave", "Pr(>Chisq)"] < 0.05) {
      emm_size_wave <- emmeans(modelo_size, specs = pairwise ~ wave, adjust = "tukey")
      write.csv(as.data.frame(emm_size_wave$contrasts), paste0(output_dir, "/table_contrasts_wave_size.csv"), row.names = FALSE)
    }
    if (anova_size_type2["grupo_analise_gena", "Pr(>Chisq)"] < 0.05) {
      cat("Main effect of Group is significant. Performing Contrasts...\n")
      emm_size_traj <- emmeans(modelo_size, specs = pairwise ~ grupo_analise_gena, adjust = "tukey")
      print(emm_size_traj$contrasts)
      write.csv(as.data.frame(emm_size_traj$contrasts), paste0(output_dir, "/table_contrasts_trajectory_size.csv"), row.names = FALSE)
    }
  }
  
  cat("\n--- 1.5. VIF and ICC (Size) ---\n")
  modelo_size_vif <- glmmTMB(tamanho_mean_average ~ wave + grupo_analise_gena + (1 | subjectid), 
                             data = nanosight_intersect_ev_pequena_long, family = gaussian(), dispformula = ~ grupo_analise_gena)
  vif_df_size <- as.data.frame(check_collinearity(modelo_size_vif))
  print(vif_df_size); write.csv(vif_df_size, paste0(output_dir, "/table_vif_size.csv"), row.names = FALSE)
  
  # ICC com proteção contra erro
  tryCatch({
    icc_size <- icc(modelo_size)
    cat("\nICC:\n"); print(icc_size)
  }, error = function(e) { cat("\nCould not calculate ICC (likely singularity/zero variance in random effects).\n") })
}

# ====================================================================
# MODEL 2: EV CONCENTRATION (CORRIGIDO PARA EVITAR ERRO DE NULL)
# ====================================================================
cat("\n\n====================================================================\n"); cat("  ANALYSIS 2: EV CONCENTRATION (concentracao_real)\n"); cat("====================================================================\n\n")

modelo_conc <- NULL

# Tentativa 1: Modelo Completo (Misto)
cat("Attempt 1: Full Negative Binomial GLMM...\n")
modelo_conc <- tryCatch({
  glmmTMB(concentracao_real ~ wave * grupo_analise_gena + (1 | subjectid), 
          data = nanosight_intersect_ev_pequena_long, family = nbinom2(link = "log"))
}, error = function(e) { cat("Attempt 1 failed.\n"); return(NULL) })

# Tentativa 2: Efeitos Principais (Misto) - Se o 1 falhou
if (is.null(modelo_conc)) {
  cat("Attempt 2: Main effects only (GLMM)...\n")
  modelo_conc <- tryCatch({
    glmmTMB(concentracao_real ~ wave + grupo_analise_gena + (1 | subjectid), 
            data = nanosight_intersect_ev_pequena_long, family = nbinom1(link = "log"))
  }, error = function(e) { cat("Attempt 2 failed.\n"); return(NULL) })
}

# Tentativa 3: GLM Simples (Sem Random Effects) - Se o 2 falhou
if (is.null(modelo_conc)) {
  cat("Attempt 3: Simple GLM (No random effects)...\n")
  modelo_conc <- tryCatch({
    glmmTMB(concentracao_real ~ wave + grupo_analise_gena, 
            data = nanosight_intersect_ev_pequena_long, family = nbinom1(link = "log"))
  }, error = function(e) { cat("Attempt 3 failed. Cannot fit concentration model.\n"); return(NULL) })
}

# SÓ EXECUTA O RESTO SE TIVERMOS UM MODELO VÁLIDO
if (!is.null(modelo_conc)) {
  resultados_finais_aic_bic <- rbind(resultados_finais_aic_bic, data.frame(Variable = "Concentration", Model = "Negative Binomial", AIC = AIC(modelo_conc), BIC = BIC(modelo_conc)))
  
  cat("\n--- 2.2. ANOVA Table (Type II) for Concentration ---\n")
  anova_conc <- Anova(modelo_conc, type = "II")
  df_anova_conc <- as.data.frame(anova_conc)
  df_anova_conc$Significance <- ifelse(df_anova_conc$`Pr(>Chisq)` < alfa, "Significant", "Not Significant")
  print(df_anova_conc); write.csv(df_anova_conc, paste0(output_dir, "/table_anova_concentration.csv"))
  
  cat("\n--- 2.3. Pairwise Comparisons (Concentration) ---\n")
  # Wave
  if (df_anova_conc["wave", "Pr(>Chisq)"] < 0.05) {
    cat("Wave is significant.\n")
    emm_conc_wave <- emmeans(modelo_conc, specs = pairwise ~ wave, adjust = "tukey")
    write.csv(as.data.frame(emm_conc_wave$contrasts), paste0(output_dir, "/table_contrasts_wave_concentration.csv"), row.names = FALSE)
  }
  # Group
  if (df_anova_conc["grupo_analise_gena", "Pr(>Chisq)"] < 0.05) {
    cat("Group is significant.\n")
    emm_conc_traj <- emmeans(modelo_conc, specs = pairwise ~ grupo_analise_gena, adjust = "tukey")
    print(emm_conc_traj$contrasts) # IMPRIMIR PARA VER NO LOG
    write.csv(as.data.frame(emm_conc_traj$contrasts), paste0(output_dir, "/table_contrasts_trajectory_concentration.csv"), row.names = FALSE)
  }
  
  cat("\n--- 2.5. DHARMa Diagnostics (Concentration) ---\n")
  residuos_conc <- simulateResiduals(fittedModel = modelo_conc)
  tabela_diag_conc <- criar_tabela_diagnostico_dharma_qq(residuos_conc)
  png(paste0(output_dir, "/plot_diagnostics_concentration.png"), width = 3000, height = 1000, res = 300)
  par(mfrow = c(1, 2)); plotQQunif(residuos_conc); plotResiduals(residuos_conc); dev.off()
  print(tabela_diag_conc); write.csv(tabela_diag_conc, paste0(output_dir, "/table_diagnostics_dharma_concentration.csv"), row.names = FALSE)
  
} else {
  cat("CRITICAL ERROR: Could not fit any model for Concentration.\n")
}

# ======================================================================
# MODEL 3: SMALL EVs PERCENTAGES
# ======================================================================
cat("\n\n======================================================================\n"); cat("  ANALYSIS 3: PERCENTAGE OF SMALL EVs (EV_pequenas_porcentagem)\n"); cat("======================================================================\n\n")

cat("--- 3.1. Fitting Beta GLMM ---\n")
modelo_perc <- tryCatch({
  glmmTMB(percentage_prop ~ wave * grupo_analise_gena + (1 | subjectid), 
          data = nanosight_intersect_ev_pequena_long, 
          family = beta_family(link = "logit"),
          dispformula = ~ grupo_analise_gena)
}, error = function(e) {
  cat("Beta model failed. Trying without random effects.\n")
  glmmTMB(percentage_prop ~ wave * grupo_analise_gena, 
          data = nanosight_intersect_ev_pequena_long, 
          family = beta_family(link = "logit"),
          dispformula = ~ grupo_analise_gena)
})

if(!is.null(modelo_perc)){
  resultados_finais_aic_bic <- rbind(resultados_finais_aic_bic, data.frame(Variable = "Percentage", Model = "Beta (GLMM)", AIC = AIC(modelo_perc), BIC = BIC(modelo_perc)))
  
  cat("\n--- 3.2. ANOVA Table (Type III) for Percentage ---\n")
  anova_perc_type3 <- Anova(modelo_perc, type = "III")
  print(anova_perc_type3); write.csv(as.data.frame(anova_perc_type3), paste0(output_dir, "/table_anova_percentage.csv"))
  
  cat("\n--- 3.3. Contrasts Percentage ---\n")
  p_interacao_perc <- as.data.frame(anova_perc_type3)["wave:grupo_analise_gena", "Pr(>Chisq)"]
  
  if (!is.na(p_interacao_perc) && p_interacao_perc < 0.05) {
    cat("Interaction significant.\n")
    emm_perc_p1 <- emmeans(modelo_perc, specs = pairwise ~ grupo_analise_gena | wave, adjust = "tukey")
    write.csv(as.data.frame(emm_perc_p1$contrasts), paste0(output_dir, "/table_contrasts_perspective1_percentage.csv"), row.names = FALSE)
  } else {
    cat("Interaction NOT significant. Main effects (Type II).\n")
    anova_perc_type2 <- Anova(modelo_perc, type = "II")
    print(anova_perc_type2)
    
    if (as.data.frame(anova_perc_type2)["grupo_analise_gena", "Pr(>Chisq)"] < 0.05) {
      emm_perc_main <- emmeans(modelo_perc, specs = pairwise ~ grupo_analise_gena, adjust = "tukey")
      print(emm_perc_main$contrasts)
      write.csv(as.data.frame(emm_perc_main$contrasts), paste0(output_dir, "/table_contrasts_group_percentage.csv"), row.names = FALSE)
    }
  }
  
  cat("\n--- 3.5. DHARMa Residual Diagnostics (Percentage) ---\n")
  residuos_perc <- simulateResiduals(fittedModel = modelo_perc)
  png(paste0(output_dir, "/plot_diagnostics_percentage.png"), width = 3000, height = 1000, res = 300)
  par(mfrow = c(1, 2)); plotQQunif(residuos_perc); plotResiduals(residuos_perc); dev.off()
}

# SUMMARY
cat("\n\n=========================================================\n"); cat("  SUMMARY TABLE\n"); cat("=========================================================\n\n")
print(resultados_finais_aic_bic)
write.csv(resultados_finais_aic_bic, paste0(output_dir, "/table_summary_aic_bic_final.csv"), row.names = FALSE)
