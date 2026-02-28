setwd ("C:/Users/Belle/Documents/Belle - Nanosight/dados_raw")

library(dplyr)

###Criar coluna da trajetoria considerando controle original
nanosight_plus_sampleinfo <- read.csv("nanosight_plus_sampleinfo.csv")
nanosight_plus_sampleinfo$grupo_analise_gena <- NA #criar coluna vazia
##Atribuir grupos baseados em 'traj_gena'
# Normalizando para maiúscula inicial para ficar bonito no gráfico
nanosight_plus_sampleinfo$grupo_analise_gena[nanosight_plus_sampleinfo$traj_gena == "Incident"] <- "Incident"
nanosight_plus_sampleinfo$grupo_analise_gena[nanosight_plus_sampleinfo$traj_gena == "Remitted"] <- "Remitted"
nanosight_plus_sampleinfo$grupo_analise_gena[nanosight_plus_sampleinfo$traj_gena == "Persistent"] <- "Persistent"
# Atribuir grupo Control baseado na coluna 'Trajetoria' (SOBRESCREVENDO qualquer valor anterior se for o caso)
nanosight_plus_sampleinfo$grupo_analise_gena[nanosight_plus_sampleinfo$Trajetoria == "Control"] <- "Control"
# Filtrar NAs (remover quem não se encaixa em nenhum dos 4 grupos)
n_antes <- nrow(nanosight_plus_sampleinfo)
nanosight_plus_sampleinfo_gena_comremitted <- nanosight_plus_sampleinfo[!is.na(nanosight_plus_sampleinfo$grupo_analise_gena), ]
n_depois <- nrow(nanosight_plus_sampleinfo)
print(table(nanosight_plus_sampleinfo_gena_comremitted$grupo_analise_gena))
nanosight_plus_sampleinfo_gena_comremitted <- nanosight_plus_sampleinfo_gena_comremitted %>%
  filter(grupo_analise_gena != "Persistent")

##Retirando outliers
nanosight_plus_sampleinfo_gena_comremitted$zscore_mode <- scale(nanosight_plus_sampleinfo_gena_comremitted$tamanho_mode_average)
nanosight_plus_sampleinfo_gena_comremitted$zscore_porcentagem <- scale(nanosight_plus_sampleinfo_gena_comremitted$EV_pequenas_porcentagem)
nanosight_plus_sampleinfo_gena_comremitted$zscore_concentracao <- scale(nanosight_plus_sampleinfo_gena_comremitted$concentracao_real)
nanosight_plus_sampleinfo_gena_comremitted$zscore_conc_nan <- scale(nanosight_plus_sampleinfo_gena_comremitted$concentracao_average)
nanosight_plus_sampleinfo_gena_comremitted_sem_outliers_mode <- nanosight_plus_sampleinfo_gena_comremitted[nanosight_plus_sampleinfo_gena_comremitted$zscore_mode >= -3.0 & nanosight_plus_sampleinfo_gena_comremitted$zscore_mode <= 3.0, ]
nanosight_plus_sampleinfo_gena_comremitted_sem_outliers_porcentagem <- nanosight_plus_sampleinfo_gena_comremitted[nanosight_plus_sampleinfo_gena_comremitted$zscore_porcentagem >= -3.0 & nanosight_plus_sampleinfo_gena_comremitted$zscore_porcentagem <= 3.0, ]
nanosight_plus_sampleinfo_gena_comremitted_sem_outliers_concentracao <- nanosight_plus_sampleinfo_gena_comremitted[nanosight_plus_sampleinfo_gena_comremitted$zscore_concentracao >= -3.0 & nanosight_plus_sampleinfo_gena_comremitted$zscore_concentracao <= 3.0, ]
nanosight_plus_sampleinfo_gena_comremitted_sem_outliers_conc_nan <- nanosight_plus_sampleinfo_gena_comremitted[nanosight_plus_sampleinfo_gena_comremitted$zscore_conc_nan >= -3.0 & nanosight_plus_sampleinfo_gena_comremitted$zscore_conc_nan <= 3.0, ]
nanosight_plus_sampleinfo_gena_comremitted_sem_outliers_mode_concentracao <- semi_join(nanosight_plus_sampleinfo_gena_comremitted_sem_outliers_mode, nanosight_plus_sampleinfo_gena_comremitted_sem_outliers_concentracao, by = "id_sample")
nanosight_plus_sampleinfo_gena_comremitted_sem_outliers_porcentagem_concentracao <- semi_join(nanosight_plus_sampleinfo_gena_comremitted_sem_outliers_porcentagem, nanosight_plus_sampleinfo_gena_comremitted_sem_outliers_concentracao, by = ("id_sample"))
nanosight_plus_sampleinfo_gena_comremitted_sem_outliers_mode_porcentagem <- semi_join(nanosight_plus_sampleinfo_gena_comremitted_sem_outliers_mode, nanosight_plus_sampleinfo_gena_comremitted_sem_outliers_porcentagem, by = "id_sample")
nanosight_intersect_gena_comremitted_mode <- semi_join(nanosight_plus_sampleinfo_gena_comremitted_sem_outliers_mode_concentracao, nanosight_plus_sampleinfo_gena_comremitted_sem_outliers_porcentagem, by = ("id_sample"))
nanosight_outliers_gena_comremitted_mode <- anti_join(nanosight_plus_sampleinfo_gena_comremitted, nanosight_intersect_gena_comremitted_mode, by = "id_sample")

##Tabela só com w1 e só w2
nanosight_w1_gena_comremitted_mode <- subset(nanosight_intersect_gena_comremitted_mode, wave == "t1")
nanosight_w2_gena_comremitted_mode <- subset(nanosight_intersect_gena_comremitted_mode, wave == "t2")
nanosight_intersect_pares_gena_comremitted_mode <- semi_join(nanosight_w1_gena_comremitted_mode, nanosight_w2_gena_comremitted_mode, by = ("subjectid"))


