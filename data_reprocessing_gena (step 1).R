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
nanosight_plus_sampleinfo_gena <- nanosight_plus_sampleinfo[!is.na(nanosight_plus_sampleinfo$grupo_analise_gena), ]
n_depois <- nrow(nanosight_plus_sampleinfo)
print(table(nanosight_plus_sampleinfo_gena$grupo_analise_gena))
nanosight_plus_sampleinfo_gena <- nanosight_plus_sampleinfo_gena %>%
  filter(grupo_analise_gena != "Persistent")
nanosight_plus_sampleinfo_gena <- nanosight_plus_sampleinfo_gena %>%
  filter(grupo_analise_gena != "Remitted")
print(table(nanosight_plus_sampleinfo_gena$grupo_analise_gena))
##Retirando outliers
nanosight_plus_sampleinfo_gena$zscore_mean <- scale(nanosight_plus_sampleinfo_gena$tamanho_mean_average)
nanosight_plus_sampleinfo_gena$zscore_porcentagem <- scale(nanosight_plus_sampleinfo_gena$EV_pequenas_porcentagem)
nanosight_plus_sampleinfo_gena$zscore_concentracao <- scale(nanosight_plus_sampleinfo_gena$concentracao_real)
nanosight_plus_sampleinfo_gena$zscore_conc_nan <- scale(nanosight_plus_sampleinfo_gena$concentracao_average)
nanosight_plus_sampleinfo_gena_sem_outliers_mean <- nanosight_plus_sampleinfo_gena[nanosight_plus_sampleinfo_gena$zscore_mean >= -3.0 & nanosight_plus_sampleinfo_gena$zscore_mean <= 3.0, ]
nanosight_plus_sampleinfo_gena_sem_outliers_porcentagem <- nanosight_plus_sampleinfo_gena[nanosight_plus_sampleinfo_gena$zscore_porcentagem >= -3.0 & nanosight_plus_sampleinfo_gena$zscore_porcentagem <= 3.0, ]
nanosight_plus_sampleinfo_gena_sem_outliers_concentracao <- nanosight_plus_sampleinfo_gena[nanosight_plus_sampleinfo_gena$zscore_concentracao >= -3.0 & nanosight_plus_sampleinfo_gena$zscore_concentracao <= 3.0, ]
nanosight_plus_sampleinfo_gena_sem_outliers_conc_nan <- nanosight_plus_sampleinfo_gena[nanosight_plus_sampleinfo_gena$zscore_conc_nan >= -3.0 & nanosight_plus_sampleinfo_gena$zscore_conc_nan <= 3.0, ]
nanosight_plus_sampleinfo_gena_sem_outliers_mean_concentracao <- semi_join(nanosight_plus_sampleinfo_gena_sem_outliers_mean, nanosight_plus_sampleinfo_gena_sem_outliers_concentracao, by = "id_sample")
nanosight_plus_sampleinfo_gena_sem_outliers_porcentagem_concentracao <- semi_join(nanosight_plus_sampleinfo_gena_sem_outliers_porcentagem, nanosight_plus_sampleinfo_gena_sem_outliers_concentracao, by = ("id_sample"))
nanosight_plus_sampleinfo_gena_sem_outliers_mean_porcentagem <- semi_join(nanosight_plus_sampleinfo_gena_sem_outliers_mean, nanosight_plus_sampleinfo_gena_sem_outliers_porcentagem, by = "id_sample")
nanosight_intersect_gena <- semi_join(nanosight_plus_sampleinfo_gena_sem_outliers_mean_concentracao, nanosight_plus_sampleinfo_gena_sem_outliers_porcentagem, by = ("id_sample"))
nanosight_outliers_gena <- anti_join(nanosight_plus_sampleinfo_gena, nanosight_intersect_gena, by = "id_sample")

##Tabela só com w1 e só w2
nanosight_w1_gena <- subset(nanosight_intersect_gena, wave == "t1")
nanosight_w2_gena <- subset(nanosight_intersect_gena, wave == "t2")
nanosight_intersect_pares_gena <- semi_join(nanosight_w1_gena, nanosight_w2_gena, by = ("subjectid"))


