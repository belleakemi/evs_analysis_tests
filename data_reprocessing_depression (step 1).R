setwd ("C:/Users/Belle/Documents/Belle - Nanosight/dados_raw")
nanosight_plus_sampleinfo <- read.csv("nanosight_plus_sampleinfo.csv")
library(dplyr)

## Atribuir grupos baseados em 'traj_dep' (Assumindo que traj_dcmadep se refere a traj_dep)
# 1.Atribuir grupos Incident, Remitted e Persistent baseado na coluna traj_dep
nanosight_plus_sampleinfo$grupo_analise_dep[nanosight_plus_sampleinfo$traj_dep == "Incident"] <- "Incident"
nanosight_plus_sampleinfo$grupo_analise_dep[nanosight_plus_sampleinfo$traj_dep == "Remitted"] <- "Remitted"
nanosight_plus_sampleinfo$grupo_analise_dep[nanosight_plus_sampleinfo$traj_dep == "Persistent"] <- "Persistent"
# 2. Atribuir grupo Control baseado na coluna 'Trajetoria' (SOBRESCREVENDO qualquer valor anterior se for o caso)
nanosight_plus_sampleinfo$grupo_analise_dep[nanosight_plus_sampleinfo$Trajetoria == "Control"] <- "Control"
# 3. Filtrar NAs (remover quem não se encaixa em nenhum dos 4 grupos)
n_antes <- nrow(nanosight_plus_sampleinfo)
nanosight_plus_sampleinfo_dep <- nanosight_plus_sampleinfo[!is.na(nanosight_plus_sampleinfo$grupo_analise_dep), ]
n_depois <- nrow(nanosight_plus_sampleinfo_dep)

cat(paste("Filtragem realizada. Linhas antes:", n_antes, "| Linhas depois:", n_depois, "\n"))
cat("Contagem por grupo:\n")
print(table(nanosight_plus_sampleinfo_dep$grupo_analise_dep))

##Retirando outliers
nanosight_plus_sampleinfo_dep$zscore_mean <- scale(nanosight_plus_sampleinfo_dep$tamanho_mean_average)
nanosight_plus_sampleinfo_dep$zscore_porcentagem <- scale(nanosight_plus_sampleinfo_dep$EV_pequenas_porcentagem)
nanosight_plus_sampleinfo_dep$zscore_concentracao <- scale(nanosight_plus_sampleinfo_dep$concentracao_real)
nanosight_plus_sampleinfo_dep$zscore_conc_nan <- scale(nanosight_plus_sampleinfo_dep$concentracao_average)
nanosight_plus_sampleinfo_dep_sem_outliers_mean <- nanosight_plus_sampleinfo_dep[nanosight_plus_sampleinfo_dep$zscore_mean >= -3.0 & nanosight_plus_sampleinfo_dep$zscore_mean <= 3.0, ]
nanosight_plus_sampleinfo_dep_sem_outliers_porcentagem <- nanosight_plus_sampleinfo_dep[nanosight_plus_sampleinfo_dep$zscore_porcentagem >= -3.0 & nanosight_plus_sampleinfo_dep$zscore_porcentagem <= 3.0, ]
nanosight_plus_sampleinfo_dep_sem_outliers_concentracao <- nanosight_plus_sampleinfo_dep[nanosight_plus_sampleinfo_dep$zscore_concentracao >= -3.0 & nanosight_plus_sampleinfo_dep$zscore_concentracao <= 3.0, ]
nanosight_plus_sampleinfo_dep_sem_outliers_conc_nan <- nanosight_plus_sampleinfo_dep[nanosight_plus_sampleinfo_dep$zscore_conc_nan >= -3.0 & nanosight_plus_sampleinfo_dep$zscore_conc_nan <= 3.0, ]
nanosight_plus_sampleinfo_dep_sem_outliers_mean_concentracao <- semi_join(nanosight_plus_sampleinfo_dep_sem_outliers_mean, nanosight_plus_sampleinfo_dep_sem_outliers_concentracao, by = "id_sample")
nanosight_plus_sampleinfo_dep_sem_outliers_porcentagem_concentracao <- semi_join(nanosight_plus_sampleinfo_dep_sem_outliers_porcentagem, nanosight_plus_sampleinfo_dep_sem_outliers_concentracao, by = ("id_sample"))
nanosight_plus_sampleinfo_dep_sem_outliers_mean_porcentagem <- semi_join(nanosight_plus_sampleinfo_dep_sem_outliers_mean, nanosight_plus_sampleinfo_dep_sem_outliers_porcentagem, by = "id_sample")
nanosight_intersect_dep <- semi_join(nanosight_plus_sampleinfo_dep_sem_outliers_mean_concentracao, nanosight_plus_sampleinfo_dep_sem_outliers_porcentagem, by = ("id_sample"))
nanosight_outliers_dep <- anti_join(nanosight_plus_sampleinfo_dep, nanosight_intersect_dep, by = "id_sample")

##Tabela só com w1 e só w2
nanosight_w1_dep <- subset(nanosight_intersect_dep, wave == "t1")
nanosight_w2_dep <- subset(nanosight_intersect_dep, wave == "t2")
nanosight_intersect_pares_dep <- semi_join(nanosight_w1_dep, nanosight_w2_dep, by = ("subjectid")) #54 pares + 4w1 + 7w2 


