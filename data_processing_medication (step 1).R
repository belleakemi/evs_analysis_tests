setwd ("C:/Users/Belle/Documents/Belle - Nanosight/dados_raw")
bhrc_treatment_subjectid <- read.csv("bhrc_treatment_subjectid.csv")
nanosight_intersect_dep <- read.csv("nanosight_intersect_dep.csv")

##Selecionar as colunas da nanosight_intersect_dep que importam pq estava dando erro no join
nanosight_intersect_dep_reduzido <- nanosight_intersect_dep %>%
  select(id_sample, concentracao_real, tamanho_mean_average, EV_pequenas_porcentagem, subjectid, wave, sex, bage, site, Trajetoria, dcany, dcmadep, dcgena, dcanyhk, grupo_analise_dep)

library(dplyr)

#Juntando as tabelas
bhrc_treatment_subjectid <- bhrc_treatment_subjectid %>%
  mutate(redcap_event_name = case_when(
    redcap_event_name == "wave0_arm_1" ~ "t0",
    redcap_event_name == "wave1_arm_1" ~ "t1",
    redcap_event_name == "wave2_arm_1" ~ "t2"
  ))
bhrc_treatment_subjectid <- bhrc_treatment_subjectid %>% rename(wave = redcap_event_name) # Renomear 'antigo' para 'novo'
nanosight_treatment_subjectid <- left_join(nanosight_intersect_dep_reduzido, bhrc_treatment_subjectid, by = c("subjectid", "wave"))

#Filtrar somente os pacientes que fizeram uso de medicaçao durante a vida
nanosight_treatment_filtrada <- nanosight_treatment_subjectid %>% filter(med_lifetime == 1)
