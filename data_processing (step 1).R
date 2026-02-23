#Importar bibliotecas e scripts necessarios
library(tidyr)
library(dplyr)
library(readr)
library(broom)

#Definir pasta de trabalho e importar scripts necessários
setwd("C:/Users/Belle/Documents/Belle - Nanosight")
source("gera_pdf_qq_csv_norm.R")
source("analise_nanosight_utils.R")

#Trazer tabelas de interesse para o R
meuwd <- setwd ("C:/Users/Belle/Documents/Belle - Nanosight/dados_raw")
isolar_exp_sum <- list.files(pattern="*ExperimentSummary.csv") 
importar_exp_sum <- lapply(isolar_exp_sum, read.csv, header = FALSE, col.names = c(paste0("V", 1:6))) #importa todas as expsum corretamente ao nomear as colunas

#Extrair dados de interesse de ExpSum
##Para ter o dado do ID
isolar_exp_sum2 <- gsub (" W1", "_w1", isolar_exp_sum)
isolar_exp_sum2 <- gsub (" W2", "_w2", isolar_exp_sum2)
id_sample <- separate(as.data.frame(isolar_exp_sum2), 1, into=c("ids", "resto"), sep=" ") [,1]
id_sample <- extract_standardized_id(id_sample) #padroniza as ids para como o inpd usa
##Para ter o dado do quanto diluiu (assumindo que está no campo diluent, e não dilution factor, como a Jess fez)
diluicao <- as.data.frame(sapply(importar_exp_sum, index_lapply_sapply, 11, 2)) #se estiver no campo diluent (dados Jessica)
diluicao[diluicao == "" | is.na(diluicao)] <- deparse(quote(1/100)) #se nao estiver em campo nenhum e n tiver precisado diluir mais (dados padrao Belle)
diluicao[diluicao == "PBS"] <- "1/100"
diluicao[diluicao == "10X diluido" | diluicao == "10x diluida" | diluicao == "10x" | diluicao == "900ul PBS + 100ul amostra"] <- "1/1000"
##Para ter a concentração da amostra antes de diluir (concentração real)
concentracao_average <- as.numeric(sapply(importar_exp_sum, index_lapply_sapply, 44, 5)) #obter média da concentração
separar_fracao_diluicao <- separate(diluicao, 1, into=c("numerador_diluicao", "denominador_diluicao"), sep="/") #gera uma tabela com o numerador em uma coluna e denominador em outra para que seja possível ter um dado do tipo numeric a partir da fração de diluição anotada (usando só as.numeric, retorna NA, pq a barra confunde o R)
separar_fracao_diluicao <- separar_fracao_diluicao %>%
  mutate(numerador_diluicao = ifelse(row_number() == 175, "1", numerador_diluicao)) %>%
  mutate(denominador_diluicao = ifelse(row_number() == 175, "10", denominador_diluicao)) #uma amostra que estava desconfigurada
numerador <- as.numeric(separar_fracao_diluicao[,1]) #isolar o numerador da tabela gerada para poder transformar em numeric
denominador <- as.numeric(separar_fracao_diluicao[,2]) #isolar o denominador da tabela gerada para poder transformar em numeric
concentracao_real <- concentracao_average/(numerador/denominador) #obtém de fato a concentração real
##Para obter o restante dos dados
tamanho_mean_average <- as.numeric(sapply(importar_exp_sum, index_lapply_sapply, 73, 5)) #obter average da mean do tamanho
tamanho_mean_average_alt <- as.numeric(sapply(importar_exp_sum, index_lapply_sapply, 77, 5)) #obter average da mean do tamanho
na_indexes_tamanho_mean <- which(is.na(tamanho_mean_average))
tamanho_mean_average_alt <- tamanho_mean_average_alt[na_indexes_tamanho_mean]
tamanho_mean_average[is.na(tamanho_mean_average)] <- tamanho_mean_average_alt
tamanho_mode_average <- as.numeric(sapply(importar_exp_sum, index_lapply_sapply, 74,5)) #obter average da moda do tamanho
tamanho_mode_average_alt <- as.numeric(sapply(importar_exp_sum, index_lapply_sapply, 78,5)) #obter average da moda do tamanho
na_indexes_tamanho_mode <- which(is.na(tamanho_mode_average))
tamanho_mode_average_alt <- tamanho_mode_average_alt[na_indexes_tamanho_mode]
tamanho_mode_average[is.na(tamanho_mode_average)] <- tamanho_mode_average_alt
##Porcentagem de vesículas pequenas e grandes (corte 128.5 nm)
concentration_EV_pequenas <- as.numeric(sapply(importar_exp_sum, function(x) {sum(as.numeric(x[85:213, 5]))}))
concentration_EV_pequenas_alt <- as.numeric(sapply(importar_exp_sum, function(x) {sum(as.numeric(x[90:218, 5]))}))
na_indexes_concentration_EV_pequenas <- which(is.na(concentration_EV_pequenas))
concentration_EV_pequenas_alt <- concentration_EV_pequenas_alt[na_indexes_concentration_EV_pequenas]
concentration_EV_pequenas[is.na(concentration_EV_pequenas)] <- concentration_EV_pequenas_alt
concentration_EV_grandes <- as.numeric(sapply(importar_exp_sum, function(x) {sum(as.numeric(x[214:1081, 5]))}))
concentration_EV_grandes_alt <- as.numeric(sapply(importar_exp_sum, function(x) {sum(as.numeric(x[219:1086, 5]))}))
na_indexes_concentration_EV_grandes <- which(is.na(concentration_EV_grandes))
concentration_EV_grandes_alt <- concentration_EV_grandes_alt[na_indexes_concentration_EV_grandes]
concentration_EV_grandes[is.na(concentration_EV_grandes)] <- concentration_EV_grandes_alt
EV_pequenas_porcentagem <- (concentration_EV_pequenas/(concentration_EV_grandes + concentration_EV_pequenas)) * 100
##Batch Nanosight
batch_nan <- sapply(importar_exp_sum, index_lapply_sapply, 8, 2)
batch_nan <- as.data.frame(batch_nan)
batch_nan <- sapply(strsplit(batch_nan$batch_nan, " "), tail, 1)
batch_nan <- as.data.frame(batch_nan)
datas_unicas <- duplicated(batch_nan$batch_nan)
datas_unicas <- subset(batch_nan, !datas_unicas) #para ver quantos grupos de datas (batch)
datas_unicas <- datas_unicas[,1]
numeros_batch <- seq(1, 22)
replacement <- setNames(numeros_batch, datas_unicas)
batch_nan$batch_nan <- replacement[batch_nan$batch_nan]
batch_nan <- batch_nan %>%
  mutate(batch_nan = case_when(
    batch_nan %in% c(1, 13) ~ "B",  # Agrupar 1 e 13 para B
    batch_nan == 2 ~ "J",           # Agrupar 2 para J
    batch_nan %in% c(4, 5, 9, 14, 20, 21) ~ "I1",
    batch_nan %in% c(3, 6, 7, 8, 10, 11, 12, 15, 16, 17, 18, 19, 22) ~ "I2",
    TRUE ~ as.character(batch_nan)  # Manter os outros números como estão
  )) #Jess, Belle e ICESP
##Juntando todos os dados obtidos do nanosight em uma tabela
tabela_tudo_nanosight <- data.frame(id_sample, concentracao_average, diluicao, concentracao_real, tamanho_mean_average, EV_pequenas_porcentagem, batch_nan, tamanho_mode_average)
colnames(tabela_tudo_nanosight)[3] <- "diluicao"
duplicatas <- duplicated(tabela_tudo_nanosight[, c("id_sample")], fromLast = TRUE) #a 2a leitura é a correta
nanosight_sem_duplicatas <- subset(tabela_tudo_nanosight, !duplicatas)
##Juntar dados do nanosight com informações da amostra
sample_information <- read.csv("sampleinformation_liriel320_atualizada.csv") #chamando/lendo/importando infos da amostra
sample_information$id_sample <- paste(sample_information$subjectid, sample_information$wave, sep = "_")
nanosight_plus_sampleinfo <- inner_join(nanosight_sem_duplicatas, sample_information, by="id_sample")
nanosight_plus_sampleinfo$Grupo <- paste(nanosight_plus_sampleinfo$Trajetoria, nanosight_plus_sampleinfo$wave, sep = "_")
sem_correspondencia <- anti_join(nanosight_sem_duplicatas, sample_information, by = "id_sample")
write.csv(nanosight_plus_sampleinfo, "nanosight_plus_sampleinfo_all.csv", row.names = FALSE)
##Amostras sem análise no Nanosight
sample_info_sem_nanosight <- anti_join(sample_information, nanosight_sem_duplicatas, by="id_sample")

#Retirar coluna que nao interessa (bage)
nanosight_plus_sampleinfo <- nanosight_plus_sampleinfo[,-12]
#Para analise de Levene as colunas de transtornos não podem ser numeros, portanto (F=False T=True)
nanosight_plus_sampleinfo <- nanosight_plus_sampleinfo %>%
  mutate(
    dcany = ifelse(dcany == 0, "F", "T"),
    dcmadep = ifelse(dcmadep == 0, "F", "T"),
    dcanyanx = ifelse(dcanyanx == 0, "F", "T"),
    dcgena = ifelse(dcgena == 0, "F", "T"),
    dcanyhk = ifelse(dcanyhk == 0, "F", "T"),
    dcpsych = ifelse(dcpsych == 0, "F", "T"),
    dcmania = ifelse(dcmania == 0, "F", "T"),
    dcptsd = ifelse(dcptsd == 0, "F", "T")
  )
#Nomeando as trajetórias para deixar menos confuso
nanosight_plus_sampleinfo <- nanosight_plus_sampleinfo %>%
  mutate(Trajetoria = case_when(
    Trajetoria == "A" ~ "Control",
    Trajetoria == "B" ~ "Incident",
    Trajetoria == "C" ~ "Remitted",
    Trajetoria == "D" ~ "Persistent"
  ))
#Alterando wave para time point
nanosight_plus_sampleinfo <- nanosight_plus_sampleinfo %>%
  mutate(wave = case_when(
    wave == "w1" ~ "t1",
    wave == "w2" ~ "t2"
  ))
##Retirando outliers
nanosight_plus_sampleinfo$zscore_mean <- scale(nanosight_plus_sampleinfo$tamanho_mean_average)
nanosight_plus_sampleinfo$zscore_porcentagem <- scale(nanosight_plus_sampleinfo$EV_pequenas_porcentagem)
nanosight_plus_sampleinfo$zscore_concentracao <- scale(nanosight_plus_sampleinfo$concentracao_real)
nanosight_plus_sampleinfo$zscore_conc_nan <- scale(nanosight_plus_sampleinfo$concentracao_average)
nanosight_plus_sampleinfo_sem_outliers_mean <- nanosight_plus_sampleinfo[nanosight_plus_sampleinfo$zscore_mean >= -3.0 & nanosight_plus_sampleinfo$zscore_mean <= 3.0, ]
nanosight_plus_sampleinfo_sem_outliers_porcentagem <- nanosight_plus_sampleinfo[nanosight_plus_sampleinfo$zscore_porcentagem >= -3.0 & nanosight_plus_sampleinfo$zscore_porcentagem <= 3.0, ]
nanosight_plus_sampleinfo_sem_outliers_concentracao <- nanosight_plus_sampleinfo[nanosight_plus_sampleinfo$zscore_concentracao >= -3.0 & nanosight_plus_sampleinfo$zscore_concentracao <= 3.0, ]
nanosight_plus_sampleinfo_sem_outliers_conc_nan <- nanosight_plus_sampleinfo[nanosight_plus_sampleinfo$zscore_conc_nan >= -3.0 & nanosight_plus_sampleinfo$zscore_conc_nan <= 3.0, ]
nanosight_plus_sampleinfo_sem_outliers_mean_concentracao <- semi_join(nanosight_plus_sampleinfo_sem_outliers_mean, nanosight_plus_sampleinfo_sem_outliers_concentracao, by = "id_sample")
nanosight_plus_sampleinfo_sem_outliers_porcentagem_concentracao <- semi_join(nanosight_plus_sampleinfo_sem_outliers_porcentagem, nanosight_plus_sampleinfo_sem_outliers_concentracao, by = ("id_sample"))
nanosight_plus_sampleinfo_sem_outliers_mean_porcentagem <- semi_join(nanosight_plus_sampleinfo_sem_outliers_mean, nanosight_plus_sampleinfo_sem_outliers_porcentagem, by = "id_sample")
nanosight_intersect <- semi_join(nanosight_plus_sampleinfo_sem_outliers_mean_concentracao, nanosight_plus_sampleinfo_sem_outliers_porcentagem, by = ("id_sample"))
nanosight_intersect_conc_nan <- semi_join(nanosight_plus_sampleinfo_sem_outliers_mean_porcentagem, nanosight_plus_sampleinfo_sem_outliers_conc_nan, by = ("id_sample"))

##Tabela só com w1 e só w2
nanosight_w1 <- subset(nanosight_plus_sampleinfo, wave == "t1")
nanosight_w2 <- subset(nanosight_plus_sampleinfo, wave == "t2")
nanosight_w1_semoutliers_mean <- subset(nanosight_plus_sampleinfo_sem_outliers_mean, wave == "t1")
nanosight_w2_semoutliers_mean <- subset(nanosight_plus_sampleinfo_sem_outliers_mean, wave == "t2")
nanosight_intersect_w1 <- subset(nanosight_intersect, wave == "t1")
nanosight_intersect_w2 <- subset(nanosight_intersect, wave == "t2")
nanosight_intersect_conc_nan_w1 <- subset(nanosight_intersect_conc_nan, wave == "t1")
nanosight_intersect_conc_nan_w2 <- subset(nanosight_intersect_conc_nan, wave == "t2")

#Tabela só com amostras que tem pares na w1 e w2
nanosight_intersect_semoutliers_mean_w1 <- semi_join(nanosight_w1_semoutliers_mean, nanosight_w2_semoutliers_mean, by = ("subjectid"))
nanosight_intersect_semoutliers_mean_w2 <- semi_join(nanosight_w2_semoutliers_mean, nanosight_w1_semoutliers_mean, by = ("subjectid"))

#Tabela de cada batch_nan
nanosight_batch_b <- subset(nanosight_intersect_conc_nan, batch_nan == "B")
nanosight_batch_i1 <- subset(nanosight_intersect_conc_nan, batch_nan == "I1")
nanosight_batch_i2 <- subset(nanosight_intersect_conc_nan, batch_nan == "I2")
nanosight_batch_j <- subset(nanosight_intersect_conc_nan, batch_nan == "J")
nanosight_batch_b_w1 <- subset(nanosight_batch_b, wave == "t1")
nanosight_batch_b_w2 <- subset(nanosight_batch_b, wave == "t2")
nanosight_batch_i1_w1 <- subset(nanosight_batch_i1, wave == "t1")
nanosight_batch_i1_w2 <- subset(nanosight_batch_i1, wave == "t2")
nanosight_batch_i2_w1 <- subset(nanosight_batch_i2, wave == "t1")
nanosight_batch_i2_w2 <- subset(nanosight_batch_i2, wave == "t2")
nanosight_batch_j_w1 <- subset(nanosight_batch_j, wave == "t1")
nanosight_batch_j_w2 <- subset(nanosight_batch_j, wave == "t2")

#Tabela só com grupo Persistente
nanosight_intersect_conc_nan_persistence <- subset(nanosight_intersect_conc_nan, Trajetoria == "Persistent")
nanosight_intersect_persistence <- subset(nanosight_intersect, Trajetoria == "Persistent")

#Tabela sem o batch_nan_i1
nanosight_batch_b_i2_j <- subset(nanosight_plus_sampleinfo, batch_nan != "I1")

#tabela considerando o batch I como um todo
nanosight_batch_i <- subset(nanosight_plus_sampleinfo, batch_nan %in% c("I1", "I2"))

##Outliers
nanosight_outliers <- anti_join(nanosight_plus_sampleinfo, nanosight_intersect, by = "id_sample")

##Amostras sem pares
linhas_sem_pares <- !duplicated(nanosight_plus_sampleinfo$subjectid) & !duplicated(nanosight_plus_sampleinfo$subjectid, fromLast = TRUE)
nanosight_sem_pares <- nanosight_plus_sampleinfo[linhas_sem_pares, ] #perde 7
linhas_sem_pares_intersect <- !duplicated(nanosight_intersect$subjectid) & !duplicated(nanosight_intersect$subjectid, fromLast = TRUE)
nanosight_sem_pares_intersect <- nanosight_intersect[linhas_sem_pares_intersect, ] #perde 7

##Amostras com pares
nanosight_plus_sampleinfo_pares <- anti_join(nanosight_plus_sampleinfo, nanosight_sem_pares, by = "id_sample")
nanosight_intersect_pares <- anti_join (nanosight_intersect, nanosight_sem_pares_intersect, by = "id_sample")

#Tabela considerando somente pares sem outliers para tamanho
nanosight_pares_semoutliers_tamanho <- semi_join(nanosight_plus_sampleinfo_sem_outliers_mean, nanosight_plus_sampleinfo_pares,by = "id_sample")
nanosight_pares_semoutliers_porcentagem <- semi_join(nanosight_plus_sampleinfo_sem_outliers_porcentagem, nanosight_plus_sampleinfo_pares,by = "id_sample")
nanosight_pares_semoutliers_concentracao <- semi_join(nanosight_plus_sampleinfo_sem_outliers_concentracao, nanosight_plus_sampleinfo_pares,by = "id_sample")

#Contagem de indivíduos com diferentes transtornos:
##w1
variaveis_w1 <- nanosight_w1 %>% 
  select(dcany, dcmadep, dcanyanx, dcgena, dcanyhk, dcpsych, dcptsd)
variaveis_w1_combinacao <- as.data.frame(rowSums(variaveis_w1 == "T"))
combinacoes_w1_frequencias <- variaveis_w1 %>%
  group_by(across(everything())) %>%
  summarise(Frequencia = n(), .groups = "drop")
##w2
variaveis_w2 <- nanosight_w2 %>% 
  select(dcany, dcmadep, dcanyanx, dcgena, dcanyhk, dcpsych, dcptsd)
variaveis_w2_combinacao <- as.data.frame(rowSums(variaveis_w2 == "T"))
combinacoes_w2_frequencias <- variaveis_w2 %>%
  group_by(across(everything())) %>%
  summarise(Frequencia = n(), .groups = "drop")
##Incident
nanosight_intersect_Incident <- subset(nanosight_intersect, Trajetoria == "Incident")
variaveis_incident <- nanosight_intersect %>% 
  select(dcany, dcmadep, dcanyanx, dcgena, dcanyhk, dcpsych, dcptsd)
variaveis_incident_combinacao <- as.data.frame(rowSums(variaveis_incident == "T"))
combinacoes_incident_frequencias <- variaveis_incident %>%
  group_by(across(everything())) %>%
  summarise(Frequencia = n(), .groups = "drop")

