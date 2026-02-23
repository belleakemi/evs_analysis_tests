setwd ("C:/Users/Belle/Documents/Belle - Nanosight/dados_raw")
nanosight_intersect_dep <- read.csv("nanosight_intersect_dep.csv")

# --- 1. PREPARAÇÃO DO AMBIENTE ---
# Função para instalar pacotes faltantes automaticamente
pacotes_necessarios <- c("lme4", "lmerTest", "sjPlot", "performance", "see", "patchwork", "ggplot2")
novos_pacotes <- pacotes_necessarios[!(pacotes_necessarios %in% installed.packages()[,"Package"])]
if(length(novos_pacotes)) install.packages(novos_pacotes)

# Carregar bibliotecas
library(lme4)
library(lmerTest)
library(sjPlot)
library(performance)
library(see)
library(ggplot2)

# --- 2. CARREGAR E PREPARAR DADOS ---
# Certifique-se de que o objeto 'nanosight_intersect_dep' já está carregado no R.
# Se não estiver, descomente a linha abaixo e coloque o caminho do seu arquivo:
# nanosight_intersect_dep <- read.csv("seu_arquivo.csv")

# Garantir que fatores estão corretos
nanosight_intersect_dep$sex <- as.factor(nanosight_intersect_dep$sex)
nanosight_intersect_dep$site <- as.factor(nanosight_intersect_dep$site)
nanosight_intersect_dep$subjectid <- as.factor(nanosight_intersect_dep$subjectid)

# --- 3. RODAR OS MODELOS OTIMIZADOS ---

# MODELO A: Tamanho Médio (Size)
# Decisão: GLMM Gamma (Log Link) devido à assimetria positiva dos dados
model_size_final <- glmer(tamanho_mean_average ~ sex + site + bage + (1 | subjectid), 
                          family = Gamma(link = "log"), 
                          data = nanosight_intersect_dep,
                          control = glmerControl(optimizer = "bobyqa", 
                                                 optCtrl = list(maxfun = 2e5))) #control porque deu problema de convergencia

# MODELO B: Concentração Real (Concentration)
# Decisão: LMM com Log10 para corrigir ordem de magnitude e normalizar resíduos
model_conc_final <- lmer(log10(concentracao_real) ~ sex + site + bage + (1 | subjectid), 
                         data = nanosight_intersect_dep)

# MODELO C: Porcentagem de EVs Pequenas (Small EVs %)
# Decisão: GLMM Gamma (Log Link) para lidar com a distribuição de porcentagem distorcida
model_peq_final  <- glmer(EV_pequenas_porcentagem ~ sex + site + bage + (1 | subjectid), 
                          family = Gamma(link = "log"), 
                          data = nanosight_intersect_dep)

# --- 4. GERAR PROVAS DAS PREMISSAS (DIAGNÓSTICOS) ---
message("Gerando diagnósticos visuais em PDF... (Isso evita erros de janela)")

# Salva os gráficos direto em PDF de alta resolução (Suplementar)
pdf("Diagnosticos_Visuais_Suplementar_dep.pdf", width = 10, height = 14)

# Página 1: Diagnóstico de Tamanho
p1 <- check_model(model_size_final)
plot(p1)
mtext("Fig S1. Diagnóstico do Modelo: Tamanho Médio (Gamma GLMM)", side = 3, line = 2, cex = 1.5)

# Página 2: Diagnóstico de Concentração
p2 <- check_model(model_conc_final)
plot(p2)
mtext("Fig S2. Diagnóstico do Modelo: Concentração (Log-LMM)", side = 3, line = 2, cex = 1.5)

# Página 3: Diagnóstico de Porcentagem
p3 <- check_model(model_peq_final)
plot(p3)
mtext("Fig S3. Diagnóstico do Modelo: % EVs Pequenas (Gamma GLMM)", side = 3, line = 2, cex = 1.5)

dev.off() # Fecha o PDF

# Salva os testes numéricos em TXT para consulta
capture.output({
  cat("=== RELATÓRIO DE VALIDAÇÃO ESTATÍSTICA ===\n\n")
  
  cat("--- 1. MODELO DE TAMANHO (GAMMA) ---\n")
  print(check_collinearity(model_size_final))
  print(model_performance(model_size_final))
  cat("\n")
  
  cat("--- 2. MODELO DE CONCENTRAÇÃO (LOG-NORMAL) ---\n")
  print(check_normality(model_conc_final)) # Deve passar (p > 0.05) ou chegar perto
  print(check_heteroscedasticity(model_conc_final))
  print(check_collinearity(model_conc_final))
  cat("\n")
  
  cat("--- 3. MODELO DE PORCENTAGEM (GAMMA) ---\n")
  print(check_collinearity(model_peq_final))
  print(model_performance(model_peq_final))
  
}, file = "Relatorio_Estatistico_Premissas.txt")

# --- 5. GERAR TABELA FINAL DE PUBLICAÇÃO ---
message("Gerando tabela final em Word...")

# Rótulos técnicos em Inglês
labels_preditores <- c(
  `(Intercept)` = "Intercept",
  `sexF` = "Sex (Female)",
  `siteSP` = "Site (São Paulo)",
  `age` = "Age (Years)"
)

# Tabela Vertical Condensada (Clean)
tab_model(
  model_size_final, model_conc_final, model_peq_final,
  
  dv.labels = c("Mean Size (nm)", "Conc. (log10)", "Small EVs (%)"),
  pred.labels = labels_preditores,
  
  # Configuração Vertical/Condensada
  collapse.ci = TRUE,      # Junta Estimate e CI
  show.p = TRUE,           # Mostra P-valor
  p.style = "scientific",  # Formato científico para p muito pequeno
  digits = 3,              # 3 casas decimais
  linebreak = TRUE,        # Quebra linha para economizar espaço
  
  title = "Table 1. Longitudinal mixed-effects analysis of demographic and procedural effects on EV characteristics",
  file = "Tabela_Final_Resultados.doc"
)

message("PROCESSO CONCLUÍDO COM SUCESSO!")
message("Verifique na sua pasta os arquivos: Tabela_Final_Resultados.doc, Diagnosticos.pdf e Relatorio.txt")
