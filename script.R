# ============================================================
# 1. PREPARACAO DO AMBIENTE E IMPORTACAO DOS DADOS
# ============================================================

# 1.1. Pacotes principais
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(BiocManager)

# 1.2. Verificar compatibilidade entre os pacotes
if (packageVersion("GeomxTools") < "2.1" &
    packageVersion("GeoMxWorkflows") >= "1.0.1") {
  stop("GeomxTools e GeoMxWorkflows possuem versoes incompatíveis.")
}

if (packageVersion("GeomxTools") > "2.1" &
    packageVersion("GeoMxWorkflows") <= "1.0.1") {
  stop("GeomxTools e GeoMxWorkflows possuem versoes incompatíveis.")
}

# 1.3. Registrar as versoes utilizadas
packageVersion("NanoStringNCTools")
packageVersion("GeomxTools")
packageVersion("GeoMxWorkflows")
BiocManager::version()
R.version.string

# 1.4. Definir diretorio de trabalho
datadir <- "C:/Users/isabe/Downloads/burns_teste_certo"

# 1.5. Criar diretorio para resultados
results_dir <- file.path(datadir, "analysis_results")

dir.create(
  results_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# 1.6. Localizar arquivos DCC e PKC
DCCFiles <- list.files(
  datadir,
  pattern = "\\.dcc$",
  full.names = TRUE
)

PKCFiles <- list.files(
  datadir,
  pattern = "\\.pkc$",
  full.names = TRUE
)

# Conferir arquivos encontrados
length(DCCFiles)
length(PKCFiles)

# 1.7. Arquivo de anotacao das amostras
SampleAnnotationFile <- file.path(
  datadir,
  "annotation_final_GSE263897.xlsx"
)

stopifnot(file.exists(SampleAnnotationFile))

# 1.8. Importar os dados GeoMx
demoData <- readNanoStringGeoMxSet(
  dccFiles = DCCFiles,
  pkcFiles = PKCFiles,
  phenoDataFile = SampleAnnotationFile,
  phenoDataSheet = "Template",
  phenoDataDccColName = "DCC"
)

# 1.9. Inspecionar objeto importado
demoData
dim(demoData)

head(pData(demoData))
head(fData(demoData))


# ============================================================
# 2. INSPECAO DOS MODULOS E ESTRUTURA DOS DADOS
# ============================================================

library(knitr)
library(dplyr)
library(networkD3)

# 2.1. Conferir os PKCs associados ao objeto
pkcs <- annotation(demoData)

modules <- gsub(
  "\\.pkc$",
  "",
  basename(pkcs)
)

kable(
  data.frame(
    PKCs = pkcs,
    Modules = modules
  )
)

# 2.2. Conferir as variaveis de anotacao
colnames(pData(demoData))

# Conferir os niveis das principais variaveis
table(pData(demoData)$sampleID, useNA = "ifany")
table(pData(demoData)$Tissue, useNA = "ifany")
table(pData(demoData)$Slide, useNA = "ifany")
table(pData(demoData)$Cell_Type, useNA = "ifany")

# 2.3. Diagrama Sankey da estrutura amostral
# Paciente -> Tissue -> Slide -> Cell Type

sankeyCols <- c("source", "target", "value")

link1 <- count(pData(demoData), sampleID, Tissue)
link2 <- count(pData(demoData), Tissue, Slide)
link3 <- count(pData(demoData), Slide, Cell_Type)

colnames(link1) <- sankeyCols
colnames(link2) <- sankeyCols
colnames(link3) <- sankeyCols

links <- rbind(link1, link2, link3)

nodes <- data.frame(
  name = unique(c(links$source, links$target))
)

links$source <- as.integer(
  match(links$source, nodes$name) - 1
)

links$target <- as.integer(
  match(links$target, nodes$name) - 1
)

sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",
  units = "AOIs",
  fontSize = 13,
  nodeWidth = 25,
  nodePadding = 15
)


# ============================================================
# 3. SEGMENT QC - CONTROLE DE QUALIDADE DOS SEGMENTOS
# ============================================================

library(ggplot2)

# 3.1. Ajustar a estrutura das contagens
demoData <- shiftCountsOne(
  demoData,
  useDALogic = TRUE
)

# 3.2. Definir parametros de QC
QC_params <- list(
  minSegmentReads = 1000,
  percentTrimmed = 80,
  percentStitched = 80,
  percentAligned = 75,
  percentSaturation = 40,
  minNegativeCount = 1,
  maxNTCCount = 9000,
  minNuclei = 20,
  minArea = 1000
)

# 3.3. Calcular flags de QC
demoData <- setSegmentQCFlags(
  demoData,
  qcCutoffs = QC_params
)

QCResults <- protocolData(demoData)[["QCFlags"]]

# 3.4. Resumir flags de qualidade
flag_columns <- colnames(QCResults)

QC_Summary <- data.frame(
  Pass = colSums(!QCResults[, flag_columns, drop = FALSE]),
  Warning = colSums(QCResults[, flag_columns, drop = FALSE])
)

# Classificar cada segmento
QCResults$QCStatus <- ifelse(
  rowSums(QCResults[, flag_columns, drop = FALSE]) == 0,
  "PASS",
  "WARNING"
)

# Resumo do numero de segmentos por status
QC_Status_Summary <- table(QCResults$QCStatus)

print(QC_Summary)
print(QC_Status_Summary)

# 3.5. Funcao para visualizar distribuicoes dos parametros de QC
QC_histogram <- function(
    assay_data = NULL,
    annotation = NULL,
    fill_by = NULL,
    thr = NULL,
    scale_trans = NULL
) {
  
  plt <- ggplot(
    assay_data,
    aes_string(
      x = paste0("unlist(`", annotation, "`)"),
      fill = fill_by
    )
  ) +
    geom_histogram(bins = 50) +
    geom_vline(
      xintercept = thr,
      linetype = "dashed",
      color = "black"
    ) +
    theme_bw() +
    guides(fill = "none") +
    facet_wrap(
      as.formula(paste("~", fill_by)),
      nrow = 4
    ) +
    scale_fill_manual(
      values = c(
        "Epithelium" = "#4A0E4E",
        "Macrophages" = "#FF66B2",
        "Stroma" = "#00F5D4"
      )
    ) +
    labs(
      x = annotation,
      y = "Segments, #",
      title = annotation
    )
  
  if (!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  
  return(plt)
}
# 3.6. Visualizar parametros de sequenciamento
col_by <- "Cell_Type"

QC_histogram(
  sData(demoData), "Trimmed (%)", col_by, 80
)

QC_histogram(
  sData(demoData), "Stitched (%)", col_by, 80
)

QC_histogram(
  sData(demoData), "Aligned (%)", col_by, 75
)

QC_histogram(
  sData(demoData), "Saturated (%)", col_by, 40
) +
  labs(
    title = "Sequencing Saturation (%)",
    x = "Sequencing Saturation (%)"
  )

# 3.7. Calcular medias geometricas dos controles negativos
negativeGeoMeans <- esBy(
  negativeControlSubset(demoData),
  GROUP = "Module",
  FUN = function(x) {
    assayDataApply(
      x,
      MARGIN = 2,
      FUN = ngeoMean,
      elt = "exprs"
    )
  }
)

protocolData(demoData)[["NegGeoMean"]] <- negativeGeoMeans

# Calcular desvios geometricos dos controles negativos.
# Necessario para a estimativa de LOQ utilizada adiante.
negativeGeoSDs <- esBy(
  negativeControlSubset(demoData),
  GROUP = "Module",
  FUN = function(x) {
    assayDataApply(
      x,
      MARGIN = 2,
      FUN = ngeoSD,
      elt = "exprs"
    )
  }
)

protocolData(demoData)[["NegGeoSD"]] <- negativeGeoSDs

# 3.8. Preparar anotacoes para visualizar os controles negativos
negCols <- paste0("NegGeoMean_", modules)

pData(demoData)[, negCols] <-
  protocolData(demoData)[["NegGeoMean"]]

for (ann in negCols) {
  print(
    QC_histogram(
      pData(demoData),
      ann,
      col_by,
      2,
      scale_trans = "log10"
    )
  )
}

# Remover colunas auxiliares antes da agregacao
pData(demoData) <- pData(demoData)[
  ,
  !colnames(pData(demoData)) %in% negCols,
  drop = FALSE
]

# 3.9. Tabela resumo de QC
kable(
  QC_Summary,
  caption = "QC Flags for Each Segment"
)

# Observacao:
# O arquivo de anotacao atual nao possui a variavel NTC.
# Executar esta verificacao somente se NTC existir.
if ("NTC" %in% colnames(sData(demoData))) {
  print(
    kable(
      table(NTC_Count = sData(demoData)$NTC),
      col.names = c("NTC Count", "# of Segments")
    )
  )
}

# 3.10. Filtrar segmentos aprovados
keep_segments <- QCResults$QCStatus == "PASS"

demoData <- demoData[, keep_segments]

# Conferir dimensoes apos o filtro
dim(demoData)
table(pData(demoData)$Cell_Type, useNA = "ifany")
table(pData(demoData)$Tissue, useNA = "ifany")


# ============================================================
# 4. PROBE QC - CONTROLE DE QUALIDADE DAS SONDAS
# ============================================================

# 4.1. Calcular flags de QC das sondas
demoData <- setBioProbeQCFlags(
  demoData,
  qcCutoffs = list(
    minProbeRatio = 0.1,
    percentFailGrubbs = 20
  ),
  removeLocalOutliers = TRUE
)

ProbeQCResults <- fData(demoData)[["QCFlags"]]

# 4.2. Resumir resultados do Probe QC
qc_df <- data.frame(
  Passed = sum(
    rowSums(ProbeQCResults[, -1, drop = FALSE]) == 0
  ),
  Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
  Local = sum(
    rowSums(ProbeQCResults[, -2:-1, drop = FALSE]) > 0 &
      !ProbeQCResults$GlobalGrubbsOutlier
  )
)

print(qc_df)

# 4.3. Excluir sondas com baixa razao ou outliers globais
ProbeQCPassed <- subset(
  demoData,
  fData(demoData)[["QCFlags"]][, "LowProbeRatio"] == FALSE &
    fData(demoData)[["QCFlags"]][, "GlobalGrubbsOutlier"] == FALSE
)

demoData <- ProbeQCPassed

# Conferir dimensoes
dim(demoData)

# Numero de alvos unicos antes da agregacao
length(unique(fData(demoData)[["TargetName"]]))

# 4.4. Agregar contagens por alvo molecular
target_demoData <- aggregateCounts(demoData)

# Conferir dimensoes apos agregacao
dim(target_demoData)

# Inspecionar as primeiras contagens
exprs(target_demoData)[1:5, 1:2]

grep(
  "^NegGeo(Mean|SD)_",
  colnames(pData(target_demoData)),
  value = TRUE
)
# ============================================================
# 5. LOQ E DETECCAO DE GENES
# ============================================================

library(scales)

# 5.1. Definir parametros de LOQ
cutoff <- 2
minLOQ <- 2

# 5.2. Calcular LOQ por modulo e segmento
LOQ <- data.frame(
  row.names = colnames(target_demoData)
)

for (module in modules) {
  
  vars <- paste0(
    c("NegGeoMean_", "NegGeoSD_"),
    module
  )
  
  if (all(vars %in% colnames(pData(target_demoData)))) {
    
    LOQ[, module] <- pmax(
      minLOQ,
      pData(target_demoData)[, vars[1]] *
        pData(target_demoData)[, vars[2]]^cutoff
    )
    
  } else {
    
    warning(
      "Dados de controles negativos ausentes para o modulo: ",
      module
    )
  }
}

# Conferir modulos para os quais o LOQ foi calculado
colnames(LOQ)
head(LOQ)

# 5.3. Armazenar LOQ no objeto
pData(target_demoData)$LOQ <- LOQ

# 5.4. Criar matriz logica de deteccao acima do LOQ
# Linhas = genes; colunas = segmentos

LOQ_Mat <- NULL

for (module in modules) {
  
  ind <- fData(target_demoData)$Module == module
  
  if (!any(ind)) {
    next
  }
  
  Mat_i <- t(
    esApply(
      target_demoData[ind, ],
      MARGIN = 1,
      FUN = function(x) {
        x > LOQ[, module]
      }
    )
  )
  
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}

# 5.5. Calcular número e proporção de genes endógenos detectados por segmento

# LOQ_Mat pode estar ordenada por módulo; alinhar CodeClass pelo nome da linha.
code_class_by_row <- setNames(
  as.character(fData(target_demoData)$CodeClass),
  rownames(fData(target_demoData))
)

if (is.null(rownames(LOQ_Mat))) {
  stop("LOQ_Mat precisa ter nomes de genes nas linhas para alinhar CodeClass.")
}

target_loq_rows <- code_class_by_row[rownames(LOQ_Mat)] == "Endogenous"
target_loq_rows[is.na(target_loq_rows)] <- FALSE

if (!any(target_loq_rows)) {
  stop(
    "Nenhuma linha Endogenous foi identificada. ",
    "Confira os valores de fData(target_demoData)$CodeClass."
  )
}

LOQ_Mat_targets <- LOQ_Mat[target_loq_rows, , drop = FALSE]

pData(target_demoData)$GenesDetected <-
  colSums(LOQ_Mat_targets, na.rm = TRUE)

pData(target_demoData)$GeneDetectionRate <-
  pData(target_demoData)$GenesDetected /
  nrow(LOQ_Mat_targets)

# 5.6. Categorizar segmentos pela taxa de deteccao
pData(target_demoData)$DetectionThreshold <- cut(
  pData(target_demoData)$GeneDetectionRate,
  breaks = c(0, 0.01, 0.05, 0.10, 0.15, 1),
  labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"),
  include.lowest = TRUE
)

# 5.7. Visualizar distribuicao da taxa de deteccao

ggplot(
  pData(target_demoData),
  aes(x = DetectionThreshold)
) +
  geom_bar(
    aes(fill = Cell_Type),
    width = 0.9
  ) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    vjust = -0.45,
    size = 5
  ) +
  theme_bw(base_size = 14) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  scale_fill_manual(
    values = c(
      "Epithelium" = "#4A0E4E",
      "Macrophages" = "#FF66B2",
      "Stroma" = "#00F5D4"
    )
  ) +
  labs(
    x = "Gene Detection Rate",
    y = "Segments, #",
    fill = "Segment Type"
  ) +
  theme(
    axis.title.x = element_text(size = 17),
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    plot.margin = margin(t = 12, r = 12, b = 12, l = 12)
  )

kable(
  table(
    pData(target_demoData)$DetectionThreshold,
    pData(target_demoData)$Cell_Type
  )
)
# 5.8. Filtrar segmentos com pelo menos 10% dos genes detectados
keep_detection <- pData(target_demoData)$GeneDetectionRate >= 0.10

target_demoData <- target_demoData[, keep_detection]

# Atualizar matriz LOQ para os segmentos mantidos
LOQ_Mat <- LOQ_Mat[
  ,
  colnames(target_demoData),
  drop = FALSE
]

dim(target_demoData)

# 5.9. Calcular taxa de detecção por gene e alinhar pelo nome

detected_by_gene <- rowSums(LOQ_Mat, na.rm = TRUE)
gene_index <- match(
  rownames(fData(target_demoData)),
  names(detected_by_gene)
)

if (anyNA(gene_index)) {
  stop("Não foi possível alinhar todos os genes de fData com LOQ_Mat.")
}

fData(target_demoData)$DetectedSegments <-
  detected_by_gene[gene_index]

fData(target_demoData)$DetectionRate <-
  fData(target_demoData)$DetectedSegments /
  ncol(target_demoData)
# 5.10. Definir genes de interesse biologico
goi <- c(
  # Sinalização Hormonal
  "ESR1", "ESR2", "PGR", "CYP19A1", "HSD17B2",
  
  # Inflamação e Quimiocinas
  "IL6", "IL1B", "TNF", "CCL2", "CXCL12", "CXCR4",
  
  # Angiogênese
  "VEGFA", "ANGPT2", "KDR",
  
  # Invasão e Remodelamento (MMPs / Cadherinas)
  "MMP2", "MMP9", "TIMP1", "TIMP2", "CDH1", "CDH2",
  
  # Fibrose, Sobrevivência e Imunidade
  "TGFB1", "BCL2", "PTPRC", "EPCAM"
)
# Verificar quais genes estao presentes
data.frame(
  Gene = goi,
  Present = goi %in% fData(target_demoData)$TargetName
)

# 5.11. Tabela de deteccao dos genes de interesse
# Pressupoe que TargetName identifica unicamente cada linha.
# Se nao for unico, usar uma correspondencia validada por ID.

goi_idx <- match(goi, fData(target_demoData)$TargetName)
goi_present <- !is.na(goi_idx)

goi_df <- data.frame(
  Gene = goi[goi_present],
  Number = fData(target_demoData)$DetectedSegments[
    goi_idx[goi_present]
  ],
  DetectionRate = percent(
    fData(target_demoData)$DetectionRate[
      goi_idx[goi_present]
    ]
  )
)

print(goi_df)

# 5.12. Proporcao de genes detectados em diferentes limiares
plot_detect <- data.frame(
  Freq = c(1, 5, 10, 20, 30, 50)
)

# Contar somente genes endógenos, excluindo controles.
is_endogenous <- as.character(
  fData(target_demoData)$CodeClass
) == "Endogenous"

target_detection_rates <-
  fData(target_demoData)$DetectionRate[is_endogenous]

plot_detect$Number <- vapply(
  c(0.01, 0.05, 0.10, 0.20, 0.30, 0.50),
  function(x) {
    sum(target_detection_rates >= x, na.rm = TRUE)
  },
  numeric(1)
)

plot_detect$Rate <-
  plot_detect$Number / length(target_detection_rates)
ggplot(
  plot_detect,
  aes(x = as.factor(Freq), y = Rate, fill = Rate)
) +
  geom_col(width = 0.9) +  # barras mais estreitas
  geom_text(
    aes(label = formatC(
      Number,
      format = "d",
      big.mark = ","
    )),
    vjust = -0.45,          # números acima das barras
    size = 5.5,              # números maiores
    color = "black"
  ) +
  scale_fill_gradient2(
    low = "#00F5D4",
    mid = "#FF66B2",
    high = "#4A0E4E",
    midpoint = 0.65,
    limits = c(0, 1),
    labels = scales::percent,
    name = "Rate"
  ) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = c(0, 0.06))
  ) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    plot.margin = margin(t = 12, r = 12, b = 12, l = 12)
  ) +
  labs(
    x = "% of Segments",
    y = "Genes Detected, % of Panel > LOQ"
  )
# 5.13. Manter genes endógenos detectados em pelo menos 10% dos segmentos.
# Controles são preservados no objeto para as etapas de diagnóstico.

is_endogenous <- as.character(
  fData(target_demoData)$CodeClass
) == "Endogenous"

is_control <- !is_endogenous
is_control[is.na(is_control)] <- TRUE

keep_genes <- is_control |
  (
    is_endogenous &
      fData(target_demoData)$DetectionRate >= 0.10
  )

target_demoData <- target_demoData[keep_genes, ]

dim(target_demoData)

# Atualizar genes de interesse presentes após o filtro.
goi <- goi[
  goi %in% fData(target_demoData)$TargetName
]

# ============================================================
# 6. DIAGNOSTICO E NORMALIZACAO
# ============================================================

library(reshape2)
library(cowplot)

# 6.1. Comparar Q3 e controles negativos
ann_of_interest <- "Cell_Type"

# Identificar as sondas negativas que foram mantidas no objeto.
is_negative_probe <- as.character(
  fData(target_demoData)$CodeClass
) == "Negative"

neg_probes <- rownames(fData(target_demoData))[is_negative_probe]

# Manter apenas sondas negativas presentes na matriz de expressão.
neg_probes <- intersect(
  neg_probes,
  rownames(target_demoData)
)

if (length(neg_probes) == 0L) {
  stop("Nenhuma sonda negativa foi encontrada no objeto.")
}

Stat_data <- data.frame(
  Segment = colnames(exprs(target_demoData)),
  Annotation = pData(target_demoData)[, ann_of_interest],
  Q3 = apply(
    exprs(target_demoData),
    2,
    quantile,
    probs = 0.75,
    na.rm = TRUE
  ),
  NegProbe = colMeans(
    exprs(target_demoData)[neg_probes, , drop = FALSE],
    na.rm = TRUE
  ),
  row.names = colnames(exprs(target_demoData))
)

Stat_data_m <- melt(
  Stat_data,
  id.vars = c("Segment", "Annotation"),
  measure.vars = c("Q3", "NegProbe"),
  variable.name = "Statistic",
  value.name = "Value"
)

# 6.2. Distribuições de Q3 e controles negativos

plt1 <- ggplot(
  Stat_data_m,
  aes(x = Value, fill = Statistic)
) +
  geom_histogram(bins = 40) +
  scale_x_continuous(trans = "log2") +
  facet_wrap(~Annotation, nrow = 1) +
  scale_fill_brewer(palette = 1, type = "qual") +
  labs(
    x = expression(Log[2] * "(Counts)"),
    y = "Segments, #"
  ) +
  theme_bw() +
  theme(
    # Mantém quadrada a área de plotagem de cada faceta
    aspect.ratio = 1,
    
    # Fundo branco e sem contorno nos títulos das facetas
    strip.background = element_rect(
      fill = "white",
      colour = NA
    ),
    strip.text = element_text(size = 12)
  )

print(plt1)

ggsave(
  file.path(results_dir, "Distributions_Q3_Negative_Controls.png"),
  plot = plt1,
  width = 12,
  height = 5,
  dpi = 300
)

# 6.3. Q3 versus media dos controles negativos
plt2 <- ggplot(
  Stat_data,
  aes(x = NegProbe, y = Q3, color = Annotation)
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "darkgray"
  ) +
  geom_point(size=2.5) +
  theme_bw() +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous(trans = "log2") +
  scale_color_manual(
    values = c(
      "Epithelium" = "#4A0E4E",
      "Macrophages" = "#FF66B2",
      "Stroma" = "#00F5D4"
    )
  ) +
  theme(aspect.ratio = 1) +
  labs(
    x = "Negative probe mean",
    y = "Q3"
  )

plt2

# 6.4. Razao entre Q3 e controles negativos
plt3 <- ggplot(
  Stat_data,
  aes(
    x = NegProbe,
    y = Q3 / NegProbe,
    color = Annotation
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "darkgray",
    linewidth = 1.2
  ) +
  geom_point(size = 3) +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous(trans = "log2", breaks = c(1, 2, 4, 8)) +
  scale_color_manual(
    values = c(
      "Epithelium" = "#4A0E4E",
      "Macrophages" = "#FF66B2",
      "Stroma" = "#00F5D4"
    )
  ) +
  labs(
    x = "Negative Probe GeoMean, Counts",
    y = "Q3 / Negative Probe Value, Counts",
    color = "Segment Type"
  ) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14)
  )

plt3

# ============================================================
# 6.5. Normalizacao por quantil Q3
# ============================================================

target_demoData <- normalize(
  target_demoData,
  norm_method = "quant",
  desiredQuantile = 0.75,
  toElt = "q_norm"
)


# ============================================================
# 6.6. Normalizacao pelo fundo dos controles negativos
# ============================================================

target_demoData <- normalize(
  target_demoData,
  norm_method = "neg",
  fromElt = "exprs",
  toElt = "neg_norm"
)


# ============================================================
# 6.7. Comparar contagens brutas e normalizadas
# ============================================================

n_show <- min(10, ncol(target_demoData))

# Raw Counts
boxplot(
  exprs(target_demoData)[, seq_len(n_show), drop = FALSE],
  col = "#9EDAE5",
  main = "Raw Counts",
  log = "y",
  names = seq_len(n_show),
  xlab = "Segment",
  ylab = "Counts, Raw"
)

# Q3 Normalized Counts
boxplot(
  assayDataElement(
    target_demoData[, seq_len(n_show)],
    elt = "q_norm"
  ),
  col = "#2CA02C",
  main = "Q3 Norm Counts",
  log = "y",
  names = seq_len(n_show),
  xlab = "Segment",
  ylab = "Counts, Q3 Normalized"
)

# Negative-Control Normalized Counts
boxplot(
  assayDataElement(
    target_demoData[, seq_len(n_show)],
    elt = "neg_norm"
  ),
  col = "#FF7F0E",
  main = "Neg Norm Counts",
  log = "y",
  names = seq_len(n_show),
  xlab = "Segment",
  ylab = "Counts, Neg. Normalized"
)



# ============================================================
# 7. ANALISE EXPLORATORIA DOS DADOS
# ============================================================

library(umap)
library(Rtsne)
library(pheatmap)

# 7.1. Definir matriz de expressao normalizada em log2
q_norm <- assayDataElement(
  target_demoData,
  elt = "q_norm"
)

log_q <- log2(q_norm)

# Conferir dimensoes
dim(log_q)

# 7.2. UMAP
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42

umap_out <- umap(
  t(log_q),
  config = custom_umap
)

pData(target_demoData)[, c("UMAP1", "UMAP2")] <-
  umap_out$layout[, 1:2]

ggplot(
  pData(target_demoData),
  aes(
    x = UMAP1,
    y = UMAP2,
    color = Cell_Type,
    shape = Tissue
  )
) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c(
      "Epithelium" = "#4A0E4E",
      "Macrophages" = "#FF66B2",
      "Stroma" = "#00F5D4"
    )
  ) +
  theme_bw() +
  labs(
    title = "UMAP - Q3 Normalized Expression",
    color = "Cell Type",
    shape = "Tissue"
  )

# 7.3. t-SNE
set.seed(42)

# Perplexity deve ser menor que o numero de segmentos.
# O valor deve ser ajustado se a amostra for pequena.
n_segments <- ncol(target_demoData)

if (n_segments > 3) {
  
  perplexity_value <- min(
    floor(n_segments * 0.15),
    floor((n_segments - 1) / 3)
  )
  
  perplexity_value <- max(1, perplexity_value)
  
  tsne_out <- Rtsne(
    t(log_q),
    perplexity = perplexity_value
  )
  
  pData(target_demoData)[, c("tSNE1", "tSNE2")] <-
    tsne_out$Y[, 1:2]
  
  print(
    ggplot(
      pData(target_demoData),
      aes(
        x = tSNE1,
        y = tSNE2,
        color = Cell_Type,
        shape = Tissue
      )
    ) +
      geom_point(size = 3) +
      scale_color_manual(
        values = c(
          "Epithelium" = "#4A0E4E",
          "Macrophages" = "#FF66B2",
          "Stroma" = "#00F5D4"
        )
      ) +
      theme_bw() +
      labs(
        title = "t-SNE - Q3 Normalized Expression",
        color = "Cell Type",
        shape = "Tissue"
      )
  )
}

# 7.4. Armazenar transformacao log2 no objeto
assayDataElement(
  object = target_demoData,
  elt = "log_q"
) <- assayDataApply(
  target_demoData,
  MARGIN = 2,
  FUN = log,
  base = 2,
  elt = "q_norm"
)

# 7.5. Calcular coeficiente de variacao dos genes
calc_CV <- function(x) {
  sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
}

CV_dat <- assayDataApply(
  target_demoData,
  elt = "log_q",
  MARGIN = 1,
  FUN = calc_CV
)

# Genes com maior CV
sort(CV_dat, decreasing = TRUE)[1:5]

# 7.6. Selecionar genes no quintil superior de variabilidade
GOI_variable <- names(CV_dat)[
  CV_dat > quantile(CV_dat, 0.80, na.rm = TRUE)
]


# 7.7. Heatmap dos 50 genes mais variáveis

# Selecionar os 50 genes com maior variância em log2(Q3).
# Limita o número de linhas para facilitar a leitura.
gene_var <- apply(log_q, 1, var, na.rm = TRUE)
gene_var <- gene_var[is.finite(gene_var)]

n_genes_heatmap <- min(50L, length(gene_var))
genes_heatmap <- names(sort(gene_var, decreasing = TRUE))[seq_len(n_genes_heatmap)]

mat <- assayDataElement(
  target_demoData[genes_heatmap, ],
  elt = "log_q"
)

mat <- as.matrix(mat)

# Remover genes com valores não finitos ou variância zero.
valid_genes <- apply(mat, 1, function(x) {
  all(is.finite(x)) && stats::var(x) > 0
})
mat <- mat[valid_genes, , drop = FALSE]

# Padronizar cada gene para mostrar expressão relativa entre AOIs.
mat_scaled <- t(scale(t(mat)))

# Limitar valores extremos para melhorar o contraste das cores.
mat_scaled[mat_scaled > 2] <- 2
mat_scaled[mat_scaled < -2] <- -2

# Manter apenas anotações relevantes e alinhar à ordem das colunas.
annotation <- as.data.frame(
  pData(target_demoData)[
    colnames(mat_scaled),
    c("Cell_Type", "Tissue"),
    drop = FALSE
  ]
)

# Garantir a ordem dos níveis nas anotações.
annotation$Cell_Type <- factor(
  annotation$Cell_Type,
  levels = c("Epithelium", "Stroma", "Macrophages")
)
annotation$Tissue <- factor(
  annotation$Tissue,
  levels = c("eutopic endometrium", "endometriotic lesion")
)

annotation_colors <- list(
  Cell_Type = c(
    "Epithelium" = "#54204E",
    "Stroma" = "#278C87",
    "Macrophages" = "#D96C91"
  ),
  Tissue = c(
    "eutopic endometrium" = "#B9A4D6",
    "endometriotic lesion" = "#D98978"
  )
)

# Paleta divergente: azul = menor expressão relativa,
# claro = média, vermelho = maior expressão relativa.
heatmap_colors <- colorRampPalette(
  c("#2878A5", "#F7F5F0", "#A83246")
)(100)

pheatmap::pheatmap(
  mat_scaled,
  scale = "none",
  color = heatmap_colors,
  breaks = seq(-2, 2, length.out = 101),
  border_color = NA,
  
  # Agrupamento por correlação.
  clustering_method = "average",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  
  # Mostrar nomes dos genes, mas manter AOIs sem rótulos.
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize = 9,
  fontsize_row = 7,
  
  annotation_col = annotation,
  annotation_colors = annotation_colors,
  annotation_names_col = TRUE,
  
  main = paste("Top", nrow(mat_scaled), "most variable genes"),
  angle_col = 45
)

# ============================================================
# 8. DOT PLOT DE MARCADORES POR TIPO CELULAR
# ============================================================

marker_sets <- list(
  Epithelium = c(
    "PAX8", "KRT8", "MSX1", "MUC1", "CLDN4",
    "CLDN7", "FOXJ1", "PAEP", "ESR1", "PGR"
  ),
  Stroma = c(
    "MME", "IGF1", "MMP11", "HAND2", "DCN", "IL1B"
  ),
  Macrophages = c(
    "LBP", "CD68", "CD4", "LYZ", "CD163", "MSR1"
  )
)
marker_table <- data.frame(
  Gene = unlist(marker_sets, use.names = FALSE),
  MarkerClass = rep(names(marker_sets), lengths(marker_sets)),
  stringsAsFactors = FALSE
)

marker_table <- marker_table[!duplicated(marker_table$Gene), ]


# ------------------------------------------------------------
# 8.1. Verificar disponibilidade dos genes marcadores
# ------------------------------------------------------------

missing_markers <- setdiff(
  marker_table$Gene,
  rownames(log_q)
)

if (length(missing_markers) > 0L) {
  message(
    "Marcadores ausentes em log_q: ",
    paste(missing_markers, collapse = ", ")
  )
}

genes_to_check <- c("MME", "CD8A", "CD8B")

for (gene in genes_to_check) {
  message(
    gene,
    if (gene %in% rownames(log_q)) {
      " está presente em log_q."
    } else {
      " NÃO está presente em log_q."
    }
  )
}

# Manter apenas os marcadores presentes na matriz de expressão
marker_table <- marker_table[
  marker_table$Gene %in% rownames(log_q),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 8.2. Alinhar a matriz de expressão e a matriz de detecção
# ------------------------------------------------------------

loq_mat <- LOQ_Mat

# Identificar genes e AOIs comuns às duas matrizes
common_genes <- intersect(
  rownames(log_q),
  rownames(loq_mat)
)

common_aois <- intersect(
  colnames(log_q),
  colnames(loq_mat)
)

# Criar matrizes com dimensões e ordenação correspondentes
log_q_dot <- log_q[
  common_genes,
  common_aois,
  drop = FALSE
]

loq_mat <- loq_mat[
  common_genes,
  common_aois,
  drop = FALSE
]

# Alinhar os metadados à ordem das colunas de expressão
meta <- pData(target_demoData)

meta <- meta[
  match(colnames(log_q_dot), rownames(meta)),
  ,
  drop = FALSE
]

rownames(meta) <- colnames(log_q_dot)

# Verificar se o alinhamento foi realizado corretamente
stopifnot(
  identical(rownames(meta), colnames(log_q_dot)),
  identical(rownames(log_q_dot), rownames(loq_mat)),
  identical(colnames(log_q_dot), colnames(loq_mat))
)

if (anyNA(meta$Cell_Type)) {
  stop(
    "Há AOIs sem correspondência nos metadados. ",
    "Verifique os nomes das amostras."
  )
}

# Garantir que os marcadores também estejam nas matrizes alinhadas
marker_table <- marker_table[
  marker_table$Gene %in% common_genes,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 8.3. Calcular expressão média e percentual acima do LOQ
# ------------------------------------------------------------

if (nrow(marker_table) > 0L) 
  
  dot_rows <- list()
  k <- 1L
  
  for (gene in marker_table$Gene) {
    
    for (cell_type in names(marker_sets)) {
      
      # Selecionar AOIs do tipo celular correspondente
      idx <- meta$Cell_Type == cell_type
      
      if (!any(idx, na.rm = TRUE)) {
        next
      }
      
      # Percentual de AOIs com expressão acima do LOQ
      detected <- if (gene %in% rownames(loq_mat)) {
        mean(
          loq_mat[gene, idx],
          na.rm = TRUE
        ) * 100
      } else {
        NA_real_
      }
      
      # Expressão média normalizada em log2
      mean_expr <- mean(
        log_q_dot[gene, idx],
        na.rm = TRUE
      )
      
      dot_rows[[k]] <- data.frame(
        Gene = gene,
        MarkerClass = marker_table$MarkerClass[
          match(gene, marker_table$Gene)
        ],
        Cell_Type = cell_type,
        MeanLog2Q3 = mean_expr,
        PercentExpressed = detected
      )
      
      k <- k + 1L
    }
  }
  
  marker_summary <- do.call(rbind, dot_rows)
  
  
  # ----------------------------------------------------------
  # 8.4. Padronizar a expressão relativa por gene
  # ----------------------------------------------------------
  
  marker_summary$RelativeExpression <- ave(
    marker_summary$MeanLog2Q3,
    marker_summary$Gene,
    FUN = function(x) {
      sx <- sd(x, na.rm = TRUE)
      
      if (is.na(sx) || sx == 0) {
        rep(0, length(x))
      } else {
        as.numeric(scale(x))
      }
    }
  )
  
  
  # ----------------------------------------------------------
  # 8.5. Definir a ordem dos fatores
  # ----------------------------------------------------------
  
  marker_summary$Gene <- factor(
    marker_summary$Gene,
    levels = marker_table$Gene
  )
  
  marker_summary$MarkerClass <- factor(
    marker_summary$MarkerClass,
    levels = names(marker_sets)
  )
  
  marker_summary$Cell_Type <- factor(
    marker_summary$Cell_Type,
    levels = c("Epithelium", "Stroma", "Macrophages")
  )
  
  
  # ----------------------------------------------------------
  # 8.6. Gerar o Dot Plot
  # ----------------------------------------------------------
  
  marker_summary$PercentExpressed01 <-
    marker_summary$PercentExpressed / 100
  
  p_marker <- ggplot(
    marker_summary,
    aes(
      x = Gene,
      y = Cell_Type,
      size = PercentExpressed01,
      color = RelativeExpression
    )
  ) +
    geom_point(alpha = 0.95) +
    
    # Escala de tamanho: proporção expressa de 0,25 a 1,00.
    # Valores abaixo de 0,25 não são mostrados.
    scale_size_continuous(
      name = "Percent Expressed",
      limits = c(0.25, 1),
      breaks = c(0.25, 0.50, 0.75, 1.00),
      labels = c("0.25", "0.50", "0.75", "1.00"),
      range = c(2, 6)  # reduz a diferença visual entre os tamanhos
    ) +
    scale_color_gradient(
      name = "Relative Expression",
      low = "#FF66B2",   # fraca em -1
      high = "#4A0E4E",  # forte em 1
      limits = c(-1, 1),
      breaks = c(-1, -0.5, 0, 0.5, 1),
      oob = scales::squish,
    ) +
    # Deixar os rótulos dos genes inclinados.
    scale_x_discrete(
      guide = guide_axis(angle = 45)
    ) +
    
    labs(
      x = NULL,
      y = NULL
    ) +
    
    # Fundo branco, sem grade ou moldura; somente os eixos x e y.
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        color = "black"
      ),
      axis.text.y = element_text(color = "black"),
      axis.line.x = element_line(color = "black", linewidth = 1.1),
      axis.line.y = element_line(color = "black", linewidth = 1.1),
      axis.ticks = element_line(color = "black"),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )
  
  print(p_marker)
  
  ggsave(
    filename = file.path(results_dir, "DotPlot_Markers.png"),
    plot = p_marker,
    width = 12,
    height = 6,
    dpi = 300
  )
  
  write.csv(
    marker_summary,
    file = file.path(results_dir, "DotPlot_Markers_Summary.csv"),
    row.names = FALSE
  )
  
# ============================================================
# PCA POR TIPO CELULAR: LESÃO ENDOMETRIÓTICA VS ENDOMÉTRIO EUTÓPICO
# ============================================================

# Metadados das AOIs
metadata <- as.data.frame(pData(target_demoData))

# Conferir se a matriz contém genes nas linhas e AOIs nas colunas
stopifnot(!is.null(rownames(log_q)))
stopifnot(!is.null(colnames(log_q)))

# Garantir correspondência entre metadados e colunas da matriz
if (!all(colnames(log_q) %in% rownames(metadata))) {
  stop("Nem todas as AOIs de log_q3 foram encontradas nos metadados.")
}

metadata <- metadata[colnames(log_q), , drop = FALSE]

# Definir fatores
metadata$Cell_Type <- as.character(metadata$Cell_Type)
metadata$Tissue <- factor(
  as.character(metadata$Tissue),
  levels = c("eutopic endometrium", "endometriotic lesion")
)
metadata$sampleID <- factor(metadata$sampleID)

# Criar pasta de resultados
pca_dir <- file.path(results_dir, "PCA_by_cell_type")
dir.create(pca_dir, showWarnings = FALSE, recursive = TRUE)

# Tipos celulares a analisar
cell_types_to_test <- c("Epithelium", "Stroma", "Macrophages")

# Armazenar objetos PCA e scores
pca_objects <- list()
pca_scores <- list()

for (cell_type in cell_types_to_test) 
  
  cat("\n-----------------------------------\n")
  cat("PCA para:", cell_type, "\n")
  
  # Selecionar AOIs do tipo celular
  idx <- !is.na(metadata$Cell_Type) &
    metadata$Cell_Type == cell_type &
    !is.na(metadata$Tissue)
  
  meta_ct <- metadata[idx, , drop = FALSE]
  
  # Confirmar representação dos dois tecidos
  cat("AOIs por tecido:\n")
  print(table(meta_ct$Tissue, useNA = "ifany"))
  
  if (nrow(meta_ct) < 4L ||
      nlevels(droplevels(meta_ct$Tissue)) < 2L) {
    warning("Dados insuficientes para PCA em ", cell_type)
    next
  }
  
  # Extrair matriz: genes nas linhas e AOIs nas colunas
  mat <- as.matrix(log_q[, rownames(meta_ct), drop = FALSE])
  
  # Remover genes sem variabilidade
  valid_genes <- apply(mat, 1, function(x) {
    all(is.finite(x)) && is.finite(stats::var(x)) &&
      stats::var(x) > 0
  })
  
  mat <- mat[valid_genes, , drop = FALSE]
  
  if (nrow(mat) < 3L) {
    warning("Genes variáveis insuficientes para PCA em ", cell_type)
    next
  }
  
  # PCAtools espera genes nas linhas e AOIs nas colunas.
  # metadata deve ter uma linha por AOI, na mesma ordem das colunas de mat.
  meta_pca <- meta_ct[
    match(colnames(mat), rownames(meta_ct)),
    c("Tissue", "sampleID"),
    drop = FALSE
  ]
  
  rownames(meta_pca) <- colnames(mat)
  
  stopifnot(
    identical(colnames(mat), rownames(meta_pca)),
    all(is.finite(mat))
  )
  
  pca_obj <- PCAtools::pca(
    mat = mat,
    metadata = meta_pca,
    removeVar = 0
  )
  pca_objects[[cell_type]] <- pca_obj
  
  # Extrair coordenadas das AOIs
  scores <- as.data.frame(pca_obj$rotated)
  scores$AOI <- rownames(scores)
  scores$Tissue <- meta_pca[scores$AOI, "Tissue"]
  scores$sampleID <- meta_pca[scores$AOI, "sampleID"]
  
  pca_scores[[cell_type]] <- scores
  
  # Salvar scores completos
  write.csv(
    scores,
    file.path(pca_dir, paste0("PCA_scores_", cell_type, ".csv")),
    row.names = FALSE
  )
  
  # Gráfico PC1 versus PC3
  if (all(c("PC1", "PC3") %in% colnames(scores))) {
    
    p_pc13 <- ggplot2::ggplot(
      scores,
      ggplot2::aes(
        x = PC1,
        y = PC3,
        color = Tissue
      )
    ) +
      ggplot2::geom_point(size = 3, alpha = 0.85) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = paste("PCA -", cell_type),
        x = paste0("PC1 (", round(pca_obj$variance[1], 1), "%)"),
        y = paste0("PC3 (", round(pca_obj$variance[3], 1), "%)"),
        color = "Tissue"
      )
    
    print(p_pc13)
    
    ggplot2::ggsave(
      file.path(pca_dir, paste0("PCA_PC1_PC3_", cell_type, ".png")),
      p_pc13,
      width = 8,
      height = 6,
      dpi = 300
    )
    
  } else {
    warning("PC3 indisponível para ", cell_type)
  }
  
  # Gráfico complementar PC1 versus PC2
  if (all(c("PC1", "PC2") %in% colnames(scores))) {
    
    p_pc12 <- ggplot2::ggplot(
      scores,
      ggplot2::aes(
        x = PC1,
        y = PC2,
        color = Tissue
      )
    ) +
      ggplot2::geom_point(size = 3, alpha = 0.85) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = paste("PCA -", cell_type),
        x = paste0("PC1 (", round(pca_obj$variance[1], 1), "%)"),
        y = paste0("PC2 (", round(pca_obj$variance[2], 1), "%)"),
        color = "Tissue"
      )
    
    print(p_pc12)
    
    ggplot2::ggsave(
      file.path(pca_dir, paste0("PCA_PC1_PC2_", cell_type, ".png")),
      p_pc12,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
# Resumo dos resultados
cat("\nPCAs concluídas para:\n")
print(names(pca_objects))


data_epi <- target_demoData[
  ,
  pData(target_demoData)$Cell_Type == "Epithelium"
]

# Testar a correlação dos PCs com as variáveis dos metadados
p_eigencor <- PCAtools::eigencorplot(
  pca_obj,
  metavars = c("Tissue", "sampleID") # Adicione outras colunas do seu metadata aqui
)
print(p_eigencor)

# ============================================================
# 9. EXPRESSAO DIFERENCIAL POR TIPO CELULAR
#    LESAO ENDOMETRIOTICA VS ENDOMETRIO EUTOPICO
# ============================================================

library(limma)
library(ggplot2)
library(dplyr)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)

# ------------------------------------------------------------
# 9.1. Preparar diretorios e parametros
# ------------------------------------------------------------

de_dir <- file.path(results_dir, "Differential_Expression")
volcano_dir <- file.path(de_dir, "Volcano")
go_dir <- file.path(de_dir, "GO_ORA")
gsea_dir <- file.path(de_dir, "GSEA")
dot_dir <- file.path(de_dir, "DotPlots")

for (d in c(de_dir, volcano_dir, go_dir, gsea_dir, dot_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cell_types <- c("Epithelium", "Stroma", "Macrophages")

tissue_reference <- "eutopic endometrium"
tissue_comparison <- "endometriotic lesion"

fdr_cutoff <- 0.05
logfc_cutoff <- 1

# Confirmar que os objetos principais existem
stopifnot(
  exists("target_demoData"),
  exists("log_q"),
  exists("results_dir")
)

# Excluir controles e manter apenas genes endógenos no teste diferencial.
feature_data <- as.data.frame(fData(target_demoData))

is_endogenous <- as.character(feature_data$CodeClass) == "Endogenous"
is_endogenous[is.na(is_endogenous)] <- FALSE

endogenous_rows <- rownames(feature_data)[is_endogenous]
missing_rows <- setdiff(endogenous_rows, rownames(log_q))

if (length(missing_rows) > 0L) {
  stop(
    "Há genes endógenos de fData ausentes em log_q. Exemplos: ",
    paste(head(missing_rows), collapse = ", ")
  )
}

expr_all <- as.matrix(
  log_q[endogenous_rows, , drop = FALSE]
)

# Metadados e alinhamento
meta_de <- as.data.frame(pData(target_demoData))

meta_de <- meta_de[
  match(colnames(expr_all), rownames(meta_de)),
  ,
  drop = FALSE
]

rownames(meta_de) <- colnames(expr_all)

stopifnot(
  identical(rownames(meta_de), colnames(expr_all)),
  all(c("sampleID", "Tissue", "Cell_Type") %in% colnames(meta_de))
)

# Conferir valores ausentes e escala
if (any(!is.finite(expr_all))) {
  stop("A matriz log_q contém valores não finitos. Investigue antes da análise.")
}

if (any(!is.finite(assayDataElement(target_demoData, elt = "q_norm"))) ||
    any(assayDataElement(target_demoData, elt = "q_norm") <= 0)) {
  stop("q_norm contém valores não positivos ou não finitos; confira a transformação log2.")
}

# Garantir os niveis e a direcao do contraste
meta_de$Tissue <- factor(
  as.character(meta_de$Tissue),
  levels = c(tissue_reference, tissue_comparison)
)

meta_de$Cell_Type <- as.character(meta_de$Cell_Type)
meta_de$sampleID <- factor(meta_de$sampleID)

# Verificar distribuicao amostral
print(with(meta_de, table(Cell_Type, Tissue)))
print(with(meta_de, table(sampleID, Tissue)))

# Salvar metadados usados
write.csv(
  meta_de,
  file.path(de_dir, "Metadata_used.csv"),
  row.names = TRUE
)


# ------------------------------------------------------------
# 9.2. Funcao para executar expressao diferencial
# ------------------------------------------------------------
# Modelo:
#   design: Tissue
#   block: paciente (sampleID)
#   duplicateCorrelation: considera correlacao entre AOIs
#   do mesmo paciente.

# ------------------------------------------------------------
# 9.2. Função para executar expressão diferencial
# ------------------------------------------------------------
# expressão ~ paciente + tecido
#
# O paciente entra como covariável fixa. A comparação é feita
# entre lesão e endométrio eutópico dentro dos pacientes pareados.

de_results_by_cell <- list()
designs_by_cell <- list()

for (ct in cell_types) {
  
  message("\nAnalisando tipo celular: ", ct)
  
  idx <- !is.na(meta_de$Cell_Type) &
    meta_de$Cell_Type == ct &
    !is.na(meta_de$Tissue) &
    !is.na(meta_de$sampleID)
  
  meta_ct <- droplevels(meta_de[idx, , drop = FALSE])
  
  meta_ct$Tissue <- factor(
    as.character(meta_ct$Tissue),
    levels = c(tissue_reference, tissue_comparison)
  )
  meta_ct$sampleID <- factor(meta_ct$sampleID)
  
  # Manter pacientes com AOIs dos dois tecidos.
  pair_table <- table(meta_ct$sampleID, meta_ct$Tissue)
  complete_patients <- rownames(pair_table)[
    rowSums(pair_table > 0) == 2L
  ]
  
  if (length(complete_patients) == 0L) {
    warning(
      "Não há pacientes com os dois tecidos para ",
      ct, "; análise omitida."
    )
    next
  }
  
  meta_ct <- droplevels(
    meta_ct[meta_ct$sampleID %in% complete_patients, , drop = FALSE]
  )
  
  mat_ct <- expr_all[
    ,
    rownames(meta_ct),
    drop = FALSE
  ]
  
  # Remover genes sem variação.
  gene_var <- apply(mat_ct, 1, var, na.rm = TRUE)
  keep_gene <- is.finite(gene_var) & gene_var > 0
  mat_ct <- mat_ct[keep_gene, , drop = FALSE]
  
  if (nrow(mat_ct) == 0L) {
    warning("Não há genes variáveis para ", ct, ".")
    next
  }
  
  if (length(unique(meta_ct$Tissue)) != 2L) {
    warning("Os dois tecidos não estão presentes para ", ct, ".")
    next
  }
  
  print(table(meta_ct$sampleID, meta_ct$Tissue))
  
  # Paciente como covariável fixa, conforme descrito no artigo.
  design <- model.matrix(
    ~ sampleID + Tissue,
    data = meta_ct
  )
  
  if (qr(design)$rank < ncol(design)) {
    warning("Desenho não estimável para ", ct, ".")
    next
  }
  
  fit <- limma::lmFit(
    mat_ct,
    design = design
  )
  
  # eBayes sem robust=TRUE para usar a configuração padrão do limma.
  fit <- limma::eBayes(fit)
  
  coef_tissue <- "Tissueendometriotic lesion"
  
  if (!coef_tissue %in% colnames(fit$coefficients)) {
    warning("Coeficiente de tecido ausente para ", ct, ".")
    next
  }
  
  res <- limma::topTable(
    fit,
    coef = coef_tissue,
    number = Inf,
    adjust.method = "BH",
    sort.by = "none"
  )
  
  res$Gene <- rownames(res)
  res$Cell_Type <- ct
  res$Contrast <- paste(
    tissue_comparison, "vs", tissue_reference
  )
  
  # O paper define DEG por FDR < 0,05.
  # Up/Down são definidos pelo sinal do logFC.
  res$DE_status <- ifelse(
    res$adj.P.Val < fdr_cutoff & res$logFC > 0,
    "Up in lesion",
    ifelse(
      res$adj.P.Val < fdr_cutoff & res$logFC < 0,
      "Down in lesion",
      "Not significant"
    )
  )
  
  de_results_by_cell[[ct]] <- list(
    results = res,
    matrix = mat_ct,
    metadata = meta_ct,
    design = design,
    fit = fit,
    coefficient = coef_tissue
  )
  
  designs_by_cell[[ct]] <- design
  
  write.csv(
    res,
    file.path(de_dir, paste0("DE_all_genes_", ct, ".csv")),
    row.names = FALSE
  )
  
  cat("\nTipo celular:", ct, "\n")
  cat("AOIs:", nrow(meta_ct), "\n")
  cat(
    "Pacientes pareados:",
    nlevels(droplevels(meta_ct$sampleID)),
    "\n"
  )
  cat("Genes testados:", nrow(res), "\n")
  cat(
    "Genes com FDR < 0.05:",
    sum(res$adj.P.Val < fdr_cutoff, na.rm = TRUE),
    "\n"
  )
}

if (length(de_results_by_cell) == 0L) {
  stop(
    "Nenhum modelo diferencial foi concluído. ",
    "Confira os avisos e metadados."
  )
}

# ------------------------------------------------------------
# 9.3. Consolidar os resultados
# ------------------------------------------------------------

de_all <- bind_rows(
  lapply(de_results_by_cell, function(x) x$results)
)

write.csv(
  de_all,
  file.path(de_dir, "DE_all_cell_types.csv"),
  row.names = FALSE
)

de_summary <- de_all %>%
  group_by(Cell_Type) %>%
  summarise(
    Genes_testados = n(),
    Genes_p_nominal = sum(P.Value < 0.05, na.rm = TRUE),
    Genes_FDR_005 = sum(adj.P.Val < 0.05, na.rm = TRUE),
    Genes_FDR_005_absLogFC1 = sum(
      adj.P.Val < 0.05 & abs(logFC) >= logfc_cutoff,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(de_summary)

write.csv(
  de_summary,
  file.path(de_dir, "DE_summary.csv"),
  row.names = FALSE
)


# ============================================================
# 9.4–9.8. VOLCANO, GO ORA E GSEA
# ============================================================

library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)
library(dplyr)
library(ggplot2)

fdr_cutoff <- 0.05

# O artigo define DEG por FDR < 0,05.
# Não aplique um corte de logFC à seleção de genes do GO,
# a menos que queira fazer uma análise adicional.
de_dir <- file.path(results_dir, "Differential_Expression")
volcano_dir <- file.path(de_dir, "Volcano")
go_dir <- file.path(de_dir, "GO_ORA")
gsea_dir <- file.path(de_dir, "GSEA")
dot_dir <- file.path(de_dir, "DotPlots")

for (d in c(de_dir, volcano_dir, go_dir, gsea_dir, dot_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ------------------------------------------------------------
# Helpers para mapeamento de IDs e validação do ranking
# ------------------------------------------------------------

# Converter símbolos para Entrez IDs e guardar:
# - pares de mapeamento;
# - símbolos sem mapeamento;
# - IDs ambíguos ou duplicados.
map_symbols_to_entrez <- function(symbols, output_prefix = NULL) {
  symbols <- unique(trimws(as.character(symbols)))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  
  if (length(symbols) == 0L) {
    stop("A lista de símbolos está vazia.")
  }
  
  raw_map <- suppressMessages(
    clusterProfiler::bitr(
      symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db,
      drop = FALSE
    )
  )
  
  raw_map <- as.data.frame(raw_map)
  raw_map$SYMBOL <- as.character(raw_map$SYMBOL)
  raw_map$ENTREZID <- as.character(raw_map$ENTREZID)
  
  unmapped <- setdiff(symbols, raw_map$SYMBOL[!is.na(raw_map$ENTREZID)])
  
  # Um símbolo que aponta para vários Entrez IDs é ambíguo.
  n_ids_per_symbol <- raw_map |>
    filter(!is.na(ENTREZID)) |>
    distinct(SYMBOL, ENTREZID) |>
    count(SYMBOL, name = "n_entrez")
  
  ambiguous_symbols <- n_ids_per_symbol |>
    filter(n_entrez > 1L) |>
    pull(SYMBOL)
  
  # Para evitar duplicar genes no ranking ou no universo, usar apenas
  # mapeamentos 1:1 símbolo–Entrez.
  one_to_one <- raw_map |>
    filter(!is.na(ENTREZID)) |>
    distinct(SYMBOL, ENTREZID) |>
    add_count(SYMBOL, name = "n_per_symbol") |>
    add_count(ENTREZID, name = "n_per_entrez") |>
    filter(n_per_symbol == 1L, n_per_entrez == 1L) |>
    select(SYMBOL, ENTREZID)
  
  if (!is.null(output_prefix)) {
    write.csv(
      raw_map,
      paste0(output_prefix, "_all_mapping_rows.csv"),
      row.names = FALSE
    )
    write.csv(
      data.frame(SYMBOL = unmapped),
      paste0(output_prefix, "_unmapped_symbols.csv"),
      row.names = FALSE
    )
    write.csv(
      data.frame(SYMBOL = ambiguous_symbols),
      paste0(output_prefix, "_ambiguous_symbols.csv"),
      row.names = FALSE
    )
    write.csv(
      one_to_one,
      paste0(output_prefix, "_one_to_one_mapping.csv"),
      row.names = FALSE
    )
  }
  
  message(
    "Símbolos: ", length(symbols),
    "; mapeados 1:1: ", nrow(one_to_one),
    "; sem mapeamento: ", length(unmapped),
    "; ambíguos: ", length(ambiguous_symbols)
  )
  
  list(
    mapping = one_to_one,
    unmapped = unmapped,
    ambiguous = ambiguous_symbols,
    raw = raw_map
  )
}

# ============================================================
# 9.4. VOLCANO PLOT COM ENHANCEDVOLCANO
# ============================================================

for (ct in names(de_results_by_cell)) {
  
  res <- de_results_by_cell[[ct]]$results
  
  stopifnot(all(c("Gene", "logFC", "adj.P.Val") %in% names(res)))
  
  res <- res |>
    dplyr::filter(
      !is.na(Gene),
      is.finite(logFC),
      is.finite(adj.P.Val)
    )
  
  # Genes significativos por FDR, ordenados do menor para o maior.
  top20 <- res |>
    dplyr::filter(adj.P.Val < fdr_cutoff) |>
    dplyr::arrange(adj.P.Val, dplyr::desc(abs(logFC))) |>
    dplyr::slice_head(n = 20)
  
  labels_top20 <- unique(top20$Gene)
  
  # Cinza para não significativos e salmão para FDR < 0,05.
  point_colors <- ifelse(
    res$adj.P.Val < fdr_cutoff,
    "#D99180",
    "#C9C7C7"
  )
  
  # Os nomes controlam os rótulos da legenda do EnhancedVolcano.
  names(point_colors) <- ifelse(
    res$adj.P.Val < fdr_cutoff,
    "FDR P",
    "Not sig."
  )
  
  p_volcano <- EnhancedVolcano::EnhancedVolcano(
    res,
    lab = res$Gene,
    x = "logFC",
    y = "adj.P.Val",
    
    # Rotular até 20 genes significativos.
    selectLab = labels_top20,
    labFace = "italic",
    boxedLabels = FALSE,
    drawConnectors = TRUE,
    widthConnectors = 0.35,
    max.overlaps = Inf,
    
    # Linhas de referência semelhantes às do exemplo.
    pCutoff = fdr_cutoff,
    FCcutoff = 0,
    hline = -log10(fdr_cutoff),
    vline = 0,
    cutoffLineType = "dashed",
    cutoffLineCol = "black",
    cutoffLineWidth = 0.5,
    
    # Cores e legenda.
    colCustom = point_colors,
    legendLabels = c("Not sig.", "FDR P"),
    legendPosition = "top",
    legendLabSize = 12,
    legendIconSize = 3,
    
    # Títulos e aparência.
    title = NULL,
    subtitle = NULL,
    xlab = expression(Log[2] ~ "(fold change)"),
    ylab = expression(-Log[10] ~ "(FDR)"),
    pointSize = 2,
    labSize = 4,
    colAlpha = 0.75
  ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        colour = "grey85",
        linewidth = 0.5
      ),
      axis.title = ggplot2::element_text(size = 16),
      axis.text = ggplot2::element_text(size = 13),
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank()
    )
  
  print(p_volcano)
  
  ggplot2::ggsave(
    file.path(volcano_dir, paste0("Volcano_", ct, ".png")),
    p_volcano,
    width = 8,
    height = 8,
    dpi = 300
  )
  
  ggplot2::ggsave(
    file.path(volcano_dir, paste0("Volcano_", ct, ".pdf")),
    p_volcano,
    width = 8,
    height = 8
  )
}
# ============================================================
# 9.5. FUNÇÃO DE MAPEAMENTO SYMBOL PARA ENTREZ ID
# ============================================================

map_symbols_to_entrez <- function(symbols, output_prefix = NULL) {
  
  symbols <- unique(trimws(as.character(symbols)))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  
  if (length(symbols) == 0L) {
    stop("A lista de símbolos está vazia.")
  }
  
  raw_map <- suppressMessages(
    clusterProfiler::bitr(
      symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db,
      drop = FALSE
    )
  )
  
  raw_map <- as.data.frame(raw_map)
  raw_map$SYMBOL <- as.character(raw_map$SYMBOL)
  raw_map$ENTREZID <- as.character(raw_map$ENTREZID)
  
  unmapped <- setdiff(
    symbols,
    raw_map$SYMBOL[!is.na(raw_map$ENTREZID)]
  )
  
  # Identificar símbolos que apontam para mais de um Entrez ID.
  n_ids_per_symbol <- raw_map |>
    dplyr::filter(!is.na(ENTREZID)) |>
    dplyr::distinct(SYMBOL, ENTREZID) |>
    dplyr::count(SYMBOL, name = "n_entrez")
  
  ambiguous_symbols <- n_ids_per_symbol |>
    dplyr::filter(n_entrez > 1L) |>
    dplyr::pull(SYMBOL)
  
  # Manter somente relações 1:1 para evitar IDs duplicados
  # no universo do GO e no ranking do GSEA.
  one_to_one <- raw_map |>
    dplyr::filter(!is.na(ENTREZID)) |>
    dplyr::distinct(SYMBOL, ENTREZID) |>
    dplyr::add_count(SYMBOL, name = "n_per_symbol") |>
    dplyr::add_count(ENTREZID, name = "n_per_entrez") |>
    dplyr::filter(n_per_symbol == 1L, n_per_entrez == 1L) |>
    dplyr::select(SYMBOL, ENTREZID)
  
  if (!is.null(output_prefix)) {
    write.csv(
      raw_map,
      paste0(output_prefix, "_all_mapping_rows.csv"),
      row.names = FALSE
    )
    
    write.csv(
      data.frame(SYMBOL = unmapped),
      paste0(output_prefix, "_unmapped_symbols.csv"),
      row.names = FALSE
    )
    
    write.csv(
      data.frame(SYMBOL = ambiguous_symbols),
      paste0(output_prefix, "_ambiguous_symbols.csv"),
      row.names = FALSE
    )
    
    write.csv(
      one_to_one,
      paste0(output_prefix, "_one_to_one_mapping.csv"),
      row.names = FALSE
    )
  }
  
  message(
    "Símbolos: ", length(symbols),
    "; mapeados 1:1: ", nrow(one_to_one),
    "; sem mapeamento: ", length(unmapped),
    "; ambíguos: ", length(ambiguous_symbols)
  )
  
  list(
    mapping = one_to_one,
    unmapped = unmapped,
    ambiguous = ambiguous_symbols,
    raw = raw_map
  )
}


# ============================================================
# 9.6. GO BIOLOGICAL PROCESS POR ORA
# ============================================================
# Genes significativos: FDR < 0,05, separados por direção.
# Universo: todos os genes testados pelo limma que foram
# mapeados para Entrez ID.

go_results <- list()

for (ct in names(de_results_by_cell)) {
  
  res <- de_results_by_cell[[ct]]$results
  
  map_info <- map_symbols_to_entrez(
    res$Gene,
    output_prefix = file.path(de_dir, paste0("GeneMap_", ct))
  )
  
  gene_map_ct <- map_info$mapping
  
  mapped_res <- res |>
    dplyr::inner_join(
      gene_map_ct,
      by = c("Gene" = "SYMBOL")
    )
  
  # Universo de genes testados e mapeados neste tipo celular.
  universe_ids <- unique(mapped_res$ENTREZID)
  
  # Separar genes significativos pelo sinal do logFC.
  up_ids <- unique(
    mapped_res$ENTREZID[
      mapped_res$adj.P.Val < fdr_cutoff &
        mapped_res$logFC > 0
    ]
  )
  
  down_ids <- unique(
    mapped_res$ENTREZID[
      mapped_res$adj.P.Val < fdr_cutoff &
        mapped_res$logFC < 0
    ]
  )
  
  for (direction in c("Up", "Down")) {
    
    gene_ids <- if (direction == "Up") up_ids else down_ids
    
    if (length(gene_ids) < 3L) {
      message(
        "GO ORA omitida para ", ct, " / ", direction,
        ": menos de 3 Entrez IDs significativos."
      )
      next
    }
    
    ego <- tryCatch(
      clusterProfiler::enrichGO(
        gene = gene_ids,
        universe = universe_ids,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.20,
        readable = TRUE
      ),
      error = function(e) {
        warning(
          "enrichGO falhou em ", ct, " / ", direction, ": ",
          conditionMessage(e)
        )
        NULL
      }
    )
    
    if (is.null(ego)) {
      next
    }
    
    go_table <- as.data.frame(ego)
    go_results[[paste(ct, direction, sep = "_")]] <- ego
    
    write.csv(
      go_table,
      file.path(go_dir, paste0("GO_BP_", ct, "_", direction, ".csv")),
      row.names = FALSE
    )
    
    if (nrow(go_table) == 0L) {
      message(
        "Nenhum termo GO passou os filtros para ",
        ct, " / ", direction
      )
      next
    }
    
    # Selecionar os 10 termos com menor p.adjust
    go_plot_data <- go_table |>
      dplyr::filter(
        is.finite(p.adjust),
        is.finite(Count)
      ) |>
      dplyr::arrange(p.adjust) |>
      dplyr::slice_head(n = 10) |>
      dplyr::mutate(
        GeneRatioNumeric = vapply(
          strsplit(as.character(GeneRatio), "/", fixed = TRUE),
          function(x) as.numeric(x[1]) / as.numeric(x[2]),
          numeric(1)
        ),
        Description = sub("^GO:", "", Description),
        Description = stringr::str_wrap(Description, width = 32)
      )
    
    # Deixar o termo mais significativo no topo
    go_plot_data$Description <- factor(
      go_plot_data$Description,
      levels = rev(go_plot_data$Description)
    )
    p_go <- ggplot2::ggplot(
      go_plot_data,
      ggplot2::aes(
        x = Count,
        y = Description,
        size = GeneRatioNumeric,
        color = p.adjust
      )
    ) +
      ggplot2::geom_point(alpha = 0.95) +
      ggplot2::scale_size_continuous(
        name = "GeneRatio",
        range = c(3, 9)
      ) +
      ggplot2::scale_color_gradientn(
        name = "p.adjust",
        colours = c(
          "#292477",
          "#54278F",
          "#9E238B",
          "#D83A75",
          "#F2E600"
        ),
        guide = ggplot2::guide_colorbar(reverse = TRUE)
      ) +
      ggplot2::scale_y_discrete(
        labels = function(x) stringr::str_wrap(x, width = 38)
      ) +
      ggplot2::labs(
        title = paste("GO Biological Process:", ct, "-", direction),
        x = "Count",
        y = NULL
      ) +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(
          colour = "grey85",
          linewidth = 0.5
        ),
        axis.title = ggplot2::element_text(size = 14),
        axis.text.x = ggplot2::element_text(size = 11),
        axis.text.y = ggplot2::element_text(size = 11),
        plot.title = ggplot2::element_text(
          size = 14,
          face = "bold"
        ),
        legend.title = ggplot2::element_text(size = 11),
        legend.text = ggplot2::element_text(size = 10)
      )
    
    print(p_go)
    
    ggplot2::ggsave(
      file.path(dot_dir, paste0("DotPlot_GO_BP_", ct, "_", direction, ".png")),
      p_go,
      width = 8,
      height = 8,
      dpi = 300
    )
  } # fecha loop de direction
}   # fecha loop de ct


# ============================================================
# 9.7. GSEA HALLMARK COM LOGFC RANQUEADO
# ============================================================
# O ranking usa o logFC de todos os genes testados mapeados,
# não somente os genes significativos.

# Compatibilidade entre versões do msigdbr.
hallmark_df <- tryCatch(
  msigdbr::msigdbr(
    species = "Homo sapiens",
    collection = "H"
  ),
  error = function(e) {
    msigdbr::msigdbr(
      species = "Homo sapiens",
      category = "H"
    )
  }
)

# Compatibilidade com nomes de coluna de versões diferentes.
if (!"entrez_gene" %in% names(hallmark_df)) {
  if ("ncbi_gene" %in% names(hallmark_df)) {
    hallmark_df$entrez_gene <- hallmark_df$ncbi_gene
  } else {
    stop("Não foi encontrada uma coluna de Entrez ID em msigdbr.")
  }
}

hallmark_term2gene <- hallmark_df |>
  dplyr::transmute(
    term = as.character(gs_name),
    gene = as.character(entrez_gene)
  ) |>
  dplyr::filter(
    !is.na(term),
    !is.na(gene),
    nzchar(gene)
  ) |>
  dplyr::distinct()

gsea_results_by_cell <- list()

for (ct in names(de_results_by_cell)) {
  
  res <- de_results_by_cell[[ct]]$results
  
  map_info <- map_symbols_to_entrez(
    res$Gene,
    output_prefix = file.path(de_dir, paste0("GeneMap_GSEA_", ct))
  )
  
  gene_map_ct <- map_info$mapping
  
  rank_df <- res |>
    dplyr::transmute(
      Gene = as.character(Gene),
      logFC = as.numeric(logFC)
    ) |>
    dplyr::filter(
      !is.na(Gene),
      is.finite(logFC)
    ) |>
    dplyr::inner_join(
      gene_map_ct,
      by = c("Gene" = "SYMBOL")
    ) |>
    dplyr::distinct(ENTREZID, .keep_all = TRUE)
  
  gene_list <- rank_df$logFC
  names(gene_list) <- rank_df$ENTREZID
  
  # Ordenar em ordem decrescente; desempatar por Entrez ID.
  ord <- order(-gene_list, names(gene_list))
  gene_list <- gene_list[ord]
  
  stopifnot(
    !anyDuplicated(names(gene_list)),
    all(is.finite(gene_list)),
    all(diff(gene_list) <= 0)
  )
  
  if (length(gene_list) < 100L) {
    warning("Menos de 100 genes mapeados para GSEA em ", ct)
    next
  }
  
  gsea_res <- tryCatch(
    clusterProfiler::GSEA(
      geneList = gene_list,
      TERM2GENE = hallmark_term2gene,
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      verbose = FALSE,
      seed = TRUE
    ),
    error = function(e) {
      warning("GSEA falhou em ", ct, ": ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(gsea_res)) {
    next
  }
  
  gsea_results_by_cell[[ct]] <- gsea_res
  gsea_table <- as.data.frame(gsea_res)
  
  write.csv(
    gsea_table,
    file.path(gsea_dir, paste0("GSEA_Hallmark_", ct, ".csv")),
    row.names = FALSE
  )
  
  if (nrow(gsea_table) == 0L) {
    message("Nenhum conjunto GSEA retornado para ", ct)
    next
  }
  
  # Contar genes da leading edge para dimensionar os pontos.
  if ("core_enrichment" %in% colnames(gsea_table)) {
    gsea_table$Count <- lengths(
      strsplit(
        as.character(gsea_table$core_enrichment),
        "/",
        fixed = TRUE
      )
    )
  } else if ("setSize" %in% colnames(gsea_table)) {
    gsea_table$Count <- gsea_table$setSize
  } else {
    stop("Não foi possível determinar a contagem de genes do resultado GSEA.")
  }
  
  gsea_plot_data <- gsea_table |>
    dplyr::filter(
      is.finite(NES),
      is.finite(p.adjust),
      p.adjust < 0.05
    ) |>
    slice_head(n = 10) |>
    dplyr::arrange(p.adjust) |>
    dplyr::arrange(p.adjust) |>
    dplyr::slice_head(n = 20) |>
    dplyr::mutate(
      Description = sub("^HALLMARK_", "", Description),
      Description = gsub("_", " ", Description),
      Description = stringr::str_wrap(Description, width = 38)
    )
  
  # Exibir os termos mais significativos no topo.
  gsea_plot_data$Description <- factor(
    gsea_plot_data$Description,
    levels = rev(gsea_plot_data$Description)
  )
  
  p_gsea <- ggplot2::ggplot(
    gsea_plot_data,
    ggplot2::aes(
      x = NES,
      y = Description,
      size = Count,
      color = p.adjust
    )
  ) +
    ggplot2::geom_point(alpha = 0.95) +
    ggplot2::scale_size_continuous(
      name = "Count",
      range = c(3, 9)
    ) +
    ggplot2::scale_color_gradientn(
      name = "p.adjust",
      colours = c(
        "#292477",
        "#54278F",
        "#9E238B",
        "#D83A75",
        "#F2E600"
      ),
      guide = ggplot2::guide_colorbar(reverse = TRUE)
    ) +
    ggplot2::labs(
      title = paste("GSEA Hallmark:", ct),
      x = "NES",
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        colour = "grey85",
        linewidth = 0.5
      ),
      axis.title = ggplot2::element_text(size = 14),
      axis.text.x = ggplot2::element_text(size = 11),
      axis.text.y = ggplot2::element_text(size = 11),
      plot.title = ggplot2::element_text(
        size = 14,
        face = "bold"
      ),
      legend.title = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 10)
    )
  
  print(p_gsea)
  
  ggplot2::ggsave(
    file.path(dot_dir, paste0("DotPlot_GSEA_Hallmark_", ct, ".png")),
    p_gsea,
    width = 8,
    height = 8,
    dpi = 300
  )
} 


# ============================================================
# 9.8. RESUMOS COMBINADOS DE GO E GSEA
# ============================================================

if (length(go_results) > 0L) {
  
  go_summary <- dplyr::bind_rows(
    lapply(names(go_results), function(nm) {
      tab <- as.data.frame(go_results[[nm]])
      
      if (nrow(tab) == 0L) {
        return(NULL)
      }
      
      tab$Analysis <- nm
      tab
    })
  )
  
  if (nrow(go_summary) > 0L) {
    write.csv(
      go_summary,
      file.path(go_dir, "GO_BP_all_cell_types.csv"),
      row.names = FALSE
    )
  }
}

if (length(gsea_results_by_cell) > 0L) {
  
  gsea_summary <- dplyr::bind_rows(
    lapply(names(gsea_results_by_cell), function(ct) {
      tab <- as.data.frame(gsea_results_by_cell[[ct]])
      
      if (nrow(tab) == 0L) {
        return(NULL)
      }
      
      tab$Cell_Type <- ct
      tab
    })
  )
  
  if (nrow(gsea_summary) > 0L) {
    write.csv(
      gsea_summary,
      file.path(gsea_dir, "GSEA_Hallmark_all_cell_types.csv"),
      row.names = FALSE
    )
  }
}
