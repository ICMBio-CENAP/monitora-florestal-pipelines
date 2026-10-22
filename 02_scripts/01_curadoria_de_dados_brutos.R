
# carregar pacotes
library(tidyverse)
library(readxl)
library(janitor)
library(stringi)
library(taxize)


# ler dados brutos
# dados brutos baixados de
# https://www.gov.br/icmbio/pt-br/assuntos/monitoramento/conteudo/dados/dados_florestal_monitora/view
data_file <- fs::dir_ls("01_entrada", regexp = "florestal.xlsx")
mydata <- read_excel(data_file, sheet = "mastoaves") |>
  print()

# renomear colunas
mydata <- mydata |>
  rename(nome_ea = estacao) |>
  select(cnuc, uc, ea, nome_ea, esforco, data, ano,
         classe, ordem, familia, genero, especie, taxon_validado,
         nivel_taxon, n_animais, cont_tot_par, distancia, plaqueta) |>
  print()

#View(mydata)


# ajustar formato de data, simplificar nomes de EAs, atualizar genero
mydata <- mydata |>
  dplyr::mutate(
    data = as.Date(data, "%d/%m/%Y"),
    ano = year(data)) |>
  mutate(esforco = as.numeric(esforco)) |>
  mutate(nome_ea = str_to_lower(nome_ea),
         nome_ea = str_replace_all(nome_ea, " ", "_")) |>
  mutate(nome_ea = stri_trans_general(nome_ea, id = "Latin-ASCII")) |>
  mutate(genero = word(taxon_validado, 1),
         especie = word(taxon_validado, 2)) |>
  mutate(classe = str_to_title(classe),
         ordem = str_to_title(ordem),
         familia = str_to_title(familia),
         genero = str_to_title(genero),
         especie = str_to_lower(especie)) |>
  #mutate(especie = replace_na(especie, "")) |> 
  select(-c(taxon_validado, nivel_taxon)) |>
  print()

# usar codigo cnuc correto (codigo completo)
mydata <- mydata %>%
  mutate(
    # converter para numerico para retirar espacos se houver, e entao formatar
    cnuc = sprintf("0000.00.%04d", as.numeric(cnuc))
  ) |>
  print()

# checar potenciais erros em esforco, datas etc
mydata |>
  ggplot(aes(x = esforco)) +
  geom_histogram()

mydata |>
  ggplot(aes(x = data, y = uc)) +
  geom_point(shape = "|", size = 3, color = "darkorchid") + # "|" creates a nice dash line
  theme_minimal() +
  labs(x = "Tempo", y = "UC")

# checar/simplificar nomes de eas
mydata |>
  distinct(uc,ea,nome_ea) |>
  arrange() |>
  print(n=Inf)
# nomes estao OK!

#----- corrigir taxonomia

# criar nova coluna com termo de busca limpo
# bases de dados funcionam melhor com formato "Genero especie"
mydata <- mydata |> 
  mutate(especie = str_remove_all(especie, '["\']'),
         search_term = case_when(genero == familia | is.na(genero) ~ NA_character_,
                                 especie == "sp." | is.na(especie)  ~ genero,
                                 TRUE ~ paste(genero, especie))) |>
  print()

# obter noves validos unicos para a busca
termos_de_busca_unicos <- mydata |> 
  filter(!is.na(search_term)) |> 
  distinct(search_term) |>
  print()


#----- buscar classificacao completa na base de dados global (GBIF)
# ATENCAO! Nao precisa fazer novamente para dados ate 2025
# pular o trcho abaixo e ir diretamente para # ler taxonomia_limpa_lookup!

###############################
###############################
#----- inicio do trecho a pular

# 'return_id = FALSE' retorna diretamente o dataframe com a taxonomia atual
classificacao_de_taxon <- classification(
  termos_de_busca_unicos$search_term, 
  db = "gbif", 
  return_id = FALSE
)
classificacao_de_taxon # checar

# transforma objeto classificacao_de_taxa em dataframe
taxonomia_limpa_long <- classificacao_de_taxon |> 
  purrr::map_dfr(~ as.data.frame(.x), .id = "search_term") |>
  print()

# coloca objeto em formato wide
taxonomia_limpa_lookup <- taxonomia_limpa_long |> 
  filter(rank %in% c("class", "order", "family", "genus", "species")) |> 
  tidyr::pivot_wider(
    id_cols = search_term, 
    names_from = rank, 
    values_from = name
  ) |> 
  print()

# reorganizar colunas
taxonomia_limpa_lookup <- taxonomia_limpa_lookup |>
arrange(class, order, family, genus, species) |>
  relocate(search_term, .after = species) |>
  print()
  
# criar nova coluna com nome de especie a ser utilizado
# substituir "especies" com nome de taxonomia superior pelo genero
# manter especies com nome original
# (assumindo que validacao pelos CNPCs esta correta)
taxonomia_limpa_lookup <- taxonomia_limpa_lookup |>
  rowwise() |>
  mutate( species_new = if_else(
    search_term %in% c(class, order, family),
    genus,
    search_term)) |>
  ungroup() |> # sempre desagrupar depois de usar rowwise()
  print(n=Inf)

# pequenos ajustes
taxonomia_limpa_lookup <- taxonomia_limpa_lookup |>
  mutate(species_new = case_when(species_new == "Cheracebus sp" ~ "Cheracebus",
                                 species_new == "Callicebus hoffmannsi" ~ "Plecturocebus hoffmannsi",
                                 species_new == "Coendou (sphiggurus)" ~ "Coendou spinosus",
                                 .default = species_new)) |>
  mutate(genus = word(species_new, 1),
         species = word(species_new, 2)) |>
  distinct(family, genus, species, .keep_all = TRUE) |>
  print(n=Inf)

# checar se ha duplicacoes
taxonomia_limpa_lookup |> 
  get_dupes(search_term)

# remover duplicacoes
taxonomia_limpa_lookup <- taxonomia_limpa_lookup |>
  distinct(class, order, family, genus, species, search_term, .keep_all = TRUE) |>
  print()

# finalizar taxonomia_limpa_lookup
taxonomia_limpa_lookup <- taxonomia_limpa_lookup |>
  select(class, order, family, genus, species, search_term) |>
  print()
  

# salvar para nao precisar consultar GBIF de novo
write_csv(taxonomia_limpa_lookup, "01_entrada/taxonomia_limpa_lookup.csv")

#----- fim do trecho a pular
############################
############################


# ler taxonomia_limpa_lookup
taxonomia_limpa_lookup <- read_csv("01_entrada/taxonomia_limpa_lookup.csv") |>
  print()

# adiciona taxonomia limpa ao mydata
mydata_standardized <- mydata|>
  # remover colunas antigas para evitar conflitos de sufixo (.x, .y)
  select(-classe, -ordem, -familia, -genero, -especie) |>
  left_join(taxonomia_limpa_lookup, by = "search_term") |>
  # renomear as colunas do lookup para portugues
  rename(
    classe = class,
    ordem = order,
    familia = family,
    genero = genus,
    especie = species
  ) |>
  relocate(classe, ordem, familia, genero, especie, .after = ano) |>
  print()


# checagem manual de especies
mydata_standardized |>
  #filter(classe == "Mammalia") |>
  distinct(classe, ordem, familia, genero, especie, search_term) |>
  arrange(classe, ordem, familia, genero, especie) |>
  print(n=Inf)
# alguns erros persistiram:
# primates, rodentia, psophia, callicebus hoffmansi, carnivora

# corrigir erros acima manualmente :(
mydata_standardized <- mydata_standardized |>
  mutate(classe = case_when(
    search_term %in% c("Carnivora", "Callicebus hoffmannsi", "Primates", "Rodentia") ~ "Mammalia",
    search_term %in% c("Psophia ") ~ "Aves",
    .default = classe
    )) |>
  mutate(ordem = case_when(
    search_term %in% c("Carnivora") ~ "Carnivora",
    search_term %in% c("Callicebus hoffmannsi", "Primates") ~ "Primates",
    search_term %in% c("Rodentia") ~ "Rodentia",
    search_term %in% c("Psophia ") ~ "Gruiformes",
    .default = ordem
  )) |>
  mutate(familia = case_when(
    search_term %in% c("Carnivora", "Primates", "Rodentia") ~ NA,
    search_term %in% c("Callicebus hoffmannsi") ~ "Pitheciidae",
    search_term %in% c("Psophia ") ~ "Psophiidae",
    .default = familia
  )) |>
  # especie antes de genero por conta de Callicebus/Plecturocebus
  mutate(especie = case_when(
    search_term %in% c("Carnivora", "Primates", "Rodentia", "Psophia ") ~ NA,
    search_term %in% c("Callicebus hoffmannsi") ~ "hoffmannsi",
    .default = especie
  )) |>
  mutate(genero = case_when(
    search_term %in% c("Carnivora", "Primates", "Rodentia", "Psophia ") ~ NA,
    search_term %in% c("Callicebus hoffmannsi") ~ "Plecturocebus",
    .default = genero
  )) |>
  print()

# checar novamente
mydata_standardized |>
  #filter(classe == "Mammalia") |>
  distinct(classe, ordem, familia, genero, especie, search_term) |>
  arrange(classe, ordem, familia, genero, especie) |>
  print(n=Inf)

# finalmente, nao precisamos mais da coluna search_term
mydata_standardized  <- mydata_standardized |>
  select(-search_term) |>
  print()


##----- preparar os dados em formato ucs, amostragens, registros, especies

#----- criar objeto ucs
ucs <- mydata_standardized |>
  distinct(cnuc, uc) |>
  mutate(uc = to_snake_case(uc)) |>
  mutate(uc = stri_trans_general(uc, "Latin-ASCII")) |>
  print()


# adicionar bioma IBGE ao objeto ucs
load("01_entrada/tabela_ucs.RData")
uc_table <- uc_table |>
  # ajustar codigo cnuc
  mutate( cnuc = sprintf("0000.00.%04d", as.numeric(cnuc))) |>
  mutate(sigla = tolower(sigla)) |>
  mutate(bioma = to_snake_case(biomaIBGE)) |>
  mutate(bioma = stri_trans_general(bioma, "Latin-ASCII")) |>
  select(cnuc, sigla, bioma) |>
  print()

ucs <- ucs |>
  left_join(uc_table, by="cnuc") |>
  mutate(categoria = case_when(sigla %in% c("esec", "parna", "rebio") ~ "protecao_integral",
                               .default = "uso_sustentavel")) |>
  relocate(bioma, .after = last_col()) |>
  drop_na() |>
  print()

# checar
sort(unique(ucs$uc))
sort(unique(ucs$sigla))
sort(unique(ucs$categoria))
sort(unique(ucs$bioma))


#----- criar objeto amostragens
amostragens <- mydata_standardized |>
  distinct(cnuc, ea, nome_ea, ano, data, .keep_all = TRUE) |>
  select(cnuc, ea, nome_ea, ano, data, esforco) |>
  arrange(cnuc, ea, nome_ea, ano, data, esforco) |>
  # se esforco for NA, completar com esforco medio na EA em questao
  group_by(cnuc, ea) |>
  mutate(esforco = case_when(is.na(esforco) ~ mean(esforco, na.rm = TRUE),
                             .default = esforco)) |>
  ungroup() |>
  mutate(esforco = round(esforco, 2)) |>
  drop_na() |>
  print()


# criar objeto registros
registros <- mydata_standardized |>
  select(cnuc, ea, ano, data,
         classe, ordem, familia, genero, especie,
         n_animais, cont_tot_par, distancia, plaqueta) |>
  # remover duplicacoes
  distinct(cnuc, ea, ano, data,
           genero,especie,
           n_animais, cont_tot_par, distancia, plaqueta,
           .keep_all = TRUE) |>
  #remover NAs de genero em diante
  drop_na(genero) |>
  # remover linhas so de NAs
  filter_all(any_vars(!is.na(.))) |>
  print()


#----- criar objeto especies
especies <- mydata_standardized |>
  distinct(classe, ordem, familia, genero, especie) |>
  arrange(classe, ordem, familia, genero, especie) |>
  print()


# adicionar categoria de ameaca
# ler lista spp ameacadas
redlist <- read_csv("01_entrada/portaria_148_2022_mma_fauna_ameacada.csv") |>
  print()

especies <- especies |>
  left_join(redlist |>
              select(genero, especie, categoria),
            by = c("genero", "especie")) |>
  mutate(categoria = case_when(categoria == "SIM" ~ NA,
                               .default = categoria)) |>
  rename(categoria = categoria) |>
  drop_na(genero) |>
  print()

especies |>
  print(n=Inf)


monitora <- list(ucs = ucs,
                 amostragens = amostragens,
                 registros = registros,
                 especies = especies)

# salvar para uso posterior
saveRDS(monitora, "01_entrada/monitora.rds")


###########################################################################
##### trecho abaixo usar somente para validacao de novos conjuntos de dados
##########################################################################
#----- Checagem adicional para especies do CENAP
# usar somente na fase de validacao

#----- Verificação 1:
# spp vs UC
# exibir sempre que houver mais de 2 spp do mesmo genero na mesma UC
# somente especies alvo do CENAP
mydata_standardized |>
  filter(ordem %in% c("Artiodactyla", "Carnivora", "Cetartiodactyla",
                      "Didelphimorphia", "Lagomorpha",
                      "Perissodactyla", "Rodentia")) |>
  distinct(uc, genero, genero, especie) |>
  group_by(uc, genero) |>
  mutate(nsp = n_distinct(genero, especie)) |>
  filter(nsp > 1) |>
  arrange(uc, genero, especie) |>
  print(n=Inf)

# para verificar:
# Sciurus NA e Sciurus aestuans em diversas UCs
# e.g. Flona do Jamari Sciurus NA e Sciurus aestuans
mydata_standardized |>
  filter(uc == "Flona do Jamari", genero == "Sciurus") |>
  group_by(genero, especie) |>
  count()


#----- Verificação 2:
# quais espécies apareceram pela primeira vez na UC em 2025
# (pode ser um novo registro, ou um erro)

mydata_standardized |>
  filter(ordem %in% c("Artiodactyla", "Carnivora", "Cetartiodactyla",
                      "Didelphimorphia", "Lagomorpha",
                      "Perissodactyla", "Rodentia")) |>
  group_by(uc, genero, especie) |>
  summarize(primeiro_ano = min(ano, na.rm = TRUE), .groups = "drop") |>
  filter(primeiro_ano == 2025) |>
  print(n=Inf)

# para verificar
# Dasyprocta prymnolopha na Esec Jari
# Dasyproctidae no Cabo Orange
# Myoprocta no Viruá
# Dasyprocta fuliginosa na Resex Ituxi
# Dasyprocta sp. em Verde Para Sempre
# Dasyprocta leporina no Rio Cajari
# Myoprocta acouchy no Rio Cajari
# Dasyproctidae no Rio Cautário

# Esec Jari provavelmente Dasyprocta leporina

mydata |>
  filter(uc == "Parna do Cabo Orange", familia == "Dasyproctidae") |>
  group_by(genero, especie) |>
  count()

mydata |>
  filter(uc == "Parna do Viruá", familia == "Dasyproctidae") |>
  group_by(genero, especie, ano) |>
  count()

mydata |>
  filter(uc == "Resex Ituxi", familia == "Dasyproctidae") |>
  group_by(genero, especie) |>
  count()

# Verde Para Sempre = Dasyprocta iacki

mydata |>
  filter(uc == "Resex do Rio Cajari", familia == "Dasyproctidae") |>
  group_by(genero, especie) |>
  count()

mydata |>
  filter(uc == "Resex do Rio Cautário", familia == "Dasyproctidae") |>
  group_by(genero, especie) |>
  count()


#----- Verificação 3:
# atualizar taxonomias
# e.g. Subulo gouazoubira, Dicotyles, 

mydata |>
  filter(grepl("gouazoubira", genero, especie) ) |>
  distinct(genero, especie) |>
  print(n=Inf)
  
mydata |>
  filter(grepl("tajacu", genero, especie) ) |>
  distinct(genero, especie) |>
  print(n=Inf)


