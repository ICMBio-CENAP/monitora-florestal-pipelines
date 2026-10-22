# Explorar os dados do programa Monitora

# limpar ambiente de trabalho
rm(list = ls())

# carregar pacotes
library(tidyverse)


# ler e ajustar dados
monitora <- readRDS("01_entrada/monitora.rds") |>
  print()

# copiar objetos da lista monitora para o ambiente de trabalho
list2env(monitora, envir = .GlobalEnv) 

# checar
ucs
amostragens
registros
especies


#----- ucs

# quantas ucs, categorias, biomas
nrow(ucs)
ucs |>
  group_by(sigla) |>
  count()
ucs |>
  group_by(categoria) |>
  count()
ucs |>
  group_by(bioma) |>
  count()


#----- amostragens

amostragens |>
  ggplot(aes(x = data, y = cnuc)) +
  geom_point(shape = "|", size = 3, color = "darkorchid") +
  theme_minimal() +
  labs(x = "", y = "UC")

# esforco
amostragens |>
  ggplot(aes(x = esforco)) +
  geom_histogram() +
  theme_minimal() +
  labs(x = "Esforço (km)")


#----- registros

# quantas ordens familias, generos, especies na base de dados
registros |>
  group_by(classe) |>
  summarise(ordens = n_distinct(ordem),
            familias = n_distinct(familia),
            generos = n_distinct(genero),
            especies = n_distinct(genero, especie),
            populacoes = n_distinct(cnuc, genero, especie))

registros |>
  filter(classe == "Aves") |>
  group_by(classe, ordem, familia) |>
  summarise(n_ucs = n_distinct(cnuc),
            n_registros = n()) |>
  arrange(classe, ordem, familia) |>
  print(n=Inf)

registros |>
  filter(classe == "Mammalia") |>
  group_by(classe, ordem, familia) |>
  summarise(n_ucs = n_distinct(cnuc),
            n_registros = n()) |>
  arrange(classe, ordem, familia) |>
  print(n=Inf)

# quais as cinco especies de aves e mamiferos mais registradas
registros |>
  group_by(classe, genero, especie) |>
  count() |>
  group_by(classe) |>
  slice_max(order_by = n, n = 5, with_ties = FALSE) |>
  ungroup() |>
  print()

# idem, mas excluindo especie == NA
registros |>
  drop_na(especie) |>
  group_by(classe, genero, especie) |>
  count() |>
  group_by(classe) |>
  slice_max(order_by = n, n = 5, with_ties = FALSE) |>
  ungroup() |>
  print()
  

# numero total de registros por especie
registros |>
  group_by(genero, especie) |>
  count() |>
  ungroup() |>
  mutate(nome_cientifico = paste(genero, especie)) |> 
  mutate(nome_cientifico = fct_reorder(nome_cientifico, n, .desc = TRUE)) |> 
  ggplot(aes(x = nome_cientifico, y = n)) +
  geom_col() +
  theme_minimal() +
  scale_x_discrete(labels = NULL) +
  labs(x = "Espécie (ranking)", y = "Número de Registros")


# numero total de registros por populacao (especie + uc)
registros |>
  group_by(cnuc, genero, especie) |>
  count() |>
  ungroup() |>
  mutate(pop = paste(cnuc,genero, especie)) |> 
  mutate(pop = fct_reorder(pop, n, .desc = TRUE)) |> 
  ggplot(aes(x = pop, y = n)) +
  geom_col() +
  theme_minimal() +
  scale_x_discrete(labels = NULL) +
  labs(x = "População (ranking)", y = "Número de Registros")


# numero de ucs em que cada especie foi registrada
registros |>
  group_by(genero, especie) |>
  summarise(n_ucs = n_distinct(cnuc)) |>
  ungroup() |>
  mutate(nome_cientifico = paste(genero, especie)) |> 
  mutate(nome_cientifico = fct_reorder(nome_cientifico, n_ucs, .desc = TRUE)) |> 
  ggplot(aes(x = nome_cientifico, y = n_ucs)) +
  geom_col() +
  theme_minimal() +
  scale_x_discrete(labels = NULL) +
  labs(x = "Espécie (ranking)", y = "Número de UCs")


# percentagem de amostragens em que especie foi registrada
registros |>
  #drop_na(especie) |>
  group_by(cnuc, genero, especie) |>
  summarise(n_reps = n_distinct(ea, ano, data))|>
  ungroup() |>
  left_join(amostragens |>
              group_by(cnuc) |>
              count() |>
              ungroup(),
            by = "cnuc") |>
  mutate(percent_reps = n_reps/n) |>
  print() |>
  ggplot(aes(x = percent_reps)) +
  geom_histogram() +
  theme_minimal() +
  labs(x = "Percentagem das amostragens")



