# Preparar dados para analises

# limpar ambiente de trabalho
rm(list = ls())

# carregar pacotes
library(tidyverse)


# ler e ajustar dados
monitora <- readRDS("01_entrada/monitora.rds") |>
  print()


ucs <- monitora$ucs |>
  #filter(bioma == "amazonia") |>
  print()


# para o modelo N-mix, precisamos somente de cnuc, ea, ano, rep
amostragens <- monitora$amostragens |>
  #filter(cnuc %in% ucs$cnuc) |> # NB! somente ucs da amazonia
  group_by(cnuc, ea, ano) |>
  arrange(data) |>
  mutate(rep = dense_rank(data)) |>
  ungroup() |>
  select(cnuc, ea, ano, data, rep) |>
  print()


registros <- monitora$registros |>
  left_join(amostragens,
            by = c("cnuc","ea","ano", "data")) |>
  relocate(rep, .before = classe) |>
  select(cnuc,ea,ano,rep,
         classe,ordem,familia,genero,especie) |>
  print()


especies <- monitora$especies |>
  print()

rm(monitora)


# manter somente 2014-2025
amostragens <- amostragens |>
  filter(ano <= 2025) |>
  print()

registros <- registros |>
  filter(ano <= 2025) |>
  print()


# manter somente UCs com >= 3 anos de dados
pelo_menos_3_anos <- amostragens |>
  group_by(cnuc) |>
  summarise(n_anos = n_distinct(ano)) |> 
  filter(n_anos > 2) |>
  pull(cnuc) |>
  print()

ucs <- ucs |>
  filter(cnuc %in% pelo_menos_3_anos) |>
  print()

amostragens <- amostragens |>
  filter(cnuc %in% pelo_menos_3_anos) |>
  print()

registros <- registros |>
  filter(cnuc %in% pelo_menos_3_anos) |>
  print()


# checar se ha EAs amostradas num unico ano
amostragens |> 
  arrange(cnuc, ea) |>
  distinct(cnuc, ea, ano) |>
  group_by(cnuc, ea) |>
  count() |>
  arrange(n) |>
  print()
# foram amostrados num unico ano:
# pa 0000.00.0123 ea 4
# pa 0000.00.0143 ea 3
# pa 0000.00.0258 ea 3
# remover?

# manter somente as ea com >= 2 anos de dados
ea_tempo_minimo <- amostragens |>
  group_by(cnuc, ea) |>
  summarise(n_anos = n_distinct(ano)) |> 
  filter(n_anos > 2) |>
  mutate(uc_ea = paste(cnuc,ea, sep = "_")) |>
  print()
# manter somente amostras que estao no objeto ea_tempo_minimo
amostragens <- amostragens |>
  filter(paste(cnuc, ea, sep = "_") %in% ea_tempo_minimo$uc_ea) |>
  print()
# manter somente ucs que estao no objeto amostragens atualizado acima
ucs <- ucs |>
  filter(cnuc %in% amostragens$cnuc) |>
  print()
# manter somente registros oriundos de eas que estao no objeto ea_tempo_minimo
# and in the updated samples object
registros <- registros |>
  filter(paste(cnuc, ea) %in% paste(amostragens$cnuc, amostragens$ea)) |>
  print()

# checar novamente se ha EAs amostradas em < 3 anos
amostragens |> arrange(cnuc, ea) |>
  distinct(cnuc, ea, ano) |>
  group_by(cnuc, ea) |>
  count() |>
  arrange(n) |>
  print()



#----- preparar dados para obter contagens por especie por ano por amostra (rep)


# checar terra do meio
amostragens |> 
  filter(cnuc == "0000.00.0047") |>
  distinct(ea, ano, data) |>
  group_by(ea, ano) |>
  count()
# sampling is too unequal
# maybe we should divide terra do meio into "two" different ucs


# algumas eas tem > 10 dias de amostragem em um unico ano
# isso pode inflar as dimensoes do y-array e aumentar o tempo de computacao
# vamos definir um teto
amostragens |>
  group_by(cnuc,ea,ano) |>
  summarise(nreps = n_distinct(data)) |>
  arrange(cnuc, ea, ano) |>
  group_by(nreps) |>
  count() |>
  arrange(desc(nreps)) |>
  print()
# vamos remover reps > #10 para reduzir tempo de computacao
amostragens <- amostragens |>
  filter(rep < 11) |>
  print()

registros <- registros |>
  filter(rep < 11) |>
  print()


#----- definir especies/populacoes-alvo

# alguns criterios a considerar:
# numero de registros: populacoes com poucos registros nao tem sinal populacional
# vamos definir "pouco" arbitrariamente como: registrada em < 15% das aostragens
# dados muito esparsos podem levar a problemas de identificabilidade p x N

# primeiro, checar visualmente
registros |>
  #drop_na(especie) |>
  group_by(cnuc, genero, especie) |>
  summarise(n_reps = n_distinct(ea, ano, rep))|>
  ungroup() |>
  left_join(amostragens |>
              group_by(cnuc) |>
              count() |>
              ungroup(),
            by = "cnuc") |>
  mutate(percent_reps = n_reps/n) |>
  print() |>
  ggplot(aes(x = percent_reps)) +
  geom_histogram()

# checar quantas populacoes resistem ao criterio
registros |>
  drop_na(especie) |>
  group_by(cnuc, genero, especie) |>
  summarise(n_reps = n_distinct(ea, ano, rep))|>
  ungroup() |>
  left_join(amostragens |>
              group_by(cnuc) |>
              count() |>
              ungroup(),
            by = "cnuc") |>
  mutate(percent_reps = n_reps/n) |>
  filter(percent_reps >= 0.15) |>
  print()

# se agruparmos tinamidae, seria possivel aumentar numero de pops
# poderiam ser incluidas dezenas de populacoes de Crypturellus sp:
registros |>
  filter(familia == 'Tinamidae') |>
  group_by(genero, especie) |>
  summarise(n = n(),
            nucs = n_distinct(cnuc))

# criar Crypturellus
registros <- registros |>
  mutate(especie = case_when(genero == "Crypturellus" ~ "sp",
                             .default = especie)) |>
  print()

# aplicar o filtro de 15%
registros <- registros |>
  drop_na(especie) |>
  group_by(cnuc, genero, especie) |>
  mutate(n_reps = n_distinct(ea, ano, rep))|>
  ungroup() |>
  left_join(amostragens |>
              group_by(cnuc) |>
              count() |>
              ungroup(),
            by = "cnuc") |>
  mutate(percent_reps = n_reps/n) |>
  filter(percent_reps >= 0.15) |>
  select(-n_reps, -n, -percent_reps) |>
  print()

# atualizar especies mantendo somente as que permaneceram em registros apos filtro
especies <-  especies |>
  filter(paste(genero, especie) %in% paste(registros$genero, registros$especie)) |>
  print()


# criar coluna pop com id_unico para cada populacao
registros <- registros |>
  arrange(cnuc, classe, ordem, familia, genero, especie) |>
  group_by(cnuc) |>
  mutate(
    # calcular rank e formato como string de dois digitos
    pop_number = sprintf("%02d", dense_rank(paste(genero, especie)))
  ) |>
  mutate(pop_id = paste(cnuc, pop_number, sep = "_")) |>
  select(-pop_number) |>
  #mutate(pop_number = as.numeric(pop_number)) |>
  ungroup() |>
  print()
  


#----- criar objeto y que sera input para modelo N-mixture
# primeiramente, criar raw_counts com contagens incluindo zeros

raw_counts <- amostragens |>
  select(cnuc,ea,ano,rep) |>
  full_join(registros |> 
              mutate(n = 1),
            relationship = "many-to-many") |>
  select(cnuc,ea,ano,rep,pop_id,n) |>
  mutate(n = case_when(is.na(n) ~ 0,
                       .default = n)) |>
  # obter contagens por pop por rep
  group_by(cnuc, pop_id, ano, ea, rep, .drop=FALSE) |>
  summarise(count = sum(n, na.rm=TRUE)) |>
  ungroup() |>
  # agrupar por uc, e dentro do agrupamento
  # cruzar todas as possiveis combinacoes de pops com
  # `(ea, ano e rep)` ja existentes nos dados
  # e preencher n = 0 para obter contagens = zero
  group_by(cnuc) |>
  complete(pop_id,
           nesting(ano, ea, rep),
           fill = list(count = 0)) |>
  arrange(cnuc, pop_id, ano, rep, ea, count) |>
  ungroup() |>
  # o primeiro full_join nesta secao criou algumas pop_id == NA
  # elas prepresentam dias de amostragem sem nenhum registro das pops alvo
  # elas precisavam ser mantidas ate usarmos a funcao complete para
  # atribuir o zero para as opos alvo nos dias sem contagem
  # agora podemos remover:
  drop_na(pop_id) |>
  # finalmente, vamos adicionar coluna pop para representar id da populacao DENTRO da UC
  # ao contrario de pop_id, pop pode ter valores repetidos (reinicia contagem para cada uc)
  group_by(cnuc) |>
  mutate(pop = dense_rank(pop_id)) |>
  relocate(pop, .after = pop_id) |>
  ungroup() |>
  print()



#----- preparar os dados em long format para input no modelo JAGS

nobs <- nrow(raw_counts)
y <- raw_counts |> pull(count)
pop_id <- raw_counts |> pull(pop_id)
pop <- raw_counts |> pull(pop)
npop <- max(raw_counts$pop)
uc <- raw_counts |> pull(cnuc)
nuc <- length(unique(uc))
ea <- raw_counts |> pull(ea)
nea <- length(unique(ea))
ano <- raw_counts |> pull(ano)
nano <- length(unique(ano))
rep <- raw_counts |> pull(rep)
nrep <- length(unique(rep))

# criar indices de uc (sera usado mais tarde para associar com a posterior)
uc_indexes <- raw_counts |>
  left_join(ucs |> distinct(cnuc, uc)) |>
  group_by(cnuc, uc) |>
  summarise(nea = n_distinct(ea),
            nano = n_distinct(ano),
            inicio = min(ano),
            fim = max(ano),
            npops = n_distinct(pop_id)) |>
  ungroup() |>
  data.frame() |>
  print()

# criar indices de populacao/especie
pop_indexes <- raw_counts |>
  left_join(registros |> distinct(genero, especie, pop_id) |> drop_na(especie)) |>
  distinct(cnuc,pop_id, pop, genero, especie) |>
  arrange(cnuc,pop_id) |>
  print()

# checar quantas populacoes de cada especie
pop_indexes |>
  group_by(genero, especie) |>
  count() |>
  arrange(genero, especie) |>
  print(n=Inf)


data_long <- list(y = y,
                  nobs = nobs,
                  pop = pop,
                  npop = npop,
                  uc = uc,
                  nuc = nuc,
                  ea = ea,
                  nea = nea,
                  ano = ano,
                  nano = nano,
                  rep = rep,
                  nrep = nrep,
                  uc_indexes = uc_indexes,
                  pop_indexes = pop_indexes)

saveRDS(data_long, "01_entrada/data_long.rds")

