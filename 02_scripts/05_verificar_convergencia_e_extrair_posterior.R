# Verificar convergencia do modelo N-mixture TSS

# NB! por enquanto, somente n, p e a1
# falta incluir no script extração de lambda, e a0

# limpar ambiente de trabalho
rm(list = ls())

# carregar pacotes
library(tidyverse)
library(jagsUI)


# carregar dados e resultados
res <- readRDS("03_saida/res_tss_multiregion.rds")
data_long <- readRDS("01_entrada/data_long.rds")
# copiar objetos da lista monitora para o ambiente de trabalho
list2env(data_long, envir = .GlobalEnv)
rm(data_long)

######################################
#----- parte 1: verificar convergencia
######################################


# o numero de eas e de populacoes diferiu entre ucs
# portanto, antes de verificar convergencia ou extrair posterior,
# temos que remover parametros "dummy"
# e.g., estimativas de abundancia para eas e populacoes que nao existem em dada uc
# para isso, vamos criar um array de indicadores lógicos
# (ou matriz de correspondência) a partir dos dados de entrada originais
# e vamos aplica-lo como mascara ao longo da dimensão da simulação mcmc
# as matrizes dos parametros tem dimensoes diferentes, portanto precisam de mascaras diferentes


#----- mascara para o parametro n

# extrair combinacoes validas de pa, pop, ea, ano (como numeros inteiros)
combos_validos <-  tibble(pa = dense_rank(uc),
                          pop = dense_rank(pop),
                          site = dense_rank(ea),
                          year = dense_rank(ano-2013)) |>
  distinct() |>
  as.matrix() # critico para o proximo passo
combos_validos # checar


# construir mascara 4d com combos validos
# definir dimensoes
param_dims <- dim(res$sims.list$n)[2:5]
# criar template preenchido com FALSE
mask_4d <- array(FALSE, dim = param_dims)
# usar combos_validos para transformar combinacoes validas em TRUE
mask_4d[combos_validos] <- TRUE

# aplicar mascara acima aos resultados mcmc
# Make a copy of your posterior array to keep the original safe
npost <- res$sims.list$n


# verificar
# checar estrutura
str(npost)
# verificar quantos NAs nas posterior
sum(is.na(npost[1, , , , ]))
# replicar a mascara 4D ao longo das simulacoes mcmc
for(i in 1:res$mcmc.info$n.samples) {
  npost[i, , , , ][!mask_4d] <- NA
}

# checar
npost[1,1,1,,] # terra do meio
npost[1,36,1,,] # cazumba, pop 1
npost[1,36,2,,] # cazumba, pop 2
npost[1,36,2,,] # cazumba, pop 2
npost[1,29,2,,] # gurupi, pop 2
npost[1,29,2,,] # gurupi, pop 2


# convergencia parametro n
str(res$sims.list$n)
dim(res$Rhat$n)
# aplicar mascara 4d
res$Rhat$n[!mask_4d] <- NA
hist(res$Rhat$n)
abline(v = 1.1, lty = "dashed", lwd = 2, col = "red")

# quais nao convergiram
dim(res$Rhat$n)

# remover da posterior as populacoes com qualquer n que nao convergiu
# isto, e, se nao convergiu em um ano, remover a populacao por completo
# # por fazer!
# logical_array <- res$Rhat$n <= 1.1
# logical_array
# substituir nomes de dimensoes no objeto Rhat
dimnames(res$Rhat$n) <- list(
  cnuc   = uc_indexes$cnuc,
  pop  = 1:dim(res$Rhat$n)[2],
  ea = 1:dim(res$Rhat$n)[3],
  ano = min(uc_indexes$inicio):max(uc_indexes$fim))
# converter o array 4D para data frame
rhat_n_df <- as.data.frame.table(res$Rhat$n, responseName = "Rhat") |>
  mutate(
    pa   = as.character(cnuc),
    pop  = as.integer(pop),
    site = as.character(ea),
    year = as.character(ano),
    Rhat = as.numeric(Rhat)
  ) |>
  drop_na(Rhat)
# sumarizar numero de N que nao convergiram por combinacao de cnuc e pop
unconverged_summary <- rhat_n_df |>
  filter(Rhat > 1.1) |>
  group_by(cnuc, pop, ea, year) |>
  summarise(
    n_unconverged = n(),
    .groups = "drop") |>
  arrange(desc(n_unconverged)) |>
  left_join(uc_indexes |>
              select(cnuc, uc)) |>
  left_join(pop_indexes |>
              select(cnuc, pop, genero, especie),
            by = c("cnuc", "pop")) |>
  select(cnuc, uc, ea, genero, especie, pop, year, n_unconverged) |>
  print(n=Inf)

# salvar para uso posterior
write_csv(unconverged_summary, "03_saida/unconverged_summary.csv")

#----- mascara para o parametro p

# extrair combinacoes validas de pa, pop, ea, ano (como numeros inteiros)
combos_validos <-  tibble(pa = dense_rank(uc),
                          pop = dense_rank(pop),
                          year = dense_rank(ano-2013)) |>
  distinct() |>
  as.matrix() # critico para o proximo passo
combos_validos # checar


# construir mascara 4d com combos validos
# definir dimensoes
param_dims <- dim(res$sims.list$p)[2:4]
# criar template preenchido com FALSE
mask_3d <- array(FALSE, dim = param_dims)
# usar combos_validos para transformar combinacoes validas em TRUE
mask_3d[combos_validos] <- TRUE

# aplicar mascara acima aos resultados mcmc
# Make a copy of your posterior array to keep the original safe
ppost <- res$sims.list$p


# verificar
# checar estrutura
str(ppost)
# verificar quantos NAs nas posterior
sum(is.na(ppost[1, , , ]))
# replicar a mascara 4D ao longo das simulacoes mcmc
for(i in 1:res$mcmc.info$n.samples) {
  ppost[i, ,  , ][!mask_3d] <- NA
}

# checar
ppost[1,1,,] # terra do meio
ppost[1,36,,] # cazumba
ppost[1,29,,] # gurupi


# convergencia parametro n
str(res$sims.list$p)
dim(res$Rhat$p)
# aplicar mascara 4d
res$Rhat$p[!mask_3d] <- NA
hist(res$Rhat$p)
abline(v = 1.1, lty = "dashed", lwd = 2, col = "red")



#----- mascara para o parametro a1 (taxa de crescimento)
dim(res$sims.list$a1)

# extrair combinacoes validas de pa, pop, ea, ano (como numeros inteiros)
combos_validos <-  tibble(pa = dense_rank(uc),
                          pop = dense_rank(pop)) |>
  distinct() |>
  as.matrix() # critico para o proximo passo
combos_validos # checar


# construir mascara 4d com combos validos
# definir dimensoes
param_dims <- dim(res$sims.list$a1)[2:3]
# criar template preenchido com FALSE
mask_2d <- array(FALSE, dim = param_dims)
# usar combos_validos para transformar combinacoes validas em TRUE
mask_2d[combos_validos] <- TRUE

# aplicar mascara acima aos resultados mcmc
# Make a copy of your posterior array to keep the original safe
a1post <- res$sims.list$a1

# verificar
# checar estrutura
str(a1post)
# See how many NAs were introduced per simulation slice
sum(is.na(a1post[1, , ]))
# replicar a mascara 2D ao longo das simulacoes mcmc
for(i in 1:res$mcmc.info$n.samples) {
  a1post[i, ,][!mask_2d] <- NA
}

# checar
a1post[1,,]


# convergencia parametro a1
str(res$sims.list$a1)
dim(res$Rhat$a1)
# aplicar mascara 2d
res$Rhat$a1[!mask_2d] <- NA
hist(res$Rhat$a1)
abline(v = 1.1, lty = "dashed", lwd = 2, col = "red")



######################################
#----- parte 2: extrair posterior
######################################

# como ja salvamos objetos npost, a1post, ppost
# podemos remover objeto res para salvar espaco
rm(res)


#----- npost

# colocar npost em formato long
# e remover NAs
dimnames(npost) <- list(
  sim  = 1:dim(npost)[1],
  cnuc   = uc_indexes$cnuc,
  pop  = 1:dim(npost)[3],
  ea = 1:dim(npost)[4],
  ano = min(uc_indexes$inicio):max(uc_indexes$fim))

# converter para tibble long-format
npost <- npost |> 
  as.data.frame.table(responseName = "n") |> 
  as_tibble() |>
  mutate(across(c(sim, pop, ea, ano, n), \(x) as.numeric(as.character(x)))) |>
  drop_na() |>
  print()

# calcular medias anuais por populacao para economizar espaco
# i.e., ignorar estimativas por ea
npost <- npost |>
  group_by(sim, cnuc, pop, ano) |>
  summarise(n = mean(n)) |>
  ungroup() |>
  print()

# adicionar nomes ucs
npost <- npost |>
  left_join(uc_indexes |>
              select(cnuc, uc)) |>
  relocate(uc, .after=cnuc) |>
  print()

# adicionar nome das especies
npost <- npost |>
  left_join(pop_indexes |>
              mutate(especie = paste(genero, especie, sep = " ")) |>
              select(cnuc, pop, especie),
            by = c("cnuc", "pop")) |>
  relocate(especie, .after = uc) |>
  select(-pop) |>
  print()

# salvar
saveRDS(npost, "03_saida/npost.rds")

# remover para liberar espaco
rm(npost)


#----- ppost


# colocar ppost em formato long
# e remover NAs
dimnames(ppost) <- list(
  sim  = 1:dim(ppost)[1],
  cnuc   = uc_indexes$cnuc,
  pop  = 1:dim(ppost)[3],
  ano = min(uc_indexes$inicio):max(uc_indexes$fim))

# converter para tibble long-format
ppost <- ppost |> 
  as.data.frame.table(responseName = "p") |> 
  as_tibble() |>
  mutate(across(c(sim, pop, ano, p), \(x) as.numeric(as.character(x)))) |>
  drop_na() |>
  print()

# adicionar nomes ucs
ppost <- ppost |>
  left_join(uc_indexes |>
              select(cnuc, uc)) |>
  relocate(uc, .after=cnuc) |>
  print()

# adicionar nome das especies
ppost <- ppost |>
  left_join(pop_indexes |>
              mutate(especie = paste(genero, especie, sep = " ")) |>
              select(cnuc, pop, especie),
            by = c("cnuc", "pop")) |>
  relocate(especie, .after = uc) |>
  select(-pop) |>
  
  print()

# salvar
saveRDS(ppost, "03_saida/ppost.rds")

# remover para liberar espaco
rm(ppost)



#----- a1post


# colocar a1post em formato long
# e remover NAs
dimnames(a1post) <- list(
  sim  = 1:dim(a1post)[1],
  cnuc   = uc_indexes$cnuc,
  pop  = 1:dim(a1post)[3])

# converter para tibble long-format
a1post <- a1post |> 
  as.data.frame.table(responseName = "a1") |> 
  as_tibble() |>
  mutate(across(c(sim, pop), \(x) as.numeric(as.character(x)))) |>
  drop_na() |>
  print()

# adicionar nomes ucs
a1post <- a1post |>
  left_join(uc_indexes |>
              select(cnuc, uc)) |>
  relocate(uc, .after=cnuc) |>
  print()

# adicionar nome das especies
a1post <- a1post |>
  left_join(pop_indexes |>
              mutate(especie = paste(genero, especie, sep = " ")) |>
              select(cnuc, pop, especie),
            by = c("cnuc", "pop")) |>
  relocate(especie, .after = uc) |>
  select(-pop) |>
  print()

# salvar
saveRDS(a1post, "03_saida/a1post.rds")

# remover para liberar espaco
rm(a1post)
