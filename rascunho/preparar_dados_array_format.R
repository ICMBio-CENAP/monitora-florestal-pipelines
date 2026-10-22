
# codigo para gerar dados de input do JAGS no formato de array multidimensional

#----- from long format to array

complete_counts <- raw_counts |>
  select(-species) |>
  complete(ea,rep,ano,pop,
           nesting(uc),
           fill = list(count = NA)) |>
  arrange(uc,pop,ano,rep,ea) |>
  print()

y_array <- with(complete_counts,
                array(data = count,
                      dim=c(length(unique(ea)),
                            length(unique(rep)),
                            length(unique(ano)),
                            length(unique(pop)),
                            length(unique(uc))),
                      dimnames = list(unique(ea),
                                      unique(rep),
                                      unique(ano),
                                      unique(pop),
                                      unique(uc))))

str(y_array)


#----- check if array is ok


str(y_array)
# uc = TdM, species = Aburria cujubi
raw_counts |>
  filter(uc==47, pop==1, ano == 2014) |>
  pivot_wider(names_from = rep, values_from = count) |>
  select(-c(uc,species,pop,ea,ano)) |>
  print()
y_array[,,1,1,1]

# uc = cazumba-iracema
raw_counts |> filter(uc==232) |> distinct(species, .keep_all = TRUE)
dimnames(y_array)
raw_counts |>
  filter(uc==232, pop==9, ano == 2024) |>
  pivot_wider(names_from = rep, values_from = count) |>
  select(-c(uc,species,pop,ea,ano)) |>
  print()
y_array[,,11,9,20]

raw_counts |>
  filter(uc==232, pop==6, ano == 2020) |>
  pivot_wider(names_from = rep, values_from = count) |>
  select(-c(uc,species,pop,ea,ano)) |>
  print()
y_array[,,7,6,20]


# uc = uatuma
raw_counts |> filter(uc==213) |> distinct(species, .keep_all = TRUE)
dimnames(y_array)
raw_counts |>
  filter(uc==213, pop==1, ano == 2019) |>
  pivot_wider(names_from = rep, values_from = count) |>
  select(-c(uc,species,pop,ea,ano)) |>
  print()
y_array[,,6,1,18]

raw_counts |>
  filter(uc==213, pop==5, ano == 2022) |>
  pivot_wider(names_from = rep, values_from = count) |>
  select(-c(uc,species,pop,ea,ano)) |>
  print()
y_array[,,9,5,18]

# check again for a ano for Guerlinguetus for a ano with records for this species
raw_counts |>
  filter(uc==213, pop==5, ano == 2014) |>
  pivot_wider(names_from = rep, values_from = count) |>
  select(-c(uc,species,pop,ea,ano)) |>
  print()
y_array[,,1,5,18]

raw_counts |>
  filter(uc==213, pop==4, ano == 2022) |>
  pivot_wider(names_from = rep, values_from = count) |>
  select(-c(uc,species,pop,ea,ano)) |>
  print()
y_array[,,9,4,18]

# hooray! It seems that long and wide format are in agreement


#----- preucre data wide format for JAGS model (ucrt 2)


data_wide <- list(y = y_array,
                  npop = dim(y_array)[4],
                  nuc = dim(y_array)[5],
                  nea = dim(y_array)[1],
                  nano = dim(y_array)[3],
                  nrep = dim(y_array)[2],
                  indexes = indexes,
                  pop_indexes = pop_indexes)

saveRDS(data_wide, here("data", "data_wide.rds"))

