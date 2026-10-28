# Explorar resultados do modelo N-mixture TSS

# limpar ambiente de trabalho
rm(list = ls())

# carregar pacotes
library(tidyverse)
#library(gridExtra)
#library(janitor)


# ler dados ---------------------------------------------------
data_long <- readRDS("01_entrada/data_long.rds")
monitora <- readRDS("01_entrada/monitora.rds")
npost <- readRDS("03_saida/npost.rds")
rpost <- readRDS("03_saida/a1post.rds")
ppost <- readRDS("03_saida/ppost.rds")



# explorar deteccao p -----------------------------------------------------

ppost |>
  group_by(uc,especie,ano) |>
  summarise(p_mean = mean(p)) |>
  ggplot(aes(x = p_mean)) +
  geom_histogram(bins = 20, 
                 color = "black", fill = "lightblue",
                 linewidth = 0.2, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Prob. detecção", y = "Freq") +
  theme(axis.text = element_text(size = 12)) +
  theme(
    axis.title.x = element_text(size = 14), # Adjust size for X-axis title
    axis.title.y = element_text(size = 14)  # Adjust size for Y-axis title
  )

ggsave(plot = last_plot(),
       scale = 4,
       units = "cm",
       width = 6,
       height = 4,
       dpi = 120,
       "03_saida/p_histograma.jpeg")

# sumarios de deteccao
p_summary <- ppost |>
  group_by(uc,especie) |>
  summarise(p_mean = mean(p)) |>
  ungroup() |>
  summarise(pmean = round(mean(p_mean),2),
            pmin = round(min(p_mean),2),
            pmax = round(max(p_mean),2)) |>
  print()



# explorar abundancia relativa n ------------------------------------------------

# grafico de distribuicao de abundancias
npost |>
  group_by(uc,especie) |>
  summarise(mean_n = mean(n)) |>
  ungroup() |>
  ggplot(aes(x = mean_n)) +
  geom_histogram(#bins = 30,
                 color = "black", fill = "lightblue",
                 linewidth = 0.2, alpha = 0.5) +
  #theme_bw() +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Abundancia media", y = "Freq") +
  theme(axis.text = element_text(size = 12)) +
  theme(
    axis.title.x = element_text(size = 14), # Adjust size for X-axis title
    axis.title.y = element_text(size = 14)  # Adjust size for Y-axis title
  )


# sumario de abundancias relativas
n_summary <- npost |>
  group_by(uc,especie) |>
  summarise(n_mean = mean(n)) |>
  ungroup() |>
  summarise(nmean = round(mean(n_mean),2),
            nmin = round(min(n_mean),2),
            nmax = round(max(n_mean),2)) |>
  print()


# checar n suspeitos (muito altos)
npost |>
  group_by(uc,especie) |>
  summarise(mean_n = mean(n)) |>
  ungroup() |>
  arrange(desc(mean_n))


# excluir n suspeitos (se for o caso)
# npost <- npost |>
#   filter(! (pa == 189 & species == "Cebus olivaceus")) |>
#   print()
# rpost <- rpost |>
#   filter(! (pa == 189 & species == "Cebus olivaceus")) |>
#   print()
# ppost <- ppost |>
#   filter(! (pa == 189 & species == "Cebus olivaceus")) |>
#   print()


# grafico tendencias todas as pops
npost |>
  mutate(grp = dense_rank(paste(uc, especie))) |>
  summarise(.by = c(grp, ano),
            n_mean = mean (n)) |>
  arrange(ano) |>
  group_by(grp) |>
  mutate(n_std = n_mean/first(n_mean)) |>
  ggplot(aes(x = ano, y = n_std)) +
  geom_smooth(aes(group = grp),
              se=FALSE,
              color = "gray80", linewidth = 0.1, alpha = 0.1) +
  geom_smooth(method = "loess", se = FALSE,
              color = "black", linewidth = 0.4) +
  theme_bw() +
  ylim(-0.7, NA) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "", y = "Abundância relativa (padronizada)")


# salvar figura
ggsave(plot = last_plot(),
       scale = 3,
       units = "cm",
       width = 6,
       height = 4,
       dpi = 120,
       "03_saida/all_trends.jpeg")

# grafico tendencias por uc media + barras de erros
npost |>
  filter(cnuc == "0000.00.0232") |>
  group_by(especie, ano) |>
  summarise(mean = mean(n, na.rm=TRUE),
            lower = quantile(n, prob = 0.1, na.rm=TRUE),
            upper = quantile(n, prob = 0.9, na.rm=TRUE)) |>
  ungroup() |>
  ggplot(aes(x = ano, y = mean)) +
  geom_line(aes(x = ano, y = mean), linewidth = 0.3, alpha = 0.8) +
  geom_point(aes(x = ano, y = mean),  shape = 21,
             color = "white", stroke = 2,
             fill = "black", size = 2) +
  geom_linerange(aes(x = ano, y = mean,
                     ymin=lower, ymax=upper),
                 size = 0.2, alpha = 0.5) +
  #coord_cartesian(xlim = c(2014, 2025), expand = FALSE) +
  scale_x_continuous(breaks = scales::breaks_pretty(n=5)) +
  ylim(0,NA) +
  xlab("") +
  ylab("Abundância relativa") +
  #theme_classic() +
  theme_bw() +
  theme(axis.text = element_text(size = 12)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(~especie, scales = "free", ncol = 3) +
  theme(strip.text = element_text(size = 14))

# salvar figura
ggsave(plot = last_plot(),
       scale = 4,
       units = "cm",
       width = 8,
       height = 6,
       dpi = 120,
       "03_saida/cazumba_trends_linerange.jpeg")


# grafico tendencias por uc ribbon
npost |>
  filter(cnuc == "0000.00.0232") |>
  group_by(especie, ano) |>
  summarise(mean = mean(n, na.rm=TRUE),
            lower = quantile(n, prob = 0.1, na.rm=TRUE),
            upper = quantile(n, prob = 0.9, na.rm=TRUE)) |>
  ungroup() |>
  ggplot(aes(x = ano, y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.9) +
  geom_line(aes(x = ano, y = mean), linewidth = 0.2, alpha = 0.5) +
  geom_point(aes(x = ano, y = mean),  shape = 21,
             color = "lightblue", stroke = 3,
             fill = "black", size = 1.5) +
  coord_cartesian(xlim = c(2014, 2025), expand = FALSE) +
  scale_x_continuous(breaks = scales::breaks_pretty(n=5)) +
  ylim(0,NA) +
  xlab("") +
  ylab("Abundância relativa") +
  #theme_classic() +
  theme_bw() +
  theme(axis.text = element_text(size = 12)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(~especie, scales = "free", ncol = 3) +
  theme(strip.text = element_text(size = 14))



# explorar lambda ----------------------------------------

#lambdapost |>
#  group_by(uc,especie) |>
#  summarise(mean_lambda = mean(lambda)) |>
#  ungroup() |>
#  ggplot(aes(x = mean_lambda)) +
#  geom_histogram(bins = 20,
#                 color = "black", fill = "lightblue",
#                 alpha = 0.6) +
#  theme_bw() +
#  theme(panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank()) +
#  labs(x = "lambda médio", y = "Freq") +
#  xlim(0,100) # NB! tem um valor muito alto!


# explorar tendencia temporal r (=a1) -----------------------------------------------------------

# tabela r
tabela_r <- rpost |>
  rename(r = a1) |>
  group_by(uc,especie) |>
  summarise(r_mean = round(mean(r),2),
            lower = round(quantile(r, prob = 0.025),2),
            upper = round(quantile(r, prob = 0.975),2),
            prob_increase = round(length(r[r>=0])/n(),2),
            prob_decrease = round(length(r[r<=0])/n(),2)) |>
  ungroup() |>
  print()

# adicionar info sobre especie
tabela_r <- tabela_r |>
  left_join(monitora$especies |>
              rename(cat_ameaca = categoria) |>
              mutate(especie = paste(genero, especie)) |>
              select(classe, especie, cat_ameaca)) |>
  #relocate(classe, cat_ameaca, .after = especie) |>
  print()

# adicionar info sobre uc
tabela_r <- tabela_r |>
  left_join(monitora$ucs |>
              rename(categoria_uc = categoria) |>
              select(-cnuc)) |>
  #relocate(sigla, categoria_uc, bioma, .after = uc) |>
  print()


# adiciona abundancia inicial
tabela_r <- tabela_r |>
  left_join(npost |>
              group_by(uc, especie, ano) |>
              summarise(mean_n = mean(n)) |>
              group_by(uc, especie) |>
              summarise(n_inicial = first(mean_n)) |>
              ungroup(),
            by = c("uc", "especie")) |>
  print()


# verificar se ha duplicacoes
tabela_r |> 
  get_dupes()


# quem tem > 90% de estar em declinio  
tabela_r |>
  filter(upper < 0) |>
  arrange(desc(r_mean)) |>
  print(n=Inf)

# quem tem > 90% de estar em aumento  
tabela_r |>
  filter(lower >= 0) |>
  arrange(desc(r_mean)) |>
  print(n=Inf)

# quais UCs com declinios  
tabela_r |>
  filter(upper < 0) |>
  group_by(uc) |>
  count() |>
  print(n=Inf)

# quais UCs com aumentos  
tabela_r |>
  filter(lower >= 0) |>
  group_by(uc) |>
  count() |>
  print(n=Inf)

# porcentagem pops de cada UC com declinios  
tabela_r |>
  mutate(dec = case_when(upper <= 0 ~ "declinio",
                         .default = "no")) |>
  group_by(uc) |>
  summarise(n_pop = n(),
            n_dec = sum(dec == "declinio"),
            prop_dec = round(sum(dec == "declinio")/n(),1)) |>
  arrange(desc(prop_dec)) |>
  filter(n_dec > 0) |>
  print(n=Inf)



# plotar histograma do r
tabela_r |>
  mutate(trend = case_when(prob_decrease > 0.95 ~ "declinio",
                           prob_increase > 0.95 ~ "aumento",
                           .default = "estavel")) |>
  mutate(trend = factor(trend, levels = c("aumento", "estavel", "declinio"))) |>
  ggplot() +
  geom_histogram(aes(x = r_mean, fill = trend),
                 #bins = 20,
                 alpha = 0.8) +
  scale_fill_manual(values = c("aumento" = "blue", "estavel" = "darkgrey", "declinio" = "red")) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.1, alpha = 0.5) +
  theme_bw() +
  labs(x = expression(alpha[1] ~ "(tendência temporal)"),
       y = "Freq") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = c(0.2, 0.7)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 14)) +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14))

ggsave(plot = last_plot(),
       scale = 4,
       units = "cm",
       width = 5,
       height = 5,
       dpi = 120,
       "03_saida/tendencia_histograma.jpeg")



# plotar ranking das tendencias
tabela_r |>
  mutate(trend = case_when(prob_decrease > 0.95 ~ "declinio",
                           prob_increase > 0.95 ~ "aumento",
                           .default = "estavel")) |>
  mutate(trend = factor(trend, levels = c("aumento", "estavel", "declinio"))) |>
  arrange(desc(r_mean)) |>
  rowid_to_column() |>
  ggplot(aes(x = rowid, y = r_mean, color = trend)) +
  geom_hline(yintercept = 0, linewidth = 0.1, alpha = 0.3) +
  geom_vline(xintercept = nrow(tabela_r)/2, linewidth = 0.1, alpha = 0.3) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                  linewidth = 0.4, alpha = 1) +
  scale_color_manual(values = c("aumento" = "blue", "estavel" = "darkgrey", "declinio" = "red")) +
  labs(x = "População (ranking)",
       y = expression(alpha[1] ~ "(tendência temporal)")) +
  theme_bw() +
  theme(legend.position="none") +
  #ylim(-0.8,0.8) +
  #scale_x_continuous(limits = c(1,NA))
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size = 12)) +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14))

ggsave(plot = last_plot(),
       scale = 4,
       units = "cm",
       width = 6,
       height = 4,
       dpi = 120,
       "03_saida/tendencia_rank_plot.jpeg")


# plotar tendencias de populacoes em declinio
decliningPops <- tabela_r |>
  filter(prob_decrease > 0.9) |>
  distinct(uc,especie) |>
  mutate(decliningPops = paste(uc, especie)) |>
  pull(decliningPops)

npost |>
  filter(paste(uc, especie) %in% decliningPops) |>
  group_by(especie, ano) |>
  summarise(mean = mean(n, na.rm=TRUE),
            lower = quantile(n, prob = 0.1, na.rm=TRUE),
            upper = quantile(n, prob = 0.9, na.rm=TRUE)) |>
  ungroup() |>
  ggplot(aes(x = ano, y = mean)) +
  geom_line(aes(x = ano, y = mean), size = 0.3, alpha = 0.8) +
  geom_point(aes(x = ano, y = mean),  shape = 21,
             color = "white", stroke = 2,
             fill = "black", size = 2) +
  geom_linerange(aes(x = ano, y = mean,
                     ymin=lower, ymax=upper),
                 size = 0.2, alpha = 0.5) +
  #coord_cartesian(xlim = c(2014, 2025), expand = FALSE) +
  scale_x_continuous(breaks = scales::breaks_pretty(n=5)) +
  ylim(0,NA) +
  xlab("") +
  ylab("Abundância relativa") +
  #theme_classic() +
  theme_bw() +
  theme(axis.text = element_text(size = 12)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(~especie, scales = "free", ncol = 3) +
  theme(strip.text = element_text(size = 14))

# salvar figura
ggsave(plot = last_plot(),
       scale = 4,
       units = "cm",
       width = 8,
       height = 6,
       dpi = 120,
       "03_saida/pops_em_declinio.jpeg")



# explorar fatores potencialmente afetando tendencias (r = a1) ----------------------------------------


# relacao entre tendencia e abundancia inicial
tabela_r |>
  filter(n_inicial < 30) |> # retira outliers
  ggplot(aes(x = n_inicial, y = r_mean, color = especie)) +
  geom_point(size = 2, color = "grey30", alpha = 0.5) +
  #geom_count(color = "grey30", alpha = 0.5) +
  geom_smooth(method = "loess",
              span = 1.5,
              se = FALSE,
              color = "black", size = 0.2, alpha = 0.7) +
  theme_bw() +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14, angle = 0, vjust = 0.5)) +
  labs(x = "Abundância inicial", y = expression(alpha[1]))


# r e categoria de ameaca
tabela_r |>
  mutate(redlist = case_when(is.na(cat_ameaca) ~ "NA",
                             .default = cat_ameaca)) |>
  mutate(redlist = factor(redlist, levels = c("CR", "EN", "VU", "NA"))) |>
  ggplot(aes(x = redlist, y = r_mean)) +
  #geom_violin() +
  geom_boxplot(aes(fill = redlist)) +
  #geom_jitter(width = 0.02, size = 1, alpha = 0.4) +
  geom_hline(yintercept = mean(tabela_r$r_mean), linetype = "dashed", alpha = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Grau de ameaça", y = "Tendência") +
  guides(fill = guide_legend(title = "Categoria"), scale = "none")


# r e categoria de ameaca histogram splitted
background_data <- tabela_r |> select(-cat_ameaca)
tabela_r |>
  mutate(redlist = case_when(is.na(cat_ameaca) ~ "Nao-ameacada",
                             .default = "Ameacada")) |>
  ggplot(aes(x = r_mean, fill = redlist)) +
  geom_histogram(data = background_data, aes(x = r_mean),
                 fill = "grey", alpha = 0.5) +
  geom_histogram(colour = "black", linewidth = 0.2) +
  guides(fill = FALSE) +
  theme_bw() +
  facet_wrap(~ redlist) +
  theme(strip.text = element_text(size = 14)) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  labs(x = expression(beta[t] ~ "(Tendência)"), y = "Count")



# r e classe histogram splitted
background_data <- tabela_r |> select(-classe)
tabela_r |>
  ggplot(aes(x = r_mean, fill = classe)) +
  geom_histogram(data = background_data, aes(x = r_mean),
                 fill = "grey", alpha = 0.5) +
  geom_histogram(colour = "black", linewidth = 0.2) +
  guides(fill = FALSE, scale = "none") +
  theme_bw() +
  facet_wrap(~ classe) +
  theme(strip.text = element_text(size = 14)) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  labs(x = expression(alpha[t] ~ "(Tendência)"), y = "Count")



# r e categoria de uc histogram splitted
background_data <- tabela_r |> select(-categoria_uc)
tabela_r |>
  ggplot(aes(x = r_mean, fill = categoria_uc)) +
  geom_histogram(data = background_data, aes(x = r_mean),
                 fill = "grey", alpha = 0.5) +
  geom_histogram(colour = "black", linewidth = 0.2) +
  guides(fill = FALSE, scale = "none") +
  theme_bw() +
  facet_wrap(~ categoria_uc) +
  theme(strip.text = element_text(size = 14)) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  labs(x = expression(beta[t] ~ "(Temporal trend)"), y = "Count")


