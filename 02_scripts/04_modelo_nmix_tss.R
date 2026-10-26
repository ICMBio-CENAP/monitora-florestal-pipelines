# Programa Monitora
# Analise de dados
# Modelo multi-regiao, multi-especies, N-mixture
# tss (time for space substitution)
# Elildo Carvalho Jr @ICMBio/CENAP 21 out 2026

rm(list = ls())

#----- carregar pacotes
library(tidyverse)
library(jagsUI)
#library(abind)


# carregar dados previamente preparados
data_wide <- readRDS("01_entrada/data_long.rds")


# juntar dados para entrada no jags
jags_data <- list(y = data_wide$y,
                  pa = dense_rank(data_wide$uc),
                  pop = data_wide$pop,
                  site = data_wide$ea,
                  year = data_wide$ano-2013,
                  nobs = data_wide$nobs,
                  npa = data_wide$nuc,
                  npop = data_wide$npop,
                  nsite = data_wide$nea,
                  nyear = data_wide$nano
                  )
str(jags_data)


# escrever modelo bayesiano em linguagem jags
sink("02_scripts/tss.jags")
cat("
model {

##----- global-level priors

# abundance
mu.a0.global ~ dnorm(0,0.01)
mu.a1.global ~ dnorm(0,0.01)
sigma.a0.global ~ dunif(0,10)
sigma.a1.global ~ dunif(0,10)
tau.a0.global <- pow(sigma.a0.global,-2)
tau.a1.global <- pow(sigma.a1.global,-2)


# detection
mu.b0.global ~ dnorm(0, 0.01)
sigma.b0.global ~ dunif(0,10)
tau.b0.global <- pow(sigma.b0.global,-2)

#----- PA-level priors
    for (pa in 1:npa) {
    
      # abundance
      mu.a0.pa[pa] ~ dnorm(mu.a0.global,tau.a0.global)
      mu.a1.pa[pa] ~ dnorm(mu.a1.global,tau.a1.global)
      sigma.a0.pa[pa] ~ dunif(0,10)
      sigma.a1.pa[pa] ~ dunif(0,10)
      tau.a0.pa[pa] <- pow(sigma.a0.pa[pa],-2)
      tau.a1.pa[pa] <- pow(sigma.a1.pa[pa],-2)

      # detection
      mu.b0.pa[pa] ~ dnorm(mu.b0.global,tau.b0.global)
      sigma.b0.pa[pa] ~ dunif(0,10)
      tau.b0.pa[pa] <- pow(sigma.b0.pa[pa],-2)

      #----- population-level priors
      for (pop in 1:npop) {
      
      # inclusion prob
      #omega[pa,pop] ~ dunif(0,1)
      
      # abundance
      a0[pa,pop] ~ dnorm(mu.a0.pa[pa],tau.a0.pa[pa])
      a1[pa,pop] ~ dnorm(mu.a1.pa[pa],tau.a1.pa[pa])
      
      # zero-inflation
      #for (site in 1:nsite) {
      #  w[pa,pop,site] ~ dbern(omega[pa,pop])
      #}#site

      # detection
      b0[pa,pop] ~ dnorm(mu.b0.pa[pa],tau.b0.pa[pa])
      
      # random year effects on detection
      tau.year[pa,pop] <- 1 / (sd.year[pa,pop]*sd.year[pa,pop])
      sd.year[pa,pop] ~ dunif(0, 1)

      }#pop
    }#pa


##----- likelihood

# ecological submodel
for (pa in 1:npa) {
  for(pop in 1:npop){
      #for(site in 1:nsite){
        for(year in 1:nyear){
      
        log(lambda[pa,pop,year]) <- a0[pa,pop] + a1[pa,pop]*(year[year]-1) # minus 1 to set 1st to zero
        
        for(site in 1:nsite){
          n[pa,pop,site,year] ~ dpois(lambda[pa,pop,year])
        }#site

   }#year
  }#pop
}#pa


# observation submodel
for (pa in 1:npa) {
  for(pop in 1:npop){
    for(year in 1:nyear){

    eps_year[pa,pop,year] ~ dnorm(0, tau.year[pa,pop])
    logit(p[pa,pop,year]) <- b0[pa,pop] + eps_year[pa,pop,year]
    
    }#year
  }#pop
}#pa 

# observation process
for(i in 1:nobs){
  
  # observation model
  y[i] ~ dbin( p[pa[i],pop[i],year[i]], n[pa[i],pop[i],site[i],year[i]] )

}#i


}#model",fill = TRUE)
sink()



# parametros a monitorar
params <- c('lambda', 'n',  'p', 'a0', 'a1')


# inicializar abundancia n usando a maior contagem observada em cada ea (site)
# n[pa,pop,site,year] # dimensoes do objeto n)
y_max <- tibble(pa = data_wide$uc,
                pop = data_wide$pop,
                site = data_wide$ea,
                year = data_wide$ano,
                y = data_wide$y) |>
  group_by(pa, pop, site, year) |>
  summarise(y = max(y)) |>
  ungroup() |>
  complete(site,year,pop,
           nesting(pa),
           fill = list(y = 0)) |>
  arrange(pa,pop,site,year) |>
  print()

# transformar em array com dimensoes [pa,pop,site,year]
# usar funcao xtabs, que cria tabelas de contingencia ou matrizes esparsas usando formulas
# sintaxe do xtabs: y ~ dim1 + dim2 + dim3 + dim4
y_array <- xtabs(y ~ pa + pop + site + year, data = y_max)
dim(y_array) # checar
y_array[1,1,,] # checar tdm, pop 1
y_array[1,2,,] # tdm, pop 2
y_array[36,1,,] # cazumba, pop 1
y_array[36,2,,] # cazumba, pop 2


inits <- function() list(n = y_array)
inits()


# parametros mcmc
ni <- 500000
nt <- 3000
nb <- 200000
nc <- 3


# chamar o jags
system.time(mod <- jagsUI::jags(data = jags_data,
                                inits = inits,
                                parameters.to.save = params,
                                model.file = "02_scripts/tss.jags",
                                n.chains = nc,
                                n.thin = nt,
                                n.iter = ni,
                                n.burnin = nb,
                                parallel = TRUE))

# salvar
saveRDS(mod, "03_saida/res_tss_multiregion.rds")


