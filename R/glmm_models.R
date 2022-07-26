library(tidyverse)
library(here)
library(glmmTMB)
library(bbmle) # for AICtab
library(lme4)
options(scipen=999)

# load data ---------------------------------------------------------------
path <- here("data", "foxtall.csv")
df <- readr::read_csv(path, col_types = "_ff_fffi_dddd", na = "NaN") %>% 
  dplyr::mutate(genotype = fct_relevel(genotype, "WT", "HET"))


# Make DV specific datasets  ----------------------------------------------
df_sr <- df %>% dplyr::select(1:7) %>% 
  dplyr::filter(!is.na(surfacerighting))

df_ng <- df %>% dplyr::select(1:6, 8) %>% 
  dplyr::filter(!is.na(negativegeotaxis))

df_ca <- df %>% dplyr::select(1:6, 9) %>% 
  dplyr::filter(!is.na(cliffaversion))

df_fg <- df %>% dplyr::select(1:6, 10) %>% 
  dplyr::filter(!is.na(forelimbgrasp))


# Surface righting models -------------------------------------------------
sr_fit1 <- glmmTMB(surfacerighting ~ genotype * sex * pnd + (pnd|litter/aididx),
                  family = Gamma(link = "log"), data = df_sr)

summary(sr_fit1)

sr_fit2 <- glmmTMB(surfacerighting ~ genotype * sex * pnd + (pnd|litter/aididx),
                   family = Gamma(link = "log"), 
                   data = (df_sr %>% dplyr::mutate(sex = fct_relevel(sex, "Male"))))

summary(sr_fit2)

# best fitting model
sr_fit3 <- glmmTMB(surfacerighting ~ genotype * sex * pnd +  
                     (pnd|aididx) + (1|experimenter),
                   family = Gamma(link = "log"), data = df_sr)

summary(sr_fit3)

AICctab(sr_fit1, sr_fit2, sr_fit3)

# Negative geotaxis models ------------------------------------------------

# models with 'pnd' in random slopes did not converge. There is a possible 
# experimenter influence on the results

ng_fit1 <- glmmTMB(negativegeotaxis ~ genotype * sex * pnd + (1|aididx),
                   family = Gamma(link = "log"), data = df_ng)

summary(ng_fit1)

# best fitting model
ng_fit2 <- glmmTMB(negativegeotaxis ~ genotype * sex * pnd + (1|aididx) + (1|experimenter),
                   family = Gamma(link = "log"), data = df_ng)

summary(ng_fit2)

AICctab(ng_fit1, ng_fit2)


# Cliff aversion models ---------------------------------------------------
ca_fit1 <- glmmTMB(cliffaversion ~ genotype * sex * pnd + (pnd|aididx),
                   family = Gamma(link = "log"), data = df_ca)

summary(ca_fit1)

# best fitting model
ca_fit2 <- glmmTMB(cliffaversion ~ genotype * sex * pnd + (pnd|aididx) + (1|experimenter),
                   family = Gamma(link = "log"), data = df_ca)

summary(ca_fit2)

AICctab(ca_fit1, ca_fit2)

g_ca <- ggplot(df_ca, aes(x = pnd, y = cliffaversion, colour = sex, shape = genotype)) + 
                 geom_point() + geom_jitter() +
  stat_summary(fun.data = "mean_cl_boot", colour = "black")
g_ca


# Forelimb grasp models ---------------------------------------------------
# best fitting model
fg_fit1 <- glmmTMB(forelimbgrasp ~ genotype * sex * pnd + (1|aididx),
                   family = Gamma(link = "log"), data = df_fg)

summary(fg_fit1)

fg_fit2 <- glmmTMB(forelimbgrasp ~ genotype * sex * pnd + (1|aididx),
                   family = gaussian(link = "log"), data = df_fg)

summary(fg_fit2)

fg_fit3 <- glmmTMB(forelimbgrasp ~ (sex + genotype) * pnd  + (1|aididx),
                   family = Gamma(link = "log"), data = df_fg)

summary(fg_fit3)

AICctab(fg_fit1, fg_fit2, fg_fit3)

g_fg <- ggplot(df_fg, aes(x = pnd, y = forelimbgrasp, colour = sex, shape = genotype)) + 
  geom_point() + geom_jitter() +
  stat_summary(fun.data = "mean_cl_boot", colour = "black")
g_fg
