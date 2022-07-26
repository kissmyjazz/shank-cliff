library(tidyverse)
library(here)
library(glmmTMB)
library(bbmle) # for AICtab

# load data ---------------------------------------------------------------
path <- here("data", "foxtall.csv")
df <- readr::read_csv(path, col_types = "_ff_fffi_dddd", na = "NaN")


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
sr_fit <- glmmTMB(surfacerighting ~ genotype * sex * pnd + (pnd|litter/aididx),
                  family = Gamma(link = "log"), data = df_sr)

summary(sr_fit)
