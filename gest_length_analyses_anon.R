
# March 13, 2022
# code to produce figures and final results tables from anonymized 
# publication data

library(lmerTest)
library(visreg)
library(MuMIn)
library(dplyr)
library(ggplot2)
library(car)
library(sjPlot)
library(lubridate)
library(ggpubr)

all_gestations = readRDS("gestation_data.rds")

multipar_gests = 
  all_gestations %>%
  filter(twin == 0, firstborn == 0)

# ----------------- gestation length analysis -------------------------------------------------------

# some quick exploratory analyses

m0_matage = lmer(gest_length ~ scale(mat_age) + I(scale(mat_age)^2) + (1|mat_id_anon), data = multipar_gests)
visreg(fit = m0_matage, xvar = "mat_age", xlab = "maternal age", ylab = "predicted gest length")

m0_tri3 = lmer(gest_length ~ wet_szn_3rd_trim + (1|mat_id_anon), data = multipar_gests) 
visreg(fit = m0_tri3, xvar = "wet_szn_3rd_trim")                 

m0_tri1 = lmer(gest_length ~ wet_szn_1st_trim + (1|mat_id_anon), data = multipar_gests) 
visreg(fit = m0_tri1, xvar = "wet_szn_1st_trim")   

gest_length_vif = lm(gest_length ~  
               mat_age_z
             + I(mat_age_z^2)
             + birth_order
             + wet_szn_3rd_trim
             + wet_szn_1st_trim
             + scale(igi)
             + last_inf_died
             + inf_sex
             + provisioned, 
             data = multipar_gests)
vif_res = vif(gest_length_vif)
vif_res

plot(multipar_gests$mat_age_z ~ multipar_gests$birth_order, 
     main = "Maternal age vs Birth Order",
     xlab = "Birth order", ylab = "Maternal age")

plot(multipar_gests$wet_szn_3rd_trim ~ multipar_gests$wet_szn_1st_trim, 
     main = "First vs third trimester in wet season")


# fit model (with REML):

gest_length_model1 = lmer(gest_length ~ 
                    mat_age_z + 
                    I(mat_age_z^2) + 
                    birth_order + 
                    wet_szn_3rd_trim + 
                    wet_szn_1st_trim +
                    scale(igi) +
                    last_inf_died + 
                    inf_sex + 
                    provisioned +
                    (1|mat_id_anon), 
                  data = multipar_gests, REML = T, na.action = na.fail)
summary(gest_length_model1)

# fit model (no REML for model comparisons):

gest_length_model2 = lmer(gest_length ~ 
                    mat_age_z + 
                    I(mat_age_z^2) + 
                    birth_order + 
                    wet_szn_3rd_trim + 
                    wet_szn_1st_trim +
                    scale(igi) +
                    last_inf_died + 
                    inf_sex + 
                    provisioned +
                    (1|mat_id_anon), 
                  data = multipar_gests, REML = F, na.action = na.fail)
summary(gest_length_model2)

# ranking models by AICc using Maximum Likelihood:

getAllTerms(gest_length_model2)

all_m_gest_length = dredge(gest_length_model2, rank = AICc, trace = 2, 
                    subset = (!(mat_age_z && birth_order) & 
                                !(`I(mat_age_z^2)` && !mat_age_z) &
                                !(wet_szn_3rd_trim && wet_szn_1st_trim)# &
                              # !(last_inf_died && `scale(igi)`)
                    ) #, 
                    # REML = FALSE 
)

all_m_gest_length = subset(all_m_gest_length, !nested(.))

best_m_gest_length = get.models(all_m_gest_length, subset = delta<6)

model.sel(best_m_gest_length)

model.avg(best_m_gest_length)

gest_length_mod_sm = lmer(gest_length ~ scale(mat_age) + I(scale(mat_age)^2) + igi + (1|mat_id_anon), 
                  data = multipar_gests, REML = T, na.action = na.fail)

r.squaredGLMM(gest_length_mod_sm)

# double check that the high R2c value not driven by females with single 
# offspring in sample: 

multipar_gests %>%
  count(mat_id_anon) %>%
  pull(n) %>%
  table()

temp_gest_length_mod_sm = 
  multipar_gests %>%
  group_by(mat_id_anon) %>%
  mutate(nobs = n()) %>%
  as.data.frame() %>%
  filter(nobs > 1) %>%
  lmer(gest_length ~ scale(mat_age) + I(scale(mat_age)^2) + igi + (1|mat_id_anon), 
       data = ., REML = T, na.action = na.fail) 

summary(temp_gest_length_mod_sm)

r.squaredGLMM(temp_gest_length_mod_sm)


theme_jtf = function (base_size = 11, base_family = "") {
  theme_linedraw() %+replace% 
    theme(
      panel.grid.major  = element_line(color = "gray94", size = 0.2),
      panel.grid.minor  = element_line(color = "white")
    )
}

theme_set(theme_jtf())


gest_length_mod_sm = lmer(gest_length ~ scale(mat_age) + I(scale(mat_age)^2) + igi + (1|mat_id_anon), 
                  data = multipar_gests, REML = T, na.action = na.fail)


plot_model(gest_length_mod_sm, type = "pred", terms = "mat_age [all]", show.data = T) +
  xlab("Maternal age (years)") +
  ylab("Gestation length (days)") +
  ggtitle(label = NULL) + 
  ylim(205, 255)


gest_length_mod_sm = lmer(gest_length ~ mat_age_z + I(mat_age_z^2) + igi + (1|mat_id_anon),
                  data = multipar_gests, REML = T, na.action = na.fail)

plot_model(gest_length_mod_sm, type = "pred", terms = "igi [all]", show.data = T) +
  xlab("Inter-gestational interval (years)") +
  ylab("Gestation length (days)") +
  ggtitle(label = NULL) +
  ylim(205, 255)


# ------------- double check about female rank -----------------

head(multipar_gests)

sum(!is.na(multipar_gests$mat_rank_concept))

nrow(multipar_gests)

multipar_gests3 =
  multipar_gests %>% 
  mutate(igiz = scale(igi),
         mat_age_z2 = mat_age_z^2) %>%
  filter(!is.na(mat_rank_concept))

nrow(multipar_gests3)

gest_length_model3 = lmer(gest_length ~ 
                    mat_age_z + 
                    mat_age_z2 + 
                    birth_order + 
                    wet_szn_3rd_trim +
                    # wet_szn_1st_trim +
                    igiz +
                    last_inf_died + 
                    inf_sex + 
                    # provisioned +
                    scale(mat_rank_concept) +
                    (1|mat_id_anon), 
                  data = multipar_gests3, REML = F, na.action = na.fail)
summary(gest_length_model3)

# Since we had some trouble with convergence doing a full model comparison, we'll  
# just compare the best model for gest length to that model plus maternal rank 
# at conception

dim(multipar_gests3)

gest_mod_sm = lmer(gest_length ~ 
                     mat_age_z 
                   + mat_age_z2
                   + igiz
                   + (1|mat_id_anon),
                   data = multipar_gests3, REML = F, na.action = na.fail)

AICc(gest_mod_sm)

gest_mod_sm2 = lmer(gest_length ~ 
                      mat_age_z 
                    + mat_age_z2
                    + igiz
                    + scale(mat_rank_concept)
                    + (1|mat_id_anon),
                    data = multipar_gests3, REML = F, na.action = na.fail)

summary(gest_mod_sm2)

gest_mod_sm3 = lmer(gest_length ~ 
                      mat_age_z 
                    + mat_age_z2
                    # + igiz
                    + scale(mat_rank_concept)
                    + (1|mat_id_anon),
                    data = multipar_gests3, REML = F, na.action = na.fail)

summary(gest_mod_sm3)

model.sel(gest_mod_sm, gest_mod_sm2, gest_mod_sm3) 

# =============================== Survival analysis =============================

# March 13, 2022
# clean code to produce figures and final results tables

library(survival)
library(coxme)
library(dplyr)
library(lubridate)
library(lmerTest)
library(MuMIn)
library(ggplot2)
library(survminer)
library(sjPlot)

all_gestations = readRDS("gestation_data.rds")

survdat = all_gestations


# this analysis does not exclude victims of infanticide/kidnapping
# or twins.  Found in supplementary materials and Table S3

nrow(survdat) # 53

cox.zph(coxph(Surv(time = lifespan, 
                   event = departtype == "D") ~ 
                inf_sex + 
                scale(mat_age) + 
                I(scale(mat_age)^2) +
                scale(gest_length) + 
                I(scale(gest_length)^2) +
                firstborn +
                orphan + 
                frailty(mat_id_anon), 
              data = survdat,
              na.action = na.fail))

survdat %>% 
  filter(orphan == 1)
# 3 offspring, all die before age 2, easiest to just leave them out

survdat = 
  survdat %>%
  filter(orphan == 0)

nrow(survdat)

res_cox = coxph(Surv(time = lifespan, 
                     event = departtype == "D") ~ 
                  inf_sex + 
                  scale(mat_age) + 
                  I(scale(mat_age)^2) +
                  scale(gest_length) + 
                  I(scale(gest_length)^2) +
                  firstborn + 
                  frailty(mat_id_anon), 
                data = survdat,
                na.action = na.fail)

cox.zph(res_cox)

ggcoxzph(cox.zph(res_cox))

# ok, we'll stratify the by infant sex to ensure data fit PH assumptions:

dim(survdat)

table(survdat$inf_sex)

# having a stratum with one individual shouldn't influence estimates, 
# effective sample size is 49 anyway, so we'll just remove the one 
# inf_sex == 0 (i.e. unknown sex) individual

survdat =
  survdat %>% 
  filter(inf_sex != 0)

res_cox = coxph(Surv(time = lifespan, 
                     event = departtype == "D") ~ 
                  strata(inf_sex) + 
                  scale(mat_age) + 
                  I(scale(mat_age)^2) +
                  scale(gest_length) + 
                  I(scale(gest_length)^2) +
                  firstborn + 
                  twin +
                  frailty(mat_id_anon), 
                data = survdat,
                na.action = na.fail)

cox.zph(res_cox)

plot(cox.zph(res_cox))

ggcoxzph(cox.zph(res_cox))

# OK looks good

# --------------- full data survival analysis ---------------------

# to make model comparisons simpler we'll define a new 
# variable for squared maternal age and squared gest length

survdat2 = 
  survdat %>%
  mutate(mat_age_z = scale(mat_age),
         mat_age_z2 = mat_age_z^2,
         gest_length_z = scale(gest_length),
         gest_length_z2 = gest_length_z^2) 

# note that departtype is coded as "D" for death, "O" for still 
# alive at latest timepoint, and "P" for permanent disappearance.
# we're considering the two "P" females to be right-censored 
# because they disappeared at a time when they could plausibly have 
# emmigrated.  See manuscript for more details.

surv_m_full2 = coxme(Surv(time = lifespan, 
                       event = departtype == "D") ~ 
                    strata(inf_sex) +
                    mat_age_z + 
                    mat_age_z2 +
                    gest_length_z + 
                    gest_length_z2 +
                    firstborn +
                    twin +
                    (1|mat_id_anon), 
                  data = survdat2,
                  na.action = na.fail)
summary(surv_m_full2)

MuMIn::getAllTerms(surv_m_full2)

all_surv_m_full = dredge(
  surv_m_full2, rank = AICc, trace = T, 
  subset = (!(mat_age_z2 && !mat_age_z) & 
              !(gest_length_z2 && !gest_length_z))
)

all_surv_m_full = subset(all_surv_m_full, !nested(.))

best_surv_m_full = get.models(all_surv_m_full, subset = delta < 6)

model.sel(best_surv_m_full) 

model.avg(best_surv_m_full)

# ------------------- Main survival analysis --------------------

survdat = all_gestations

survdat = 
  survdat %>% 
  filter(unnat_death == 0, twin == 0)

# this analysis excludes victims of infanticide/kidnapping
# and twins.  Produces results found in Table 3

nrow(survdat) # 42

cox.zph(coxph(Surv(time = lifespan, 
                   event = departtype == "D") ~ 
                inf_sex + 
                scale(mat_age) + 
                I(scale(mat_age)^2) +
                scale(gest_length) + 
                I(scale(gest_length)^2) +
                firstborn +
                orphan + 
                frailty(mat_id_anon), 
              data = survdat,
              na.action = na.fail))

survdat %>% 
  filter(orphan == 1)
# 3 offspring, all die before age 2, we will leave them out

survdat = 
  survdat %>%
  filter(orphan == 0)

nrow(survdat)

surv_mod_init = coxph(Surv(time = lifespan, 
                     event = departtype == "D") ~ 
                  inf_sex + 
                  scale(mat_age) + 
                  I(scale(mat_age)^2) +
                  scale(gest_length) + 
                  I(scale(gest_length)^2) +
                  firstborn + 
                  frailty(mat_id_anon), 
                data = survdat,
                na.action = na.fail)

cox.zph(surv_mod_init)

ggcoxzph(cox.zph(surv_mod_init))

# ok, we'll stratify the by infant sex to ensure data fit PH assumptions:

dim(survdat)

table(survdat$inf_sex)

# having a stratum with one individual shouldn't influence estimates, 
# effective sample size is 49 anyway, so we'll just remove the one 
# inf_sex == 0 (i.e. unknown sex) individual

survdat =
  survdat %>% 
  filter(inf_sex != 0)

surv_mod_init = coxph(Surv(time = lifespan, 
                     event = departtype == "D") ~ 
                  strata(inf_sex) + 
                  scale(mat_age) + 
                  I(scale(mat_age)^2) +
                  scale(gest_length) + 
                  I(scale(gest_length)^2) +
                  firstborn + 
                  frailty(mat_id_anon), 
                data = survdat,
                na.action = na.fail)

cox.zph(surv_mod_init)

plot(cox.zph(surv_mod_init))

ggcoxzph(cox.zph(surv_mod_init))

# ------------------- run survival models (main analysis) -------------------------

# to simplify model comparison we'll create new scaled age/ gestation variables
survdat2 =
  survdat %>%
  mutate(mat_age_z = scale(mat_age),
         mat_age_z2 = mat_age_z^2,
         gest_length_z = scale(gest_length),
         gest_length_z2 = gest_length_z^2) 

surv_m_full2 = coxme(Surv(time = lifespan, 
                       event = departtype == "D") ~ 
                    strata(inf_sex) + 
                    mat_age_z + 
                    mat_age_z2 +
                    gest_length_z + 
                    gest_length_z2 +
                    firstborn +
                    (1|mat_id_anon), 
                  data = survdat2,
                  na.action = na.fail)
summary(surv_m_full2)

all_surv_m_main = dredge(surv_m_full2, rank = AICc, trace = 2, 
                      subset = (!(mat_age_z2 && !mat_age_z) & 
                                  !(gest_length_z2 && !gest_length_z))
)

all_surv_m_main = subset(all_surv_m_main, !nested(.))

best_surv_m_main = get.models(all_surv_m_main, subset = delta < 6)

model.sel(best_surv_m_main) 

model.avg(best_surv_m_main)

# ------------- plot lifetime survival (main survival analysis) -------------------

best_surv_main_mod = coxph(Surv(time = lifespan, 
                           event = departtype == "D") ~ 
                        strata(inf_sex) + 
                        gest_length + 
                        firstborn, 
                      data = survdat2,
                      na.action = na.fail)

summary(best_surv_main_mod)

survdat2 %>% count(inf_sex, departtype) # only 2 of 15 females died, 14 of 23 males died

fb_df = data.frame(firstborn = c(1, 0),
                   gest_length = mean(survdat2$gest_length),
                   inf_sex = 1)

fit = survfit(best_surv_main_mod, newdata = fb_df)

ggsurvplot(fit = fit, data = fb_df, conf.int = T, censor = FALSE,
           main = "Survival by infant firstborn status", 
           xlab = "Offspring age (years)",
           ggtheme = theme_jtf(),
           legend.title = "Firstborn status:",
           legend.labs = c("Yes", "No")) +
  labs(tag = "A")

gest_qs = c(218, 237)

gest_df = data.frame(firstborn = 0,
                     gest_length = unname(gest_qs),
                     inf_sex = 1)

fit = survfit(best_surv_main_mod, newdata = gest_df)

ggsurvplot(fit = fit, data = gest_df, conf.int = T, censor = FALSE,
           main = "Survival by gestation length", 
           xlab = "Offspring age (years)",
           # palette = "RdYlBu",
           ggtheme = theme_jtf(),
           legend.title = "Gestation length:",
           legend.labs = c("218 days", "237 days")) + 
  labs(tag = "C")

sex_df = data.frame(firstborn = 0,
                    gest_length = mean(survdat2$gest_length),
                    inf_sex = c(1, -1))

fit = survfit(best_surv_main_mod, newdata = sex_df)

ggsurvplot(fit = fit, data = sex_df, conf.int = T, censor = TRUE,
           main = "Survival by infant firstborn status", 
           xlab = "Offspring age (years)",
           ggtheme = theme_jtf(),
           legend.title = "Offspring sex:",
           legend.labs = c("Male", "Female")) + 
  labs(tag = "B")




#  ----------------- only offspring with known maternal rank at conception ----------------

survdat3 = 
  survdat2 %>%
  filter(!is.na(mat_rank_concept))

# now 28 offspring

survdat3 %>%
  count(inf_sex, right_cens)

table(survdat3$right_cens)

ranksurv2 = 
  coxme(Surv(time = lifespan, 
             event = departtype == "D") ~ 
          strata(inf_sex)
        + mat_age_z 
        + mat_age_z2
        + mat_rank_concept*gest_length_z
        + gest_length_z2
        +firstborn
        + (1|mat_id_anon)
        , 
        data = survdat3,
        na.action = na.fail)

all_mod_rank = dredge(ranksurv2, rank = AICc, trace = 2, 
                      subset = (!(mat_age_z2 && !mat_age_z) & 
                                  !(gest_length_z2 && !gest_length_z))
)

all_mod_rank = subset(all_mod_rank, !nested(.))

best_mod_rank = get.models(all_mod_rank, subset = delta < 6)

model.sel(best_mod_rank) 

model.avg(best_mod_rank)

# ------------- survival by firstborn status in bigger sample -----------------

fb2 = readRDS("big_firstborn_data_anon.rds")

fitfb = coxme(Surv(time = longev, 
                   event = (departtype == "D")) ~ 
                sex + 
                firstborn + 
                (1|momid_anon), 
              data = fb2,
              na.action = na.fail)

summary(fitfb)

# plot:

fitfb2 = coxph(Surv(time = longev, 
                    event = (departtype == "D")) ~ 
                 sex + 
                 firstborn, 
               data = fb2,
               na.action = na.fail)
summary(fitfb2)

fb_df = data.frame(firstborn = c(TRUE, FALSE),
                   sex = 0)
# fb_df

fit = survfit(fitfb2, newdata = fb_df)

a = ggsurvplot(fit = fit, data = fb_df, conf.int = T, censor = FALSE,
               main = "Survival by infant firstborn status", 
               xlab = "Offspring age (years)",
               # palette = "Spectral",
               ggtheme = theme_jtf(),
               legend.title = "Firstborn:",
               legend.labs = c("True", "False")) +
  labs(tag = "A")

sex_df = data.frame(firstborn = FALSE,
                    sex = c(1, -1))

fit = survfit(fitfb2, newdata = sex_df)

b = ggsurvplot(fit = fit, data = sex_df, conf.int = T, censor = FALSE,
               main = "Survival by infant firstborn status", 
               xlab = "Offspring age (years)",
               # palette = "Spectral",
               ggtheme = theme_jtf(),
               legend.title = "Offspring sex:",
               legend.labs = c("Male", "Female")) +
  labs(tag = "B")

arrange_ggsurvplots(list(a, b), ncol = 2, nrow = 1, title = "Predicted Infant Survival")
