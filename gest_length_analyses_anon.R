
# Sept 21, 2021
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

all_gests = readRDS("gestation_data.rds")

multipar = 
  all_gests %>%
  filter(twin == 0, firstborn == 0)

# ----------------- gestation length analysis -------------------------------------------------------

# some quick exploratory analyses

m0_matage = lmer(gest_length ~ scale(mat_age) + I(scale(mat_age)^2) + (1|mat_id_anon), data = multipar)
visreg(fit = m0_matage, xvar = "mat_age", xlab = "maternal age", ylab = "predicted gest length")

m0_tri3 = lmer(gest_length ~ wet_szn_3rd_trim + (1|mat_id_anon), data = multipar) 
visreg(fit = m0_tri3, xvar = "wet_szn_3rd_trim")                 

m0_tri1 = lmer(gest_length ~ wet_szn_1st_trim + (1|mat_id_anon), data = multipar) 
visreg(fit = m0_tri1, xvar = "wet_szn_1st_trim")   

igi_vif = lm(gest_length ~  
               mat_age_z
             + I(mat_age_z^2)
             + birth_order
             + wet_szn_3rd_trim
             + wet_szn_1st_trim
             + scale(igi)
             + last_inf_died
             + inf_sex
             + provisioned, 
             data = multipar)
vif_res = vif(igi_vif)
vif_res

plot(multipar$mat_age_z ~ multipar$birth_order, 
     main = "Maternal age vs Birth Order",
     xlab = "Birth order", ylab = "Maternal age")

plot(multipar$wet_szn_3rd_trim ~ multipar$wet_szn_1st_trim, 
     main = "First vs third trimester in wet season")


# fit model (with REML):

igi_model1 = lmer(gest_length ~ 
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
                  data = multipar, REML = T, na.action = na.fail)
summary(igi_model1)

# fit model (no REML for model comparisons):

igi_model2 = lmer(gest_length ~ 
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
                  data = multipar, REML = F, na.action = na.fail)
summary(igi_model2)

# ranking models by AICc using Maximum Likelihood:

getAllTerms(igi_model2)

all_m_gest = dredge(igi_model2, rank = AICc, trace = 2, 
                    subset = (!(mat_age_z && birth_order) & 
                                !(`I(mat_age_z^2)` && !mat_age_z) &
                                !(wet_szn_3rd_trim && wet_szn_1st_trim)# &
                              # !(last_inf_died && `scale(igi)`)
                    ) #, 
                    # REML = FALSE 
)

all_m_gest = subset(all_m_gest, !nested(.))

best_m_gest = get.models(all_m_gest, subset = delta<6)

model.sel(best_m_gest)

model.avg(best_m_gest)

igi_mod_sm = lmer(gest_length ~ scale(mat_age) + I(scale(mat_age)^2) + igi + (1|mat_id_anon), 
                  data = multipar, REML = T, na.action = na.fail)

r.squaredGLMM(igi_mod_sm)

# double check that the high R2c value not driven by females with single 
# offspring in sample: 

multipar %>%
  count(mat_id_anon) %>%
  pull(n) %>%
  table()

temp_igi_mod_sm = 
  multipar %>%
  group_by(mat_id_anon) %>%
  mutate(nobs = n()) %>%
  as.data.frame() %>%
  filter(nobs > 1) %>%
  lmer(gest_length ~ scale(mat_age) + I(scale(mat_age)^2) + igi + (1|mat_id_anon), 
       data = ., REML = T, na.action = na.fail) 

summary(temp_igi_mod_sm)

r.squaredGLMM(temp_igi_mod_sm)


theme_jtf = function (base_size = 11, base_family = "") {
  theme_linedraw() %+replace% 
    theme(
      panel.grid.major  = element_line(color = "gray94", size = 0.2),
      panel.grid.minor  = element_line(color = "white")
    )
}

theme_set(theme_jtf())


igi_mod_sm = lmer(gest_length ~ scale(mat_age) + I(scale(mat_age)^2) + igi + (1|mat_id_anon), 
                  data = multipar, REML = T, na.action = na.fail)


plot_model(igi_mod_sm, type = "pred", terms = "mat_age [all]") +
  xlab("Maternal age (years)") +
  ylab("Gestation length (days)") +
  ggtitle(label = NULL) + 
  ylim(205, 255)

# plot_model(igi_mod_sm, type = "pred", terms = "igi [all]") +
#   xlab("Inter-gestational interval (years)") +
#   ylab("Gestation length (days)")

igi_mod_sm = lmer(gest_length ~ mat_age_z + I(mat_age_z^2) + igi + (1|mat_id_anon),
                  data = multipar, REML = T, na.action = na.fail)

plot_model(igi_mod_sm, type = "pred", terms = "igi [all]") +
  xlab("Inter-gestational interval (years)") +
  ylab("Gestation length (days)") +
  ggtitle(label = NULL) +
  ylim(205, 255)


# ------------- double check about female rank -----------------

head(multipar)

sum(!is.na(multipar$mat_rank_concept))

nrow(multipar)

multipar3 =
  multipar %>% 
  mutate(igiz = scale(igi),
         mat_age_z2 = mat_age_z^2) %>%
  filter(!is.na(mat_rank_concept))

nrow(multipar3)

igi_model3 = lmer(gest_length ~ 
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
                  data = multipar3, REML = F, na.action = na.fail)
summary(igi_model3)

# Since we had some trouble with convergence doing a full model comparison, we'll  
# just compare the best model for gest length to that model plus maternal rank 
# at conception

dim(multipar3)

gest_mod_sm = lmer(gest_length ~ 
                     mat_age_z 
                   + mat_age_z2
                   + igiz
                   + (1|mat_id_anon),
                   data = multipar3, REML = F, na.action = na.fail)

AICc(gest_mod_sm)

gest_mod_sm2 = lmer(gest_length ~ 
                      mat_age_z 
                    + mat_age_z2
                    + igiz
                    + scale(mat_rank_concept)
                    + (1|mat_id_anon),
                    data = multipar3, REML = F, na.action = na.fail)

summary(gest_mod_sm2)

model.sel(gest_mod_sm, gest_mod_sm2, gest_mod_sm3) 

# =============================== Survival analysis =============================

# Sept 21, 2021
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

all_gests = readRDS("gestation_data.rds")

surv2 = all_gests


nrow(surv2) # 53

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
              data = surv2,
              na.action = na.fail))

surv2 %>% 
  filter(orphan == 1)
# 3 offspring, all die before age 2, easiest to just leave them out

surv2 = 
  surv2 %>%
  filter(orphan == 0)

nrow(surv2)

res_cox = coxph(Surv(time = lifespan, 
                     event = departtype == "D") ~ 
                  inf_sex + 
                  scale(mat_age) + 
                  I(scale(mat_age)^2) +
                  scale(gest_length) + 
                  I(scale(gest_length)^2) +
                  firstborn + 
                  frailty(mat_id_anon), 
                data = surv2,
                na.action = na.fail)

cox.zph(res_cox)

ggcoxzph(cox.zph(res_cox))

# ok, we'll stratify the by infant sex to ensure data fit PH assumptions:

dim(surv2)

table(surv2$inf_sex)

# having a stratum with one individual shouldn't influence estimates, 
# effective sample size is 49 anyway, so we'll just remove the one 
# inf_sex == 0 (i.e. unknown sex) individual

surv2 =
  surv2 %>% 
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
                data = surv2,
                na.action = na.fail)

cox.zph(res_cox)

plot(cox.zph(res_cox))

ggcoxzph(cox.zph(res_cox))

# OK looks good

# --------------- full data survival analysis ---------------------

# to make model comparisons simpler we'll define a new 
# variable for squared maternal age and squared gest length

surv4 = 
  surv2 %>%
  mutate(mat_age_z = scale(mat_age),
         mat_age_z2 = mat_age_z^2,
         gest_length_z = scale(gest_length),
         gest_length_z2 = gest_length_z^2) 

# note that departtype is coded as "D" for death, "O" for still 
# alive at latest timepoint, and "P" for permanent disappearance.
# we're considering the two "P" females to be right-censored 
# because they disappeared at a time when they could plausibly have 
# emmigrated.  See manuscript for more details.

cox_full2 = coxme(Surv(time = lifespan, 
                       event = departtype == "D") ~ 
                    strata(inf_sex) +
                    mat_age_z + 
                    mat_age_z2 +
                    gest_length_z + 
                    gest_length_z2 +
                    firstborn +
                    twin +
                    (1|mat_id_anon), 
                  data = surv4,
                  na.action = na.fail)
summary(cox_full2)

MuMIn::getAllTerms(cox_full2)

all_mod_full = dredge(
  cox_full2, rank = AICc, trace = T, 
  subset = (!(mat_age_z2 && !mat_age_z) & 
              !(gest_length_z2 && !gest_length_z))
)

all_mod_full = subset(all_mod_full, !nested(.))

best_mod_full = get.models(all_mod_full, subset = delta < 6)

model.sel(best_mod_full) 

model.avg(best_mod_full)

# ---------------------- plot best model lifetime survival full dataset --------------

# best model

best_mod_full = coxph(Surv(time = lifespan, 
                           event = departtype == "D") ~ 
                        strata(inf_sex) + 
                        gest_length + 
                        mat_age_z +
                        I(mat_age_z^2), 
                      data = surv4,
                      na.action = na.fail)
summary(best_mod_full)

# ok thing to remember is that only 4 females died
surv4 %>% count(inf_sex, right_cens)

age_qs = quantile(surv4$mat_age_z, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

age_df = data.frame(mat_age_z = unname(age_qs),
                    gest_length = mean(surv4$gest_length),
                    inf_sex = 1)

fit = survfit(best_mod_full, newdata = age_df)

ggsurvplot(fit = fit, data = age_df, conf.int = F, censor = FALSE,
               main = "Survival by maternal age",
               xlab = "Offspring age (years)",
               ggtheme = theme_jtf(),
               legend.title = "Mom age:",
               legend.labs = names(age_qs)) + 
  labs(tag = "A")

gest_qs = quantile(surv4$gest_length, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

gest_df = data.frame(mat_age_z = 0,
                     gest_length = unname(gest_qs),
                     inf_sex = 1)

fit = survfit(best_mod_full, newdata = gest_df)

ggsurvplot(fit = fit, data = gest_df, conf.int = F, censor = FALSE,
               main = "Survival by gestation length", 
               xlab = "Offspring age (years)",
               # palette = "RdYlBu",
               ggtheme = theme_jtf(),
               legend.title = "Gest length:",
               legend.labs = names(gest_qs)) + # or maybe just do high vs low...
  labs(tag = "C")

sex_df = data.frame(mat_age_z = 0,
                    gest_length = mean(surv4$gest_length),
                    inf_sex = c(1, -1))

fit = survfit(best_mod_full, newdata = sex_df)

ggsurvplot(fit = fit, data = sex_df, conf.int = T, censor = T,
               main = "Survival by infant firstborn status", 
               xlab = "Offspring age (years)",
               ggtheme = theme_jtf(),
               legend.title = "Offspring sex:",
               legend.labs = c("Male", "Female")) + 
  labs(tag = "B")

# ------------------- remove twin births -------------------------

surv4 = 
  surv4 %>%
  filter(twin == 0)

cox_full2 = coxme(Surv(time = lifespan, 
                       event = departtype == "D") ~ 
                    strata(inf_sex) + 
                    mat_age_z + 
                    mat_age_z2 +
                    gest_length_z + 
                    gest_length_z2 +
                    firstborn +
                    # twin +
                    (1|mat_id_anon), 
                  data = surv4,
                  na.action = na.fail)
summary(cox_full2)

all_mod_full = dredge(cox_full2, rank = AICc, trace = 2, 
                      subset = (!(mat_age_z2 && !mat_age_z) & 
                                  !(gest_length_z2 && !gest_length_z))
)

all_mod_full = subset(all_mod_full, !nested(.))

best_mod_full = get.models(all_mod_full, subset = delta < 6)

model.sel(best_mod_full) 

model.avg(best_mod_full)

# ------------- plot lifetime survival no twin births -------------------

best_mod_full = coxph(Surv(time = lifespan, 
                           event = departtype == "D") ~ 
                        strata(inf_sex) + 
                        gest_length + 
                        firstborn, 
                      data = surv4,
                      na.action = na.fail)

summary(best_mod_full)

surv4 %>% count(inf_sex, departtype) # only 4 of 17 females died, 17 of 26 males died

fb_df = data.frame(firstborn = c(1, 0),
                   gest_length = mean(surv4$gest_length),
                   inf_sex = 1)

fit = survfit(best_mod_full, newdata = fb_df)

ggsurvplot(fit = fit, data = fb_df, conf.int = T, censor = FALSE,
               main = "Survival by infant firstborn status", 
               xlab = "Offspring age (years)",
               ggtheme = theme_jtf(),
               legend.title = "Firstborn status:",
               legend.labs = c("Yes", "No")) +
  labs(tag = "A")

gest_qs = quantile(surv4$gest_length, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

gest_df = data.frame(firstborn = 0,
                     gest_length = unname(gest_qs),
                     inf_sex = 1)

fit = survfit(best_mod_full, newdata = gest_df)

ggsurvplot(fit = fit, data = gest_df, conf.int = F, censor = FALSE,
               main = "Survival by gestation length", 
               xlab = "Offspring age (years)",
               ggtheme = theme_jtf(),
               legend.title = "Gest length:",
               legend.labs = names(gest_qs)) + 
  labs(tag = "C")

sex_df = data.frame(firstborn = 0,
                    gest_length = mean(surv4$gest_length),
                    inf_sex = c(1, -1))

fit = survfit(best_mod_full, newdata = sex_df)

ggsurvplot(fit = fit, data = sex_df, conf.int = T, censor = TRUE,
               main = "Survival by infant firstborn status", 
               xlab = "Offspring age (years)",
               ggtheme = theme_jtf(),
               legend.title = "Offspring sex:",
               legend.labs = c("Male", "Female")) + 
  labs(tag = "B")


# ------------------- KK sample with female Elo at siring ---------------------

# include twins again:

surv5 = 
  surv2 %>%
  mutate(mat_age_z = scale(mat_age),
         mat_age_z2 = mat_age_z^2,
         gest_length_z = scale(gest_length),
         gest_length_z2 = gest_length_z^2) %>%
  filter(!is.na(mat_rank_concept))

dim(surv5) # 38 rows

table(surv5$right_cens)

surv5 %>% 
  count(inf_sex, right_cens) # only 3/15 females died, 19/23 males died

ranksurv2 = 
  coxme(Surv(time = lifespan, 
             event = departtype == "D") ~ 
          strata(inf_sex)
        + mat_age_z 
        + mat_age_z2
        + mat_rank_concept*gest_length_z
        + gest_length_z2
        + twin
        +firstborn
        + (1|mat_id_anon)
        , 
        data = surv5,
        na.action = na.fail)

all_mod_rank = dredge(ranksurv2, rank = AICc, trace = 2, 
                      subset = (!(mat_age_z2 && !mat_age_z) & 
                                  !(gest_length_z2 && !gest_length_z))
)

all_mod_rank = subset(all_mod_rank, !nested(.))

best_mod_rank = get.models(all_mod_rank, subset = delta < 6)

model.sel(best_mod_rank) 

model.avg(best_mod_rank)

#  ----------------- exclude twins from KK Elo sample ----------------

surv5 = 
  surv5 %>%
  filter(twin == 0)

# now 32 offspring

surv5 %>%
  count(inf_sex, right_cens)

table(surv5$right_cens)

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
        data = surv5,
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
# sex_df

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
