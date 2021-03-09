#=======================#
# Individual Viral load #
#=======================#

# source previous datasets and functions
source("00_wgtSurvey.R")
source("01_createVariable.R")

# load packages
library(gee)
library(geepack)
library(tidyverse)
library(ggplot2)


#-----------------------#
# Table 2               #
# Trajectory plot       #
#-----------------------#
vl.wide <- c %>% select(c(studyid, VL1, VL2, VL3, VL4, VL5)) 
vl.long = gather(vl.wide, visit, vl, VL1:VL5, factor_key = TRUE)
vl.long$visit = substr(vl.long$visit,3,4) 
vl.long$vl2 = ifelse(vl.long$vl == "LDL", 200, vl.long$vl)
vl.long$vl2 = as.numeric(vl.long$vl2)
vl.long$vl_fin = ifelse(vl.long$vl2 < 400, 200, vl.long$vl2)
# Convert to longitudinal data
vl.long2 <- vl.long %>% filter(!is.na(vl_fin)) %>% arrange(studyid, desc(visit)) %>% 
  group_by(studyid) %>% mutate(visit_fin = seq(1:n())) %>% ungroup()
# Dataset for trajectory of viral load from NASCOP databbase
vl.long_final <- vl.long2 %>% mutate(vl_log = log(vl_fin), studyid = as.character(studyid))


#-----------------------#
# Table 4               #
#-----------------------#
vl <- c %>% select(c(studyid, VL1, VL2,  VL3, VL4, VL5))
dt <- c %>% select(c(studyid, Date1, Date2,  Date3, Date4, Date5))
vl2 = gather(vl, num, vl, VL1:VL5, factor_key = TRUE) %>% 
  mutate(studyid = as.numeric(studyid), num = as.numeric(num)) %>%
  mutate(vl = as.numeric(ifelse(vl == "LDL", 200, vl)))
dt2 = gather(dt, num, dt, Date1:Date5, factor_key = TRUE) %>% mutate(studyid = as.numeric(studyid), num = as.numeric(num))
colnames(dt2); colnames(vl2)
vl3 <- vl2 %>% left_join(dt2, by = c("studyid", "num")) %>% filter(!is.na(vl)) %>%
  arrange(studyid, dt) %>% 
  group_by(studyid) %>% mutate(visit = seq(1:n())) %>% ungroup() %>% select(-c("num"))

e <- vl3 %>% full_join(c) %>% mutate(vl_fin = ifelse(is.na(vl), new_lastvl, vl)) %>%
  filter(!is.na(vl_fin)) # only keep patients with viral load results
e$t4_vl_bin = ifelse(e$vl_fin >= 1000, 1, 0) #surpressed vs unsurpressed

# HIV viral load
e$t4_vl_cat = ifelse(e$vl_fin >= 1000, ">=1000", "<1000") #surpressed vs unsurpressed

# Timing of viral load result
e$t4_vl_preg <-  case_when(#is.na(e$vl_fin) | c$nascop_vl_known != "Known" ~ NA_character_, #no VL result
  is.na(e$dt) ~ "Unknown", #has VL result but not date
  e$dt < e$tmp_datepregstart ~ "Before pregnancy", 
  e$dt >= e$tmp_datepregstart & e$dt < e$deliverydate ~ "During pregnancy",
  e$dt >= e$deliverydate ~ "After pregnancy")

# Timing of unsuppressed viral load result
e$t4_vl_preg_unsup <- ifelse(e$t4_vl_cat == ">=1000", e$t4_vl_preg, NA)

# Timing of suppressed viral load result
e$t4_vl_preg_sup <- ifelse(e$t4_vl_cat == "<1000", e$t4_vl_preg, NA)

# HIV viral load result before pregnancy
e$t4_beforepreg <- ifelse(e$t4_vl_preg == "Before pregnancy" & e$t4_vl_cat == ">=1000", ">=1000", 
                          ifelse(e$t4_vl_preg == "Before pregnancy" & e$t4_vl_cat == "<1000", "<1000", NA))
# HIV viral load result during pregnancy
e$t4_durpreg <- ifelse(e$t4_vl_preg == "During pregnancy" & e$t4_vl_cat == ">=1000", ">=1000", 
                       ifelse(e$t4_vl_preg == "During pregnancy" & e$t4_vl_cat == "<1000", "<1000", NA))
# HIV viral load result after pregnancy
e$t4_afterpreg <- ifelse(e$t4_vl_preg == "After pregnancy" & e$t4_vl_cat == ">=1000", ">=1000", 
                         ifelse(e$t4_vl_preg == "After pregnancy" & e$t4_vl_cat == "<1000", "<1000", NA))


# Timing of unsuppressed viral load result
e_addl <- e %>% group_by(studyid, fac_name, weight) %>% summarize(t4_unsup_count = sum(vl_fin >= 1000)) %>% ungroup()
sum(e_addl$t4_unsup_count) # 86 
e_addl$t4_multi_unsup = ifelse(e_addl$t4_unsup_count == 0, "none", 
                               ifelse(e_addl$t4_unsup_count == 1, "one", 
                               ifelse(e_addl$t4_unsup_count == 2, "two", "overtwo")))

# Add dummy continuous variable to use the macro des
e_addl$cont_dummy <- runif(nrow(e_addl), 0, 10)

# Summarize results
new_var=c(colnames(e)[grepl('t4_.',colnames(e))])
t4_pt = e %>% select(c(new_var, "fac_name", "weight")) #Delivery Mode and Location
t4_pt_wgt <- svydesign(id=~fac_name, weights=~weight, data = t4_pt)
t4_pt_unwgt <- t4_pt %>% select(-c("fac_name", "weight"))
t4.all <- des(t4_pt_wgt, t4_pt_unwgt, "all", 2)

t4_addl_pt = e_addl %>% select(c("t4_multi_unsup", "cont_dummy","fac_name", "weight")) #Delivery Mode and Location
t4_addl_pt_wgt <- svydesign(id=~fac_name, weights=~weight, data = t4_addl_pt)
t4_addl_pt_unwgt <- t4_addl_pt %>% select(-c("fac_name", "weight"))
t4_addl.all <- des(t4_addl_pt_wgt, t4_addl_pt_unwgt, "all", 2)

t4_all_fin <- rbind(t4.all, t4_addl.all)

t4_summary_final <- t4_all_fin[, c("var", "n", "n.miss", "wgt.mean", "wgt.sd", "wgt.lower", "wgt.upper", "unwgt.mean")]
rownames(t4_summary_final) <- NULL

# Export dataset
file_name = paste("~/Google Drive/UW Seattle/Christine McGrath/T4a.csv", sep = "")
write.csv(t4_summary_final, file_name, row.names = FALSE)


#-----------------------#
# Table 4               #
# Trajectory plot       #
#-----------------------#
unsup_pt <- unique(e[e$t4_vl_cat == ">=1000",]$studyid) # 63 patients ever unsuppressed

# Only include patients with more than 1 vl (at least 1 of them is unsuppressed)
vl4 <- e[e$studyid %in% unsup_pt,] %>% filter(!is.na(dt)) %>% 
  mutate(log_vl = log(vl_fin), studyid = as.character(studyid)) %>% 
  group_by(studyid) %>% 
  mutate(count = max(visit)) %>% ungroup() %>% 
  filter(count != 1)
vl5 <-  vl4 %>% arrange(studyid, dt) %>% group_by(studyid, t4_vl_preg) %>% mutate(preg_visit = seq(1:n()))
vl5$preg_visit_fin = paste(vl5$t4_vl_preg, vl5$preg_visit)

unsup_pt2 <- unique(vl5$studyid) #51 patients
vl5$group <- case_when(vl5$studyid %in% unsup_pt2[1:10] ~ "Group 1",
                       vl5$studyid %in% unsup_pt2[11:20] ~ "Group 2",
                       vl5$studyid %in% unsup_pt2[21:30] ~ "Group 3",
                       vl5$studyid %in% unsup_pt2[31:40] ~ "Group 4",
                       vl5$studyid %in% unsup_pt2[41:51] ~ "Group 5")
vl5$preg_visit_fin <- factor(vl5$preg_visit_fin, levels = c("Before pregnancy 1", "Before pregnancy 2", "Before pregnancy 3",
                                                            "During pregnancy 1", "During pregnancy 2",
                                                            "After pregnancy 1", "After pregnancy 2", "After pregnancy 3", "After pregnancy 4"))
ggplot(vl5, aes(x = preg_visit_fin, y = log_vl, group = studyid, color = studyid)) + 
  geom_hline(yintercept = log(1000), color = "grey") + 
  facet_wrap(~group, ncol = 1) + labs(x = "Time", y = "Log of viral load") + 
  geom_line() + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1)) +
  ggtitle("Trajectory of viral load \nfor participants with unsuppressed results")




#-----------------------#
# Table 4               #
# Poison GEE            #
#-----------------------#

# Marital status: Married/cohabiting vs Never married/divorced/separated/widowed
e$tt_marital = ifelse(is.na(e$new_marital), NA, ifelse(e$new_marital == "Married/cohabiting", "married", "single"))
e$tt_marital <- relevel(as.factor(e$tt_marital), ref = "single")

# Employment: Salaried/Self-employed vs housewide/unemployed
e$tt_employ = ifelse(is.na(e$new_employ), NA, ifelse(e$new_employ == "Salaried/Self-employed","employed", "unemployed"))
e$tt_employ <- relevel(as.factor(e$tt_employ), ref = "unemployed")

# Highest level of education: completed secondary or higher vs none/completed primary
e$tt_educat = ifelse(is.na(e$new_educat), NA, 
                     ifelse(grepl('none|primary', e$new_educat, ignore.case = TRUE), "primary", "secondary"))

# ART regimen: TDF + 3TC + EFV (first-line regimen) vs other regimens
e$tt_martregimen = ifelse(is.na(e$new_martregimen) | e$new_martregimen == "Unknown", NA, as.character(e$new_martregimen))

# Time on ART: > 1 year vs <= 1 year 
e$tt_timeonartcat = ifelse(is.na(e$new_timeonartcat) | e$new_timeonartcat == "Unknown", NA, 
                           ifelse(grepl('5', e$new_timeonartcat), ">1year", "<=1year"))
e$tt_timeonartcat <- relevel(as.factor(e$tt_timeonartcat), ref = "<=1year")

# Crowding (>=3 people per room)
e$tt_crowd <- relevel(as.factor(e$new_crowd), ref = "Not crowding")
# Partner HIV status: positive vs negative
e$tt_ptnrhivstatus = ifelse(is.na(e$new_ptnrhivstatus) | e$new_ptnrhivstatus == "Unknown", NA, 
                            ifelse(grepl('Negative', e$new_ptnrhivstatus), 'negative', 'positive'))
# Pregnancies: Primigravida vs multigravida
e$tt_pregnum <- relevel(as.factor(e$new_pregnum), ref = "Primigravida")

# Timing of HIV diagnosis: Before pregnancy (reference)
e$tt2_diagtimepoint = ifelse(is.na(e$new_diagtimepoint), NA, 
                             ifelse(grepl('labo', e$new_diagtimepoint), 'Postpartum', as.character(e$new_diagtimepoint)))
e$tt2_diagtimepoint <- factor(e$tt2_diagtimepoint, levels = c("Prior to preg", "Pregnancy", "Postpartum"))

new_var = c(colnames(e)[grepl('tt_.',colnames(e))], "new_mage", "new_numanc", "new_hivdisclose")

# Run Poisson GEE
# Binary variable
geeglm_bin = lapply(1:length(new_var), function(i) {
  f = as.formula(paste("t4_vl_bin ~", new_var[i]))
  tmp <- e %>% filter(!is.na(eval(as.name(paste(new_var[i])))))
  mod <- geeglm(formula = f, 
               data = tmp, id=studyid, family=binomial(link="logit"), corstr="independence")
  # Defualt 'san.se' is the usual robust estimate
  s =  summary(mod)
  print(s)
  c = s$coefficients
  c('var' = rownames(c)[2], 'n' = nrow(tmp), 'coef_b' = c[2, 1], "se_b" = c[2, 2], "p_b" = c[2,4])
})
tmp <- e %>% filter(!is.na(tt2_diagtimepoint))
mod <- geeglm(formula = t4_vl_bin ~ tt2_diagtimepoint, 
             data = tmp, 
             id=studyid, family=binomial(link="logit"), corstr="independence")
# Defualt 'san.se' is the usual robust estimate
s =  summary(mod)

c <- s$coefficients[2:3, c("Estimate", "Std.err", "Pr(>|W|)")]
time = cbind('var' = rownames(c), 'n' = nrow(tmp), "coef_b" = c$"Estimate", "se_b" = c$"Std.err", "p_b" = c$"Pr(>|W|)")
#time
b_comb = as.data.frame(rbind(do.call(rbind, geeglm_bin), time))

geeglm_poi = lapply(1:length(new_var), function(i) {
  f = as.formula(paste("t4_vl_bin ~", new_var[i]))
  tmp <- e %>% filter(!is.na(eval(as.name(paste(new_var[i])))))
  mod <- geeglm(formula = f, 
               data =tmp, id=studyid, family=poisson, corstr="independence")
  # Defualt 'san.se' is the usual robust estimate
  s =  summary(mod)
  c = s$coefficients
  c('var' = rownames(c)[2], 'coef_p' = c[2, 1], "se_p" = c[2, 2], "p_p" = c[2,4])
})
mod <- geeglm(formula = t4_vl_bin ~ tt2_diagtimepoint, 
             data = e %>% filter(!is.na(tt2_diagtimepoint)), 
             id=studyid, family=poisson, corstr="independence")
s =  summary(mod)

# Orgnaize the output
c <- s$coefficients[2:3, c("Estimate", "Std.err", "Pr(>|W|)")]
time = cbind('var' = rownames(c), "coef_p" = c$"Estimate", "se_p" = c$"Std.err", "p_p" = c$"Pr(>|W|)")
p_comb = as.data.frame(rbind(do.call(rbind, geeglm_poi), time))
b_comb$exp_coef_b = exp(as.numeric(as.character(b_comb$coef_b)))
b_comb$exp_lower_b = exp(as.numeric(as.character(b_comb$coef_b)) - 1.96 * as.numeric(as.character(b_comb$se_b)))
b_comb$exp_upper_b = exp(as.numeric(as.character(b_comb$coef_b)) + 1.96 * as.numeric(as.character(b_comb$se_b)))
p_comb$exp_coef_p = exp(as.numeric(as.character(p_comb$coef_p)))
p_comb$exp_lower_p = exp(as.numeric(as.character(p_comb$coef_p)) - 1.96 * as.numeric(as.character(p_comb$se_p)))
p_comb$exp_upper_p = exp(as.numeric(as.character(p_comb$coef_p)) + 1.96 * as.numeric(as.character(p_comb$se_p)))
t4b_final = p_comb %>% left_join(b_comb, by = "var")

# Export
file_name = paste("~/Google Drive/UW Seattle/Christine McGrath/T4b.csv", sep = "")
write.csv(t4b_final, file_name, row.names = FALSE)
