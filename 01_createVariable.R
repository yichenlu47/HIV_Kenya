#=======================#
# Create variables      #
#=======================#

# load pacakge and data
library(readstata13) #read.dta13
library(dplyr)
library(survey)
library(readxl)
library(lubridate)
d = read.dta13("OB+_fac_data_or_analysis.dta")

#-----------------------#
# Table 1               #
#-----------------------#

# Maternal HIV status
d$new2_mhivstatus = ifelse(d$mhivstatus == 1, "Positive",  
                           ifelse(d$mhivstatus == 0, "Negative","Unknown"))

# Maternal age, mean (SD)
d$new_mage <- d$mage
d$new_mage_cat = ifelse(is.na(d$mage), NA, 
                        ifelse(d$mage <= 19, "15-19", 
                               ifelse(d$mage <= 24, "20-24", 
                                      ifelse(d$mage <= 29, "25-29", 
                                             ifelse(d$mage <= 34, "30-34", "35+")))))

# Marital status
d$new_marital = ifelse(grepl('div|widow', d$fac_marital, ignore.case = TRUE), "Divorced/separated/widowed",
                       ifelse(grepl('currently married|come', d$fac_marital, ignore.case = TRUE), "Married/cohabiting", 
                              as.character(d$fac_marital)))

# Employment
d$new_employ = ifelse(grepl('Salar|Self', d$fac_employ), "Salaried/Self-employed", 
                      ifelse(grepl('No', d$fac_employ), NA, as.character(d$fac_employ)))


# Highest level of education
d$new_educat <- as.character(d$fac_educat)

# Crowding (>=3 people per room)
d$tmp = d$fac_numpeoplehh/d$fac_room
d$tmp2 = ifelse(is.na(d$tmp) | d$tmp == 0 | is.infinite(d$tmp), NA, d$tmp)
d$new_crowd = ifelse(is.na(d$tmp2), NA, ifelse(d$tmp2 >= 3, "Crowding", "Not crowding"))

# Pregnancies
d$new_pregnum = ifelse(d$fac_pregnant != "Yes" & d$fac_pregnum == 1, "Primigravida", "Multigravida")

# Number of living children, mean (SD)
d$new_numchild <- ifelse(is.na(d$fac_numchild), 1, 
                         ifelse(d$fac_numchild == 0, NA, d$fac_numchild))

# 4 or more ANC visits
d$tmp_numanc = ifelse(grepl('Don', d$fac_numancunk) | d$fac_numanc == 999, NA, d$fac_numanc)
d$new_numanc = ifelse(d$tmp_numanc >= 4, "4 or more ANC visits", "0-3 ANC visits")

# Mean gestational age at 1st ANC, weeks, mean (SD)
d$new_ancgestage <- ifelse(d$tmp_numanc == "0" | is.na(d$tmp_numanc), NA, d$fac_ancgestage)

# Vaginal delivery
d$new_deliverymode = ifelse(grepl('vaginal', d$fac_deliverymode, ignore.case = TRUE), "Vaginal delivery", 
                            as.character(d$fac_deliverymode))


# Preterm delivery (<37 weeks)
d$new_pretermbirth = ifelse(grepl('MCH|Know', d$fac_pretermbirth1, ignore.case = TRUE), NA, 
                            as.character(d$fac_pretermbirth1))

# Non-facility delivery
d$new_deliverylocation = 
  ifelse(grepl('facility', d$fac_deliverylocation, ignore.case = TRUE), as.character(d$fac_deliverylocation), 
         ifelse(is.na(d$fac_deliverylocation), NA, "Non-facility delivery"))

# Infant age, mean (SD)
d$new_infagemo <- d$infagemo
d$new_infagemo_cat = ifelse(is.na(d$infagemo), NA, 
                            ifelse(d$infagemo <= 4, "0-4", 
                                   ifelse(d$infagemo <= 9, "5-9", 
                                          ifelse(d$infagemo <= 13, "10-13", 
                                                 ifelse(d$infagemo <= 24, "14-24", "25+")))))

# Infant sex
d$new_infsex <- as.character(d$fac_infsex)

# Currently breastfeeding
d$new_breastfed <- as.character(d$fac_breastfed)



#-----------------------#
# Table 2               #
#-----------------------#

# load NASCOP dataset and combine with the main dataset
nascop = read_xlsx("Viral Load_deidentified_ylu.xlsx")
c = d %>% filter(mhivstatus == 1) %>% left_join(nascop, by = "studyid")

#Timing of HIV diagnosis
c$new_diagtimepoint = ifelse(c$diagtimepoint == 1, "Pregnancy",
                             ifelse(c$diagtimepoint == 2, "Postpartum",
                                    ifelse(c$diagtimepoint == 3, "labor and deliv",
                                           ifelse(c$diagtimepoint == 4, "Prior to preg", NA))))

# Timing of ART initiation
c$new_artstartpreg = ifelse(c$studyid == "x", "Before last pregnancy", as.character(c$fac_artstartpreg))

# Mother on Antiretroviral Therapy (ART)
c$new_onart = ifelse(grepl('know', c$fac_onart), NA, as.character(c$fac_onart))

# ART regimen
c$new_martregimen = ifelse(is.na(c$martregimen) | c$fac_onart != "Yes", NA, 
                           ifelse(c$martregimen == 1, "T3E", 
                                  ifelse(c$martregimen == 88, "Unknown", "Other regimens")))

# Time on ART
c$new_timeonartcat = ifelse(c$timeonartcat == 1, "<=6mo",
                            ifelse(c$timeonartcat == 2, ">6mo and <=1year",
                                   ifelse(c$timeonartcat == 3, ">1year & <=5year",
                                          ifelse(c$timeonartcat == 4, ">=5 year",
                                                 ifelse(c$timeonartcat == 88 | c$timeonartcat == 999, "Unknown", NA)))))

# CD4 cell count KNOWN
table(c$fac_cd4result, useNA = "always") #225 Received results
table(c[c$fac_cd4result == "Yes",]$fac_cd4countunk, useNA = "always") # 135 out of 225 select unknown
table(c$fac_cd4count, useNA = "always") #183 participants with 0 cd4 count (all unknown if received results)
c$new_cd4count_known = ifelse(is.na(c$fac_cd4result) | c$fac_cd4result == "No", NA,
                              ifelse(grepl("Don't know", c$fac_cd4countunk) | c$fac_cd4count == 0, "Unknown", "Known"))  #based on Eddy's email, those with cd4 = 0 as unknown
table(c$new_cd4count_known, useNA = "always")  

# Last CD4 count, cells/mm3, mean (SD)
c$new_cd4count = ifelse(c$new_cd4count_known == "Known", c$fac_cd4count, NA)
table(c$new_cd4count, useNA = "always")

c$new_cd4countc = ifelse(is.na(c$new_cd4count), NA, ifelse(c$new_cd4count>= 250, ">=250", "<250"))
table(c$new_cd4countc)

# HIV Viral Load KNOWN
c$new_vl_known = ifelse(is.na(c$fac_vlresult) | c$fac_vlresult == "No", NA, #not received or unsure excluded
                        ifelse(grepl("know", c$fac_lastvlunk) | c$fac_lastvl == 0, "Unknown", "Known"))



# Last HIV viral load, mean (SD)
c$new_lastvl = ifelse(c$new_vl_known == "Known", c$fac_lastvl, NA)
table(c$new_lastvl, useNA = "always") #56 vl

c$new_lastvlc = ifelse(is.na(c$new_lastvl), NA, ifelse(c$new_lastvl>= 1000, ">=1000", "<1000"))

# Last HIV viral load
c$new_lastvlc2 = ifelse(is.na(c$new_lastvl), NA, ifelse(c$new_lastvl>= 400, ">=400", "<400"))

# Partner HIV status
c$tmp_ptnrhivtest = ifelse(c$fac_ptnrhivtest == "Yes" | c$fac_recptnrhivtest == "Yes", "Yes", NA)

c$new_ptnrhivstatus = ifelse(grepl('Positive', c$fac_recptnrhivstatus) | grepl('Positive', c$fac_ptnrhivstatus), "Positive",
                             ifelse(grepl('Negative', c$fac_recptnrhivstatus) | grepl('Negative', c$fac_ptnrhivstatus), "Negative", 
                                    ifelse(grepl('know', c$fac_recptnrhivstatus) | grepl('know', c$fac_ptnrhivstatus), "Unknown", NA)))
# Disclosed HIV status to male partner
c$new_hivdisclose = as.character(c$fac_hivdisclose)

#-----------------------#
# Table 2               #
# NASCOP                #
#-----------------------#

# HIV viral load KNOWN 
c$tmp_lastvl <- c$"VL1"
table(c$tmp_lastvl, useNA = "always") #no 0, 49 missing
table(c$new_vl_known, useNA = "always")
table(c[is.na(c$tmp_lastvl),]$new_vl_known, useNA = "always") #44 missing fron NASCOP - 5 known, 29 unknown and 15 not received results from main
c$nascop_vl_known <- ifelse(is.na(c$tmp_lastvl) & is.na(c$new_vl_known), NA, #if missing from NASCOP and never received result from main
                            ifelse(is.na(c$tmp_lastvl) & c$new_vl_known == "Unknown", "Unknown", "Known"))
table(c$nascop_vl_known, useNA = "always")

# Last HIV viral load
c$tmp_lastvl2 <- ifelse(c$nascop_vl_known == "Unknown" | is.na(c$nascop_vl_known), NA,
                        ifelse(is.na(c$tmp_lastvl), c$new_lastvl, as.character(c$tmp_lastvl))) #if known but not from nascop, then from main
table(c$tmp_lastvl2, useNA = "always")

c$nascop_lastvl <- as.numeric(ifelse(is.na(c$tmp_lastvl2), NA, 
                                     ifelse(c$tmp_lastvl2 == "LDL", 200, c$tmp_lastvl2)))
table(c$nascop_lastvl, useNA = "always") #final nascop number (cont)
c$nascop_lastvlc2 <- ifelse(is.na(c$nascop_lastvl), NA, 
                            ifelse(c$nascop_lastvl < 1000, "<1000", ">=1000"))

# Last HIV viral load
c$nascop_lastvlc <- ifelse(is.na(c$nascop_lastvl), NA, 
                           ifelse(c$nascop_lastvl < 400, "<400", ">=400"))

table(c$nascop_lastvlc, useNA = "always"); table(c$nascop_lastvlc2, useNA = "always")


#-----------
# IF HIV Viral Load KNOWN, first HIV Viral Load (continuous + categorical)
#-----------
c$tmp_vl2 <- c$"VL2"; c$tmp_vl3 <- c$"VL3"; c$tmp_vl4 <- c$"VL4"; c$tmp_vl5 <- c$"VL5"

sum(!is.na(c$tmp_vl2));sum(!is.na(c$tmp_vl3)); sum(!is.na(c$tmp_vl4)); sum(!is.na(c$tmp_vl5));
c$tmp_firstvlc <- ifelse(is.na(c$nascop_lastvl), NA,
                         ifelse(!is.na(c$tmp_vl5), "vl5", 
                                ifelse(!is.na(c$tmp_vl4), "vl4", 
                                       ifelse(!is.na(c$tmp_vl3),"vl3",
                                              ifelse(!is.na(c$tmp_vl2), "vl2", "vl1")))))
c$tmp_firstvl <- ifelse(is.na(c$nascop_lastvl), NA,
                        ifelse(!is.na(c$tmp_vl5), c$tmp_vl5, 
                               ifelse(!is.na(c$tmp_vl4), c$tmp_vl4, 
                                      ifelse(!is.na(c$tmp_vl3), c$tmp_vl3,
                                             ifelse(!is.na(c$tmp_vl2), c$tmp_vl2, c$nascop_lastvl)))))
c$nascop_firstvl <- as.numeric(ifelse(is.na(c$tmp_firstvl), NA, 
                                      ifelse(c$tmp_firstvl == "LDL", 200, c$tmp_firstvl)))

c$nascop_firstvlc <- ifelse(is.na(c$nascop_firstvl), NA, 
                            ifelse(c$nascop_firstvl < 400, "<400", ">=400"))
c$nascop_firstvlc2 <- ifelse(is.na(c$nascop_firstvl), NA, 
                             ifelse(c$nascop_firstvl < 1000, "<1000", ">=1000"))

# Number of viral load (VL) results, mean (SD)
c$nascop_multi <- case_when(c$nascop_vl_known == "Known" & is.na(c$tmp_vl2) ~ 1,
                            c$nascop_vl_known == "Known" & !is.na(c$tmp_vl2) & is.na(c$tmp_vl3) ~ 2,
                            c$nascop_vl_known == "Known" & !is.na(c$tmp_vl3) & is.na(c$tmp_vl4) ~ 3,
                            c$nascop_vl_known == "Known" & !is.na(c$tmp_vl4) & is.na(c$tmp_vl5) ~ 4,
                            c$nascop_vl_known == "Known" & !is.na(c$tmp_vl5) ~ 5)
c$nascop_multic <-case_when(c$nascop_vl_known == "Known" & !is.na(c$tmp_vl2) & !is.na(c$tmp_vl3) ~ "over2", 
                            c$nascop_vl_known == "Known" & !is.na(c$tmp_vl2) & is.na(c$tmp_vl3) ~ "two", 
                            c$nascop_vl_known == "Known" & is.na(c$tmp_vl2) ~ "single")

c$tmp_lastvl_date <- c$"Date1"
table(c$tmp_lastvl_date , useNA = "always")  #55 missing
table(c$datepregstart, useNA = "always") #87 misisng, 2 in 2020

c$tmp_datepregstart <- ifelse(is.na(c$datepregstart) | year(c$datepregstart) == 2020 | c$deliverydate <= c$datepregstart, #2 patients with delivery date before pregnancy date
                              c$deliverydate - 40 * 7 + 1, c$datepregstart) # calculate pregnancy start date if missing
table(c$tmp_datepregstart, useNA = "always")
class(c$tmp_datepregstart) <- "Date"
table(c$tmp_datepregstart, useNA = "always"); table(c$deliverydate, useNA = "always")

# Timing of last HIV viral load result
c$nascop_lastvl_preg <-  case_when(is.na(c$nascop_vl_known) | c$nascop_vl_known != "Known" ~ NA_character_, #no VL result
                                   is.na(c$tmp_lastvl_date) ~ "Unknown", #has VL result but not date
                                   c$tmp_lastvl_date < c$tmp_datepregstart ~ "Before pregnancy", 
                                   c$tmp_lastvl_date >= c$tmp_datepregstart & c$tmp_lastvl_date < c$deliverydate ~ "During pregnancy",
                                   c$tmp_lastvl_date >= c$deliverydate ~ "After pregnancy")
table(c$nascop_lastvl_preg, useNA = "always")   #n = 12 unknown here means has vl result but no date 

# Last HIV viral load result before pregnancy
c$nascop_lastvl_beforepreg <- ifelse(c$nascop_lastvl_preg == "Before pregnancy" & c$nascop_lastvl >= 1000, ">=1000", 
                                     ifelse(c$nascop_lastvl_preg == "Before pregnancy" & c$nascop_lastvl < 1000, "<1000", NA))
# Last HIV viral load result during pregnancy
c$nascop_lastvl_durpreg <- ifelse(c$nascop_lastvl_preg == "During pregnancy" & c$nascop_lastvl >= 1000, ">=1000", 
                                  ifelse(c$nascop_lastvl_preg == "During pregnancy" & c$nascop_lastvl < 1000, "<1000", NA))
# Last HIV viral load result after pregnancy
c$nascop_lastvl_afterpreg <- ifelse(c$nascop_lastvl_preg == "After pregnancy" & c$nascop_lastvl >= 1000, ">=1000", 
                                    ifelse(c$nascop_lastvl_preg == "After pregnancy" & c$nascop_lastvl < 1000, "<1000", NA))
table(c$nascop_lastvl_beforepreg, useNA = "always"); table(c$nascop_lastvl_durpreg, useNA = "always"); table(c$nascop_lastvl_afterpreg, useNA = "always")


# Timing of first HIV viral load result
c$tmp_vl2_date <- c$"Date2"; table(c$tmp_vl2_date)
c$tmp_vl3_date <- c$"Date3"; table(c$tmp_vl3_date)
c$tmp_vl4_date <- c$"Date4"; table(c$tmp_vl4_date)
c$tmp_vl5_date <- c$"Date5"; table(c$tmp_vl5_date)

c$tmp_firstvl_date <- ifelse(is.na(c$tmp_lastvl_date), NA,
                             ifelse(!is.na(c$tmp_vl5_date), as.Date(c$tmp_vl5_date), 
                                    ifelse(!is.na(c$tmp_vl4_date), as.Date(c$tmp_vl4_date), 
                                           ifelse(!is.na(c$tmp_vl3_date), as.Date(c$tmp_vl3_date),
                                                  ifelse(!is.na(c$tmp_vl2_date), as.Date(c$tmp_vl2_date), as.Date(c$tmp_lastvl_date))))))

class(c$tmp_firstvl_date) <- "Date"; table(c$tmp_firstvl_date)


c$nascop_firstvl_preg <-  case_when(is.na(c$nascop_vl_known) | c$nascop_vl_known != "Known" ~ NA_character_, #no VL result
                                    is.na(c$tmp_firstvl_date) ~ "Unknown", #has VL result but not date
                                    c$tmp_firstvl_date < c$tmp_datepregstart ~ "Before pregnancy", 
                                    c$tmp_firstvl_date >= c$tmp_datepregstart & c$tmp_firstvl_date < c$deliverydate ~ "During pregnancy",
                                    c$tmp_firstvl_date >= c$deliverydate ~ "After pregnancy")
table(c$nascop_firstvl_preg, useNA = "always")   #n = 12 unknown here means has vl result but no date 

# First HIV viral load result before pregnancy
c$nascop_firstvl_beforepreg <- ifelse(c$nascop_firstvl_preg == "Before pregnancy" & c$nascop_firstvl >= 1000, ">=1000", 
                                      ifelse(c$nascop_firstvl_preg == "Before pregnancy" & c$nascop_firstvl < 1000, "<1000", NA))
# First HIV viral load result during pregnancy
c$nascop_firstvl_durpreg <- ifelse(c$nascop_firstvl_preg == "During pregnancy" & c$nascop_firstvl >= 1000, ">=1000", 
                                   ifelse(c$nascop_firstvl_preg == "During pregnancy" & c$nascop_firstvl < 1000, "<1000", NA))

# First HIV viral load result during pregnancy
c$nascop_firstvl_afterpreg <- ifelse(c$nascop_firstvl_preg == "After pregnancy" & c$nascop_firstvl >= 1000, ">=1000", 
                                     ifelse(c$nascop_firstvl_preg == "After pregnancy" & c$nascop_firstvl < 1000, "<1000", NA))
table(c$nascop_firstvl_beforepreg, useNA = "always"); table(c$nascop_firstvl_durpreg, useNA = "always"); table(c$nascop_firstvl_afterpreg, useNA = "always")



# Timing of last HIV viral load result during pregnancy
c$nascop_lastvl_durpreg_c<-  case_when(c$nascop_lastvl_preg == "During pregnancy" & c$tmp_lastvl_date < c$tmp_datepregstart + 11 * 7 ~ "1tri", 
                                       c$nascop_lastvl_preg == "During pregnancy" & c$tmp_lastvl_date >= c$tmp_datepregstart + 11 * 7 & c$tmp_lastvl_date < c$tmp_datepregstart + 23 * 7 ~ "2tri", 
                                       c$nascop_lastvl_preg == "During pregnancy" & c$tmp_lastvl_date >= c$tmp_datepregstart + 23 * 7 ~ "3tri")

c$nascop_lastvl_durpreg1tri <- case_when(c$nascop_lastvl_durpreg_c == "1tri" & c$nascop_lastvl >= 1000 ~ ">=1000", 
                                         c$nascop_lastvl_durpreg_c == "1tri" & c$nascop_lastvl < 1000 ~ "<1000")
c$nascop_lastvl_durpreg2tri <- case_when(c$nascop_lastvl_durpreg_c == "2tri" & c$nascop_lastvl >= 1000 ~ ">=1000", 
                                         c$nascop_lastvl_durpreg_c == "2tri" & c$nascop_lastvl < 1000 ~ "<1000")
c$nascop_lastvl_durpreg2tri <- factor(c$nascop_lastvl_durpreg2tri, levels = c("<1000", ">=1000"))
c$nascop_lastvl_durpreg3tri <- case_when(c$nascop_lastvl_durpreg_c == "3tri" & c$nascop_lastvl >= 1000 ~ ">=1000", 
                                         c$nascop_lastvl_durpreg_c == "3tri" & c$nascop_lastvl < 1000 ~ "<1000")
table(c$nascop_lastvl_durpreg1tri, useNA = "always"); table(c$nascop_lastvl_durpreg2tri, useNA = "always"); table(c$nascop_lastvl_durpreg3tri, useNA = "always")



# Timing of first HIV viral load result during pregnancy
c$nascop_firstvl_durpreg_c <-  case_when(c$nascop_firstvl_preg == "During pregnancy" & c$tmp_firstvl_date < c$tmp_datepregstart + 11 * 7 ~ "1tri", 
                                         c$nascop_firstvl_preg == "During pregnancy" & c$tmp_firstvl_date >= c$tmp_datepregstart + 11 * 7 & c$tmp_firstvl_date < c$tmp_datepregstart + 23 * 7 ~"2tri", 
                                         c$nascop_firstvl_preg == "During pregnancy" & c$tmp_firstvl_date >= c$tmp_datepregstart + 23 * 7 ~ "3tri")


c$nascop_firstvl_durpreg1tri <- case_when(c$nascop_firstvl_durpreg_c  == "1tri" & c$nascop_firstvl >= 1000 ~ ">=1000", 
                                          c$nascop_firstvl_durpreg_c == "1tri" & c$nascop_firstvl < 1000 ~ "<1000")
c$nascop_firstvl_durpreg2tri <- case_when(c$nascop_firstvl_durpreg_c == "2tri" & c$nascop_firstvl >= 1000 ~ ">=1000", 
                                          c$nascop_firstvl_durpreg_c  == "2tri" & c$nascop_firstvl < 1000 ~ "<1000")
c$nascop_firstvl_durpreg3tri <- case_when(c$nascop_firstvl_durpreg_c  == "3tri" & c$nascop_firstvl >= 1000 ~ ">=1000", 
                                          c$nascop_firstvl_durpreg_c  == "3tri" & c$nascop_firstvl < 1000 ~ "<1000")
table(c$nascop_firstvl_durpreg1tri, useNA = "always"); table(c$nascop_firstvl_durpreg2tri, useNA = "always"); table(c$nascop_firstvl_durpreg3tri, useNA = "always")


# Last HIV viral load result after survey date
table(c$today, useNA = "always") #0 missing value
#only among those with date
c$tmp_lastvl_today <-  difftime(c$tmp_lastvl_date, c$today, units = "days") #difference in days
table(c$tmp_lastvl_today, useNA = "always") #still 55 missing

c$nascop_lastvl_survey_c = ifelse(is.na(c$tmp_lastvl_today), NA,
                                 ifelse(c$tmp_lastvl_today < 0, "vl_after_today", "vl_before_today"))

table(c$nascop_lastvl_survey_c, useNA = "always") #0 missing value

c$nascop_lastvl_before_today <- ifelse(is.na(c$nascop_lastvl_survey_c) | c$tmp_lastvl_today < 0, NA, c$tmp_lastvl_today/30.25)
c$nascop_lastvl_after_today  <- ifelse(is.na(c$nascop_lastvl_survey_c) | c$tmp_lastvl_today >= 0, NA, abs(c$tmp_lastvl_today/30.25))

# First HIV viral load result after survey date
c$tmp_firstvl_today <-  difftime(c$tmp_firstvl_date, c$today, units = "days") #difference in days
table(c$tmp_firstvl_today, useNA = "always") #still 55 missing

c$nascop_firstvl_survey_c = ifelse(is.na(c$tmp_firstvl_today), NA,
                                ifelse(c$tmp_firstvl_today < 0, "vl_after_today", "vl_before_today"))
c$nascop_firstvl_before_today <- ifelse(is.na(c$nascop_firstvl_survey_c) | c$tmp_firstvl_today < 0, NA, c$tmp_firstvl_today/30.25)
c$nascop_firstvl_after_today  <- ifelse(is.na(c$nascop_firstvl_survey_c) | c$tmp_firstvl_today >= 0, NA, abs(c$tmp_firstvl_today/30.25))
