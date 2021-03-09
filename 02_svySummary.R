#=======================#
# Summary results       #
#=======================#

# source previous datasets and functions
source("00_wgtSurvey.R")
source("01_createVariable.R")

#-----------------------#
# Table 1               #
#-----------------------#
new_var=c(colnames(d)[grepl('new_.',colnames(d))]) # extract variables

# Prepare datasets
t1_pt = d %>% select(c(new_var,  "new2_mhivstatus", "fac_name", "weight")) 

t1_pt_wgt <- svydesign(id=~fac_name, weights=~weight, data = t1_pt)
t1_pt_wgt.pos <- subset(t1_pt_wgt, new2_mhivstatus == "Positive") # add groups
t1_pt_wgt.neg <- subset(t1_pt_wgt, new2_mhivstatus == "Negative")

t1_pt_unwgt <- t1_pt %>% select(-c("fac_name", "weight"))
t1_pt_unwgt.pos <- subset(t1_pt_unwgt, new2_mhivstatus == "Positive") %>% select(-"new2_mhivstatus")
t1_pt_unwgt.neg <- subset(t1_pt_unwgt, new2_mhivstatus == "Negative") %>% select(-"new2_mhivstatus")

# Descriptive analysis
t1.all <- des(t1_pt_wgt, t1_pt_unwgt, "all", 2)
t1.pos <- des(t1_pt_wgt.pos, t1_pt_unwgt.pos, "pos", 3)
t1.neg <- des(t1_pt_wgt.neg, t1_pt_unwgt.neg, "neg", 3)

# Run tests
t1.test <- test2grp(t1_pt_wgt, t1_pt_unwgt, "all", "new2_mhivstatus")

# Combine results
t1_summary_final <- rbind(t1.all, t1.pos, t1.neg) %>% left_join(t1.test, by = c("mega.var", "group"))
t1_summary_final$group.var <- paste0(t1_summary_final$group, t1_summary_final$var)
t1_summary_final <- t1_summary_final[, c("group.var", "group", "var", "n", "n.miss", "wgt.mean", "wgt.sd", "wgt.lower", "wgt.upper", "unwgt.mean", "pvalue")]
rownames(t1_summary_final) <- NULL

#-----------------------#
# Table 2               #
#-----------------------#
new_var=c(colnames(c)[grepl('new_.',colnames(c)) | grepl('nascop_.',colnames(c))])
t2_pt = c %>% select(c(new_var,  "fac_name", "weight")) #Delivery Mode and Location
# t2_pt <- t2_pt %>% select(-"nascop_lastvl_durpreg2tri")
t2_pt_wgt <- svydesign(id=~fac_name, weights=~weight, data = t2_pt)
t2_pt_unwgt <- t2_pt %>% select(-c("fac_name", "weight"))

t2.all <- des(t2_pt_wgt, t2_pt_unwgt, "all", 2)
t2_summary_final <- t2.all[, c("var", "n", "n.miss", "wgt.mean", "wgt.sd", "wgt.lower", "wgt.upper", "unwgt.mean")]
rownames(t2_summary_final) <- NULL


#-----------------------#
# Table 3               #
#-----------------------#
new_var=c("new_mage", "new_mage_cat", "new_marital", "new_employ", "new_educat", "new_crowd",
          "new_pregnum", "new_numanc",
          "new_diagtimepoint", "new_martregimen", "new_timeonartcat",
          "new_ptnrhivstatus", "new_hivdisclose",
          "nascop_multi", "nascop_multic")
t3_pt = c %>% select(c(new_var,  "nascop_lastvlc2", "fac_name", "weight")) #Delivery Mode and Location

t3_pt_wgt <- svydesign(id=~fac_name, weights=~weight, data = t3_pt)
t3_pt_wgt.sup <- subset(t3_pt_wgt, nascop_lastvlc2 == "<1000")
t3_pt_wgt.unsup <- subset(t3_pt_wgt, nascop_lastvlc2 == ">=1000")

t3_pt_unwgt <- t3_pt %>% select(-c("fac_name", "weight"))
t3_pt_unwgt.sup <- subset(t3_pt_unwgt, nascop_lastvlc2 == "<1000") %>% select(-"nascop_lastvlc2")
t3_pt_unwgt.unsup <- subset(t3_pt_unwgt, nascop_lastvlc2 == ">=1000") %>% select(-"nascop_lastvlc2")

t3.all <- des(t3_pt_wgt, t3_pt_unwgt, "all", 2)
t3.sup <- des(t3_pt_wgt.sup, t3_pt_unwgt.sup, "sup", 3)
t3.unsup <- des(t3_pt_wgt.unsup, t3_pt_unwgt.unsup, "unsup", 3)
t3.test <- test2grp(t3_pt_wgt, t3_pt_unwgt, "all","nascop_lastvlc2" )

t3_summary_final <- rbind(t3.all, t3.sup, t3.unsup) %>% left_join(t3.test, by = c("mega.var", "group"))
t3_summary_final$group.var <- paste0(t3_summary_final$group, t3_summary_final$var)
t3_summary_final <- t3_summary_final[, c("group.var", "group", "var", "n", "n.miss", "wgt.mean", "wgt.sd", "wgt.lower", "wgt.upper", "unwgt.mean", "pvalue")]
rownames(t3_summary_final) <- NULL

# output results
file_name = paste("~/Google Drive/UW Seattle/Christine McGrath/T1raw.csv", sep = "")
write.csv(t1_summary_final, file_name, row.names = FALSE)
file_name = paste("~/Google Drive/UW Seattle/Christine McGrath/T2raw.csv", sep = "")
write.csv(t2_summary_final, file_name, row.names = FALSE)
file_name = paste("~/Google Drive/UW Seattle/Christine McGrath/T3raw.csv", sep = "")
write.csv(t3_summary_final, file_name, row.names = FALSE)


