#=======================#
# Weighted survey       #
#=======================#

# Compare the two groups
test2grp <- function(pt_wgt, pt_unwgt, g, v){ 
  num_col <- unlist(lapply(pt_unwgt, is.numeric))
  cat_col <- unlist(lapply(pt_unwgt, is.character))
  
  # categorical: weighted chi-square
  cat.test <- lapply(1:(ncol(pt_unwgt[, cat_col])-1), function(i){
    f = as.formula(paste("~", colnames(pt_unwgt[, cat_col])[i], "+", v))
    test.res = svychisq(f, pt_wgt)
    data.frame("mega.var" = colnames(pt_unwgt[, cat_col])[i], "pvalue" = test.res$p.value)
  })
  cat.test2 <- do.call(rbind, cat.test)
  
  # continuous: weighted t-test
  cont.test <- lapply(1:ncol(pt_unwgt[, num_col]), function(i){
    f = as.formula(paste(colnames(pt_unwgt[, num_col])[i], "~", v))
    test.res = svyttest(f, pt_wgt) 
    data.frame("mega.var" = colnames(pt_unwgt[, num_col])[i], "pvalue" = test.res$p.value)
  })
  cont.test2 <- do.call(rbind, cont.test )
  
  test <- rbind(cont.test2, cat.test2)
  test$group <- g # add group
  test
}

# Calculate weighted and unweighted summary statistics
des <- function(pt_wgt, pt_unwgt, g, r){ 
  l = ncol(pt_wgt) - r #number of variables to analyze, we need to remove the split variable when analyzing subgroups
  
  # weighted
  w <- lapply(1:l, function(i){
    f = as.formula(paste("~",colnames(pt_wgt)[i]))
    
    m = svymean(f, pt_wgt, na.rm = TRUE)
    c = confint(m) #95% confidence interval
    m = as.data.frame(m)
    colnames(m)[2] <- "SE"
    cbind(var = rownames(m), m, c)
  })
  ww <- as.data.frame(do.call(rbind, w))
  colnames(ww)[2:5] <- c("wgt.mean", "wgt.sd", "wgt.lower", "wgt.upper")
  
  # unweighted
  nn <- as.data.frame(t(sapply(pt_unwgt, function(x) c(n.miss = sum(is.na(x)), n = sum(!is.na(x))))))
  nn$mega.var = rownames(nn)
  
  # continuous
  num_col <- unlist(lapply(pt_unwgt, is.numeric))
  m.cont  <- as.data.frame(t(sapply(pt_unwgt[, num_col], function(x) c(unwgt.mean = mean(x, na.rm = TRUE), unwgt.sd = sd(x, na.rm = TRUE)))))
  m.cont$mega.var = m.cont$var = rownames(m.cont)
  
  # categorical
  cat_col <- pt_unwgt[, unlist(lapply(pt_unwgt, is.character))]
  m.cat.list  <- lapply(1:ncol(cat_col), function(i) {
    tmp = as.data.frame(table(cat_col[,colnames(cat_col)[i]]))
    tmp$mega.var = colnames(cat_col)[i]
    tmp
  })
  m.cat <- do.call(rbind, m.cat.list)
  
  m.cat2 = m.cat %>% mutate(unwgt.mean = Freq, var = paste(mega.var, Var1, sep = ""), unwgt.sd = NA) %>% 
    select(-c(Var1, Freq))
  
  rbind(m.cont, m.cat2) %>% left_join(nn, by = "mega.var") %>% left_join(ww, by = "var") %>% mutate(group = g)
}
