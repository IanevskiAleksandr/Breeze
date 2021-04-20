lapply(c("plotly", "scales", "parallel", "foreach", "gridExtra", "grid", "graphics", "gplots", "ggplot2", "raster", "xtable","Rcpp","dplyr"), library, character.only = !0)
lapply(c("drc", "caTools", "ggplot2", "gsubfn", "gtools", "data.table", "doSNOW","stringr","MESS"), library, character.only = !0)
options(stringsAsFactors = F)


# INPUT YOUR DATA (drug/pair doses)
dose = c(0.1, 1, 10, 100, 1000); 
inhibition_percent = c(0, 5.5, 15, 60, 99)
drug_name = "MyDrugName"
DSS_type = 2


dss<-function(ic50,slope,max,min.conc.tested,max.conc.tested,y=10,DSS.type=2,concn_scale=1e-9){
  
  a=as.numeric(unname(max)); b=as.numeric(unname(slope)); d=0; ic50 = as.numeric(unname(ic50));
  min.conc.tested = as.numeric(unname(min.conc.tested)); max.conc.tested = as.numeric(unname(max.conc.tested));
  Min.Conc<- log10(min.conc.tested*concn_scale); Max.Conc<- max.conc.tested; x2<-log10(Max.Conc*concn_scale)  
  
  if(is.na(ic50)||is.na(b)||is.na(a)||is.na(Min.Conc)||is.na(Max.Conc)){ dss<-NA } else if(isTRUE(ic50>=Max.Conc)){ dss<-0 }
  else if(isTRUE(b==0)){ dss<-0 } else{
    if(a>100){ a<-100  }
    if(isTRUE(b<0)){ b<--b  }
    c<-log10(ic50*concn_scale)
    if(a>y){
      if(y!=0){
        x1<-(c - ((log(a-y)-log(y-d))/(b*log(10))))
        if(isTRUE(x1 < Min.Conc)){x1<-Min.Conc} else if(isTRUE(x1 > x2)){x1<-x2}
      }
      else {x1<-Min.Conc}
      int_y=(((((a-d)*log(1+10^(b*(c-x2))))/(b*log(10)))+a*x2)-((((a-d)*log(1+10^(b*(c-x1))))/(b*log(10)))+a*x1)) - (y*(x2-x1))
      total_area<-(x2-Min.Conc)*(100-y)
      
      if(DSS.type==1){norm_area<-((int_y/total_area)*100)} #DSS1
      if(DSS.type==2){norm_area<-((int_y/total_area)*100)/log10(a)#DSS2 #AUC1
        if(isTRUE(norm_area > 50)){ norm_area <- 0}
      }
      if(DSS.type==3){norm_area<-((int_y/total_area)*100)*(log10(100)/log10(a))*((x2-x1)/(x2-Min.Conc))}#DSS3 #AUC5
      if(isTRUE(norm_area < 0|norm_area > 100)){ dss<-0 } else {
        dss<-round(norm_area,digits=4)}
    } else {dss<-0} 
  } 
  return (dss)
}


    #combine the data and sort by dose.
    mat_tbl <- data.frame(inhibition_percent, dose, logconc = log10(dose), 100-inhibition_percent)
    mat_tbl <- mat_tbl[order(mat_tbl[,"dose"]),]  
  
    #############    IC50
    estimate_param <- tryCatch({drm(inhibition_percent ~ logconc, data = mat_tbl, fct = LL.4(fixed = c(NA, NA, NA,NA),names = c("SLOPE","MIN","MAX","IC50")),logDose=10,control = drmc(errorm = F))}, 
                               warning=function(w){drm(inhibition_percent ~ logconc, data = mat_tbl, fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10)},
                               error=function(e){drm(inhibition_percent ~ logconc, data = mat_tbl, fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10)})
    # (extract and name coefficients)
    coef_estim <- coef(estimate_param); names(coef_estim) <- c("SLOPE","MIN","MAX","IC50")
    # see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4696819/
    coef_estim["SLOPE"] <- coef_estim["SLOPE"]*-1 
    
    # if curve decreases or IC50 is higher than max (i.e. IC50 is "outlier"), set IC50 to max conc.
    coef_estim["IC50"] <- ifelse(coef_estim["MAX"]<=coef_estim["MIN"] | coef_estim["IC50"]>max(mat_tbl$dose,na.rm=T), max(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
    # if IC50 is less than 0 set it to min. conc. and if even min. conc. < 0, then set IC50 to mean of all conc.
    coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<0,min(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
    coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<0,mean(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
    # similar to previous step but now compare log10(IC50) with log(min. conc.).
    coef_estim["IC50"] <- log10(coef_estim["IC50"])
    coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<min(mat_tbl$logconc),max(mat_tbl$logconc),coef_estim["IC50"])
    # if all inhib. < 0 set IC50 to max. log. conc !!!!! not obvious why!
    coef_estim["IC50"] <- ifelse(all(mat_tbl$inhibition_percent<0),max(mat_tbl$logconc,na.rm=T),coef_estim["IC50"])
    #(Trying to fix curves that need outlier kickout)
    coef_estim["MIN"] <- 0; coef_estim["MAX"] <- max(mat_tbl$inhibition_percent,na.rm=T)
    #(Fix off minimums) Find lowest inhibition_percent value. If it is not in (0:100), fix it whether to 0 or 99.  
    min_lower <- ifelse(min(mat_tbl$inhibition_percent,na.rm=T) > 0,min(mat_tbl$inhibition_percent,na.rm=T),0)
    min_lower <- ifelse(min_lower >= 100,99,min_lower)
    #similar to previous step but for MAX
    coef_estim["MAX"] <- ifelse(coef_estim["MAX"]>100,100,coef_estim["MAX"])
    coef_estim["MAX"] <- ifelse(coef_estim["MAX"]<0,100,coef_estim["MAX"])
    #max_lower and max_upper - lower and upper bounds for 'nl2sol' algorithm in nonlinear least-squares
    max_lower <- ifelse(max(mat_tbl$inhibition_percent,na.rm=T)>100,coef_estim["MAX"],max(mat_tbl$inhibition_percent,na.rm=T))
    max_lower <- ifelse(max_lower < 0,coef_estim["MAX"],max(mat_tbl$inhibition_percent,na.rm=T))
    max_lower <- ifelse(max_lower < 0,0,max_lower)
    max_lower <- ifelse(max_lower > 100,100,max_lower)
    #(Fix upper maximum for negative slopes)
    run_avg <- caTools::runmean(mat_tbl$inhibition_percent, 10)
    max_upper <- ifelse(any(run_avg[-nrow(mat_tbl)]>run_avg[nrow(mat_tbl)]),max(mat_tbl$inhibition_percent[run_avg>run_avg[nrow(mat_tbl)]]),coef_estim["MAX"])
    max_upper <- ifelse(any(mat_tbl$inhibition_percent > max_upper),mean(mat_tbl$inhibition_percent[mat_tbl$inhibition_percent > max_upper])+5,max_upper)
    max_upper <- ifelse(max_upper < 0,coef_estim["MAX"],max_upper)
    max_upper <- ifelse(max_upper > 100,100,max_upper) #coef_estim["MAX"]
    max_upper <- ifelse(max_lower > max_upper,coef_estim["MAX"],max_upper)
    # left it as it was, just rewritten a bit (ALEKS). not clear how values 25, 60 and 5 are chosen. 
    mean_inh_last = mean(tail(mat_tbl$inhibition_percent,2),na.rm=T)
    if(mean_inh_last < 60) {
      if(mean_inh_last > 25) coef_estim["IC50"] <- mean(mat_tbl$logconc,na.rm=T)
      else if(mean_inh_last < 25) coef_estim["IC50"] <- max(mat_tbl$logconc,na.rm=T)}
    if(mean(mat_tbl$inhibition_percent[1:3],na.rm=T)<5) coef_estim["IC50"] <- max(mat_tbl$logconc,na.rm=T)
    #add a bit of positive noise to MAX if it is the same as MIN. 
    if(unname(coef_estim["MIN"]) == unname(coef_estim["MAX"])) coef_estim["MAX"] <- coef_estim["MAX"] + 0.001
    
    #adaptive nonlinear Least-Squares algorithm NL2SOL to estimate parameters.
    nls_result_ic50_old <- function(){
      tryCatch({
        nls(inhibition_percent ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port", start=list(SLOPE=1,MIN=coef_estim["MIN"][[1]],MAX=coef_estim["MAX"][[1]],IC50=coef_estim["IC50"][[1]]), lower=list(SLOPE=0.1,MIN=0,MAX=max_lower, IC50=min(mat_tbl$logconc)), upper=list(SLOPE=2.5,MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),control=list(warnOnly=T,minFactor = 1/2048))
      }, error = function(e) {
        
        # allows higher residual sum-of-squares
        minpack.lm::nlsLM(inhibition_percent ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl,
                          start=list(SLOPE=1, MIN=coef_estim["MIN"][[1]],MAX=coef_estim["MAX"][[1]],IC50=coef_estim["IC50"][[1]]),
                          lower=c(SLOPE=0.1, MIN=0,MAX=max_lower, IC50=min(mat_tbl$logconc)),
                          upper=c(SLOPE=2.5, MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)))
      })
    } 
    
    # IC50 first
    nls_result_ic50 <- nls_result_ic50_old();
    
    # IC50 second
    nls_result_ic50_2 <- tryCatch({
      # allows higher residual sum-of-squares
      nls(inhibition_percent ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port",  start=list(SLOPE=1,MIN=coef_estim["MIN"][[1]],MAX=coef_estim["MAX"][[1]], IC50=median(mat_tbl$logconc)),lower=list(SLOPE=0.1,MIN=0,MAX=max_lower, IC50=min(mat_tbl$logconc)),upper=list(SLOPE=2.5,MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),control=list(warnOnly=T,minFactor = 1/2048))
    },warning = function(w) {
      nls_result_ic50_old()
    },error = function(e) {
      nls_result_ic50_old()
    })
    
    #element (4, 4) is zero, so the inverse cannot be computed
    nls_result_ic50 = tryCatch({summary(nls_result_ic50); nls_result_ic50},error=function(e){nls_result_ic50_2})
    
    #Calculate the standard error scores
    sumIC50 = list(summary(nls_result_ic50), summary(nls_result_ic50_2))
    
    ic50std_resid <- round(sqrt(sum((sumIC50[[1]]$residuals)^2)/(length(sumIC50[[1]]$residuals)-1)),1);
    ic50std_resid2 <- round(sqrt(sum((sumIC50[[2]]$residuals)^2)/(length(sumIC50[[2]]$residuals)-1)),1);
    
    # continue with the best
    switch_ = which.min(c(ic50std_resid, ic50std_resid2))
    nls_result_ic50 = list(nls_result_ic50, nls_result_ic50_2)[[switch_]]
    
    
    #if SLOPE <= 0.2, decrease IC50, change lower bound for SLOPE to 0.1 and repeat.
    if(coef(nls_result_ic50)["SLOPE"] <= 0.2)
    {
      if(mean_inh_last > 60)
        coef_estim["IC50"] <- min(mat_tbl$logconc,na.rm=T)
      nls_result_ic50 <- nls(inhibition_percent ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port",start=list(SLOPE=1, MIN=unname(coef_estim["MIN"]),MAX=unname(coef_estim["MAX"]),IC50=unname(coef_estim["IC50"])),lower=list(SLOPE=0.1,MIN=min_lower,MAX=max_lower,IC50=min(mat_tbl$logconc)),upper=list(SLOPE=2.5, MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),control=list(warnOnly=T,minFactor = 1/2048))
    }
    
    #Calculate the standard error scores
    sumIC50 = summary(nls_result_ic50); 
    ic50std_Error <- sumIC50$coefficients["IC50","Std. Error"]; #tec50std_Error <- sumTEC50$coefficients["TEC50","Std. Error"]
    ic50std_resid <- round(sqrt(sum((sumIC50$residuals)^2)/(length(sumIC50$residuals)-1)),1);
    max_signal <- max(mat_tbl$dose,na.rm=T); min_signal <- min(mat_tbl$dose,na.rm=T)
    
    #############################  
    #############   Final modification & STD error
    
    #prepare final data and convert IC50 back from log scale (inverse)
    coef_ic50 <- coef(nls_result_ic50)[c("IC50", "SLOPE","MAX","MIN")]; coef_ic50["IC50"] <- 10^coef_ic50["IC50"]
    #(Fix ic50 for curves in wrong direction)
    coef_ic50["IC50"] <- ifelse(coef_ic50["SLOPE"]<0,max_signal,coef_ic50["IC50"])
    #(Fix based on MAX)
    coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"]<0,max_signal,coef_ic50["IC50"])
    coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"]<10,max_signal,coef_ic50["IC50"])
    coef_ic50["MAX"] <- ifelse(coef_ic50["MAX"]<0,0,coef_ic50["MAX"])
    #(Fix over sensitive drugs)
    coef_ic50["IC50"] <- ifelse(all(c(max(mat_tbl$inhibition_percent,na.rm=T),min(mat_tbl$inhibition_percent,na.rm=T))>50),min_signal,coef_ic50["IC50"])
    
    # for ploting
    x <- seq(min(mat_tbl$logconc),max(mat_tbl$logconc), length=100)
    yic <- predict(nls_result_ic50, data.frame(logconc=x))
    auc <- MESS::auc(x,yic)
    
    ##average replicates
    mat_tblCp <- mat_tbl[, c("inhibition_percent", "dose")]
    cols_ <- colnames(mat_tblCp)[!grepl("inhibition_percent", colnames(mat_tblCp))] # columns which should be equal to average PI
    X <- as.data.table(mat_tblCp)
    mat_tblCp <- as.data.frame(X[,list(inhibition_percent = mean(inhibition_percent)),cols_], stringAsFactors = !1)
    
    
    perInh <- t(matrix(mat_tblCp[,"inhibition_percent"],dimnames=
                         list(paste0(rep("D", length(mat_tblCp[,"inhibition_percent"])), 1:length(mat_tblCp[,"inhibition_percent"])))))
    
    coef_tec50 = coef_ic50; 
    coef_tec50["IC50"] <- ifelse(coef_tec50["MAX"] > 25, coef_tec50["IC50"], max(mat_tbl$dose,na.rm=T))
    names(coef_tec50) <- c("EC50","SLOPE","MAX","MIN");
    coef_tec50["SLOPE"] = -1 * coef_tec50["SLOPE"]; # min - 0, max - 77 in ec50 it is max - 100, min - 23
    tmp = coef_tec50["MAX"]; coef_tec50["MAX"] = 100 - coef_tec50["MIN"]; coef_tec50["MIN"] = 100 - tmp; ytec <- 100 - yic;
    perViaTox <- 100 - perInh;

    
    
    ############################# 
    #############    DSS
    
    dss_score <- round(as.numeric(dss(coef_ic50["IC50"],coef_ic50["SLOPE"],coef_ic50["MAX"],min_signal,max_signal, DSS.type=as.integer(DSS_type))),1);
    coef_ic50 <- c(coef_ic50,Min.Conc.tested=min_signal,Max.Conc.tested=max_signal,IC50_std_error=ic50std_Error)
    coef_tec50 <- c(coef_tec50,Min.Conc.tested=min_signal,Max.Conc.tested=max_signal,TEC50_std_error=ic50std_Error)
    
    ####
    # Absolute IC50
    xIC50ABS <- seq(min(mat_tbl$logconc),max(mat_tbl$logconc)*15, length=5000)
    yicIC50ABS <- predict(nls_result_ic50, data.frame(logconc=xIC50ABS))
    if(all(yicIC50ABS < 50)) coef_ic50ABS= Inf else coef_ic50ABS = 10**xIC50ABS[which.min(abs(yicIC50ABS - 50))]
    ####
    
    #dataframe for IC50
    IC50_dataframe <- data.frame(DRUG_NAME=drug_name,ANALYSIS_NAME="IC50", t(as.matrix(coef_ic50)), perInh, 
                                 DSS = as.numeric(dss_score), sDSS = "", SE_of_estimate = as.numeric(ic50std_resid),AUC=auc)

    
    #round by 2 dex. all the numeric colums
    numeric_cols <- sapply(IC50_dataframe, is.numeric); IC50_dataframe[,numeric_cols] <- round(IC50_dataframe[,numeric_cols],1)
    
    # plot IC50
    #mat_tbl$inhibition_percent = xpr_tbl$inhibition_percent_percent[idx_filt]; # if we have all values < 0, they will be replaced
    #mat_tbl$viability = 100 - mat_tbl$inhibition_percent;  # we are replacing them back here.
    icpl <- ggplot2::ggplot(mat_tbl, aes(logconc, inhibition_percent)) + scale_x_continuous(breaks=mat_tbl$logconc,labels=mat_tbl$dose) +
      geom_point(color = "blue", size = 2.8) + geom_line(data = data.frame(x = x, y = yic), aes(x, yic), color="blue", size = 0.8) +
      geom_vline(xintercept = log10(coef_ic50["IC50"]), colour="grey", size = 0.8) + ggtitle(paste0(strtrim(drug_name, 15)," (dss:",dss_score,")\n"))+
      theme_bw() + labs(y = "% inhibition_percent", x = "conc(nM)")  +  ylim(-25, 125) +
      geom_text(mapping=aes(x2,y2,label = text2), data=data.frame(x2=log10(coef_ic50["IC50"])*0.95, y2=115, text2="IC50"), color="grey", parse=T) +
      theme(plot.background = element_rect(fill = "transparent",colour = NA),
            panel.background =element_rect(fill = "transparent",colour = NA), plot.title = element_text(hjust = 0.5, size = 12.5))
    
    
    print(icpl)
    print(IC50_dataframe)
