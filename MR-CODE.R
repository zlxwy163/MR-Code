library(TwoSampleMR) 


Exposure<-extract_instruments(
  outcomes='GWAS ID',
  p1 = 5e-06, #threshold
  
  clump=TRUE, 
  
  r2=0.1,  #threshold
  
  kb=100,
  
  access_token= NULL) 

dim(S_U_P_A_R) #dimension; 

View(S_U_P_A_R)

#查看结果 

write.csv(S_U_P_A_R, file="D:\\学习  工作\\孟德尔随机化讲课0919\\单变量孟德尔随机化讲课示例\\S_U_P_A_R.csv")



#剔除混杂、结局的SNP
#剔除后，重新读取数据

S_U_P_A_R <-read.csv("S_U_P_A_R.csv", header = TRUE)
dim(S_U_P_A_R) 
View(S_U_P_A_R)



SCZ<- extract_outcome_data(
  
  snps=RA$SNP,
  
  
  outcomes='ieu-b-5070',
  
  
  proxies = FALSE,
  
  
  maf_threshold = 0.01,
  
  access_token = NULL)


dim(SCZ)
View(SCZ)

write.csv(SCZ,file="SCZ.csv")

#备用计算R2，F


Mydata$EAF2 <- (1 - Mydata$eaf.exposure)
Mydata$MAF <- pmin(Mydata$eaf.exposure, Mydata$EAF2)
PVEfx <- function(BETA, MAF, SE, N){
  pve <- (2*(BETA^2)*MAF*(1 - MAF))/((2*(BETA^2)*MAF*(1 - MAF)) + ((SE^2)*2*N*MAF*(1 - MAF)))
  return(pve) 
}
Mydata$PVE <- mapply(PVEfx, Mydata$beta.exposure, Mydata$MAF, Mydata$se.exposure, N = Mydata$samplesize.exposure)
Mydata$FSTAT <- ((Mydata$samplesize.exposure - 1 - 1)/1)*(Mydata$PVE/(1 - Mydata$PVE))


Mydata
write.csv(Mydata, file="D:\\学习  工作\\孟德尔随机化讲课0919\\单变量孟德尔随机化讲课示例\\Mydata.csv")


Mydata <- harmonise_data(
  exposure_dat=S_U_P_A_R,
  outcome_dat=myocardial_infarction,
  action= 2)                                
View(Mydata)
dim(Mydata)

write.csv(Mydata,file="D:\\学习  工作\\孟德尔随机化讲课0919\\单变量孟德尔随机化讲课示例\\Mydata.csv")


Mydata<-read.csv("Mydata.csv",header = T)
res <-mr(Mydata)


#res <-mr(Mydata, method_list=c("mr_egger_regression", "mr_weighted_median","mr_ivw_mre", "mr_ivw_fe","mr_simple_mode","mr_weighted_mode" ))

res  #结果计算

write.csv(res,file="D:\\学习  工作\\孟德尔随机化讲课0919\\单变量孟德尔随机化讲课示例\\res.csv")
##########################################################################
########################################################
##########################################################



#异质性，敏感性，多效性分析


mr_scatter_plot(res,Mydata)#散点图

res_single <- mr_singlesnp(Mydata)
mr_forest_plot(res_single)#森林图




mr_heterogeneity(Mydata, method_list=c("mr_egger_regression", "mr_ivw"))  #异质性检验——IVWorMRegger
#I2=[Q-(K-1)]/Q=[5.73-(11-1)]/5.73=4.27/5.73=


pleio <- mr_pleiotropy_test(Mydata) 
pleio   #多效性检验——MR egger 
View(pleio)
write.csv(pleio, file="pleio.csv")




single <- mr_leaveoneout(Mydata)
mr_leaveoneout_plot(single)   #留一法检验敏感性

write.csv(single, file="single.csv")


mr_funnel_plot(res_single)#漏斗图
result<-generate_odds_ratios(res)#算OR值

result
write.csv(result, file="result.csv")




#MRPRESSO
library(MRPRESSO)
mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data =Mydata, NbDistribution = 1000, SignifThreshold = 0.05)