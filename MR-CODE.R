#install.packages("devtools")#

#devtools::install_github("MRCIEU/TwoSampleMR")#  



#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")




library(TwoSampleMR) 




Exposure<-extract_instruments(
  
  
  outcomes='ExposureGWAS ID',
  
  
  p1 = 5e-06, 
  
  
  clump=TRUE, 
  
  r2=0.001,  
  
  kb=10000,
  
  access_token= NULL) 

dim(Exposure) #dimension; 
View(Exposure)
write.csv(Exposure, file="Exposure.csv")


Exposure <-read.csv("Exposure.csv", header = TRUE)

##############################################################

OUTCOM<- extract_outcome_data(
  
  snps=OUTCOM$SNP,
  
  
  outcomes='OUTCOM GWAS ID',
  
  
  proxies = FALSE,
  

  maf_threshold = 0.01,
  
  access_token = NULL)


dim(OUTCOM)
View(OUTCOM)

write.csv(OUTCOM,file="OUTCOM.csv")

Mydata <- harmonise_data(
  exposure_dat=Exposure,
  outcome_dat=OUTCOM,
  action= 2)                             
View(Mydata)
dim(Mydata)

write.csv(Mydata, file="Mydata.csv")

Mydata<-read.csv("Mydata.csv",header = T)

res <-mr(Mydata)
res  

##########################################################################
########################################################
#########################################################

mr_scatter_plot(res,Mydata)

res_single <- mr_singlesnp(Mydata)
mr_forest_plot(res_single)




mr_heterogeneity(Mydata, method_list=c("mr_egger_regression", "mr_ivw"))  #


pleio <- mr_pleiotropy_test(Mydata) 
pleio   #——MR egger 
View(pleio)
write.csv(pleio, file="pleio.csv")




single <- mr_leaveoneout(Mydata)
mr_leaveoneout_plot(single)   #

write.csv(single, file="single.csv")


mr_funnel_plot(res_single)#
result<-generate_odds_ratios(res)#OR

result
write.csv(result, file="result.csv")




#MRPRESSO
library(MRPRESSO)
mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data =Mydata, NbDistribution = 1000, SignifThreshold = 0.05)


