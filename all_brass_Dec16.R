#Packages
library(MasterBayes)
library(dplyr)

#Parameters
group<-"All" #or "Batch"
Plot<-"Track" #pick the batch you want
Type<-"Part"
Vrs<-c(15:40) #variables want to calculate
trait.fun<-function(Plot, type="Full"){
  if (Plot == "B30" | Plot == "Comp") {
    if (type == "Full") {traits<-c(1:14,41:43)} else if (
      type == "Part") {traits<-c(41:43)}} else if (
        Plot == "Fwest" | Plot == "Track") {
        if (type == "Full") {traits<-c(15:43)} else if (
          type == "Part") {traits<-c(41:43)}}
  
  return(traits)
} #select traits based on plots #select traits based on plots
zscore.fun<-function(df,trait,trait.name){
  temp.df<-transmute(df, temp=(trait-mean(trait))/sd(trait))
  names(temp.df)<-paste("z",trait.name,sep=".")
  temp.df<-bind_cols(df,temp.df)
  return(temp.df)
} #create and bind z-score columns to df
included_loci<-c(1:10)
number_loci<-length(included_loci)
epsilon1<-c(0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005)	#allelic dropout
epsilon2<-c(0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005)	#stochastic (misscore) error
if (Plot == "B30" | Plot == "Comp"){
  p.list<-c(2:15) #column indices for calculating z-scores
} else if (Plot == "Track" | Plot == "Fwest"){
  p.list<-c(2:15) #column indices for calculating z-scores #need to figure out
}


#Functions
zscore<-function(X){(X-mean(X))/sd(X)} 	#function to standardize traits from Austen

#Set working directory
save_location<-"~/Users/lauraleventhal/Dropbox/ISfolder_Laura_Leventhal"
setwd("~/Dropbox/ISfolder_Laura_Leventhal/Nov_7_all_sheets")

#code
if (group == "All") {
  #script for Batch dur and total flowers, looking at the rest of the pop
  BrapaG<-read.csv(paste("BrapaG_", Plot, ".csv", sep=""))
  BrapaP<-read.csv(paste("BrapaP_", Plot, ".csv", sep=""))
  BrapaG$X<-NULL
  BrapaP$X.1<-NULL
  BrapaG$sex<-NULL
  BrapaG<-merge(BrapaP,BrapaG,by="id")
  BrapaG<-subset(BrapaG, select=-c(offspring,terr,X,Y,ft,dur,tot_flwrs,sex))

  
  if (Plot == "Fwest"){
    fwest_prob<-BrapaP[grep("2808c_FWest_bG_F|2809c_FWest_bG_F|2810c_FWest_bG_F|2811c_FWest_bG_F|2814c_FWest_bG_F",BrapaP$id),]
    
    fwest_prob$terr<-paste("Fwest_bG_F")
    fwest_1<-data.frame(do.call('rbind', strsplit(as.character(fwest_prob$id),'_',fixed=TRUE)))
    
    fwest_prob<-cbind(fwest_1,fwest_prob)
    
    fwest_prob$X2<-paste("Fwest")
    
    fwest_prob$id<-NULL
    
    fwest_prob$id<-paste(fwest_prob$X1, fwest_prob$X2,fwest_prob$X3, fwest_prob$X4, sep="_")
    fwest_prob<-fwest_prob[,-c(1:4)]
    
    BrapaP<-BrapaP[-grep("2808c_FWest_bG_F|2809c_FWest_bG_F|2810c_FWest_bG_F|2811c_FWest_bG_F|2814c_FWest_bG_F",BrapaP$id),]
    
    BrapaP<-rbind(BrapaP,fwest_prob)
    
    ######same thing for BrapaG####
    fwest_probG<-BrapaG[grep("2808c_FWest_bG_F|2809c_FWest_bG_F|2810c_FWest_bG_F|2811c_FWest_bG_F|2814c_FWest_bG_F",BrapaG$id),]
    
    fwest_1G<-data.frame(do.call('rbind', strsplit(as.character(fwest_probG$id),'_',fixed=TRUE)))
    
    fwest_probG<-cbind(fwest_1G,fwest_probG)
    
    fwest_probG$X2<-paste("Fwest")
    
    fwest_probG$id<-NULL
    
    fwest_probG$id<-paste(fwest_probG$X1, fwest_probG$X2,fwest_probG$X3, fwest_probG$X4, sep="_")
    fwest_probG<-fwest_probG[,-c(1:4)]
    
    BrapaG<-BrapaG[-grep("2808c_FWest_bG_F|2809c_FWest_bG_F|2810c_FWest_bG_F|2811c_FWest_bG_F|2814c_FWest_bG_F",BrapaG$id),]
    
    BrapaG<-rbind(BrapaG,fwest_probG)
  } #end if statement over Fwest problem#dont need because did in earlier code?
  if (Plot == "B30" | Plot == "Comp"){
    BrapaG<-BrapaG[-grep("bA2",BrapaG$timevar),] 
    BrapaP<-BrapaP[-grep("bA2",BrapaP$timevar),]
  }#end of if statement getting rid of bA2 in comp and B30
  
  if (Plot == "B30") {
    flwcount<-read.csv(paste("b30phenZerosJune7-2013.csv", sep=""))
    flwcount$id<-paste("B30",flwcount$pos, sep="_")
  } #end of if statement adding pheno data to B30
  if (Plot == "Comp") {
    flwcount<-read.csv(paste("compphenZerosJune7-2013.csv", sep=""))
    flwcount$id<-paste("Comp",flwcount$pos, sep="_")
  }#end of if statement adding pheno data to Comp
  if (Plot == "Track") {
    flwcount<-read.csv(paste("tkphenZerosJune7-2013.csv", sep=""))
    flwcount$id<-paste("Track",flwcount$pos, sep="_")
  }#end of if statement adding pheno data to Track
  if (Plot == "Fwest") {
    flwcount<-read.csv(paste("fwphenZerosJune7-2013.csv", sep=""))
    flwcount$id<-paste("Fwest",flwcount$pos, sep="_")
  }#end of if statement adding pheno data to Fwest
  
  
  BrapaP_male<-BrapaP[grep("Male", BrapaP$sex),]
  BrapaP_male<-merge(BrapaP_male,flwcount,by="id")
  
  if (Plot == "B30" | Plot == "Comp"){
    BrapaP_test<-BrapaP_male
    flwr_test<-BrapaP_male
    BrapaP_test<-select(BrapaP_male, f.178:f.241)
    flwr_test<-select(BrapaP_male, id:sex)
  }#end of if statment for mating opporunity trait for B30 and Comp
  if (Plot == "Fwest" | Plot == "Track"){
    BrapaP_test<-select(BrapaP_male, d228:d286)
    flwr_test<-select(BrapaP_male, id:sex)
  }#end of if statment for mating opporunity trait for Fwest and Track
  
  for (k in 1:nrow(BrapaP_test)){
    SimplifiedIntervals<-ifelse(BrapaP_test[k,]>0, 1, 0)
    BrapaP_test[k,]<-SimplifiedIntervals }#end loop for creation of dur and sum column
  BrapaP_male<-cbind(BrapaP_test,flwr_test)
  
  if (Plot == "B30" | Plot == "Comp"){
    BrapaP_male$batchdurA<-rowMeans(BrapaP_male[,3:4])
    BrapaP_male$batchdurB<-rowMeans(BrapaP_male[,5:6])
    BrapaP_male$batchdurB2<-rowMeans(BrapaP_male[,6:8])
    BrapaP_male$batchdurC<-rowMeans(BrapaP_male[,8:9])
    BrapaP_male$batchdurD<-rowMeans(BrapaP_male[,9:10])
    BrapaP_male$batchdurD2<-rowMeans(BrapaP_male[,10:11])
    BrapaP_male$batchdurE<-rowMeans(BrapaP_male[,11:19])
    
    BrapaP_male<-BrapaP_male[,-c(1:19)]
    BrapaP_male<-merge(BrapaP_male,flwcount,by="id")
    
    BrapaP_male$batchsumA<-rowSums(BrapaP_male[20:21]) #bA
    BrapaP_male$batchsumB<-rowSums(BrapaP_male[22:23])#bB
    BrapaP_male$batchsumB2<-rowSums(BrapaP_male[23:25])#bB2
    BrapaP_male$batchsumC<-rowSums(BrapaP_male[25:26])#bC
    BrapaP_male$batchsumD<-rowSums(BrapaP_male[26:27])#bD
    BrapaP_male$batchsumD2<-rowSums(BrapaP_male[27:28])#bD2
    BrapaP_male$batchsumE<-rowSums(BrapaP_male[28:36])#bE
    
    BrapaP_male<-BrapaP_male[,-c(2:10,18:40)]
    
  } #end of loop calculating row means and row sums
  
  if (Plot == "Fwest" | Plot == "Track"){
    BrapaP_male$batchdurG<-rowMeans(BrapaP_male[,1:3])
    BrapaP_male$batchdurH<-rowMeans(BrapaP_male[,3:5])
    BrapaP_male$batchdurI<-rowMeans(BrapaP_male[,5:7])
    BrapaP_male$batchdurJ<-rowMeans(BrapaP_male[,7:9])
    BrapaP_male$batchdurK<-rowMeans(BrapaP_male[,9:11])
    BrapaP_male$batchdurL<-rowMeans(BrapaP_male[,11:13])
    BrapaP_male$batchdurM<-rowMeans(BrapaP_male[,13:15])
    BrapaP_male$batchdurN<-rowMeans(BrapaP_male[,15:17])
    BrapaP_male$batchdurN2<-rowMeans(BrapaP_male[,17:19])
    BrapaP_male$batchdurO<-rowMeans(BrapaP_male[,19:20])
    BrapaP_male$batchdurO2<-rowMeans(BrapaP_male[,20:21])
    BrapaP_male$batchdurP<-rowMeans(BrapaP_male[,21:22])
    BrapaP_male$batchdurQ<-rowMeans(BrapaP_male[,22:23])
    
    BrapaP_male<-BrapaP_male[,-c(1:24)]
    BrapaP_male<-merge(BrapaP_male,flwcount,by="id")
    
    BrapaP_male$batchsumG<-rowSums(BrapaP_male[,25:27])
    BrapaP_male$batchsumH<-rowSums(BrapaP_male[,27:29])
    BrapaP_male$batchsumI<-rowSums(BrapaP_male[,29:31])
    BrapaP_male$batchsumJ<-rowSums(BrapaP_male[,31:33])
    BrapaP_male$batchsumK<-rowSums(BrapaP_male[,33:35])
    BrapaP_male$batchsumL<-rowSums(BrapaP_male[,35:37])
    BrapaP_male$batchsumM<-rowSums(BrapaP_male[,37:39])
    BrapaP_male$batchsumN<-rowSums(BrapaP_male[,39:41])
    BrapaP_male$batchsumN2<-rowSums(BrapaP_male[,41:43])
    BrapaP_male$batchsumO<-rowSums(BrapaP_male[,43:44])
    BrapaP_male$batchsumO2<-rowSums(BrapaP_male[,44:45])
    BrapaP_male$batchsumP<-rowSums(BrapaP_male[,45:46])
    BrapaP_male$batchsumQ<-rowSums(BrapaP_male[,46:47])
    
    BrapaP_male<-BrapaP_male[,-c(2:10,24:50)] ## need to see which ones to delete
    BrapaP<-merge(BrapaP_male,BrapaP,by="id",all.y=T)
  }#creates batchsum and batch dur proportion variables in Fwest or Track ####fix for fwest and track later
  
 
 
  
  ###z score
  if (Plot == "B30" | Plot == "Comp"){
    traits<-trait.fun(Plot, Type)
    for (t in traits){
      if (t == 1){batchdurA<-"Duration for batch A"}
      if (t == 2){batchdurB<-"Duration for batch B"}
      if (t == 3){batchdurB2<-"Duration for batch B2"}
      if (t == 4){batchdurC<-"Duration for batch C"}
      if (t == 5){batchdurD<-"Duration for batch D"}
      if (t == 6){batchdurD2<-"Duration for batch D2"}
      if (t == 7){batchdurE<-"Duration for batch E"}
      if (t == 8){batchsumA<-"Total flowers for batch A"}
      if (t == 9){batchsumB<-"Total flowers for batch B"}
      if (t == 10){batchsumB2<-"Total flowers for batch B2"}
      if (t == 11){batchsumC<-"Total flowers for batch C"}
      if (t == 12){batchsumD<-"Total flowers for batch D"}
      if (t == 13){batchsumD2<-"Total flowers for batch D2"}
      if (t == 14){batchsumE<-"Total flowers for batch E"}
      
      if (t == 15){batchdurG<-"Duration for batch G"}
      if (t == 16){batchdurH<-"Duration for batch H"}
      if (t == 17){batchdurI<-"Duration for batch I"}
      if (t == 18){batchdurJ<-"Duration for batch J"}
      if (t == 19){batchdurK<-"Duration for batch K"}
      if (t == 20){batchdurL<-"Duration for batch L"}
      if (t == 21){batchdurM<-"Duration for batch M"}
      if (t == 22){batchdurN<-"Duration for batch N"}
      if (t == 23){batchdurN2<-"Duration for batch N2"}
      if (t == 24){batchdurO<-"Duration for batch O"}
      if (t == 25){batchdurO2<-"Duration for batch O2"}
      if (t == 26){batchdurP<-"Duration for batch P"}
      if (t == 27){batchdurQ<-"Duration for batch Q"}
      if (t == 28){batchsumG<-"Total flowers for batch G"}
      if (t == 29){batchsumH<-"Total flowers for batch H"}
      if (t == 30){batchsumI<-"Total flowers for batch I"}
      if (t == 31){batchsumJ<-"Total flowers for batch J"}
      if (t == 32){batchsumK<-"Total flowers for batch K"}
      if (t == 33){batchsumL<-"Total flowers for batch L"}
      if (t == 34){batchsumM<-"Total flowers for batch M"}
      if (t == 35){batchsumN<-"Total flowers for batch N"}
      if (t == 36){batchsumN2<-"Total flowers for batch N2"}
      if (t == 37){batchsumO<-"Total flowers for batch O"}
      if (t == 38){batchsumO2<-"Total flowers for batch O2"}
      if (t == 39){batchsumP<-"Total flowers for batch P"}
      if (t == 40){batchsumQ<-"Total flowers for batch Q"}
      
      if (t == 41){tot_flwrs<-"Total flowers for a season"}
      if (t == 42){tot_flwrs<-"Duration for a season"}
      if (t == 43){var.dist<-"X and Ycoordinates of an individual"}
      
    } #end for loop
    #grep out all the males so you can z score all of their traits
    #BrapaP_male<-BrapaP[grep("Male",BrapaP$sex),]
    foo<-data.frame(do.call('rbind', strsplit(as.character(BrapaP_male$id),'_',fixed=TRUE)))
    BrapaP_male<-bind_cols(BrapaP_male, foo[2])
    #names(BrapaP_male)[dim(BrapaP_male)[2]]<-"plantnum"
    
    #Convert the traits we want to estimate selection on to z-scores 
    trait_file<-BrapaP_male
    p.list<-c(2:15)
    for (p in p.list){
      BrapaP_male<-zscore.fun(BrapaP_male,BrapaP_male[,p],paste(names(BrapaP_male)[p],sep="."))
    }
    
  }#end of loop for zscore for B30 and Comp
  
  if (Plot == "Fwest" | Plot == "Track"){
    traits<-trait.fun(Plot, Type)
    for (t in traits){
      if (t == 1){batchdurA<-"Duration for batch A"}
      if (t == 2){batchdurB<-"Duration for batch B"}
      if (t == 3){batchdurB2<-"Duration for batch B2"}
      if (t == 4){batchdurC<-"Duration for batch C"}
      if (t == 5){batchdurD<-"Duration for batch D"}
      if (t == 6){batchdurD2<-"Duration for batch D2"}
      if (t == 7){batchdurE<-"Duration for batch E"}
      if (t == 8){batchsumA<-"Total flowers for batch A"}
      if (t == 9){batchsumB<-"Total flowers for batch B"}
      if (t == 10){batchsumB2<-"Total flowers for batch B2"}
      if (t == 11){batchsumC<-"Total flowers for batch C"}
      if (t == 12){batchsumD<-"Total flowers for batch D"}
      if (t == 13){batchsumD2<-"Total flowers for batch D2"}
      if (t == 14){batchsumE<-"Total flowers for batch E"}
      
      if (t == 15){batchdurG<-"Duration for batch G"}
      if (t == 16){batchdurH<-"Duration for batch H"}
      if (t == 17){batchdurI<-"Duration for batch I"}
      if (t == 18){batchdurJ<-"Duration for batch J"}
      if (t == 19){batchdurK<-"Duration for batch K"}
      if (t == 20){batchdurL<-"Duration for batch L"}
      if (t == 21){batchdurM<-"Duration for batch M"}
      if (t == 22){batchdurN<-"Duration for batch N"}
      if (t == 23){batchdurN2<-"Duration for batch N2"}
      if (t == 24){batchdurO<-"Duration for batch O"}
      if (t == 25){batchdurO2<-"Duration for batch O2"}
      if (t == 26){batchdurP<-"Duration for batch P"}
      if (t == 27){batchdurQ<-"Duration for batch Q"}
      if (t == 28){batchsumG<-"Total flowers for batch G"}
      if (t == 29){batchsumH<-"Total flowers for batch H"}
      if (t == 30){batchsumI<-"Total flowers for batch I"}
      if (t == 31){batchsumJ<-"Total flowers for batch J"}
      if (t == 32){batchsumK<-"Total flowers for batch K"}
      if (t == 33){batchsumL<-"Total flowers for batch L"}
      if (t == 34){batchsumM<-"Total flowers for batch M"}
      if (t == 35){batchsumN<-"Total flowers for batch N"}
      if (t == 36){batchsumN2<-"Total flowers for batch N2"}
      if (t == 37){batchsumO<-"Total flowers for batch O"}
      if (t == 38){batchsumO2<-"Total flowers for batch O2"}
      if (t == 39){batchsumP<-"Total flowers for batch P"}
      if (t == 40){batchsumQ<-"Total flowers for batch Q"}
      
      if (t == 41){tot_flwrs<-"Total flowers for a season"}
      if (t == 42){dur<-"Duration for a season"}
      if (t == 43){var.dist<-"X and Ycoordinates of an individual"}
      
    } #end for loop
    #grep out all the males so you can z score all of their traits
    #BrapaP_male<-BrapaP[grep("Male",BrapaP$sex),]
    foo<-data.frame(do.call('rbind', strsplit(as.character(BrapaP_male$id),'_',fixed=TRUE)))
    BrapaP_male<-bind_cols(BrapaP_male, foo[2])
    #names(BrapaP_male)[dim(BrapaP_male)[2]]<-"plantnum"
    
    #Convert the traits we want to estimate selection on to z-scores 
    trait_file<-BrapaP_male
    p.list<-c(2:27)
    for (p in p.list){
      BrapaP_male<-zscore.fun(BrapaP_male,BrapaP_male[,p],paste(names(BrapaP_male)[p],sep="."))
    }
    
  }#end of loop for zscore for Fwest and Track
  #then
  if (Plot == "B30" | Plot == "Comp"){
    #BrapaP<-BrapaP[,-c(2:27)]
    BrapaP_male<-BrapaP_male[,-c(2:16)]
    BrapaP<-merge(BrapaP, BrapaP_male, by = "id",all.x=T)
    #for some reason, this full_join is going way over board... like it should be about 11,000 cells individuals from BrapaP and 1,500 from BrapaP_males to get 12,500. But instead i get 82,902! I am going to see if erasing duplciates fixes this
    BrapaP<-distinct(BrapaP)
    # it worked.. woohoo! but still werid...
  }
  if (Plot == "Fwest" | Plot == "Track"){
    BrapaP<-BrapaP[,-c(2:27)]
    BrapaP_male<-BrapaP_male[,-c(2:28)]
    BrapaP<-merge(BrapaP, BrapaP_male, by = "id",all.x=T)
    #for some reason, this full_join is going way over board... like it should be about 11,000 cells individuals from BrapaP and 1,500 from BrapaP_males to get 12,500. But instead i get 82,902! I am going to see if erasing duplciates fixes this
    BrapaP<-distinct(BrapaP)
    # it worked.. woohoo! but still werid...
  }

  
  BrapaP_for_MB<-BrapaP
  BrapaG_for_MB<-BrapaG
  
  ####Restrictions####
  #Restriction 0: Offspring cannot be parents.
  res.osnotpar<-expression(varPed(x="offspring", restrict=0))
  #Restriction 1: Offspring's "terr" (i.e., it's mothers ID) must match its mother's terr (i.e. id)
  res.osterr<-expression(varPed(x="terr", gender="Female", relational="OFFSPRING", restrict="=="))
  #Restriction 3: Father must have been 'active' (flowering) during the tagging interval in which offpsring was produced
  res.phen<-expression(varPed(x="timevar", gender="Male", relational="OFFSPRING", restrict="=="))	#Only parent rows where timevar == os timevar can have non-zero paternity
  
  ####variables####
  for (v in Vrs){
    if (v == 1){var1<-expression(varPed(x="z.batchdurA", gender="Male"))}
    if (v == 2){var2<-expression(varPed(x="z.batchdurB", gender="Male"))}
    if (v == 3){var3<-expression(varPed(x="z.batchdurB2", gender="Male"))}
    if (v == 4){var4<-expression(varPed(x="z.batchdurC", gender="Male"))}
    if (v == 5){var5<-expression(varPed(x="z.batchdurD", gender="Male"))}
    if (v == 6){var6<-expression(varPed(x="z.batchdurD2", gender="Male"))}
    if (v == 7){var7<-expression(varPed(x="z.batchdurE", gender="Male"))}
    if (v == 8){var8<-expression(varPed(x="z.batchsumA", gender="Male"))}
    if (v == 9){var9<-expression(varPed(x="z.batchsumB", gender="Male"))}
    if (v == 10){var10<-expression(varPed(x="z.batchsumB2", gender="Male"))}
    if (v == 11){var11<-expression(varPed(x="z.batchsumC", gender="Male"))}
    if (v == 12){var12<-expression(varPed(x="z.batchsumD", gender="Male"))}
    if (v == 13){var13<-expression(varPed(x="z.batchsumD2", gender="Male"))}
    if (v == 14){var14<-expression(varPed(x="z.batchsumE", gender="Male"))}
    if (v == 15){var15<-expression(varPed(x="z.batchsumG", gender="Male"))}
    if (v == 16){var16<-expression(varPed(x="z.batchsumH", gender="Male"))}
    if (v == 17){var17<-expression(varPed(x="z.batchsumI", gender="Male"))}
    if (v == 18){var18<-expression(varPed(x="z.batchsumJ", gender="Male"))}
    if (v == 19){var19<-expression(varPed(x="z.batchsumK", gender="Male"))}
    if (v == 20){var20<-expression(varPed(x="z.batchsumL", gender="Male"))}
    if (v == 21){var21<-expression(varPed(x="z.batchsumM", gender="Male"))}
    if (v == 22){var22<-expression(varPed(x="z.batchsumN", gender="Male"))}
    if (v == 23){var23<-expression(varPed(x="z.batchsumN2", gender="Male"))}
    if (v == 24){var24<-expression(varPed(x="z.batchsumO", gender="Male"))}
    if (v == 25){var25<-expression(varPed(x="z.batchsumO2", gender="Male"))}
    if (v == 26){var26<-expression(varPed(x="z.batchsumP", gender="Male"))}
    if (v == 27){var27<-expression(varPed(x="z.batchsumQ", gender="Male"))}
    if (v == 28){var28<-expression(varPed(x="z.batchdurG", gender="Male"))}
    if (v == 29){var29<-expression(varPed(x="z.batchdurH", gender="Male"))}
    if (v == 30){var30<-expression(varPed(x="z.batchdurI", gender="Male"))}
    if (v == 31){var31<-expression(varPed(x="z.batchdurJ", gender="Male"))}
    if (v == 32){var32<-expression(varPed(x="z.batchdurK", gender="Male"))}
    if (v == 33){var33<-expression(varPed(x="z.batchdurL", gender="Male"))}
    if (v == 34){var34<-expression(varPed(x="z.batchdurM", gender="Male"))}
    if (v == 35){var35<-expression(varPed(x="z.batchdurN", gender="Male"))}
    if (v == 36){var36<-expression(varPed(x="z.batchdurN2", gender="Male"))}
    if (v == 37){var37<-expression(varPed(x="z.batchdurO", gender="Male"))}
    if (v == 38){var38<-expression(varPed(x="z.batchdurO2", gender="Male"))}
    if (v == 39){var39<-expression(varPed(x="z.batchdurP", gender="Male"))}
    if (v == 40){var40<-expression(varPed(x="z.batchdurQ", gender="Male"))}
    
    if (v == 41){var41<-expression(varPed(x="z.tot_flwrs", gender="Male"))} 
    if (v == 42){var42<-expression(varPed(x="z.dur_flwrs", gender="Male"))} 
    if (v == 43){var.dist<-expression(varPed(x=c("X", "Y"), gender="Male", relational="MATE"))}
  } #end for loop
  Vrs.list<-lapply(ls(pattern="var*"),get)
  
  BrapaG_alleles<-BrapaG_for_MB
  BrapaG_alleles$timevar<-NULL
  alleles<-extractA(BrapaG_alleles)#estimate true allele frequencies based on largest available sample
  #only need gender in phenotype not genotype
  BrapaG_for_MB$timevar<-NULL
  
  
  ####GdP####
  GdP<-GdataPed(G=BrapaG_for_MB, categories=NULL, perlocus=T)	#need to specify 'per locus' to allow variable error rates by locus.
  
  mod_subfolder<-"data_files_3/files_generated_by_PatAnalysis_Script"
  #MODEL: Distance, start date, total flowers, duration of flowering, skew, etc. by batch
  mod_subfolder<-"data_files/Analysis_Dec_19_LL"	#This empty folder should already be present.
  
  
  
  ####PdP####
  #PdP<-if(sum(Vrs) == 105){PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var1,var2, var3,var4,var5,var6,var7,var8,var9,var10,var11,var12,var13,var14), data=BrapaP_for_MB)} else if (
 #  sum(Vrs) == 715) {PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var15,var16, var17,var18,var19,var20,var21,var22,var23,var24,var25,var26,var27,var28,var29,var30,var31,var32,var33,var34,var35,var36,var37,var38,var39,var40), data=BrapaP_for_MB)} else if (
 #    sum(Vrs) == 126){PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var41,var42,var.dist), data=BrapaP_for_MB)} else if (
   #     sum(Vrs) == 231) {PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var1,var2, var3,var4,var5,var6,var7,var8,var9,var10,var11,var12,var13,var14,var41,var42,var.dist), data=BrapaP_for_MB)} else if (
  #        sum(Vrs) == 841) {PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var15,var16, var17,var18,var19,var20,var21,var22,var23,var24,var25,var26,var27,var28,var29,var30,var31,var32,var33,var34,var35,var36,var37,var38,var39,var40,var41,var42,var.dist), data=BrapaP_for_MB)}
  
PdP<-PdataPed(formula=list(res.osnotpar, res.osterr, res.phen,var15, var16, var17, var18, var19, var20, var21, var22, var23, var24, var25, var26, var27, var28, var30, var31, var32, var33, var34, var35, var36, var37, var38, var39, var40 ), data=BrapaP_for_MB)
  #  system is computationally singular: reciprocal condition number = 2.48059e-19

 # var15, var16, var17, var18, var19, var20, var21, var22, var23, var24, var25, var26, var27, var28, var30, var31, var32, var33, var34, var35, var36, var37, var38, var39, var40
  #var1, var2, var3, var4, var5, var6, var7, var8, var9, var10, var11, var12, var13, var14
  
  ####chain####
  
  sP.chain1<-startPed(estG=F, E1=epsilon1, E2=epsilon2, A=alleles, estE1=F, estE2=F, estA=F)
  sP.chain2<-startPed(estG=F, E1=epsilon1, E2=epsilon2, A=alleles, estE1=F, estE2=F, estA=F, beta=rep(0, times=1))
  sP.chain3<-startPed(estG=F, E1=epsilon1, E2=epsilon2, A=alleles, estE1=F, estE2=F, estA=F, beta=rnorm(n=1, mean=c(-0.1, 0.1), sd=0.05))####
  ##three chains with different betas
  
  #find a tunePed(beta=??) that gives acceptance rate ~ 20% - 50% (aim for 30 to 35).  Increasing beta generally reduces acceptance rate
  tP<-tunePed(beta=0.325)
  #m.quick<-MCMCped(PdP=PdP12, GdP=GdP, sP=sP.chain1, verbose=T, mm.tol=3, nitt=5000, tP=tP)
  #autocorr(m.quick$beta)	
  pP<-priorPed(USsire=list(mu=log(5), sigma=0.5))
  #BrapaP_for_MB$population_num=NULL
  #BrapaP_for_MB$type_of_site=NULL
  
  ch1<-MCMCped(PdP=PdP, GdP=GdP, sP=sP.chain1, tP=tP,verbose=T,mm.tol=5, nitt=10000+(10000*10), thin=10*10, burnin=10000)
  ch1.dur<-MCMCped(PdP=PdP, GdP=GdP, sP=sP.chain1, tP=tP,verbose=T,mm.tol=5, nitt=10000+(10000*10), thin=10*10, burnin=10000)
  mod_subfolder<-"data_files_3/files_generated_by_PatAnalysis_Script"
  save(ch1, file=paste(save_location, mod_subfolder, "m12.ch1.Rdata", sep="/"))
  write.csv(ch1, file = "Fwest_sum_slices_ch1.csv")
  save(ch1, file=paste(save_location,mod_subfolder, sep="/"))
  save(ch1, file=paste(save_location, mod_subfolder, "ch1.Fwest.sum_slices.RData", sep="/"))
summary(ch1$beta)
quartz()
summary(ch1.dur$beta)
} else if (group == "Batch") {
  #script for doing analysis just by batch 
  BrapaG<-read.csv(paste("BrapaG_", Plot, ".csv", sep=""))
  BrapaP<-read.csv(paste("BrapaP_", Plot, ".csv", sep=""))
  BrapaG$X<-NULL
  BrapaP$X.1<-NULL
  BrapaG$sex<-NULL
  BrapaG<-merge(BrapaP,BrapaG,by="id")
  BrapaG<-subset(BrapaG, select=-c(offspring,terr,X,Y,ft,dur,tot_flwrs,sex))

}
  
  