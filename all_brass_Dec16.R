#Packages
library(MasterBayes)
library(dplyr)
library(tidyr)

#Parameters
group<-"All" #or "Batch"
Plot<-"Track" #pick the batch you want
Type<-"Part"
extra<-"Batch"
Distance<-"Not-Fancy" #Fancy-meaning special clump distances 
Vrs<-c(1:51) #variables want to calculate
trait.fun<-function(Plot, type="Full"){
  if (Plot == "B30" | Plot == "Comp") {
    if (type == "Full") {traits<-c(1:14,41:50)} else if (
      type == "Part") {traits<-c(41:50)}} else if (
        Plot == "Fwest" | Plot == "Track") {
        if (type == "Full") {traits<-c(15:50)} else if (
          type == "Part") {traits<-c(41:50)}}
  
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
  p.list<-c(2:18) #column indices for calculating z-scores
} else if (Plot == "Track" | Plot == "Fwest"){
  p.list<-c(2:33) #column indices for calculating z-scores #need to figure out
}


#Functions
zscore<-function(X){(X-mean(X))/sd(X)} 	#function to standardize traits from Austen

#Set working directory
save_location<-"~/Users/lauraleventhal/Dropbox/ISfolder_Laura_Leventhal"
setwd("~/Dropbox/ISfolder_Laura_Leventhal/All_brass_files")

#code
if (group == "All") {
  #script for Batch dur and total flowers, looking at the rest of the pop
  BrapaG<-read.csv(paste("BrapaG_", Plot, ".csv", sep=""))
  BrapaP<-read.csv(paste("BrapaP_", Plot, ".csv", sep=""))
  BrapaG$X<-NULL
  BrapaP$X.1<-NULL
  BrapaG$sex<-NULL
  BrapaG<-merge(BrapaP,BrapaG,by="id")
  if (Plot == "B30" | Plot == "Comp"){
    BrapaG<-subset(BrapaG, select=-c(offspring,terr,X,Y,ft,dur,tot_flwrs,sex))
  }
  if (Plot == "Fwest" | Plot == "Track"){
    BrapaG<-subset(BrapaG, select=-c(offspring,terr,X,Y,ft,dur_inpop,dur_total,tot_flwrs,sex))
  }

  
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
    flwr_test<-select(BrapaP_male, id:dur_total)
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
    
    BrapaP_male<-BrapaP_male[,-c(2:10,18:36,38:40)]
    BrapaP<-merge(BrapaP,BrapaP_male,by="id",all.x=T)
    BrapaP<-distinct(BrapaP)
    BrapaP<-BrapaP[,-c(7)]
    
  } #end of loop calculating row means and row sums
  
  if (Plot == "Fwest" | Plot == "Track"){
    BrapaP_male$batchdurG<-rowMeans(BrapaP_male[,1:3])
    BrapaP_male$batchdurH<-rowMeans(BrapaP_male[,3:5])
    BrapaP_male$batchdurI<-rowMeans(BrapaP_male[,5:7])
    BrapaP_male$batchdurJ<-rowMeans(BrapaP_male[,7:9])
    BrapaP_male$batchdurK<-rowMeans(BrapaP_male[,9:11])
    BrapaP_male$batchdurL<-rowMeans(BrapaP_male[,11:13])
    BrapaP_male$batchdurM<-rowMeans(BrapaP_male[,13:15])
    BrapaP_male$batchdurM2<-rowMeans(BrapaP_male[,15:17])
    BrapaP_male$batchdurN<-rowMeans(BrapaP_male[,17:19])
    BrapaP_male$batchdurN2<-rowMeans(BrapaP_male[,19:20])
    BrapaP_male$batchdurO<-rowMeans(BrapaP_male[,20:21])
    BrapaP_male$batchdurO2<-rowMeans(BrapaP_male[,21:22])
    BrapaP_male$batchdurP<-rowMeans(BrapaP_male[,22:23])
    BrapaP_male$batchdurQ<-rowMeans(BrapaP_male[,23:24])
    
    BrapaP_male<-BrapaP_male[,-c(1:24)]
    BrapaP_male<-merge(BrapaP_male,flwcount,by="id")
    
    BrapaP_male$batchsumG<-rowSums(BrapaP_male[,26:28])
    BrapaP_male$batchsumH<-rowSums(BrapaP_male[,28:30])
    BrapaP_male$batchsumI<-rowSums(BrapaP_male[,30:32])
    BrapaP_male$batchsumJ<-rowSums(BrapaP_male[,32:34])
    BrapaP_male$batchsumK<-rowSums(BrapaP_male[,34:36])
    BrapaP_male$batchsumL<-rowSums(BrapaP_male[,36:38])
    BrapaP_male$batchsumM<-rowSums(BrapaP_male[,38:40])
    BrapaP_male$batchsumM2<-rowSums(BrapaP_male[,40:42])
    BrapaP_male$batchsumN<-rowSums(BrapaP_male[,42:44])
    BrapaP_male$batchsumN2<-rowSums(BrapaP_male[,44:45])
    BrapaP_male$batchsumO<-rowSums(BrapaP_male[,45:46])
    BrapaP_male$batchsumO2<-rowSums(BrapaP_male[,46:47])
    BrapaP_male$batchsumP<-rowSums(BrapaP_male[,47:48])
    BrapaP_male$batchsumQ<-rowSums(BrapaP_male[,48:49])
    
    BrapaP_male<-BrapaP_male[,-c(2:11,25:50)] 
    BrapaP<-merge(BrapaP,BrapaP_male,by="id",all.x=T)
    BrapaP<-distinct(BrapaP)
  }#creates batchsum and batch dur proportion variables in Fwest or Track 
  
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
      if (t == 22){batchdurM<-"Duration for batch M2"}
      if (t == 23){batchdurN<-"Duration for batch N"}
      if (t == 24){batchdurN2<-"Duration for batch N2"}
      if (t == 25){batchdurO<-"Duration for batch O"}
      if (t == 26){batchdurO2<-"Duration for batch O2"}
      if (t == 27){batchdurP<-"Duration for batch P"}
      if (t == 28){batchdurQ<-"Duration for batch Q"}
      if (t == 29){batchsumG<-"Total flowers for batch G"}
      if (t == 30){batchsumH<-"Total flowers for batch H"}
      if (t == 31){batchsumI<-"Total flowers for batch I"}
      if (t == 32){batchsumJ<-"Total flowers for batch J"}
      if (t == 33){batchsumK<-"Total flowers for batch K"}
      if (t == 34){batchsumL<-"Total flowers for batch L"}
      if (t == 35){batchsumM<-"Total flowers for batch M"}
      if (t == 36){batchsumN<-"Total flowers for batch N"}
      if (t == 37){batchsumN2<-"Total flowers for batch N2"}
      if (t == 38){batchsumO<-"Total flowers for batch O"}
      if (t == 39){batchsumO2<-"Total flowers for batch O2"}
      if (t == 40){batchsumP<-"Total flowers for batch P"}
      if (t == 41){batchsumQ<-"Total flowers for batch Q"}
    
      if (t == 42){tot_flwrs<-"Total flowers for a season"}
      if (t == 43){dur<-"Duration for a season"}
      if (t == 44){sd<-"Start Date"}
      if (t == 45){sd_field<-"Start date in the field"}
      if (t == 46){var.dist<-"X and Ycoordinates of an individual"}
      
    } #end for loop
    #grep out all the males so you can z score all of their traits
    #BrapaP_male<-BrapaP[grep("Male",BrapaP$sex),]
   # foo<-data.frame(do.call('rbind', strsplit(as.character(BrapaP_male$id),'_',fixed=TRUE)))
    #BrapaP_male<-bind_cols(BrapaP_male, foo[2])
    #names(BrapaP_male)[dim(BrapaP_male)[2]]<-"plantnum"
    
    #Convert the traits we want to estimate selection on to z-scores 
    BrapaP_male2<-BrapaP[grep("Male",BrapaP$sex),]
    trait_file<-select(BrapaP_male2, id, dur, tot_flwrs, batchdurA:batchsumE)
    p.list<-c(2:18)
    for (p in p.list){
      trait_file<-zscore.fun(trait_file,trait_file[,p],paste(names(trait_file)[p],sep="."))
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
      if (t == 22){batchdurM2<-"Duration for batch M2"}
      if (t == 23){batchdurN<-"Duration for batch N"}
      if (t == 24){batchdurN2<-"Duration for batch N2"}
      if (t == 25){batchdurO<-"Duration for batch O"}
      if (t == 26){batchdurO2<-"Duration for batch O2"}
      if (t == 27){batchdurP<-"Duration for batch P"}
      if (t == 28){batchdurQ<-"Duration for batch Q"}
      if (t == 29){batchsumG<-"Total flowers for batch G"}
      if (t == 30){batchsumH<-"Total flowers for batch H"}
      if (t == 31){batchsumI<-"Total flowers for batch I"}
      if (t == 32){batchsumJ<-"Total flowers for batch J"}
      if (t == 33){batchsumK<-"Total flowers for batch K"}
      if (t == 34){batchsumL<-"Total flowers for batch L"}
      if (t == 35){batchsumM<-"Total flowers for batch M"}
      if (t == 36){batchsumM2<-"Total flowers for batch M2"}
      if (t == 37){batchsumN<-"Total flowers for batch N"}
      if (t == 38){batchsumN2<-"Total flowers for batch N2"}
      if (t == 39){batchsumO<-"Total flowers for batch O"}
      if (t == 40){batchsumO2<-"Total flowers for batch O2"}
      if (t == 41){batchsumP<-"Total flowers for batch P"}
      if (t == 42){batchsumQ<-"Total flowers for batch Q"}
      
      if (t == 43){tot_flwrs<-"Total flowers for a season"}
      if (t == 44){dur<-"Duration for a season"}
      if (t == 45){dur_inpop<-"Duration in the population"}#for fwest and track only
      if (t == 46){dur_total<-"Duration in and out of the population"}#for fwest and track only
      if (t == 47){sd<-"Start Date"}
      if (t == 48){sd_field<-"Start date in the field"}
      if (t == 49){var.dist<-"X and Ycoordinates of an individual"}
      
    } #end for loop
    #grep out all the males so you can z score all of their traits
    #BrapaP_male<-BrapaP[grep("Male",BrapaP$sex),]
    #foo<-data.frame(do.call('rbind', strsplit(as.character(BrapaP_male$id),'_',fixed=TRUE)))
   # BrapaP_male<-bind_cols(BrapaP_male, foo[2])
    #names(BrapaP_male)[dim(BrapaP_male)[2]]<-"plantnum"
    #make same as above
    #Convert the traits we want to estimate selection on to z-scores
    BrapaP_male2<-BrapaP[grep("Male",BrapaP$sex),]
    trait_file<-select(BrapaP_male2, id, tot_flwrs, dur_inpop:batchsumQ)
    p.list<-c(2:33)
    for (p in p.list){
      trait_file<-zscore.fun(trait_file,trait_file[,p],paste(names(trait_file)[p],sep="."))
    }
    
  }#end of loop for zscore for Fwest and Track
  #then
  if (Plot == "B30" | Plot == "Comp"){
    #BrapaP<-BrapaP[,-c(2:27)]
    trait_file<-trait_file[,-c(2:18)]
    BrapaP<-merge(BrapaP, trait_file, by = "id",all.x=T)
    BrapaP<-distinct(BrapaP)
    BrapaP<-BrapaP[,-c(7,8,10:24)]
    #for some reason, this full_join is going way over board... like it should be about 11,000 cells individuals from BrapaP and 1,500 from BrapaP_males to get 12,500. But instead i get 82,902! I am going to see if erasing duplciates fixes this 
    # it worked.. woohoo! but still werid...
  }
  if (Plot == "Fwest" | Plot == "Track"){
    BrapaP<-BrapaP[,-c(7:8,10:40)]
    trait_file<-trait_file[,-c(2:33)]
    BrapaP<-merge(BrapaP, trait_file, by = "id",all.x=T)
    #for some reason, this full_join is going way over board... like it should be about 11,000 cells individuals from BrapaP and 1,500 from BrapaP_males to get 12,500. But instead i get 82,902! I am going to see if erasing duplciates fixes this
    BrapaP<-distinct(BrapaP)
    # it worked.. woohoo! but still werid...
  }
if (Distance == "Fancy"){
  if (Plot == "B30" | Plot == "Track"){
    large_clump_info<-read.csv("EvenSite_LrgClump.csv")
    medium_clump_info<-read.csv("EvenSite_MedClump.csv")
    small_clump_info<-read.csv("EvenSite_SmClump.csv")
    
    male2<-BrapaP[grep("Male", BrapaP$sex),]
    male3<-data.frame(do.call('rbind', strsplit(as.character(male2$id),'_',fixed=TRUE)))
    colnames(male3)[2]<-"pat_pos"
    male4<-merge(male3,small_clump_info,by="pat_pos")
    male5<-merge(male4,medium_clump_info,by="pat_pos")
    male6<-merge(male5,large_clump_info,by="pat_pos")
    male6$id<- paste(male6$X1,male6$pat_pos, sep="_")
    male6<-male6[-c(1,2,4:9)]
    
    female1<-BrapaP[grep("Female", BrapaP$sex),]
    female1<-female1[,-c(25:27)]
    female2<-data.frame(do.call('rbind', strsplit(as.character(female1$id),'_',fixed=TRUE)))
    female3<- cbind(female2,med_clump=rep(female2$X3))
    female4<- cbind(female3,lrg_clump=rep(female3$X3))
    female5<- cbind(female4,sm_clump=rep(female4$X3))
    female5$id<-paste(female5$X1,female5$X2,female5$X3,sep="_")
    female5<-female5[,-c(1:3)]
    female1<-merge(female1,female5,by="id")
    if (Plot == "B30"){
      female1<-female1[,-c(2:24)]#may not work for B30... have to do another loop with in here to make sure it works
    }
    if (Plot == "Track"){
      female1<-female1[,-c(2:35)]#may not work for B30... have to do another loop with in here to make sure it works
    }
    
    test<-rbind(female1,male6)
    test<-distinct(test)
    
    BrapaP<-merge(BrapaP,test,by="id",all.x=T)#creates more individuals.... because same dad can have possibly 3 different letters.. this may affect something down the line... should make this a separate part of the code.
    BrapaP<-distinct(BrapaP)
    
    #need to create fake clump for offspring: "Z"
    BrapaP$clump_pos[BrapaP$offspring == 1]<-"Z"
    
  }#end of loop for applying small, med, and large restrictions for even plots
  
  if (Plot == "Comp" | Plot == "Fwest"){
    male1<-BrapaP[grep("Male",BrapaP$sex),]
    male2<-data.frame(do.call('rbind', strsplit(as.character(male1$id),'_',fixed=TRUE)))
    colnames(male2)[2]<-"pat_pos"
    male3<-data.frame(do.call('rbind', strsplit(as.character(male1$id),'_',fixed=TRUE)))
    male3$X1<-NULL
    male3<-cbind(male3,male2)
    male3$id<-paste(male3$X1,male3$X2,sep="_")
    male3<-male3[,-c(1:2)]
    male3$pat_pos<-as.numeric(as.character(male3$pat_pos))#why do you have to do it like this... when just as.numeric, the numbers change to their true integers. Why does adding as.character change that?
    male3$clump_pos[male3$pat_pos >= 1 & male3$pat_pos <= 36]<-"A"
    male3$clump_pos[male3$pat_pos >= 37 & male3$pat_pos <= 73]<-"B"
    male3$clump_pos[male3$pat_pos >= 74 & male3$pat_pos <= 110]<-"C"
    male3$clump_pos[male3$pat_pos >= 111 & male3$pat_pos <= 147]<-"D"
    male3$clump_pos[male3$pat_pos >= 148 & male3$pat_pos <= 184]<-"E"
    male3$clump_pos[male3$pat_pos >= 185 & male3$pat_pos <= 221]<-"F"
    male3$pat_pos<-NULL
    
    female1<-BrapaP[grep("Female",BrapaP$sex),]
    female2<-data.frame(do.call('rbind', strsplit(as.character(female1$id),'_',fixed=TRUE)))
    female2$clump_pos<-female2$X3
    female2$id<-paste(female2$X1,female2$X2,female2$X3,sep="_")
    female2<-female2[,-c(1:3)]
    
    clump_pos<-rbind(male3,female2)
    BrapaP<-merge(BrapaP,clump_pos,by="id",all.x=T)
    BrapaP<-distinct(BrapaP)
    
    
    #need to create fake clump for offspring: "Z"
    BrapaP$clump_pos[BrapaP$offspring == 1]<-"Z"
  }
  
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
    if (v == 42){var42<-expression(varPed(x="z.dur", gender="Male"))} #comp and b30 only
    if (v == 43){var43<-expression(varPed(x="z.dur_inpop", gender="Male"))} #fwest and track only
    if (v == 44){var44<-expression(varPed(x="z.dur_total", gender="Male"))}#fwest and track only 
    if (v == 45){var45<-expression(varPed(x="z.sd", gender="Male"))} #all
    if (v == 46){var46<-expression(varPed(x="z.sd_field", gender="Male"))}#fwest and track only 
    if (v == 47){var.dist<-expression(varPed(x=c("X", "Y"), gender="Male", relational="MATE"))}
    
    if(v == 48){var.small<-expression(varPed(x="sm_clump",gender="Male",relational="MATE"))}
    if(v == 49){var.medium<-expression(varPed(x="med_clump",gender="Male",relational="MATE"))}
    if(v == 50){var.large<-expression(varPed(x="lrg_clump",gender="Male",relational="MATE"))}
    if(v == 51){var.clump<-expression(varPed(x="clump_pos",gender="Male",relational="MATE"))}
    
  } #end for loop
  Vrs.list<-lapply(ls(pattern="var*"),get)
  
  BrapaG_alleles<-BrapaG_for_MB
  BrapaG_alleles$timevar<-NULL
  alleles<-extractA(BrapaG_alleles)#estimate true allele frequencies based on largest available sample
  #only need gender in phenotype not genotype
  BrapaG_for_MB$timevar<-NULL
  
  
  ####GdP####
  GdP<-GdataPed(G=BrapaG_for_MB, categories=NULL, perlocus=T)	#need to specify 'per locus' to allow variable error rates by locus.
  
  #mod_subfolder<-"data_files_3/files_generated_by_PatAnalysis_Script"
  #MODEL: Distance, start date, total flowers, duration of flowering, skew, etc. by batch
 # mod_subfolder<-"data_files/Analysis_Dec_19_LL"	#This empty folder should already be present.
  
  
  
  ####PdP####
  #PdP<-if(sum(Vrs) == 105){PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var1,var2, var3,var4,var5,var6,var7,var8,var9,var10,var11,var12,var13,var14), data=BrapaP_for_MB)} else if (
 #  sum(Vrs) == 715) {PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var15,var16, var17,var18,var19,var20,var21,var22,var23,var24,var25,var26,var27,var28,var29,var30,var31,var32,var33,var34,var35,var36,var37,var38,var39,var40), data=BrapaP_for_MB)} else if (
 #    sum(Vrs) == 126){PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var41,var42,var.dist), data=BrapaP_for_MB)} else if (
   #     sum(Vrs) == 231) {PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var1,var2, var3,var4,var5,var6,var7,var8,var9,var10,var11,var12,var13,var14,var41,var42,var.dist), data=BrapaP_for_MB)} else if (
  #        sum(Vrs) == 841) {PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var15,var16, var17,var18,var19,var20,var21,var22,var23,var24,var25,var26,var27,var28,var29,var30,var31,var32,var33,var34,var35,var36,var37,var38,var39,var40,var41,var42,var.dist), data=BrapaP_for_MB)}
  
#PdP<-PdataPed(formula=list(res.osnotpar, res.osterr, var.dist), data=BrapaP_for_MB)
  PdP<-PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var.dist, var41, var43,var44,var45,var46), data=BrapaP_for_MB)
  # when just one alone, it is missing covariate data
  

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
  
#  test<-as.matrix(PdP)
 # df <- data.frame(matrix(unlist(PdP), nrow=1585, byrow=T))
  #df2 <- data.frame(matrix(unlist(PdP), nrow=1585, byrow=T),stringsAsFactors=FALSE)
  
  ch1<-MCMCped(PdP=PdP, GdP=GdP, sP=sP.chain1, tP=tP,verbose=T,mm.tol=5, nitt=10000+(10000*10), thin=10*10, burnin=10000)

  ####saving output####
  Beta.out<-ch1$beta
  B<-length(Beta.out)
  Beta.df<-as.data.frame(matrix(nrow=B,ncol=1))
  for (b in 1:B){
    Beta.df[b,1]<-Beta.out[b]
  }
  write.csv(Beta.df, "Beta.df.fwest.phenotraits.csv", row.names=FALSE)
  
  summary(ch1$beta)
  
####by batch####
if (extra == "Batch"){

  if (Plot == "B30" | Plot == "Comp"){
    BrapaG_bA<-BrapaG[grep("bA",BrapaG$timevar),]
    BrapaG_bB<-BrapaG[-grep("bB2",BrapaG$timevar),]
    BrapaG_bB<-BrapaG[grep("bB",BrapaG_bB$timevar),]
    BrapaG_bB2<-BrapaG[grep("bB2",BrapaG$timevar),]
    BrapaG_bC<-BrapaG[grep("bC",BrapaG$timevar),]
    BrapaG_bD<-BrapaG[-grep("bD2",BrapaG$timevar),]
    BrapaG_bD<-BrapaG[grep("bD",BrapaG_bD$timevar),]
    BrapaG_bD2<-BrapaG[grep("bD2",BrapaG$timevar),]
    BrapaG_bE<-BrapaG[grep("bE",BrapaG$timevar),]    
    
    BrapaP_bA<-BrapaP[grep("bA",BrapaP$timevar),]
    BrapaP_bB<-BrapaP[-grep("bB2",BrapaP$timevar),]
    BrapaP_bB<-BrapaP[grep("bB",BrapaP_bB$timevar),]
    BrapaP_bB2<-BrapaP[grep("bB2",BrapaP$timevar),]
    BrapaP_bC<-BrapaP[grep("bC",BrapaP$timevar),]
    BrapaP_bD<-BrapaP[-grep("bD2",BrapaP$timevar),]
    BrapaP_bD<-BrapaP[grep("bD",BrapaP_bD$timevar),]
    BrapaP_bD2<-BrapaP[grep("bD2",BrapaP$timevar),]
    BrapaP_bE<-BrapaP[grep("bE",BrapaP$timevar),] 
    
    BrapaG_for_MB_bA<-BrapaG_bA
    BrapaG_for_MB_bB<-BrapaG_bB
    BrapaG_for_MB_bB2<-BrapaG_bB2
    BrapaG_for_MB_bC<-BrapaG_bC
    BrapaG_for_MB_bD<-BrapaG_bD
    BrapaG_for_MB_bD2<-BrapaG_bD2
    BrapaG_for_MB_bE<-BrapaG_bE
    
    BrapaP_for_MB_bA<-BrapaP_bA
    BrapaP_for_MB_bB<-BrapaP_bB
    BrapaP_for_MB_bB2<-BrapaP_bB2
    BrapaP_for_MB_bC<-BrapaP_bC
    BrapaP_for_MB_bD<-BrapaP_bD
    BrapaP_for_MB_bD2<-BrapaP_bD2
    BrapaP_for_MB_bE<-BrapaP_bE
    
    BrapaG_alleles_bA<-BrapaG_for_MB_bA
    BrapaG_alleles_bB<-BrapaG_for_MB_bB
    BrapaG_alleles_bB2<-BrapaG_for_MB_bB2
    BrapaG_alleles_bC<-BrapaG_for_MB_bC
    BrapaG_alleles_bD<-BrapaG_for_MB_bD
    BrapaG_alleles_bD2<-BrapaG_for_MB_bD2
    BrapaG_alleles_bE<-BrapaG_for_MB_bE
    
    BrapaG_alleles_bA$timevar<-NULL
    BrapaG_alleles_bB$timevar<-NULL
    BrapaG_alleles_bB2$timevar<-NULL
    BrapaG_alleles_bC$timevar<-NULL
    BrapaG_alleles_bD$timevar<-NULL
    BrapaG_alleles_bD2$timevar<-NULL
    BrapaG_alleles_bE$timevar<-NULL
    
    
    
    alleles_bA<-extractA(BrapaG_alleles_bA)
    alleles_bB<-extractA(BrapaG_alleles_bB)
    alleles_bB2<-extractA(BrapaG_alleles_bB2)
    alleles_bC<-extractA(BrapaG_alleles_bC)
    alleles_bD<-extractA(BrapaG_alleles_bD)
    alleles_bD2<-extractA(BrapaG_alleles_bD2)
    alleles_bE<-extractA(BrapaG_alleles_bE)
    #estimate true allele frequencies based on largest available sample
    #only need gender in phenotype not genotype
    
    BrapaG_for_MB_bA$timevar<-NULL
    BrapaG_for_MB_bB$timevar<-NULL
    BrapaG_for_MB_bB2$timevar<-NULL
    BrapaG_for_MB_bC$timevar<-NULL
    BrapaG_for_MB_bD$timevar<-NULL
    BrapaG_for_MB_bD2$timevar<-NULL
    BrapaG_for_MB_bE$timevar<-NULL
    
    GdP_bA<-GdataPed(G=BrapaG_for_MB_bA, categories=NULL, perlocus=T)
    GdP_bB<-GdataPed(G=BrapaG_for_MB_bB, categories=NULL, perlocus=T)
    GdP_bB2<-GdataPed(G=BrapaG_for_MB_bB2, categories=NULL, perlocus=T)
    GdP_bC<-GdataPed(G=BrapaG_for_MB_bC, categories=NULL, perlocus=T)
    GdP_bD<-GdataPed(G=BrapaG_for_MB_bD, categories=NULL, perlocus=T)
    GdP_bD2<-GdataPed(G=BrapaG_for_MB_bD2, categories=NULL, perlocus=T)
    GdP_bE<-GdataPed(G=BrapaG_for_MB_bE, categories=NULL, perlocus=T)
    
    PdP_bA<-PdataPed(formula=list(res.osnotpar, res.osterr,res.phen, var.dist), data=BrapaP_for_MB_bA)
    PdP_bB<-PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var.dist), data=BrapaP_for_MB_bB)
    PdP_bB2<-PdataPed(formula=list(res.osnotpar, res.osterr, res.phen,var.dist), data=BrapaP_for_MB_bB2)
    PdP_bC<-PdataPed(formula=list(res.osnotpar, res.osterr, res.phen,var.dist), data=BrapaP_for_MB_bC)
    PdP_bD<-PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var.dist), data=BrapaP_for_MB_bD)
    PdP_bD2<-PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var.dist), data=BrapaP_for_MB_bD2)
    PdP_bE<-PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var.dist), data=BrapaP_for_MB_bE)
    
    sP.chain1<-startPed(estG=F, E1=epsilon1, E2=epsilon2, A=alleles, estE1=F, estE2=F, estA=F)
    sP.chain2<-startPed(estG=F, E1=epsilon1, E2=epsilon2, A=alleles, estE1=F, estE2=F, estA=F, beta=rep(0, times=1))
    sP.chain3<-startPed(estG=F, E1=epsilon1, E2=epsilon2, A=alleles, estE1=F, estE2=F, estA=F, beta=rnorm(n=1, mean=c(-0.1, 0.1), sd=0.05))
    
    tP<-tunePed(beta=0.325)
    pP<-priorPed(USsire=list(mu=log(5), sigma=0.5))

        ch1_bA<-MCMCped(PdP=PdP_bA, GdP=GdP_bA, sP=sP.chain1, tP=tP,verbose=T,mm.tol=5, nitt=10000+(10000*10), thin=10*10, burnin=10000)
        ch1_bB<-MCMCped(PdP=PdP_bB, GdP=GdP_bB, sP=sP.chain1, tP=tP,verbose=T,mm.tol=5, nitt=10000+(10000*10), thin=10*10, burnin=10000)
        ch1_bB2<-MCMCped(PdP=PdP_bB2, GdP=GdP_bB2, sP=sP.chain1, tP=tP,verbose=T,mm.tol=5, nitt=10000+(10000*10), thin=10*10, burnin=10000)
        ch1_bC<-MCMCped(PdP=PdP_bC, GdP=GdP_bC, sP=sP.chain1, tP=tP,verbose=T,mm.tol=5, nitt=10000+(10000*10), thin=10*10, burnin=10000)
        ch1_bD<-MCMCped(PdP=PdP_bD, GdP=GdP_bD, sP=sP.chain1, tP=tP,verbose=T,mm.tol=5, nitt=10000+(10000*10), thin=10*10, burnin=10000)
        ch1_bD2<-MCMCped(PdP=PdP_bD2, GdP=GdP_bD2, sP=sP.chain1, tP=tP,verbose=T,mm.tol=5, nitt=10000+(10000*10), thin=10*10, burnin=10000)
        ch1_bE<-MCMCped(PdP=PdP_bE, GdP=GdP_bE, sP=sP.chain1, tP=tP,verbose=T,mm.tol=5, nitt=10000+(10000*10), thin=10*10, burnin=10000) 
        
        summary(ch1_bA$beta)
        summary(ch1_bB$beta)
        summary(ch1_bB2$beta)
        summary(ch1_bC$beta)
        summary(ch1_bD$beta)
        summary(ch1_bD2$beta)
        summary(ch1_bE$beta)
        
  }#end loop for B30 and Comp by batch 
  
  if (Plot == "Track" | Plot == "Fwest"){
    BrapaG_bG<-BrapaG[grep("bG",BrapaG$timevar),]
    BrapaG_bH<-BrapaG[grep("bH",BrapaG$timevar),]
    BrapaG_bI<-BrapaG[grep("bI",BrapaG$timevar),]
    BrapaG_bJ<-BrapaG[grep("bJ",BrapaG$timevar),]
    BrapaG_bK<-BrapaG[grep("bK",BrapaG$timevar),]
    BrapaG_bL<-BrapaG[grep("bL",BrapaG$timevar),]
    BrapaG_bM<-BrapaP[-grep("bM2",BrapaG$timevar),]
    BrapaG_bM<-BrapaG[grep("bM",BrapaG_bM$timevar),]
    BrapaG_bM2<-BrapaG[grep("bM2",BrapaG$timevar),]
    BrapaG_bN<-BrapaP[-grep("bN2",BrapaG$timevar),]
    BrapaG_bN<-BrapaG[grep("bN",BrapaG_bN$timevar),]    
    BrapaG_bN2<-BrapaG[grep("bN2",BrapaG$timevar),]
    BrapaG_bO<-BrapaP[-grep("bO2",BrapaG$timevar),]
    BrapaG_bO<-BrapaG[grep("bO",BrapaG_bO$timevar),]    
    BrapaG_bO2<-BrapaG[grep("bO2",BrapaG$timevar),]
    BrapaG_bP<-BrapaG[grep("bP",BrapaG$timevar),]  
    BrapaG_bQ<-BrapaG[grep("bQ",BrapaG$timevar),]   
    
    BrapaP_bG<-BrapaP[grep("bG",BrapaP$timevar),]
    BrapaP_bH<-BrapaP[grep("bH",BrapaP$timevar),]
    BrapaP_bI<-BrapaP[grep("bI",BrapaP$timevar),]
    BrapaP_bJ<-BrapaP[grep("bJ",BrapaP$timevar),]
    BrapaP_bK<-BrapaP[grep("bK",BrapaP$timevar),]
    BrapaP_bL<-BrapaP[grep("bL",BrapaP$timevar),]
    BrapaP_bM<-BrapaP[-grep("bM2",BrapaP$timevar),]
    BrapaP_bM<-BrapaP[grep("bM",BrapaP_bM$timevar),]
    BrapaP_bM2<-BrapaP[grep("bM2",BrapaP$timevar),]
    BrapaP_bN<-BrapaP[-grep("bN2",BrapaP$timevar),]
    BrapaP_bN<-BrapaP[grep("bN",BrapaP_bN$timevar),]    
    BrapaP_bN2<-BrapaP[grep("bN2",BrapaP$timevar),]
    BrapaP_bO<-BrapaP[-grep("bO2",BrapaP$timevar),]
    BrapaP_bO<-BrapaP[grep("bO",BrapaP_bO$timevar),]    
    BrapaP_bO2<-BrapaP[grep("bO2",BrapaP$timevar),]
    BrapaP_bP<-BrapaP[grep("bP",BrapaP$timevar),]  
    BrapaP_bQ<-BrapaP[grep("bQ",BrapaP$timevar),]   
    
    
    BrapaG_for_MB_bG<-BrapaG_bG
    BrapaG_for_MB_bH<-BrapaG_bH
    BrapaG_for_MB_bI<-BrapaG_bI
    BrapaG_for_MB_bJ<-BrapaG_bJ
    BrapaG_for_MB_bK<-BrapaG_bK
    BrapaG_for_MB_bL<-BrapaG_bL
    BrapaG_for_MB_bM<-BrapaG_bM
    BrapaG_for_MB_bM<-BrapaG_bM
    BrapaG_for_MB_bM2<-BrapaG_bM2
    BrapaG_for_MB_bN<-BrapaG_bN
    BrapaG_for_MB_bN2<-BrapaG_bN2
    BrapaG_for_MB_bO<-BrapaG_bO
    BrapaG_for_MB_bO2<-BrapaG_bO2
    BrapaG_for_MB_bP<-BrapaG_bP
    BrapaG_for_MB_bQ<-BrapaG_bQ  
    
    BrapaP_for_MB_bG<-BrapaP_bG
    BrapaP_for_MB_bH<-BrapaP_bH
    BrapaP_for_MB_bI<-BrapaP_bI
    BrapaP_for_MB_bJ<-BrapaP_bJ
    BrapaP_for_MB_bK<-BrapaP_bK
    BrapaP_for_MB_bL<-BrapaP_bL
    BrapaP_for_MB_bM<-BrapaP_bM
    BrapaP_for_MB_bM<-BrapaP_bM
    BrapaP_for_MB_bM2<-BrapaP_bM2
    BrapaP_for_MB_bN<-BrapaP_bN
    BrapaP_for_MB_bN2<-BrapaP_bN2
    BrapaP_for_MB_bO<-BrapaP_bO
    BrapaP_for_MB_bO2<-BrapaP_bO2
    BrapaP_for_MB_bP<-BrapaP_bP
    BrapaP_for_MB_bQ<-BrapaP_bQ  
    
    
    BrapaG_alleles_bG<-BrapaG_for_MB_bG
    BrapaG_alleles_bH<-BrapaG_for_MB_bH
    BrapaG_alleles_bI<-BrapaG_for_MB_bI
    BrapaG_alleles_bJ<-BrapaG_for_MB_bJ
    BrapaG_alleles_bK<-BrapaG_for_MB_bK
    BrapaG_alleles_bL<-BrapaG_for_MB_bL
    BrapaG_alleles_bM<-BrapaG_for_MB_bM
    BrapaG_alleles_bM<-BrapaG_for_MB_bM
    BrapaG_alleles_bM2<-BrapaG_for_MB_bM2
    BrapaG_alleles_bN<-BrapaG_for_MB_bN
    BrapaG_alleles_bN2<-BrapaG_for_MB_bN2
    BrapaG_alleles_bO<-BrapaG_for_MB_bO
    BrapaG_alleles_bO2<-BrapaG_for_MB_bO2
    BrapaG_alleles_bP<-BrapaG_for_MB_bP
    BrapaG_alleles_bQ<-BrapaG_for_MB_bQ 
    
    BrapaG_alleles_bG$timevar<-NULL
    BrapaG_alleles_bH$timevar<-NULL
    BrapaG_alleles_bI$timevar<-NULL
    BrapaG_alleles_bJ$timevar<-NULL
    BrapaG_alleles_bK$timevar<-NULL
    BrapaG_alleles_bL$timevar<-NULL
    BrapaG_alleles_bM$timevar<-NULL
    BrapaG_alleles_bM$timevar<-NULL
    BrapaG_alleles_bM2$timevar<-NULL
    BrapaG_alleles_bN$timevar<-NULL
    BrapaG_alleles_bN2$timevar<-NULL
    BrapaG_alleles_bO$timevar<-NULL
    BrapaG_alleles_bO2$timevar<-NULL
    BrapaG_alleles_bP$timevar<-NULL
    BrapaG_alleles_bQ$timevar<-NULL

    
    alleles_bG<-extractA(BrapaG_alleles_bG)
    alleles_bH<-extractA(BrapaG_alleles_bG)
    alleles_bI<-extractA(BrapaG_alleles_bI)
    alleles_bJ<-extractA(BrapaG_alleles_bJ)
    alleles_bK<-extractA(BrapaG_alleles_bK)
    alleles_bL<-extractA(BrapaG_alleles_bL)
    alleles_bM<-extractA(BrapaG_alleles_bM)
    alleles_bM2<-extractA(BrapaG_alleles_bM2)
    alleles_bN<-extractA(BrapaG_alleles_bN)
    alleles_bN2<-extractA(BrapaG_alleles_bN2)
    alleles_bO<-extractA(BrapaG_alleles_bO)
    alleles_bO2<-extractA(BrapaG_alleles_bO2)
    alleles_bP<-extractA(BrapaG_alleles_bP)
    alleles_bQ<-extractA(BrapaG_alleles_bQ)
    #estimate true allele frequencies based on largest available sample
    #only need gender in phenotype not genotype
    
    BrapaG_for_MB_bG$timevar<-NULL
    BrapaG_for_MB_bH$timevar<-NULL
    BrapaG_for_MB_bI$timevar<-NULL
    BrapaG_for_MB_bJ$timevar<-NULL
    BrapaG_for_MB_bK$timevar<-NULL
    BrapaG_for_MB_bL$timevar<-NULL
    BrapaG_for_MB_bM$timevar<-NULL
    BrapaG_for_MB_bM$timevar<-NULL
    BrapaG_for_MB_bM2$timevar<-NULL
    BrapaG_for_MB_bN$timevar<-NULL
    BrapaG_for_MB_bN2$timevar<-NULL
    BrapaG_for_MB_bO$timevar<-NULL
    BrapaG_for_MB_bO2$timevar<-NULL
    BrapaG_for_MB_bP$timevar<-NULL
    BrapaG_for_MB_bQ$timevar<-NULL 
    
    GdP_bG<-GdataPed(G=BrapaG_for_MB_bG, categories=NULL, perlocus=T)
    GdP_bH<-GdataPed(G=BrapaG_for_MB_bH, categories=NULL, perlocus=T)
    GdP_bI<-GdataPed(G=BrapaG_for_MB_bI, categories=NULL, perlocus=T)
    GdP_bJ<-GdataPed(G=BrapaG_for_MB_bJ, categories=NULL, perlocus=T)
    GdP_bK<-GdataPed(G=BrapaG_for_MB_bK, categories=NULL, perlocus=T)
    GdP_bL<-GdataPed(G=BrapaG_for_MB_bL, categories=NULL, perlocus=T)
    GdP_bM<-GdataPed(G=BrapaG_for_MB_bM, categories=NULL, perlocus=T)
    GdP_bM2<-GdataPed(G=BrapaG_for_MB_bM2, categories=NULL, perlocus=T)
    GdP_bN<-GdataPed(G=BrapaG_for_MB_bN, categories=NULL, perlocus=T)
    GdP_bN2<-GdataPed(G=BrapaG_for_MB_bN2, categories=NULL, perlocus=T)
    GdP_bO<-GdataPed(G=BrapaG_for_MB_bO, categories=NULL, perlocus=T)
    GdP_bO2<-GdataPed(G=BrapaG_for_MB_bO2, categories=NULL, perlocus=T)
    GdP_bP<-GdataPed(G=BrapaG_for_MB_bP, categories=NULL, perlocus=T)
    GdP_bQ<-GdataPed(G=BrapaG_for_MB_bQ, categories=NULL, perlocus=T)
    
    PdP_bG<-PdataPed(formula=list(res.osnotpar, res.osterr,res.phen, var.dist), data=BrapaP_for_MB_bG)
    PdP_bH<-PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var.dist), data=BrapaP_for_MB_bH)
    PdP_bI<-PdataPed(formula=list(res.osnotpar, res.osterr, res.phen,var.dist), data=BrapaP_for_MB_bI)
    PdP_bJ<-PdataPed(formula=list(res.osnotpar, res.osterr, res.phen,var.dist), data=BrapaP_for_MB_bJ)
    PdP_bK<-PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var.dist), data=BrapaP_for_MB_bK)
    PdP_bL<-PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var.dist), data=BrapaP_for_MB_bL)
    PdP_bM<-PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var.dist), data=BrapaP_for_MB_bM)
    PdP_bM2<-PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var.dist), data=BrapaP_for_MB_bM2)
    PdP_bN<-PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var.dist), data=BrapaP_for_MB_bN)
    PdP_bN2<-PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var.dist), data=BrapaP_for_MB_bN2)
    PdP_bO<-PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var.dist), data=BrapaP_for_MB_bO)
    PdP_bO2<-PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var.dist), data=BrapaP_for_MB_bO2)
    PdP_bP<-PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var.dist), data=BrapaP_for_MB_bP)
    PdP_bQ<-PdataPed(formula=list(res.osnotpar, res.osterr, res.phen, var.dist), data=BrapaP_for_MB_bQ)
    
    
    sP.chain1<-startPed(estG=F, E1=epsilon1, E2=epsilon2, A=alleles, estE1=F, estE2=F, estA=F)
    sP.chain2<-startPed(estG=F, E1=epsilon1, E2=epsilon2, A=alleles, estE1=F, estE2=F, estA=F, beta=rep(0, times=1))
    sP.chain3<-startPed(estG=F, E1=epsilon1, E2=epsilon2, A=alleles, estE1=F, estE2=F, estA=F, beta=rnorm(n=1, mean=c(-0.1, 0.1), sd=0.05))
    
    tP<-tunePed(beta=0.325)
    pP<-priorPed(USsire=list(mu=log(5), sigma=0.5))
    
    ch1_bG<-MCMCped(PdP=PdP_bG, GdP=GdP_bG, sP=sP.chain1, tP=tP,verbose=T,mm.tol=5, nitt=10000+(10000*10), thin=10*10, burnin=10000)
    ch1_bH<-MCMCped(PdP=PdP_bH, GdP=GdP_bH, sP=sP.chain1, tP=tP,verbose=T,mm.tol=5, nitt=10000+(10000*10), thin=10*10, burnin=10000)
    ch1_bI<-MCMCped(PdP=PdP_bI, GdP=GdP_bI, sP=sP.chain1, tP=tP,verbose=T,mm.tol=5, nitt=10000+(10000*10), thin=10*10, burnin=10000)
    ch1_bJ<-MCMCped(PdP=PdP_bJ, GdP=GdP_bJ, sP=sP.chain1, tP=tP,verbose=T,mm.tol=5, nitt=10000+(10000*10), thin=10*10, burnin=10000)
    ch1_bK<-MCMCped(PdP=PdP_bK, GdP=GdP_bK, sP=sP.chain1, tP=tP,verbose=T,mm.tol=5, nitt=10000+(10000*10), thin=10*10, burnin=10000)
    ch1_bL<-MCMCped(PdP=PdP_bL, GdP=GdP_bL, sP=sP.chain1, tP=tP,verbose=T,mm.tol=5, nitt=10000+(10000*10), thin=10*10, burnin=10000)
    ch1_bM<-MCMCped(PdP=PdP_bM, GdP=GdP_bM, sP=sP.chain1, tP=tP,verbose=T,mm.tol=5, nitt=10000+(10000*10), thin=10*10, burnin=10000) 
    ch1_bM2<-MCMCped(PdP=PdP_bM2, GdP=GdP_bM2, sP=sP.chain1, tP=tP,verbose=T,mm.tol=5, nitt=10000+(10000*10), thin=10*10, burnin=10000)
    ch1_bN<-MCMCped(PdP=PdP_bN, GdP=GdP_bN, sP=sP.chain1, tP=tP,verbose=T,mm.tol=5, nitt=10000+(10000*10), thin=10*10, burnin=10000) 
    ch1_bN2<-MCMCped(PdP=PdP_bN2, GdP=GdP_bN2, sP=sP.chain1, tP=tP,verbose=T,mm.tol=5, nitt=10000+(10000*10), thin=10*10, burnin=10000)
    ch1_bO<-MCMCped(PdP=PdP_bO, GdP=GdP_bO, sP=sP.chain1, tP=tP,verbose=T,mm.tol=5, nitt=10000+(10000*10), thin=10*10, burnin=10000) 
    ch1_bO2<-MCMCped(PdP=PdP_bO2, GdP=GdP_bO2, sP=sP.chain1, tP=tP,verbose=T,mm.tol=5, nitt=10000+(10000*10), thin=10*10, burnin=10000)
    ch1_bP<-MCMCped(PdP=PdP_bP, GdP=GdP_bP, sP=sP.chain1, tP=tP,verbose=T,mm.tol=5, nitt=10000+(10000*10), thin=10*10, burnin=10000) 
    ch1_bQ<-MCMCped(PdP=PdP_bQ, GdP=GdP_bQ, sP=sP.chain1, tP=tP,verbose=T,mm.tol=5, nitt=10000+(10000*10), thin=10*10, burnin=10000)
    
    summary(ch1_bG$beta)
    summary(ch1_bH$beta)
    summary(ch1_bI$beta)
    summary(ch1_bJ$beta)
    summary(ch1_bK$beta)
    summary(ch1_bL$beta)
    summary(ch1_bM$beta)
    summary(ch1_bM2$beta)
    summary(ch1_bN$beta)
    summary(ch1_bN2$beta)
    summary(ch1_bO$beta)
    summary(ch1_bO2$beta)
    summary(ch1_bP$beta)
    summary(ch1_bQ$beta)
    #multiple records for individuals but time var does not vary... happens for B30 and Track... must be because of small med and large clumps... because can have same individ with a different thing in each column 
}#end lop for Fwest and Track by Batch
  

    #G<-c(1:7)
    #for (g in alone_batch_G){
     # if (g == 1){BrapaG_bA<-BrapaG[grep("bA",BrapaG$timevar)]}
      #if (g == 2){BrapaG_bB<-BrapaG[-grep("bB2",BrapaG$timevar)] %>%
      #BrapaG[grep("bB",BrapaG$timevar)]}
      #if (g == 3){BrapaG_bB2<-BrapaG[grep("bB2",BrapaG$timevar)]}
      #if (g == 4){BrapaG_bC<-BrapaG[grep("bC",BrapaG$timevar)]}
      #if (g == 5){BrapaG_bD<-BrapaG[-grep("bD2",BrapaG$timevar)] %>%
      # BrapaG[grep("bD",BrapaG$timevar)]}
      #if (g == 6){BrapaG_bD2<-BrapaG[grep("bD2",BrapaG$timevar)]}
      #if (g == 7){BrapaG_bE<-BrapaG[grep("bE",BrapaG$timevar)]}
    #}
    
    #P<-c(1:7)
    #for (p in alone_batch_P){
     # if (p == 1){BrapaP_bA<-BrapaP[grep("bA",BrapaP$timevar)]}
      #if (p == 2){BrapaP_bB<-BrapaP[-grep("bB2",BrapaP$timevar)] %>%
       # BrapaP[grep("bB",BrapaP$timevar)]}
      #if (p == 3){BrapaP_bB2<-BrapaP[grep("bB2",BrapaP$timevar)]}
      ##if (p == 4){BrapaP_bC<-BrapaP[grep("bC",BrapaP$timevar)]}
      #if (p == 5){BrapaP_bD<-BrapaP[-grep("bD2",BrapaP$timevar)] %>%
       # BrapaP[grep("bD",BrapaP$timevar)]}
    #  if (p == 6){BrapaP_bD2<-BrapaP[grep("bD2",BrapaP$timevar)]}
     # if (p == 7){BrapaP_bE<-BrapaP[grep("bE",BrapaP$timevar)]}
    #}
    
   # for (G in 1:G){
   #   BrapaG_for_MB<-BrapaG
   #   Beta.df[b,1]<-Beta.out[b]
  #  }
    
  

 
 }
else if (group == "Batch") {
  #script for doing analysis just by batch 
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
  } 
 }
  
  