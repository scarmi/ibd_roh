setwd("/lustre/scratch115/projects/autozyg/BiB/code_for_popgen_paper/")
outdir="/lustre/scratch115/projects/autozyg/BiB/code_for_popgen_paper/plots/"
library(data.table)
library(scales)

min.cm=5
max.cm=30
### load a table of consanguinity rates estimated from self-reported data in the Bradford Pakistani mothers, split by homogeneous sub-group defined in FineStructure --> we just use this to determine kinship values used and then plot these
load("consanguinity_rates_estimated_naively_from_self_reported_information.homogenous_groups_from_fineSTRUCTURE.RData")

### load the empirical ROH and IBD length distributions/footprints, generated with various different algorithms with or without filtering
load(paste("empirical_ROH_and_IBD_distributions_to_plot/IBD_length_distribution_between_",min.cm,"_and_",max.cm,"cM.IBDseq.with_full_filtering.homogenous_groups_from_fineSTRUCTURE.RData",sep=""))
load(paste("empirical_ROH_and_IBD_distributions_to_plot/IBD_footprint_between_",min.cm,"_and_",max.cm,"cM.IBDseq.with_full_filtering.homogenous_groups_from_fineSTRUCTURE.RData",sep=""))
load(paste("empirical_ROH_and_IBD_distributions_to_plot/ROH_length_distribution_between_",min.cm,"_and_",max.cm,"cM..homogenous_groups_from_fineSTRUCTURE.RData",sep=""))
load(paste("empirical_ROH_and_IBD_distributions_to_plot/ROH_footprint_between_",min.cm,"_and_",max.cm,"cM.homogenous_groups_from_fineSTRUCTURE.RData",sep="") )
load(paste("empirical_ROH_and_IBD_distributions_to_plot/IBD_length_distribution_between_",min.cm,"_and_",max.cm,"cM.IBDseq.with_NO_filtering.homogenous_groups_from_fineSTRUCTURE.RData",sep=""))#prop.ibdseq.segs.per.bin.by.group2.raw
load(paste("empirical_ROH_and_IBD_distributions_to_plot/IBD_footprint_between_",min.cm,"_and_",max.cm,"cM.IBDseq.with_NO_filtering.homogenous_groups_from_fineSTRUCTURE.RData",sep=""))#footprint.ibdseq.segs.per.bin.by.group2.raw
load(paste("empirical_ROH_and_IBD_distributions_to_plot/ROH_length_distribution_between_",min.cm,"_and_",max.cm,"cM.different_ROH_callers.homogenous_groups_from_fineSTRUCTURE.RData",sep=""))#compare.prop.roh.segs.per.bin.by.group
load(paste("empirical_ROH_and_IBD_distributions_to_plot/ROH_footprint_between_",min.cm,"_and_",max.cm,"cM.different_ROH_callers.homogenous_groups_from_fineSTRUCTURE.RData",sep=""))#compare.footprint.roh.segs.per.bin.by.group


######varying max_T_exact_model
max_T_exact_model.options=c(10,20,30,40,50)

### scaling factors used for the rates vector
scaling.factors=c(1,1.1,1.5,1.55,2,2.18,2.2,2.5,3,3.1)

### function to get the filename of the output from run_code_to_generate_expected_IBD_and_ROH_distribution.flexible_options.R
get.label=function(pop.label,ne.label,max_T_exact_model,old.kinship,which.recent.consang,kinship){
  if(which.recent.consang==1){
    label=paste(pop.label,".",ne.label,".empirical_consang_rates_for_",max_T_exact_model,"_gens",sep="")
  } else {
    label=paste(pop.label,".",ne.label,".mean_kinship_",kinship,"_for_",max_T_exact_model,"_gens",sep="")
  }
  label=paste(label,".then_mean_kinship_",old.kinship,sep="")
  return(label)
}

### specify some plotting parameters
my.lwd=3
roh.lty=5
ibd.lty=2

#### effect of changing population size on ROH footprint - Supplementary Figure 17
all.my.ylim=c()
pdf(paste(outdir,"Figure_S17.effect_of_Ne_and_kinship_on_ROH_footprint.max_T_exact_model_50.Pathan_versus_JattChoudry.pdf",sep=""),height=8,width=8)
my.ltys=c(1,3)
names(my.ltys)=c("PATHAN","Jatt/Choudhary")
for(max_T_exact_model in c(50)){
  store.kinship=c()
  for(pop in c("PATHAN","Jatt/Choudhary")){
    if(pop=="PATHAN"){
      pop.label="PATHAN_in_Cluster_8"
    } else {
      pop.label="JattChoudhry_in_cluster_10"
    }
    which.recent.consang=scaling.factors[c(1)]
    rates = c(0,consang.by.bir2.robust[pop,] * which.recent.consang)
    n = length(rates)
    kinship = sum(rates * 4^(-(1:n)))
    store.kinship=c(store.kinship,kinship)
### loop over old kinship
    for(old.kinship in c(-1)){
      if(old.kinship<0){
        old.kinship=kinship
      }
      constant.ne.values=c(1e3,10e3,100e3)
      ne.labels = c(paste("constant_Ne_",constant.ne.values/1e3,"k",sep=""))
      cols=c("#b2abd2","#DA1AE5","#542788")
      names(cols)=ne.labels
      for(i in 1:length(ne.labels)){
          load(paste("output/output.",get.label(pop.label,ne.labels[i],max_T_exact_model,old.kinship,which.recent.consang,kinship),".RData",sep=""))
          save.output$ells=save.output$ells[-length(save.output$ells)]
          my.ylim=range(c(save.output$shared_within_inf))
          all.my.ylim=c(all.my.ylim,my.ylim)
            my.ylim=c(0.0002 ,0.0021761522)
          if(i==1 & which.recent.consang==1 & pop.label=="PATHAN_in_Cluster_8"){
            plot(save.output$ells,save.output$shared_within_inf,log="xy",xlab="segment length, in bins of 1cM",ylab="ROH footprint",col=cols[ne.labels[i]],lwd=my.lwd,lty=my.ltys[pop],type="l",ylim=my.ylim,
                 main=paste("ROH footprint, exact model for t< ",max_T_exact_model,sep=""))
          } else {
            lines(save.output$ells,save.output$shared_within_inf,col=cols[ne.labels[i]],lwd=my.lwd,lty=my.ltys[pop])
          }
          rm(save.output)
        }
    }
  }
}
legend("topright",c(paste("N =",c("1k","10k","100k"),sep=" "),paste("kinship = ",round(store.kinship,4),sep="") ),lwd=my.lwd,lty=c(1,1,1,my.ltys[1:2]),col=c(cols,"grey","grey"),cex=1.5)
dev.off()

constant.ne.values=c(1e3,5e3,10e3,20e3,30e3,40e3,50e3,100e3)
ne.labels = c("Ne_trajectory_from_IBDNe_for_50_gens_then_constant", "Ne_trajectory_from_IBDNe_LWR_CI_for_50_gens_then_constant", "Ne_trajectory_from_IBDNe_UPR_CI_for_50_gens_then_constant",
  paste("constant_Ne_",constant.ne.values/1e3,"k",sep=""))
names(ne.labels)[1:3]=c("lower limit of 95% CI of Ne from IBDNe","Ne trajectory from IBDNe","upper limit of 95% CI of Ne from IBDNe")

for(pop in c("PATHAN","Jatt/Choudhry")){
#####specify which population you are modelling
if(pop=="PATHAN"){
  pop1="Pathan"
#  outdir=paste(outdir,"/Pathan/",sep="")
  pop.label="PATHAN_in_Cluster_8"
  pop2="Pathan (Cluster 8)"
  which.scaling.factors0=c(7,8) #specify which scaling factors you want to use for the consanguinity rates vector rates()
} else {
  pop1="Jatt/Choudhry"
  pop.label="JattChoudhry_in_cluster_10"
  pop2="Jatt/Choudhry (Cluster 10)"
  rownames(consang.by.bir2.robust)=gsub("Jatt/Choudhary","Jatt/Choudhry",rownames(consang.by.bir2.robust))
#  outdir=paste(outdir,"/JattChoudhry/",sep="")
  which.scaling.factors0=c(4,5) #specify which scaling factors you want to use for the consanguinity rates vector rates()
}

### calculate kinship values for each of the scaling factors, starting with the empirical consanguinity rates in the consang.by.bir2.robust data frame
kinship.values=sapply(scaling.factors,function(which.recent.consang){
rates = c(0,consang.by.bir2.robust[pop,] * which.recent.consang)
n = length(rates)
kinship = sum(rates * 4^(-(1:n)))
return(kinship)
})


############################################################################################################################################
### effect of changing recent kinship - Figure 5
  my.ylims=c()
max_T_exact_model=50
### loop over old kinship
k=3
old.kinship=c(1e-4,0.05,-1)[k]
### loop over the Ne trajectory
j=1
ne.label=ne.labels[j]
pdf(paste(outdir,"Figure_5.",pop.label,".comparing_effect_of_recent_kinship.",ne.label,".different_ROH_and_IBD_calls.max_T_exact_model_",max_T_exact_model,".old_kinship_",old.kinship,".pdf",sep=""),height=7,width=7)
if(pop=="PATHAN"){
  which.scaling.factors=c(10,which.scaling.factors0,1)
} else {
  which.scaling.factors=c(8,4,5,1)
}
for(i in which.scaling.factors){
  which.recent.consang=scaling.factors[i]
    rates = c(0,consang.by.bir2.robust[pop,] * which.recent.consang)
    n = length(rates)
    kinship = sum(rates * 4^(-(1:n)))
    if(k==3){
      old.kinship=kinship
    }
    load(paste("output/output.",get.label(pop.label,ne.label,max_T_exact_model,old.kinship,which.recent.consang,kinship),".RData",sep=""))
    save.output$ells=save.output$ells[-length(save.output$ells)]
        my.ylim=range(c(save.output$shared_within_inf,save.output$shared_between_inf,footprint.ibdseq.segs.per.bin.by.group2[[pop2]][,2]), footprint.roh.segs.per.bin.by.group[[pop]][,2])
        my.ylims=c(my.ylims,my.ylim)
    if(pop=="Jatt/Choudhry"){
      my.ylim=c( 7.802103e-06 ,2.398668e-3)
    }
  cols=c("#4dac26","#f1b6da","#b8e186","#d01c8b")
  names(cols)=as.character(kinship.values[which.scaling.factors])
  this.col=cols[as.character(kinship.values[i])]

    if(i==which.scaling.factors[1]){
      plot(save.output$ells,save.output$shared_within_inf,log="xy",xlab="segment length, in bins of 1cM",ylab="footprint",col=this.col,lwd=my.lwd,lty=roh.lty,type="l",ylim=my.ylim,
           main=pop,cex.axis=1.5,cex.main=1.5,cex.lab=1.5)
    } else {
      lines(save.output$ells,save.output$shared_within_inf,col=this.col,lwd=my.lwd,lty=roh.lty)
    }
    lines(save.output$ells,save.output$shared_between_inf,col=this.col,lwd=my.lwd,lty=ibd.lty)
    rm(save.output)
}
lines(footprint.ibdseq.segs.per.bin.by.group2[[pop2]][,1],footprint.ibdseq.segs.per.bin.by.group2[[pop2]][,2],lty=ibd.lty,lwd=my.lwd,col="black")
lines(footprint.roh.segs.per.bin.by.group[[pop]][,1],footprint.roh.segs.per.bin.by.group[[pop]][,2],col="black",lwd=my.lwd,lty=roh.lty)
  
legend("bottomleft",c("IBD","ROH","observed",paste("expected with recent kinship = ", round(sort(kinship.values[which.scaling.factors],decreasing=TRUE),4),sep="")),lwd=my.lwd,
       lty=c(ibd.lty,roh.lty,rep(1,length(which.scaling.factors)+1)),
       col=c("grey","grey","black",cols[c(1,3,2,4)]),cex=1,bty="n")

dev.off()

print(range(my.ylims))



############################################################################################################################################
### effect of changing recent kinship, but also showing ROH and IBD calls with different filtering - Supplementary Figure 18
  my.ylims=c()
max_T_exact_model=50
### loop over old kinship
k=3
old.kinship=c(1e-4,0.05,-1)[k]
pdf(paste(outdir,"Figure_S18.",pop.label,".comparing_effect_of_recent_kinship.varying_Ne_trajectory.max_T_exact_model_",max_T_exact_model,".old_kinship_",old.kinship,".pdf",sep=""),height=8,width=12)
par(mfrow=c(2,3))
### loop over the Ne trajectory
for(j in 1:3){
  ne.label=ne.labels[j]
if(pop=="PATHAN"){
which.scaling.factors=c(10,which.scaling.factors0,1)
} else {
which.scaling.factors=c(8,4,5,1)
}
for(i in which.scaling.factors){
which.recent.consang=scaling.factors[i]
rates = c(0,consang.by.bir2.robust[pop,] * which.recent.consang)
n = length(rates)
kinship = sum(rates * 4^(-(1:n)))
if(k==3){
  old.kinship=kinship
}
load(paste("output/output.",get.label(pop.label,ne.label,max_T_exact_model,old.kinship,which.recent.consang,kinship),".RData",sep=""))
save.output$ells=save.output$ells[-length(save.output$ells)]
my.ylim=range(c(save.output$shared_within_inf,save.output$shared_between_inf,footprint.ibdseq.segs.per.bin.by.group2[[pop2]][,2], footprint.roh.segs.per.bin.by.group[[pop]][,2],
  footprint.ibdseq.segs.per.bin.by.group2.raw[[pop2]][,2],compare.footprint.roh.segs.per.bin.by.group[["bcftools.raw"]][[pop]][,2],compare.footprint.roh.segs.per.bin.by.group[["GARLIC.raw"]][[pop]][,2],
  compare.footprint.roh.segs.per.bin.by.group[["PLINK.raw"]][[pop]][,2]))
my.ylims=c(my.ylims,my.ylim)
cols=c("#4dac26","#f1b6da","#b8e186","#d01c8b")
names(cols)=as.character(kinship.values[which.scaling.factors])
this.col=cols[as.character(kinship.values[i])]
if(i==which.scaling.factors[[1]]){
  plot(save.output$ells,save.output$shared_within_inf,log="xy",xlab="segment length, in bins of 1cM",ylab="footprint",col=this.col,lwd=my.lwd,lty=roh.lty,type="l",ylim=my.ylim,
       main=names(ne.labels)[j],cex.axis=1.5,cex.main=1.5,cex.lab=1.5)
} else {
          lines(save.output$ells,save.output$shared_within_inf,col=this.col,lwd=my.lwd,lty=roh.lty)
        }
        lines(save.output$ells,save.output$shared_between_inf,col=this.col,lwd=my.lwd,lty=ibd.lty)
        rm(save.output)
  }
  lines(footprint.ibdseq.segs.per.bin.by.group2[[pop2]][,1],footprint.ibdseq.segs.per.bin.by.group2[[pop2]][,2],lty=ibd.lty,lwd=my.lwd,col="black")
  lines(footprint.roh.segs.per.bin.by.group[[pop]][,1],footprint.roh.segs.per.bin.by.group[[pop]][,2],col="black",lwd=my.lwd,lty=roh.lty)

  lines(footprint.ibdseq.segs.per.bin.by.group2.raw[[pop2]][,1],footprint.ibdseq.segs.per.bin.by.group2.raw[[pop2]][,2],lty=ibd.lty,lwd=my.lwd,col="blue")
  lines(compare.footprint.roh.segs.per.bin.by.group[["bcftools.raw"]][[pop]][,1],compare.footprint.roh.segs.per.bin.by.group[["bcftools.raw"]][[pop]][,2],col="blue",lwd=my.lwd,lty=roh.lty)

  lines(compare.footprint.roh.segs.per.bin.by.group[["GARLIC.raw"]][[pop]][,1],compare.footprint.roh.segs.per.bin.by.group[["GARLIC.raw"]][[pop]][,2],col="#2c7fb8",lwd=my.lwd,lty=roh.lty)
  lines(compare.footprint.roh.segs.per.bin.by.group[["PLINK.raw"]][[pop]][,1],compare.footprint.roh.segs.per.bin.by.group[["PLINK.raw"]][[pop]][,2],col="purple",lwd=my.lwd,lty=roh.lty)

}

plot(0,0,col=NA)
plot(0,0,col=NA)
    legend("bottomleft",c("IBD","ROH","observed IBD - IBDseq, filtered","observed IBD - IBSseq, unfiltered",
"observed ROH - bcftools/roh, filtered","observed ROH - bcftools/roh, unfiltered","observed ROH - GARLIC","observed ROH - PLINK",
                          "expected based on model with:",paste("kinship = ", round(sort(kinship.values[which.scaling.factors],decreasing=TRUE),3),sep="")),lwd=my.lwd,
           lty=c(ibd.lty,roh.lty,ibd.lty,ibd.lty,roh.lty,roh.lty,roh.lty,roh.lty,NA,rep(1,length(which.scaling.factors))),
           col=c("grey","grey","black","blue","black","blue","#2c7fb8","purple",NA,cols[c(1,3,2,4)]),cex=1.3,bty="n")
plot(0,0,col=NA)
dev.off()
print(range(my.ylims))


  if(pop=="PATHAN"){
### effect of changing old kinship - Supplementary Figure 19
max_T_exact_model=50
  #loop over new kinship

i=7

which.recent.consang=scaling.factors[i]
rates = c(0,consang.by.bir2.robust[pop,] * which.recent.consang)
n = length(rates)
kinship = sum(rates * 4^(-(1:n)))
cat(pop,"\t",kinship,"\n")
pdf(paste(outdir,"Figure_S19.",pop.label,".comparing_effect_of_kinship_prior_to_max_T_exact_model.varying_Ne_trajectory.max_T_exact_model_",max_T_exact_model,".recent_kinship_",kinship,".pdf",sep=""),height=7,width=7)

### loop over the Ne trajectory

j=1
ne.label=ne.labels[j]
### loop over old kinship
cols=c("#1b9e77","#d95f02","#7570b3")
for(k in 1:3){
old.kinship=c(1e-4,0.05,kinship)[k]
load(paste("output/output.",get.label(pop.label,ne.label,max_T_exact_model,old.kinship,which.recent.consang,kinship),".RData",sep=""))
save.output$ells=save.output$ells[-length(save.output$ells)]
my.ylim=range(c(save.output$shared_within_inf,save.output$shared_between_inf,footprint.ibdseq.segs.per.bin.by.group2[[pop2]][,2]), footprint.roh.segs.per.bin.by.group[[pop]][,2])
this.col=cols[k]
if(k==1){
  plot(save.output$ells,save.output$shared_within_inf,log="xy",xlab="segment length, in bins of 1cM",ylab="footprint",col=this.col,lwd=my.lwd,lty=roh.lty,type="l",ylim=my.ylim,
       main=pop,cex=0.9,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
} else {
  lines(save.output$ells,save.output$shared_within_inf,col=this.col,lwd=my.lwd,lty=roh.lty)
}
lines(save.output$ells,save.output$shared_between_inf,col=this.col,lwd=my.lwd,lty=ibd.lty)
rm(save.output)
}
lines(footprint.ibdseq.segs.per.bin.by.group2[[pop2]][,1],footprint.ibdseq.segs.per.bin.by.group2[[pop2]][,2],lty=ibd.lty,lwd=my.lwd,col="black")
lines(footprint.roh.segs.per.bin.by.group[[pop]][,1],footprint.roh.segs.per.bin.by.group[[pop]][,2],col="black",lwd=my.lwd,lty=roh.lty)

legend("bottomleft",c("IBD","ROH","observed",paste("expected with old kinship = ", c(1e-4,0.05,round(kinship,4)),sep="")),lwd=my.lwd,
       lty=c(ibd.lty,roh.lty,rep(1,length(which.recent.consang)+1)),
       col=c("grey","grey","black",cols),cex=1,bty="n")
dev.off()

}


}
