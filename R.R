library(data.table)
library(TwoSampleMR)
library(gwasglue)
library(VariantAnnotation)
library(ieugwasr)
library(vroom)
library(MRPRESSO)
library(tidyverse)
source("get_f.R")
bim_vcf <- readVcf("ieu-b-7.vcf.gz")
a <- gwasvcf_to_TwoSampleMR(vcf=bim_vcf,type="exposure")
b<-subset(a,pval.exposure<5e-8)
b$id=b$exposure
c=ieugwasr::ld_clump_local(dplyr::tibble(rsid=b$SNP, pval=b$pval.exposure, id=b$id),
                           clump_r2=0.001,clump_kb=10000,clump_p=1,
                           plink_bin="D:/R-4.3.0/library/plinkbinr/bin/plink_Windows.exe", ## 欧洲的EUR
                           bfile="E:/Rword/OSA/1kg.v3/1kg.v3/EUR")
exp_dat=b[match(c$rsid,b$SNP),]
setwd("eight")
ls=list.files(".")
case_id=unlist(str_split(ls,"[.]",simplify=T))[,1]
Allres=data.frame()
Allheterogeneity=data.frame()
AllPRESSO=data.frame()
for (i in 1:8) {
  a2 <- fread(paste0(ls[i]),header = T)
  pick=a2[match(c$rsid,a2$rsids),]
  pick$id=paste0(case_id[i])
  out_dat<- format_data(pick,
                        type = "outcome",
                        snp_col ="rsids",
                        phenotype_col = "id",
                        beta_col = "beta",
                        se_col = "sebeta",
                        eaf_col = "af_alt_cases",
                        effect_allele_col ="alt",
                        other_allele_col = "ref",
                        pval_col = "pval"
                        ,chr_col = "#chrom",
                        pos_col = "pos",
                        id_col = "id")
  dat <- harmonise_data(
    exposure_dat = exp_dat,
    outcome_dat = out_dat
  )
  dir.create(paste0(case_id[i],"_result"))
  setwd(paste0(case_id[i],"_result"))           
  res_presso <-mr_presso(BetaOutcome ="beta.outcome", 
                         BetaExposure = "beta.exposure", 
                         SdOutcome ="se.outcome", 
                         SdExposure = "se.exposure",
                         OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = dat
                         , NbDistribution = 1000,SignifThreshold = 0.05)
  write.csv(res_presso$`MR-PRESSO results`$`Global Test`$Pvalue,
            file = paste0(case_id[i], "_PRESSO.csv"))
  res <- mr(dat)
  save(res, file = paste0(case_id[i], "_result.rdata"))
  OR <-generate_odds_ratios(res)
  write.csv(OR,file = paste0(case_id[i], "_OR.csv"))
  #异质性检测
  heterogeneity <- mr_heterogeneity(dat)
  write.csv(heterogeneity,file = paste0(case_id[i], "_heterogeneity.csv"))
  #多效性检验
  pleiotropy <- mr_pleiotropy_test(dat)
  write.csv(pleiotropy,file = paste0(case_id[i], "_pleiotropy.csv"))
  Allres <- rbind(Allres,res)
  Allheterogeneity<- rbind(Allheterogeneity,heterogeneity)
  AllPRESSO <- rbind(AllPRESSO,res_presso$`MR-PRESSO results`$`Global Test`$Pvalue)
  setwd("..")
}
#散点图
mr_scatter_plot(results,dat)
#逐个剔除检验，做留一图
leaveoneout <- mr_leaveoneout(dat)
mr_leaveoneout_plot(leaveoneout)
#森林图
results_single <- mr_singlesnp(dat)
mr_forest_plot(results_single)
#漏斗图
mr_funnel_plot(results_single)


