########################         二、 symbol 转为 ENSEMBL              ######################
rm(list = ls())
options(stringsAsFactors = F)

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

getwd()


sig_genefile="hub.txt"
#读入显著基因
sig_gene=data.table::fread(sig_genefile, header = F)
genes=unique(as.vector(sig_gene$V1))
#基因格式转换，ENSEMBL转为symbol
ensembl <- bitr(genes,
                fromType = "SYMBOL",#现有的ID类型
                toType = "ENSEMBL",#需转换的ID类型
                OrgDb = "org.Hs.eg.db",drop =F)
ensembl2 <- na.omit(ensembl)
write.csv(ensembl2,file = "2DEG_ensembl_noNA.csv",row.names = F)
save(ensembl2,file ='2DEG_ensembl_noNA.rdata')



########################         三、孟德尔随机化 暴露：差异基因eqlt  结局：gwas疾病            ######################
rm(list = ls())
options(stringsAsFactors = F)
dir.create("3allgene")
library(TwoSampleMR) # 加载TwoSampleMR库
library(MRPRESSO) # 加载MRPRESSO库
library(tidyverse) # 加载tidyverse库
library(GagnonMR) # 加载GagnonMR库

load("2DEG_ensembl_noNA.rdata") # 加载2DEG_ensembl_noNA.rdata数据文件
load("all_eqtlgen_id.rdata") # 加载all_eqtlgen_id.rdata数据文件

out_id="finn-b-C3_GBM_EXALLC" # 设置结局的IEU id

genes <- unique(ensembl2$ENSEMBL) %>% sort() # 提取ensembl2中ENSEMBL列的唯一值并排序
genes <- genes[genes %in% eqtl_id$trait] # 筛选出在eqtl_id的trait列中出现过的基因

# 提取暴露的工具变量
Allres <- data.frame() # 创建一个空数据框

j=1    # 如果循环中断，就把1改为最后那个数字再次开始循环

# 开始循环，需要联网
gene <- genes[1:length(genes)] # 设置gene为genes中的所有基因
# 设置最大重试次数
max_attempts <- 5

# 创建一个空的向量来保存失败的基因
failed_genes <- c()

# 把基因列表分成多个批次
batch_size <- 100
gene_batches <- split(gene, ceiling(seq_along(gene)/batch_size))

# 对每个批次进行处理
for (batch in seq_along(gene_batches)) {
  # 获取当前批次的基因
  genes <- gene_batches[[batch]]

  # 对每个基因进行处理
  for (i in genes) {
    print(j)
    attempts <- 0  # 设置当前重试次数

    while(attempts <= max_attempts) {
      attempts <- attempts + 1
      tryCatch({
        exp <- extract_instruments(outcomes = paste0('eqtl-a-',i))
        if (is.null(exp)) {
          j <- j+1
          break
        }
        else{
          out <- extract_outcome_data(snps = exp$SNP, outcomes = out_id)
          if (is.null(out)) {
            j <- j+1
            break
          }
          else{
            dat <- TwoSampleMR::harmonise_data(exposure_dat = exp, outcome_dat = out,action = 2)
            dat <- subset(dat, dat$mr_keep == TRUE)
            if (nrow(dat) == 0) {
              j <- j+1
              break
            }
            else{
              res <- GagnonMR::primary_MR_analysis(dat = dat)
              if (res$pval>0.05 | res$pval=='NaN') {
                j <- j+1
                break
              }
              else{
                j <- j+1
                print(j)
                save(exp,out,dat,res,file = paste0('./3allgene/',i,"_mrResult.rdata"))
                Allres <- rbind(Allres,res)
              }
            }
          }
        }
      }, error = function(e) {
        if(attempts <= max_attempts) {
          print(paste("Attempt", attempts, "for gene", i, "failed. Retrying..."))
          Sys.sleep(10)  # 添加时间延迟
        } else {
          print(paste("Attempt", attempts, "for gene", i, "failed. Skipping this gene..."))
          failed_genes <- c(failed_genes, i)  # 添加失败的基因到列表
        }
      })
      if (attempts > max_attempts) {
        next  # Skip to the next gene if the maximum number of attempts is reached
      }
    }
  }
  # 保存批次结果
  save(Allres, file = paste0('./3allgene/batch_', batch, "_result.rdata"))
}

# 打印失败的基因
print(failed_genes)



#全部循环完再运行，中间勿运行其他程序
#检测重复SNP
duplicated_rows <- duplicated(Allres$id.exposure)
#删除重复SNP
Allres <- Allres[!duplicated_rows, ]


#save(Allres,file = paste0(afdir,"/3AlleQTL_mrres.rdata"))
save(Allres,file = paste0("3AlleQTL_mrres.rdata"))
#write.csv(Allres,file = paste0(afdir,"/3AlleQTL_mrres.csv"))
write.csv(Allres,file = paste0("3AlleQTL_mrres.csv"))



########################         四、 ENSEMBL转为symbol            ######################
rm(list = ls())
options(stringsAsFactors = F)
library(clusterProfiler)
library(org.Hs.eg.db)

load("3AlleQTL_mrres.rdata")

ensem <-sapply(strsplit(Allres$exposure,split = ' '),'[',1)

#基因格式转换，ENSEMBL转为symbol
genes <- bitr(ensem,
              fromType = "ENSEMBL",#现有的ID类型
              toType = "SYMBOL",#需转换的ID类型
              OrgDb = "org.Hs.eg.db",drop =F)

write.csv(genes,file = "4all_ensem_genes.csv",row.names = F)
save(genes,file ='4all_ensem_genes.rdata')



########################         五、 绘图           ######################
rm(list = ls())
options(stringsAsFactors = F)
afdir <- paste0(getwd(),"/eqtl-gwas")

library(tidyverse)
genes=read_csv('4all_ensem_genes.csv')
#load("4all_ensem_genes.rdata")
load("3AlleQTL_mrres.rdata")

Allres$ensem <- sapply(strsplit(Allres$exposure,split = ' '),'[',1)
Allres$gene <- genes$SYMBOL[match(genes$ENSEMBL,Allres$ensem)]

OR <-generate_odds_ratios(Allres)
str(OR)
#加载ggplot2包
library(ggplot2)

#绘制森林图
pdf(file = paste0(afdir,"/5.森林图.pdf"),width=10,height=8)
ggplot(data = OR, aes(y = id.exposure, x = or, xmin = or_lci95, xmax = or_uci95)) +
  geom_point() +
  geom_errorbarh(height = .1) +
  scale_y_discrete(name = "", labels = OR$gene) +
  scale_x_log10(name = "Odds Ratio") +
  labs(title = "Odds Ratio by Gene", y = "Gene") +
  geom_vline(xintercept = 1, color = "black", linetype = "dashed", alpha = .5) +
  theme_minimal()
dev.off()


#开始画图
head(Allres)
pdf(file = paste0(afdir,"/5.Dot_Alleqtl_MR.pdf"),width=10,height=8)
#ggsave(filename = 'Dot_Alleqtl_MR.pdf',width=10,height=8)
ggplot(Allres,aes(b,forcats::fct_reorder(gene, b))) +
  geom_segment(aes(xend=0, yend = gene)) +
  geom_point(aes(color=pval, size = abs(b))) +
  scale_size_continuous(range=c(3, 9)) +
  theme_bw() +
  xlab("b") +
  ylab(NULL) +
  ggtitle("All MR results")
dev.off()
#x 轴是 MR 分析结果中的 MR 估计值 b，y 轴是基因名称 gene，
#每个基因的颜色代表了 MR 分析中的 P 值 pval，点的大小则表示 MR 估计值的绝对值



