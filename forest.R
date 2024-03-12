library(forestplot)
OR=read.csv("PD_outcome.csv")
hz <- paste(round(OR$or,3),
            "(",round(OR$or_lci95,3),
            "-",round(OR$or_uci95,3),")",sep = "")


tabletext <- cbind(c(NA,"PD as outcome ",OR[,1]),###改
                   c(NA,"nSNPs",round(OR$nsnp)),
                   c(NA,"P",ifelse(OR$pval<0.001,"P < 0.001",round(OR$pval,3))),
                   c(NA,"OR(95% CI)",hz))
tabletext[tabletext=="NA(NA-NA)"] <- NA
forestplot(labeltext=tabletext, 
                      graph.pos=3,  #为Pvalue箱线图所在的位置
                      col=fpColors(box="#0020C2", lines="#3BB9FF", zero = "gray50"),
                      mean=c(NA,NA,OR$or),
                      lower=c(NA,NA,OR$or_lci95), #95%置信区间下限
                      upper=c(NA,NA,OR$or_uci95), #95%置信区间上限
                      xlab="OR",
                      boxsize=0.2,lwd.ci=2,   #箱子大小，线的宽度
                      ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
                      zero=1,lwd.zero=1,      #zero线宽 基准线的位置
                      colgap=unit(10,"mm"),    #列间隙
                      xticks = c(0.8,0.9,1,1.1,1.2), #横坐标刻度
                      lwd.xaxis=1,            #X轴线宽
                      lineheight = unit(2,"cm"),#固定行高
                      graphwidth = unit(.2,"npc"), #图在表中的宽度比例
                      cex=0.9, fn.ci_norm = fpDrawNormalCI, #误差条显示方式
                      hrzl_lines=list("2" = gpar(lwd=2, col="black"),
                                      "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
                                      "48" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"29"中数字为nrow(tabletext)+1
                      mar=unit(rep(0.5, times = 4), "mm"),#图形页边距
                      clip=c(0.8,1.2),
                      is.summary=c(F,T,T,F,F,F,F,T,F,F,F,F,T,F,F,F,F,T,F,F,F,F,T,F,F,F,F,T,F,F,F,F,T,F,F,F,F,T,F,F,F,F,T,F,F,F,F,F,F,T),
                      # #fpTxtGp函数中的cex参数设置各个组件的大小
                      txt_gp=fpTxtGp(label=gpar(cex=1),
                                     ticks=gpar(cex=1),
                                     xlab=gpar(cex = 1),
                                     title=gpar(cex = 1))
)
