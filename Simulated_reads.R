library(readr)
library(ggplot2)

gold <- read_delim("/media/Disco3/santiagor/Datos_TCruzi_Optimizacion_Mapeos/Simulated-Reads/gold_standard.tab", "\t", escape_double = FALSE, comment = "#",trim_ws = TRUE)
gold <- gold[,c("Geneid","Shortstack_Epi.bam")]
colnames(gold) <- c("id","gold")
bowtie <- read_delim("/media/Disco3/santiagor/Datos_TCruzi_Optimizacion_Mapeos/Simulated-Reads/bowtie_counts.tab", "\t", escape_double = FALSE, comment = "#",trim_ws = TRUE)
bowtie <- bowtie[,c("Geneid","bwt_simulated.sam")]
colnames(bowtie) <- c("id","bowtie")
shortstack <- read_delim("/media/Disco3/santiagor/Datos_TCruzi_Optimizacion_Mapeos/Simulated-Reads/shortstack_counts.tab", "\t", escape_double = FALSE, comment = "#",trim_ws = TRUE)
shortstack <- shortstack[,c("Geneid","ShortStack/simulated.bam")]
colnames(shortstack) <- c("id","shortstack")

temp <- merge(gold, bowtie, by=c("id"))
merged <- merge(temp, shortstack, by=c("id"))
rm(temp)


corr_bwt_spearman <- cor.test(x=merged$gold, y=merged$bowtie, method = 'spearman')
corr_shortstack_spearman <- cor.test(x=merged$gold, y=merged$shortstack, method = 'spearman')

corr_bwt_pearson <- cor.test(x=merged$gold, y=merged$bowtie, method = 'pearson')
corr_shortstack_pearson <- cor.test(x=merged$gold, y=merged$shortstack, method = 'pearson')


scatter_plot <- ggplot(merged, aes(log(gold), log(bowtie)))
scatter_plot + geom_point() + labs(title="All Counts: Gold vs Bowtie",x = "Log(Gold)", y = "Log(Bowtie)") + geom_smooth(method="lm")

scatter_plot <- ggplot(merged, aes(log(gold), log(shortstack)))
scatter_plot + geom_point() + labs(title="All Counts: Gold vs ShortStack",x = "Log(Gold)", y = "Log(ShortStack)") + geom_smooth(method="lm")




# Ts subset


Ts_list <- read_csv("/media/Disco3/santiagor/Datos_TCruzi_Optimizacion_Mapeos/Simulated-Reads/Ts_list.txt", col_names = FALSE)

Ts_counts <- subset(merged,merged$id %in% Ts_list$X1)

corr_bwt_spearman_TS <- cor.test(x=Ts_counts$gold, y=Ts_counts$bowtie, method = 'spearman')
corr_shortstack_spearman_TS <- cor.test(x=Ts_counts$gold, y=Ts_counts$shortstack, method = 'spearman')

corr_bwt_pearson_TS <- cor.test(x=Ts_counts$gold, y=Ts_counts$bowtie, method = 'pearson')
corr_shortstack_pearson_TS <- cor.test(x=Ts_counts$gold, y=Ts_counts$shortstack, method = 'pearson')

scatter_plot <- ggplot(Ts_counts, aes(log(gold), log(bowtie)))
scatter_plot + geom_point() + labs(title="TS Counts: Gold vs Bowtie",x = "Log(Gold)", y = "Log(Bowtie)") + geom_smooth(method="lm")

scatter_plot <- ggplot(Ts_counts, aes(log(gold), log(shortstack)))
scatter_plot + geom_point() + labs(title="TS Counts: Gold vs ShortStack",x = "Log(Gold)", y = "Log(ShortStack)") + geom_smooth(method="lm")



