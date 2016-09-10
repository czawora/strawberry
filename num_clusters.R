##exploring the results of WGCNA clustering


setwd("/Users/Chris/Documents/strawberry/ZCL/wgcna_out")

dirs <- read.table("dir.txt")
dirs <- as.character(dirs$V1)

config_clusters <- data.frame(config = character(), cluster_num = double(), subset = character(), norm = character(), power = double(), modsize= double(), merge=double())

for(name in dirs){

  print(name)
  splits <- strsplit(name, "_")
  #print(splits)
  merged <- splits[[1]][length(splits[[1]])]

  if(splits[[1]][1] == "hd"){
    if(merged == 1){
      clusters <- read.table(paste(name, "/merged_clusters.csv", sep=""), sep="," , header=TRUE)
    } else {
      clusters <- read.table(paste(name, "/original_clusters.csv", sep=""), sep="," , header=TRUE)
    }


    num_cluster <- length(unique(clusters$group))

    if(splits[[1]][1] == "hd"){
      config_clusters <- rbind(config_clusters, data.frame(config=name, cluster_num=num_cluster, subset = "hd_zh", norm=splits[[1]][4], power=as.numeric(splits[[1]][6]), modsize=as.numeric(splits[[1]][8]), merge = as.numeric(merged)))
    }
#     } else {
#
#       config_clusters <- rbind(config_clusters, data.frame(config=name, cluster_num=num_cluster, subset = "lcm", norm=splits[[1]][3], power=splits[[1]][5], modsize=splits[[1]][7], merge = merged))
#
#     }
  }
}

# config_clusters$norm <- as.numeric(config_clusters$norm)
# config_clusters$power <- as.numeric(config_clusters$power)
# config_clusters$modsize <- as.numeric(config_clusters$modsize)
#

write.csv(config_clusters, file = "config_clusters.csv")


###read the config_back in

#40, 60, 90, 120, 150, 180, 210
config_clusters$color = "black"
config_clusters$color[config_clusters$modsize >= 120] = "red"
config_clusters$color[config_clusters$modsize < 120] = "blue"

sub <- config_clusters[ which(config_clusters$cluster_num <= 100), ]
sub2 <- config_clusters[ which(config_clusters$cluster_num > 100), ]

png(filename = "/Users/Chris/Documents/strawberry/ZCL/under100_hist.png")
hist(sub$cluster_num, main="hisotgram of modules created under 100 clusters", xlab= "")
dev.off()

print(sum(sub$merge))

ordered <- config_clusters[order(config_clusters$cluster_num), ]

png(filename = "/Users/Chris/Documents/strawberry/ZCL/cluster_plot.png", height=2000, width=3000)
plot(ordered$cluster_num, main="cluster per configuration", type="p", xlab="", ylab="# of clusters", xaxt='n', col = ordered$color, pch=16)
axis(1, at=seq(1,nrow(ordered),by=1), labels=ordered[, "modsize"], las = 2)
dev.off()

##plot cluaster disribution historgrams for each normalization type
lrpkm_rows <- which(ordered$norm =="lrpkm")
lcpm_rows <- which(ordered$norm =="lcpm")
log_rows <- which(ordered$norm =="log")
rld_rows <- which(ordered$norm =="rld")

png(filename = "/Users/Chris/Documents/strawberry/ZCL/rld_hist.png")
hist(config_clusters[rld_rows, "cluster_num"], main="hisotgram of modules created for rld", xlab= "")
dev.off()


#hd_zh_sub_rld_power_4_modulesize_180_TOM_0_merge_0
#hd_zh_sub_lrpkm_power_1_modulesize_90_TOM_0_merge_1

#testing checking module overlap of two results



