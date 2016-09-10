
suppressMessages(library(cbcbSEQ))


setwd("/Users/Chris/Documents/strawberry/ZCL/counts/")

count_files <- c("rpkm_hd_lcm_zh.csv", "rld_hd_lcm_zh.csv", "lrpkm_hd_lcm_zh.csv",
                 "log_raw_hd_lcm_zh.csv", "lcpm_hd_lcm_zh.csv", "cpm_hd_lcm_zh.csv",
                 "hd_lcm_zh.csv")

groups <- c("hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd","hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd","hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd","hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd","hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd","hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd","hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd", "hd","hd","hd","hd","hd","lcm","lcm","lcm","lcm","lcm","lcm","lcm","lcm","lcm","lcm","zh","zh","zh","zh","zh","zh","zh","zh")

dates <- read.delim(file = "/Users/Chris/Documents/strawberry/ZCL/counts/sample_dates.txt", sep="\t")


for(read_file in count_files){

   #read_file <- "hd_lcm_zh.csv"
  print(read_file)

  input <- read.csv(read_file)

  if (read_file == "hd_lcm_zh.csv"){
    input$X.1 <- NULL
  }

  dates$combo <- dates$date#paste(dates$date, "_", dates$lane, sep = "")

  #change rownames to first colomn
  rownames(input) <- input$X
  input$X <- NULL


  #now subsetting so transpose briefly
  t_input <- as.data.frame(t(input))

  t_input$groups <- groups

  #mtch batch info with row names
  matches <- match(rownames(t_input), dates$sampleID, nomatch=nrow(dates))
  #add batch info to data
  t_input$batch <- dates[matches,"combo"]
  #rename
  #coerce to character
  t_input$batch <- as.character(t_input$batch)

  #create design df
  design <- data.frame(row.names= rownames(t_input), method = t_input$groups, date=t_input$batch)

  #make SVD
  res <- makeSVD(input)

  pcr <- pcRes(res$v, res$d, design$method, design$date)

  write.table(pcr , file = paste("/Users/Chris/Documents/strawberry/ZCL/graphs/pca/", read_file,"/pcres.txt",sep=""), sep = "\t")

}
