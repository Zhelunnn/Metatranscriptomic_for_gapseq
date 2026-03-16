############################# Running NOISeq for all MAGS ##################################

library(ggplot2)
library(NOISeq)
library(stringr)
directories<-c('030616_reaction','050416_reaction','050617_reaction','050716_reaction','061016_reaction','070217_reaction','150816_reaction')
#model_names<-unlist(str_split('bin_10 bin_16-notfinal CCGammaproteobacteria_3 bin_17 bin_21 bin_22-polished bin_6 CCAlphaproteobacteria_3 CCAlphaproteobacteria_6 CCAlphaproteobacteria_7 CCBdellovibrionales_7 CCGammaproteobacteria_9 CCLitorimonas_sp_1 CCNitrosopumilus_sp_1 CCNitrospirales_2 CCPseudohongiella_sp_1 CCRhizobiaceae_1 CCRhizobiaceae_2 CCRhizobiaceae_3 CCRobiginitomaculum_sp_1 CCSynoicihabitans_sp_1 CCThermodesulfobacteriota_1 CCVerrucomicrobiales_1',' '))
model_names<-unlist(str_split('bin_16_ctg8312 alt_Areni bin_22 CCDehalococcoidia_1'))

saturation_list<-list()
# Loop over the file names
for (model in model_names){
  Organism<-c() ## the organism name_sample date (e.g. bin_10_030616)
  Date<-c() ## store the sample date which is used as condition
  sample_date_list<-list() ## store the df containing gene expression level under this data (3 replicates)
  for(directory in directories) {
    print(model)
    print(directory)
    # List all files that end with '_reaction_reads_count.csv'
    file_names <- list.files(path=directory,pattern = "_reaction_reads_count\\.tsv$")
    target_file<-file_names[grepl(model,file_names)] ## get the target file
    
    # Read only the first 4 columns
    o_date<-gsub('_reaction_reads_count.tsv','',target_file)

    df <- read.table(paste(directory,'/',target_file,sep=''), sep = "\t", header = TRUE, colClasses = c("character", "numeric", "numeric", "numeric", rep("NULL", 8)))
    ##remove bad blast features
    df<-df[!grepl('_bad_blast_',df$Biocyc_ID),]
    
    ### check and remove if any sample (column) only contains 0
    df<-df[,!apply(df,2,function(x) all(x==0))]
    if (all(df==0)){
      cat('None expressed genes are detected in ',o_date,'\n')
      next
    }
    colnames(df) <- c('Biocyc_ID', paste(o_date, '_', seq_len(ncol(df) - 1), sep = ''))
    sample_date_list[[o_date]]<-df
    
    Organism<-c(Organism,rep(o_date,ncol(df)-1))
    d<-sub(".*_", "", o_date)
    Date<-c(Date,rep(d,ncol(df)-1))
  }
  merged_df <- Reduce(function(x, y) merge(x, y, by = "Biocyc_ID", all = TRUE), sample_date_list)
  

  rownames(merged_df)<-merged_df[,1];merged_df<-merged_df[,-1];merged_df[is.na(merged_df)] <- 0
  myfactors<-data.frame(sample=Organism,date=Date)
  # Perform the analysis with the current file's data
  
  mydata <- readData(data = merged_df, factors = myfactors) # Ensure readData and subsequent functions are correctly specified for your context
  mysaturation = dat(mydata, k = 0, ndepth = 10, type = "saturation", factor = "date") # Adjust function call as needed

  saturation_list[[model]]<-mysaturation
}

# Open a graphics device to save the plot
plot_filename <- sub("_combined\\.csv$", "", file_name)
png(paste0(plot_filename, ".png"))

explo.plot(saturation_list[["bin_6"]], toplot = 1, samples = 1:12, yleftlim = NULL, yrightlim = NULL)

explo.plot(a)

length(saturation_list[[1]]@dat$depth)
