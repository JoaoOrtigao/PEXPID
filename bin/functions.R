################################################

CREATE_CTS = function(DATA_DIR){
  
  # DATA_DIR => dir containing STAR output
  
  
}
  
################################################

IMPORT_CTS = function(FILE){
  
  # file => table containing raw counts 
  
  cts_all = read.table(FILE)
  row.names(cts_all) = cts_all$ensembl
  cts_all = cts_all[,-1]
  return(cts_all)
  
}

################################################

IMPORT_COLDATA = function(FILE){
  
  coldata = read.table(FILE,sep = ";",header = T)
  
  row.names(coldata) = coldata$Run
  
  return(coldata)
  
}

################################################

FILTER_CTS = function(cts_all,coldata,BIOPROJECT){
  
  # return a list containing coldata and cts filtered by BIOPROJECT
  
  coldata_tmp = coldata[coldata$BioProject == BIOPROJECT,]
  row.names(coldata_tmp) = coldata_tmp$Run
  
}

