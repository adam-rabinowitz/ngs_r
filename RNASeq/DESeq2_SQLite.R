require(DESeq2)
source('~/github/ngs_r/RNASeq/DESeq2_Functions.R')
setwd('/camp/stp/babs/working/rabinoa/projects/stockingerg/amina.metidji/')
dds <- createDESeqDataSet(
  'alignments/pairedEnd/combined_gene.exp.txt',
  'deseq2/pairedEnd/allSamples/IN_KOvsWT/IN_KOvsWT.comparison.txt'  
)
ddsAnalysis <- DESeq(dds)


###############################################################################
## Function to create a SQLite datbase from a DESeq2DataSet
###############################################################################
createDESeq2SQLite <- function(dbname, DESeq2DataSet, altnames=NA) {
  # Connect to database
  if (file.exists(dbname)) {
    stop('File already exists')
  }
  db <- dbConnect(SQLite(), dbname=dbname, flags=SQLITE_RWC)
  # Read in results table
  resultData <- results(DESeq2DataSet)
  # Create results table
  resultsTable <- 'CREATE TABLE results(\"GENE\" TEXT PRIMARY KEY, \"'
  resultsTable <- paste0(resultsTable, paste0(
    colnames(resultData), collapse = '" REAL, "'))
  resultsTable <- paste0(resultsTable, '\" REAL)')
  dbExecut
  # Read in counts
  counts <- as.data.frame(counts(DESeq2DataSet, normalized=F))
  # Create counts table
  countsTable <- 'CREATE TABLE counts(\"GENE\" TEXT PRIMARY KEY, \"'
  countsTable <- paste0(countsTable, paste0(
    colnames(counts), collapse = '" REAL, "'))
  countsTable <- paste0(countsTable, '\" REAL)')
  dbExecute(db, countsTable)
  # Reformat counts
  counts$GENE <- row.names(counts)
  counts <- counts[,c(ncol(counts), 1:(ncol(counts)-1))]
  print(head(counts))
  dbExecute(db, countsTable)
  dbWriteTable(db, 'counts', counts, append=T, row.names=F)
  dbDisconnect(db)
}

###############################################################################
## Add dataframe to SQLite database
###############################################################################
createTableFromDataFrame <- function(df, name, db, pkCol=1) {
  # Extract data types
  types <- apply(df, 2, function(z) {
    if (is.numeric(z)) {
      return('REAL')
    } else if (is.character(z)) {
      return('TEXT')
    } else {
      stop('Unrecognide column format')
    }
  })
  # Create command to create table
  tableCommand <- paste('CREATE TABLE', name, '(')
  for(i in 1:length(colnames(df))) {
    tableCommand <- paste0(
      tableCommand,
      colnames(df)[i],
      types[i]
    )
  }
  
  tableCommand <- c(
    tableCommand,
    pasteO(colnames(counts), )
  )
}

apply(df, 2, function(z) {
  if (is.numeric(z)) {
    return('REAL')
  } else if (is.character(z)) {
    return('TEXT')
  } else {
    stop('Unrecognide column format')
  }
})



createDESeq2SQLite('Test9.sqlite', ddsAnalysis)
  
  setwd('/camp/stp/babs/working/rabinoa/projects/stockingerg/amina.metidji/alignments/pairedEnd/')
counts <- read.table('combined_gene.exp.txt', header=T, sep='\t', row.names=1)
results <- read.table('')
db <- dbConnect(SQLite(), dbname='Test.sqlite')
dbWriteTable(db, 'counts', counts)
