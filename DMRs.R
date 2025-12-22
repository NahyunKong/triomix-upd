suppressMessages(library(DSS))
suppressMessages(require(bsseq))
suppressMessages(library(optparse))

option_list <- list(
  make_option(c("-a", "--hap1"), type="character", help="hap1 bedfile"),
  make_option(c("-b", "--hap2"), type="character", help="hap2 bedfile"),
  make_option(c("-o", "--output"), type="character", help="outfile"),
  make_option(c("-t", "--threshold"), type="numeric", help="Threshold difference to be considered DMR")
 # make_option(c("-h", "--help"), action="help", help="Show this help message and exit")
)

# Parse the options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#generate statistics (only first time)
{
  #input data
  dat1.1 = read.table(file.path(opt$hap1), header=F) #hap1
  dat2.1 = read.table(file.path(opt$hap2), header=F) #hap2
  colnames(dat1.1) <- c("chr", "pos", "N" , "X") #N: total reads, X: methylated reads
  colnames(dat2.1) <- c("chr", "pos", "N" , "X")
  
  #create BSobj
  BSobj = makeBSseqData( list(dat1.1,dat2.1),
                         c("Hap1","Hap2") ) #[1:1119212,] #total number of CpGs
  #statistical test (taks a while)
  dmlTest.sm = DMLtest(BSobj, group1=c("Hap1"), group2=c("Hap2"), 
                       smoothing=TRUE) #smoothing.span: size of smoothing window in bp. Default 500
  #ncores: Number of CPUcores in parallel computing
}

#differential loci
dmls = callDML(dmlTest.sm, p.threshold=0.001, delta=opt$threshold) #delta: define threshold of differential methylation
head(dmls)
#select region of interest
#result <- dmls[dmls$chr == "chr16" & dmls$pos >= 29638676 & dmls$pos <= 30188531, ]

#differential region
dmrs = callDMR(dmlTest.sm, p.threshold=0.01, delta=opt$threshold) #p.threshold: loci with p-values< be picked and joint to form DMRs
#minlen: Minimum length for DMR, Default 50bp
#minCG: Minimum CpG sites for DMR, Default 3
#dis.merge: When two DMRs are very close, distance< be merged into one, Default 50bp
#pct.sig: In DMR, percentage of CGs with significant p-values must >, Default 0.5. mainly correct merging of nearby DMRs
head(dmrs)

#save results
BSobj_bc2027 <- BSobj
dmls_bc2027 <- dmls
dmlTest.sm_bc2027 <- dmlTest.sm
dmrs_bc2027 <- dmrs

write.table(dmrs, opt$output, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


