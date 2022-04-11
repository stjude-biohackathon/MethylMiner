##########################################################################################################################################################
### Created by P. Garg (March 06, 2019)
##########################################################################################################################################################
# findDMR function
findDMR<-function(sign=NULL, chrCol="CpG_chrm", startCol = "CpG_beg", dmrSep = 1000){
  
  #Calculate pairwise distance of consecutive probes
  s = sign[,startCol]
  s2 = s[-1] - s[1:(length(s)-1)]
  
  #Check chromosome of consecutive probes
  c = sign[,chrCol]
  c2 = c[-1] == c[1:(length(c)-1)]
  
  #Consider seperate DMRs when two DMRs are locaated >dmrSep apart and are on different chromosomes
  dmrEnd = c(which(s2>dmrSep | c2 == F), nrow(sign))
  dmrStart = c(1,dmrEnd[1:(length(dmrEnd)-1)]+1)
  return(list(dmrStart = dmrStart,dmrEnd=dmrEnd))
}


##########################################################################################################################################################
# getDMRlist function

args=(commandArgs(TRUE))
if(length(args)==0)
{
    print("No arguments supplied. Please try again.")
    ##supply default values
   quit("no")
}else{
        args=c(args)
}

print(paste0("Input file: ", args[1]))
# print(paste0("Plot list file: ", args[2]))

input_file = args[1]
qCutMin = paste0("Quan",args[2])
qCutMax = paste0("Quan",args[3])

# if(length(args)>1){ 
#   qc = read.table(args[2], header = T, as.is = T, sep = "\t", check.names = F)
#   rownames(qc) = qc[,1]
#   if(class(qc$predictedGender)=="logical"){ qc$predictedGender[!qc$predictedGender] = "F"}
# }
args

print(paste0("Reading input file: ", input_file, " ..."))
x = as.data.frame(read.table(input_file, sep = "\t", as.is = T, header = T, check.names = F))
head(x)

###Editted by Paras Mar 15,2019
if(nrow(x) ==0){return(0); }

start = findDMR(x)$dmrStart
end = findDMR(x)$dmrEnd

dir = getwd()
cohort = rev(unlist(strsplit(dir, "/")))[2]
length(start)
dmrlist = lapply(1:length(start), function(i){ 
  dmr = x[start[i]:end[i],]
  ####added by Paras on Jan 21, 2019
  ind = unique(unlist(strsplit(dmr[,"Sign_individuals_t0.1_n3_w1k_BOTH"],",")))
#  individual = paste(unique(unlist(strsplit(dmr[,5],","))), collapse = ",")

  d = lapply(ind, function(s){ 
    hyper = which(dmr[,s] >= dmr[,qCutMin]+0.15)
    hypo = which(dmr[,s] <= dmr[,qCutMax]-0.15)

    o1 = c(); o2 = c()
    if(length(hyper) >=3){
      o1 = data.frame(Chr_DMR = unique(dmr[hyper,"CpG_chrm"]), Start_DMR = min(dmr[hyper,"CpG_beg"]), End_DMR = max(dmr[hyper,"CpG_end"]),
                 Cohort=cohort,EpiInd = s, direction = "HYPER")
    }
    if(length(hypo) >=3){
      o2 = data.frame(Chr_DMR = unique(dmr[hypo,"CpG_chrm"]), Start_DMR = min(dmr[hypo,"CpG_beg"]), End_DMR = max(dmr[hypo,"CpG_end"]),
                Cohort=cohort,EpiInd = s, direction = "HYPO")
    }
    rbind(o1, o2)
  })
  do.call(rbind, d)
})

dmrlist = do.call(rbind, dmrlist)
dim(dmrlist)
#dmrlist = do.call(rbind, dmrlist)
write.table(dmrlist, sub(".sig", ".sig_dmr.txt", input_file), row.names =F, sep = "\t", quote =F)
head(warnings())




