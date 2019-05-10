args = commandArgs(trailingOnly=TRUE)

require(gtools)

PDB_FILENAME = args
PDB = readLines(PDB_FILENAME)
ATOMS.inx = sapply(PDB, FUN=function(X) substr(X,1,4) == "ATOM")

FIND_PAIR_Command = paste("find_pair",PDB_FILENAME)
FP_Results = system(FIND_PAIR_Command,intern = T)

FP_Results_Base_Pairs = FP_Results[grep("DA]|DC]|DG]|DT]",FP_Results)]

print(FP_Results_Base_Pairs)



Base_Pairs = data.frame(
CHAIN.A = unlist(sapply(as.character(FP_Results_Base_Pairs),FUN=function(X) substr(X,31,31))),
RES.A = as.character(sapply(as.character(FP_Results_Base_Pairs),FUN=function(X) substr(X,41,42))),
RESN.A = gsub("\\.","",sapply(as.character(FP_Results_Base_Pairs),FUN=function(X) substr(X,33,36))),
CHAIN.B = unlist(sapply(as.character(FP_Results_Base_Pairs),FUN=function(X) substr(X,63,63))),
RES.B = unlist(sapply(as.character(FP_Results_Base_Pairs),FUN=function(X) substr(X,53,54))),
RESN.B = gsub("\\.","",sapply(as.character(FP_Results_Base_Pairs),FUN=function(X) substr(X,57,60)))
)
rownames(Base_Pairs) = rep(NULL,nrow(Base_Pairs))
print(Base_Pairs)

Pent.Start.Stop.inx = c()
for(B in 1:(nrow(Base_Pairs)-4))
{
  Pent.Start.Stop.inx = rbind(Pent.Start.Stop.inx, c(B:(B+4)))
}

print(Pent.Start.Stop.inx)

system(paste("mkdir ",PDB_FILENAME,"_Models/",sep=""))

apply(Pent.Start.Stop.inx,MARGIN=1, function(X){
  
  BASES.Pick = apply(Base_Pairs[X,],MARGIN=1, function(X) c(paste(X[c(2,1,3)],collapse=" "),paste(X[c(5,4,6)],collapse=" ")))
  GREP_Base_qry = paste(BASES.Pick,collapse ="|")
  
  NEW_PDB = c(
    ## Protein -> Residue name has 3 characters
    PDB[ATOMS.inx] [nchar(gsub(" ","",substr(PDB[ATOMS.inx] ,18,20)))==3],
    ## Selected Pentamer Base Pairs
    PDB[ATOMS.inx][grep(GREP_Base_qry,gsub("\\s+"," ",PDB[ATOMS.inx]))]
  )
  
  write(NEW_PDB,paste(PDB_FILENAME,"_Models/",X[1],".pdb",sep=""))
  

  
})#apply(Pent.Start.Stop.inx,MARGIN=1, function(X){


Base_Combinations = permutations(n=4,v = c("DA","DC","DG","DT"),repeats.allowed=T,r=5)
Base_Compliments_Combinations = Base_Combinations[nrow(Base_Combinations):1,]


apply(Pent.Start.Stop.inx,MARGIN=1, function(X){
  for(K in 1:nrow(Base_Combinations))
  {

    Base_Mutations = paste(c(paste(
      paste("c=",Base_Pairs$CHAIN.A[X],sep=""),
      paste("s=",Base_Pairs$RESN.A[X],sep=""),
      paste("m=",Base_Combinations[K,],sep="")
      ),
      paste(
      paste("c=",Base_Pairs$CHAIN.B[X],sep=""),
      paste("s=",Base_Pairs$RESN.B[X],sep=""),
      paste("m=",Base_Compliments_Combinations[K,],sep="")
    )),collapse = ";") #Base_Mutations = paste(c(paste(


    Seq.EXT = paste(gsub("D","",Base_Combinations[K,]),collapse="")
    system_Command = paste("mutate_bases \"",Base_Mutations,"\" ",PDB_FILENAME,"_Models/",X[1],".pdb ",PDB_FILENAME,"_Models/",X[1],"_",Seq.EXT,".pdb",sep="")
    #print(system_Command)
    system(system_Command)
  }#for(K in 1:nrow(Base_Combinations))
})#apply(Pent.Start.Stop.inx,MARGIN=1, function(X){


