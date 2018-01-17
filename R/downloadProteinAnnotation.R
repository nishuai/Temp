
downloadProteinAnnotation <- function(organism = "HUMAN",
                                      maindatadir = ".",
                                      ProtFeaturesFile = "ProtFeatures.rda",
                                      RBDmapFTfile = "RBDmapFT.rda",
                                      UniprotVersionFile = "UniprotVersion.rda") {

dir.create(maindatadir,showWarnings=FALSE)
tmpdir = file.path(maindatadir,"tmp")
dir.create(tmpdir,showWarnings=FALSE)

#if(missing(idmappingfile)) {
  idmappingfile = switch(organism,
    "ARATH" = "ARATH_3702_idmapping.dat",
    "CAEEL" = "CAEEL_6239_idmapping.dat",
    "CHICK" = "CHICK_9031_idmapping.dat",
    "DANRE" = "DANRE_7955_idmapping.dat",
#    "DICDI" = "DICDI_44689_idmapping.dat",
    "DROME" = "DROME_7227_idmapping.dat",
#    "ECOLI" = "ECOLI_83333_idmapping.dat",
    "HUMAN" = "HUMAN_9606_idmapping.dat",
    "MOUSE" = "MOUSE_10090_idmapping.dat",
    "RAT"   = "RAT_10116_idmapping.dat",
#    "SCHPO" = "SCHPO_284812_idmapping.dat",
    "YEAST" = "YEAST_559292_idmapping.dat",
    default = NA)
## }
## if (missing(fastafile)) {
  fastafile = switch(organism,
    "ARATH"   = "ARATH.fasta",
#    "BOVIN"   = "BOVIN.fasta",
    "CAEEL"   = "CAEEL.fasta",
#    "CANFA"   = "CANFA.fasta",
    "CHICK"   = "CHICK.fasta",
    "DANRE"   = "DANRE.fasta",
    "DROME"   = "DROME.fasta",
    "HUMAN"   = "HUMAN.fasta",
    "MOUSE"   = "MOUSE.fasta",
#    "PIG"   = "PIG.fasta",
    "RAT"   = "RAT.fasta",
    "YEAST"   = "YEAST.fasta",
    default = NA)
## }
ensemblDataset = switch(organism,
    "ARATH"   = "athaliana_gene_ensembl", # needs review
    "CAEEL"   = "celegans_gene_ensembl",
    "CHICK"   = "ggallus_gene_ensembl",
    "DANRE"   = "drerio_gene_ensembl",
    "DROME"   = "dmelanogaster_gene_ensembl",
    "HUMAN"   = "hsapiens_gene_ensembl",
    "MOUSE"   = "mmusculus_gene_ensembl",
    "RAT"   = "rnorvegicus_gene_ensembl",
    "YEAST"   = "scerevisiae_gene_ensembl",
    default = NA)

if (any(is.na(c(idmappingfile,fastafile)))) {
  stop('organism not supported. Try one of ARATH, CAEEL, CHICK, DANRE, DROME, HUMAN, RAT, YEAST')
}

####################################
## check for command line tools
####################################

errmsg = c()
res=system("ftp --help",intern=FALSE,ignore.stdout=TRUE,ignore.stderr=TRUE)
if (res > 1) {
  errmsg = c(errmsg, "command line ftp not found.\nPlease install a command line ftp before running downloadProteinAnnotation.\n")
}
res=system("gunzip --help",intern=FALSE,ignore.stdout=TRUE,ignore.stderr=TRUE)
if (res > 1) {
  errmsg = c(errmsg, "command line gunzip not found.\nPlease install a command line gunzip before running downloadProteinAnnotation.\n")
}
res=system("iupred --help",intern=FALSE,ignore.stdout=TRUE,ignore.stderr=TRUE)
if (res > 1) {
  errmsg = c(errmsg, "command line iupred not found.\nPlease install a command line iupred before running downloadProteinAnnotation.\n")
}
res=system("perl --help",intern=FALSE,ignore.stdout=TRUE,ignore.stderr=TRUE)
if (res > 1) {
  errmsg = c(errmsg, "command line perl not found.\nPlease install a command line perl before running downloadProteinAnnotation.\n")
}
if (length(errmsg) > 0) {
  stop(errmsg)
}

####################################
## download swissprot protein sequences and gene mapping
####################################

cat("Download idmapping file ", idmappingfile, "\n")
system(sprintf("ftp ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/%s.gz",idmappingfile))
system(sprintf("gunzip %s.gz",idmappingfile))
file.rename(idmappingfile, sprintf("%s/%s",tmpdir, idmappingfile))

cat("Download fasta file ", fastafile, "\n")
system(sprintf("ftp ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/proteomes/%s.gz",fastafile))
system(sprintf("gunzip %s.gz",fastafile))
file.rename(fastafile, sprintf("%s/%s",tmpdir, fastafile))

cat("Download release notes\n")
system("ftp ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/proteomes/relnotes.txt")
UniprotVersion = list(Uniprot = readLines("relnotes.txt"))
file.remove("relnotes.txt")

###################################################
### Read swissprot fasta file
###################################################
cat("read fasta file ",sprintf("%s/%s",tmpdir, fastafile),"\n")
ProtSeq = readAAStringSet(sprintf("%s/%s",tmpdir, fastafile))
desc = sapply(strsplit(names(ProtSeq), split="|",fixed=TRUE), function(x) { x[2] } )
names(ProtSeq) = desc
K = grepl("-", names(ProtSeq), fixed=TRUE)
ProtSeq = ProtSeq[which(!K)]
# ProtSeq = sapply(Fasta, function(x) { x$seq } )
# g = regexpr("[*]", ProtSeq)
# ProtSeq[g>0] = substr(ProtSeq[g>0],1,g[g>0]-1)
save(ProtSeq, desc, file=file.path(tmpdir,"ProtSeq.rda"))

###################################################
### Read swissprot IDmapping file
###################################################
cat("read IDmapping file ",sprintf("%s/%s",tmpdir, idmappingfile),"\n")
IDmapping = read.table(sprintf("%s/%s",tmpdir, idmappingfile), sep="\t", header=FALSE,
                       comment.char="",quote="",stringsAsFactors=FALSE)
if (organism == "YEAST") {
  IDmapping = IDmapping[IDmapping[[2]] == "CYGD",]
  IDmapping[[3]] = toupper(IDmapping[[3]])
} else {
  if (organism == "DROME") {
    IDmapping = IDmapping[IDmapping[[2]] == "FlyBase",]
  } else {
    IDmapping = IDmapping[IDmapping[[2]] == "Ensembl",]
  }
}
GeneName = tapply(IDmapping[[3]], factor(IDmapping[[1]],levels=names(ProtSeq)),
       function(x) { sort(x)[1] } )

###################################################
### Write one fasta file per protein and run iupred
###################################################
cat("\n")
for (i in 1:length(ProtSeq)) {
  cat("\rwrite fasta files seq ",i," out of ",length(ProtSeq))
  writeXStringSet(x=ProtSeq[i], filepath=sprintf("%s/tmp%d.seq",tmpdir,i))
}
cat("\n")


iupred = list()
for (i in 1:length(ProtSeq)) {
  cat("\riupred seq ",i," out of ",length(ProtSeq))
  X = system(sprintf("iupred %s/tmp%d.seq long",tmpdir,i), intern=TRUE)
  X = X[10:length(X)]
  iupred[[i]] = sapply(strsplit(X,split=" "),function(x) { as.numeric(x[length(x)]) } )
}
names(iupred) = names(ProtSeq)
save(iupred, file=file.path(tmpdir,"iupred.rda"))
cat("\n")
for (i in 1:length(ProtSeq)) {
  file.remove(sprintf("%s/tmp%d.seq",tmpdir,i))
}

###################################################
### Download Interpro domain mapping
###################################################
QueryXML = readLines(system.file("Biomart/Query.xml", package="RBDmap"))
file.copy(system.file("Biomart/webQuery.pl", package="RBDmap"), "webQuery.pl")

J = grep("PROTID",QueryXML)
cat("\n")
L = list()
Interpro = list()
Size = 500
n = names(ProtSeq)
N = ceiling(length(n) / Size)
cat("\n")
for (i in 1:N) {
  I = 1:Size + (i-1)*Size
  I = I[I %in% seq_along(n)]
  QXML = QueryXML  
  QXML[J] = sub(pattern="PROTID", paste(n[I],collapse=","), QueryXML[J])
  writeLines(QXML, "Query.xml")
  L[[i]] = system("perl webQuery.pl Query.xml",intern=TRUE)
  cat("\nDownload Interpro. Query ",i, " from ",N," length=",length(L[[i]]),"\n")
  while (length(L[[i]]) < 1) {
    cat("Reload data")
    L[[i]] = system("perl webQuery.pl Query.xml",intern=TRUE)
    cat("\nDownload Interpro. Query ",i, " from ",N," length=",length(L[[i]]),"\n")
  }
  Interpro[[i]] = read.table(text=L[[i]], sep="\t", quote="", 
                             comment.char="", stringsAsFactors=FALSE)
#  print(head(L[[i]]))
}
cat("\n")
Interpro = do.call(what=rbind, args=Interpro)
colnames(Interpro) = c(
  "proteinAccession","crc64","md5","proteinLength","methodId",
  "methodName","methodDatabaseName","posFrom","posTo",
  "matchStatus","matchScore")
K = as.integer(factor(Interpro$proteinAccession, levels=names(ProtSeq))) * 
      (max(Interpro$posTo)+1) + Interpro$posFrom
Interpro = Interpro[order(K),]
save(Interpro, file=file.path(tmpdir,"Interpro.rda"))
file.remove("Query.xml")
file.remove("webQuery.pl")

cat("\nDownloading Interpro release notes\n")
system("ftp ftp://ftp.ebi.ac.uk/pub/databases/interpro/release_notes.txt")
UniprotVersion$Interpro = readLines("release_notes.txt")
file.remove("release_notes.txt")
save(UniprotVersion, file=file.path(maindatadir,UniprotVersionFile))


###################################################
### Download Uniprot annotation of protein features
###################################################
cat("\n")
L = list()
PID = list()
UniprotFT = list()
Size = 200
n = names(ProtSeq)
N = ceiling(length(n) / Size)
for (i in 1:N) {
  I = 1:Size + (i-1)*Size
  I = I[I %in% seq_along(n)]
  L[[i]] = system(paste('wget -O- "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=uniprotkb&id=',paste(n[I],collapse="%2C"),'&format=default&style=raw&Retrieve=Retrieve"',sep=""),intern=TRUE,ignore.stderr=TRUE)
  cat("Downloading Uniprot. Query ",i, " from ",N," length=",length(L[[i]]),"\n")
  while (length(L[[i]]) < 1) {
    cat("Reload data")
    L[[i]] = system(paste('wget -O- "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=uniprotkb&id=',paste(n[I],collapse="%2C"),'&format=default&style=raw&Retrieve=Retrieve"',sep=""),intern=TRUE,ignore.stderr=TRUE)
    cat("Downloading Uniprot. Query ",i, " from ",N," length=",length(L[[i]]),"\n")
  }
  L[[i]] = c(L[[i]][substr(L[[i]],1,2) %in% c("AC","CC","FT")],"AC")

  K = which(substr(L[[i]],1,2) == "AC")
  Y = sapply(n[I], function(x) { grep(x, L[[i]][K])[1] })
  if (any(is.na(Y))) {
    cat(sum(is.na(Y))," proteins not found in UniprotKB in batch ",i,"!\n")
  }
  Z = lapply(Y[!is.na(Y)], function(x) { if (K[x+1] - K[x] > 1) { (K[x]+1):(K[x+1]-1) } else { c() } })
  PID[[i]] = rep(n[I[!is.na(Y)]], listLen(Z))
  L[[i]] = L[[i]][unlist(Z)]
}
cat("\n")
PID = do.call(what=c, args=PID)
L = do.call(what=c, args=L)
Caution = PID[which(grepl("CC   -!- CAUTION: The sequence shown here is derived from an Ensembl",L))]

I = which(substr(L,1,2) == "FT")
PID = PID[I]
L = L[I]
humsavar = NULL
if (organism == "HUMAN") {
  humsavar = system('wget -O- http://www.uniprot.org/docs/humsavar.txt',intern=TRUE,ignore.stderr=TRUE)
}
save(L, PID, humsavar,file="Ltmp.rda")
FTuniprot = parseUniprotFT(L, PID, humsavar)
cat("Downloading Uniprot. Query ",i, " from ",N," length=",length(L[[i]]),"\n")

save(FTuniprot, desc, file=file.path(tmpdir,"FTuniprot.rda"))
save(Caution, desc, file=file.path(tmpdir,"Caution.rda"))

###################################################
### Download ELM
###################################################

# X = system("wget http://elm.eu.org/elms/browse_instances.tsv?q=DOC_&taxon=homo+sapiens")
# X = system("wget http://elm.eu.org/elms/browse_instances.tsv?q=DOC_&taxon=mus+musculus")
# http://elm.eu.org/elms/browse_instances.tsv?q=P12931

###################################################
### Map iupred disorderdness features
###################################################
RR = lapply(iupred,function(x) {
  if (length(x) < 9) {
    start = c()
    end = c()
    res = matrix(0,nr=0,nc=2)
  } else {
    disordered = x > 0.4
    disordered = filter(disordered,rep(1,9)) >= 5
    disordered[1:4] = disordered[5]
    disordered[-4:0 + length(disordered)] = disordered[length(disordered)-5]
    I = which(disordered)
    start = I[which(!((I-1) %in% I))]
    end = I[which(!((I+1) %in% I))]
    res = cbind(start,end)
  }
  res
} )

Start = sapply(RR,function(x) { as.vector(x[,1]) })
Stop = sapply(RR,function(x) { as.vector(x[,2]) })
tmp = rep(names(RR),sapply(Start,length))
FTiupred = data.frame(UniProtID = tmp,
  Source = "iupred",Type="disorder",Name="disorder",
                      Start=unlist(Start),Stop=unlist(Stop),
                      ShortName="disorder",
                      Info="",
                      stringsAsFactors=FALSE)
save(FTiupred, desc, file=file.path(tmpdir,"FTiupred.rda"))


###################################################
### Map Interpro features
###################################################
FTinterpro = Interpro
FTinterpro = cbind(FTinterpro, Source="Interpro",Info="",stringsAsFactors=FALSE)
FTinterpro = FTinterpro[,c("proteinAccession","Source","methodDatabaseName","methodName","posFrom","posTo","methodName","Info")]
colnames(FTinterpro)  = c("UniProtID","Source","Type","Name","Start","Stop","ShortName","Info")
save(FTinterpro, desc, file=file.path(tmpdir,"FTinterpro.rda"))


###################################################
### Combine Uniprot, Interpro, and disorderdness protein features for RBDmap
###################################################
RBDmapFT = rbind(FTinterpro,FTuniprot,FTiupred)

tmp = sprintf("%s:::%s:::%s:::%0.7d",RBDmapFT[,1],RBDmapFT[,2],RBDmapFT[,3],RBDmapFT[,5])
RBDmapFT = RBDmapFT[order(tmp),]
row.names(RBDmapFT) = 1:nrow(RBDmapFT)

I = which(!duplicated(RBDmapFT[,1]))
J = which(!duplicated(RBDmapFT[,1],fromLast=TRUE))
Index = data.frame(RBDmapFT[I,1],From=I,To=J)
attr(RBDmapFT,which="Index") = Index
save(RBDmapFT, file=file.path(maindatadir,RBDmapFTfile))

###################################################
### Download gene names
###################################################
ensembl = useMart("ensembl", dataset = ensemblDataset)
ENSG2name <- getBM(attributes = c("ensembl_gene_id", "external_gene_id"), 
                   mart = ensembl)
ENSG2name <- ENSG2name[nchar(ENSG2name$external_gene_id) > 0,]
ENSG2name <- unique(ENSG2name)

###################################################
### List with protein features for mRNAinteractome
###################################################
ProtFeatures = list()
ProtFeatures$names = names(ProtSeq)
ProtFeatures$ProtSeq = ProtSeq
ProtFeatures$GeneName = GeneName
ProtFeatures$Symbol = paste("sp:",ProtFeatures$names,sep="")
m = match(ProtFeatures$GeneName, ENSG2name$ensembl_gene_id)
ProtFeatures$Symbol[!is.na(m)] = ENSG2name$external_gene_id[m[!is.na(m)]]
names(ProtFeatures$Symbol) = ProtFeatures$names
ProtFeatures$iupred = iupred
IP = Interpro[ Interpro$methodDatabaseName %in% c("Pfam","PfamB"), ]
ProtFeatures$domainDesc = tapply(IP$methodName,
                                 factor(IP$proteinAccession, levels=names(ProtSeq)),
                                 function(x) { x } )
ProtFeatures$domainStart= tapply(IP$posFrom,
                                 factor(IP$proteinAccession, levels=names(ProtSeq)),
                                 function(x) { x } )
ProtFeatures$domainEnd = tapply(IP$posTo,
                                 factor(IP$proteinAccession, levels=names(ProtSeq)),
                                 function(x) { x } )
ProtFeatures$Caution = ProtFeatures$names %in% Caution
names(ProtFeatures$Caution) = ProtFeatures$names
ProtFeatures$H = getComplexity(ProtSeq)
save(ProtFeatures, file=file.path(maindatadir,ProtFeaturesFile))

Index = toKmers(ProtFeatures$ProtSeq, verbose=TRUE)
save(Index, file=file.path(maindatadir,"Index.rda"))

invisible(NULL)
}

# http://www.uniprot.org/docs/humpvar.txt
# http://www.uniprot.org/docs/humsavar.txt
