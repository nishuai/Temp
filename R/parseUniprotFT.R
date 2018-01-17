
parseUniprotFT <- function(Uniprot, PID, humsavar=NULL) {
  PF = data.frame(UniProtID = PID,
                  Type = gsub(" ","", substr(Uniprot,6,13), fixed=TRUE),
                  Start = gsub(" ","",substr(Uniprot,14,20), fixed=TRUE),
                  Stop = gsub(" ","",substr(Uniprot,21,27), fixed=TRUE),
                  Name = substr(Uniprot,28,10000), Source="Uniprot", stringsAsFactors=FALSE)
  I = which((nchar(PF$Start) == 0) & (nchar(PF$Stop) == 0))
  I = I[!((I-1) %in% I)]
  while (length(I) > 0) {
    PF$Name[I-1] = paste(PF$Name[I-1], PF$Name[I],sep="")
    PF = PF[-I,]
    I = which((nchar(PF$Start) == 0) & (nchar(PF$Stop) == 0))
    I = I[!((I-1) %in% I)]
  }
  PF = PF[which(!(grepl("[?<>]", PF$Start) | grepl("[?<>]", PF$Stop))),]
  PF$Start = as.integer(PF$Start)
  PF$Stop = as.integer(PF$Stop)
  PF = PF[,c(1,6,2,5,3,4)]

  n = sapply(gregexpr("./FTId=", PF$Name), function(x) { x[1] })
  I = which((n > 0) & (PF$Type == "VARIANT"))
  PF$Info = ""
  PF$Info[I] = substr(PF$Name[I],n[I]+7,n[I]+16)
  I = which(n > 0)
  PF$Name[I] = substr(PF$Name[I],rep(1,nrow(PF)),n[I]-1)

  if(!is.null(humsavar)) {
    a = which(humsavar == "_________ __________ __________  ____________   ____________  ___________ _____________________")
    humsavar = humsavar[(a+1):length(humsavar)]
    humsavar = humsavar[1:(length(humsavar)-5)]
    var = substr(humsavar,22,31)
    dbSNP = substr(humsavar,63,73)
    varDbSNP = var[which(substr(dbSNP,1,1) != "-")]
    A = I[PF$Info[I] %in% varDbSNP]
    PFvariant = PF[A,,drop=FALSE]
    PFvariant$Info = ""
    
    df = data.frame(var = substr(humsavar,22,31),
                    Info = substr(humsavar,75,100000),stringsAsFactors=FALSE)
    df = df[df$Info != "-",]
    A = I[PF$Info[I] %in% df$var]
    PFvariantMIM = PF[A,]
    PFvariantMIM$var = PFvariantMIM$Info
    PFvariantMIM$Info = NULL
    PFvariantMIM = merge(PFvariantMIM, df)
    PFvariantMIM$var = NULL
    PFvariantMIM$Type = "VARIANTMIM"
    PF = PF[-I,]
    PF = rbind(PF, PFvariant, PFvariantMIM)
  }

  NUMBERS = as.character(0:9)
#  I=seq_len(nrow(PF))
  PF$Type = sapply(strsplit(PF$Type, split=""), function(x) { paste(x[!(x %in% NUMBERS)],collapse="") } )

  PF$Name = gsub("^\\s+|\\s+$","",PF$Name)
  PF$Name = gsub("\\s+"," ",PF$Name)
  PF$ShortName = PF$Name
  PF$ShortName = gsub(" (Probable)","",PF$ShortName,fixed=TRUE)
  PF$ShortName = gsub(" (Potential)","",PF$ShortName,fixed=TRUE)
  PF$ShortName = gsub(" (By similarity)","",PF$ShortName,fixed=TRUE)
  PF$ShortName = gsub(" (By       similarity)","",PF$ShortName,fixed=TRUE)
  PF$ShortName = gsub("Probable.","",PF$ShortName,fixed=TRUE)
  PF$ShortName = gsub("Potential.","",PF$ShortName,fixed=TRUE)
  PF$ShortName = gsub("By similarity.","",PF$ShortName,fixed=TRUE)
  g = gregexpr("[;]",text=PF$ShortName)
  p = sapply(g, function(x) { x[1] })
  p[p<0] = nchar(PF$ShortName)[p<0]+1
  PF$ShortName = substr(PF$ShortName,1,p-1)

  I = which(PF$Type %in% c("VARIANT","VARIANTMIM"))
  st = substr(PF$ShortName[I],7,10000)
  g = gregexpr(pattern="[. ]",st)
  p = sapply(g, function(x) { x[1] })
  PF$ShortName[I] = substr(PF$ShortName[I],1,5+p)

#  I = which(PF$Type == "VARIANT")
#   J = I[which(grepl(" (",PF$Name[I],fixed=TRUE) & !grepl("dbSNP",PF$Name[I],fixed=TRUE))]
#   PF$Type[J] = "VARIANTMIM"
#   Disease = rep("",nrow(PF))
#   g1 = sapply(gregexpr("[(]",PF$Name[J]), function(x){ x[1]})
#   g2 = sapply(gregexpr("[)]",PF$Name[J]), function(x){ x[length(x)]})
#   Disease[J] = substr(PF$Name[J],g1+1,g2-1)
#   g2 = sapply(gregexpr("[;]",Disease[J]), function(x){ x[1] })
#   Disease[J[g2>0]] = substr(Disease[J[g2>0]], 1, g2[g2>0]-1)
#   g1 = sapply(gregexpr("in ",Disease[J],fixed=TRUE), function(x){ x[1]})
#   Disease[J[g1>0]] = substr(Disease[J[g1>0]], g1[g1>0]+3,10000)
#   PF$Info = Disease

  g = gregexpr("[:]",text=PF$ShortName)
  p = sapply(g, function(x) { x[1] })
  n = nchar(PF$ShortName)
  K = which((p < n) & (p > 0) & (PF$Type == "MUTAGEN"))
  PF$ShortName[K] = substr(PF$ShortName[K],1,p[K]-1)

  g = gregexpr("[.]",text=PF$ShortName)
  p = sapply(g, function(x) { x[1] })
  n = nchar(PF$ShortName)
  K = which((p < n) & (p > 0) & !(PF$Type %in% c("CHAIN", "CONFLICT", "VAR_SEQ", "CARBOHYD","DOMAIN","DNA_BIND")))
  PF$ShortName[K] = substr(PF$ShortName[K],1,p[K]-1)

  g = gregexpr("[)]",text=PF$ShortName)
  p = sapply(g, function(x) { x[1] })
  n = nchar(PF$ShortName)
  K = which((p < n) & (p > 0) & (PF$Type == "CARBOHYD"))
  PF$ShortName[K] = substr(PF$ShortName[K],1,p[K])

  g = gregexpr(" (interchain",text=PF$ShortName,fixed=TRUE)
  p = sapply(g, function(x) { x[1] })
  n = nchar(PF$ShortName)
  K = which((p < n) & (p > 0) & (PF$Type == "CROSSLNK"))
  PF$ShortName[K] = substr(PF$ShortName[K],1,p[K]-1)

  g = gregexpr("[.]",text=PF$ShortName)
  p = sapply(g, function(x) { x[length(x)] })
  n = nchar(PF$ShortName)
  K = which((p == n) & (p > 0))
  PF$ShortName[K] = substr(PF$ShortName[K],1,p[K]-1)

  g = gregexpr("[0123456789 -]",text=PF$ShortName)
  p = sapply(g, function(x) { x[length(x)] })
  n = nchar(PF$ShortName)
  K = which((p == n) & (p > 0) & ( PF$Type == "REPEAT"))
  while(length(K) > 0) {
    PF$ShortName[K] = substr(PF$ShortName[K],1,p[K]-1)
    g = gregexpr("[0123456789 -]",text=PF$ShortName)
    p = sapply(g, function(x) { x[length(x)] })
    n = nchar(PF$ShortName)
    K = which((p == n) & (p > 0) & ( PF$Type == "REPEAT"))
  }

  PF$ShortName[ PF$Type == "CHAIN" ] = "CHAIN"
  PF$ShortName[ PF$Type == "CONFLICT" ] = "CONFLICT"
  PF$ShortName[ PF$Type == "VAR_SEQ" ] = "VAR_SEQ"
  PF$ShortName[ (PF$Source == "Uniprot") & (nchar(PF$ShortName) == 0) ] = PF$Type[ (PF$Source == "Uniprot") & (nchar(PF$ShortName) == 0) ]

  PF = PF[,c("UniProtID","Source","Type","Name","Start","Stop","ShortName","Info")]
  tmp = sprintf("%s:::%s:::%s:::%0.7d",PF[,1],PF[,2],PF[,3],PF[,5])
  PF = PF[order(tmp),]
  
  PF
}

