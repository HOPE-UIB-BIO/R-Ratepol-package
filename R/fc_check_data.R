fc_check_data <- function (data.source.check, proportion = F, Debug=F, Species=T, Samples=T)
{
  # check if there is a sample that do not have a pollen data and delete it
  # & 
  # check if there are any specie without pollen record and delete them
  
  kill.all <- function(data.source.check.kill, Species= Species, Samples=Samples)
  {
    if(any(rowSums(data.source.check.kill$Pollen, na.rm = T)==0) & Samples == T) # if there are some samples without pollen
    {
      data.source.check.kill$Age <- data.source.check.kill$Age[rowSums(data.source.check.kill$Pollen, na.rm = T)>0,]
      data.source.check.kill$Age.un <- data.source.check.kill$Age.un[,rowSums(data.source.check.kill$Pollen, na.rm = T)>0]
      data.source.check.kill$Pollen <- data.source.check.kill$Pollen[rowSums(data.source.check.kill$Pollen, na.rm = T)>0,]
    }
    
    if(any(colSums(data.source.check.kill$Pollen, na.rm=T)==0) & Species == T) # if there are some species without pollen
    {
      data.source.check.kill$Pollen<- data.source.check.kill$Pollen[,colSums(data.source.check.kill$Pollen, na.rm = T)>0]
    }
    
    # count the species and sampels
    data.source.check.kill$Dim.val[1] <- ncol(data.source.check.kill$Pollen)
    data.source.check.kill$Dim.val[2] <- nrow(data.source.check.kill$Pollen)
    data.source.check.kill$Dim.val[3] <- nrow(data.source.check.kill$Age)
    
    return(data.source.check.kill)
  }
  
  data.source.check <- kill.all(data.source.check, Species = Species, Samples = Samples) 
  
  if(Debug==T)
  {
    cat("", fill=T)
    cat(paste("Pollen data have",data.source.check$Dim.val[1],"species with pollen record and",
              data.source.check$Dim.val[2],"samples. Age data have",data.source.check$Dim.val[3],"samples"),fill=T)
    cat("", fill=T)
    cat(paste("Age data has values of min",min(data.source.check$Age$age),", max",max(data.source.check$Age$age),",mean",
                round(mean(data.source.check$Age$age),2),",and median",round(median(data.source.check$Age$age),2)), fill=T)
    cat("", fill=T)
    
  }
  
  # check if all values is new age are in positive values and interpolate if necesery
  if(any(data.source.check$Age$newage<0, na.rm = T))
  {
    data.source.check$Age$newage <- data.source.check$Age$newage + min(data.source.check$Age$newage)*(-1)
  }
  
  if (proportion == T)
  {
    if (Debug==T){ cat("POllen values converted to proportions", fill=T)}

    # convert the values pollen data to proportion of sum of each sample
    p.counts.row.sums <- apply(data.source.check$Pollen, 1, sum)
    data.source.check$Pollen <- as.data.frame(lapply(data.source.check$Pollen, function(x) x/p.counts.row.sums))
    data.source.check <- kill.all(data.source.check, Species = Species, Samples = Samples)
  }
  
  return(data.source.check)
}