
## Read in delimited file with sequence identifiers and trait values
traitDataFile <- tcltk::tk_choose.files(default = "", caption = "Please select one trait data file", multi=FALSE)
traitData <- read.delim(traitDataFile, header=TRUE, row.names=1, na.strings = "-999") ## -999 is TASSEL output

## Read in delimited files with sequence names and rpkm values
rpkmA.temp <- read.table ("RPKM_AT2018.txt", header=TRUE, sep= "\t", row.names=1) # can use the check.names=FALSE argument if there's "unusual" symbols in the column headers  2018-01-15

outputName=""
outputName=readline(prompt="Please enter output identifier...   ")
print("reading files...")

traitLabel = colnames(traitData)

meansCutOff = 0.4

## remove low rpkm means
print(paste0("Removing markers with a mean RPKM < ", meansCutOff))
delrpkmA = rpkmA.temp[rowMeans(rpkmA.temp) >= meansCutOff, ]

rpkmA <- t(delrpkmA)

## Delete missing varieties from rpkm files
rpkmmergeA <- merge(traitData,rpkmA, by="row.names")
rownames(rpkmmergeA) <- rpkmmergeA[,1]

## Function for apply()ing. Does a linear model between trait data and genes
performLinearRegression <- function(exp_vector) {
  model <- lm(rpkmmergeA[,2] ~ exp_vector)
  return(as.numeric(anova(model)[1,]))
}

## Do the linear models. Each gene is in a column hence apply() over
## columns. rpkmmergeA contains lines in column 1 and the trait data
## in column 2 hence these are dropped in the apply() call. The beauty
## of this way is that it makes the row names the genes in the final
## lmResults output.
lmResults <- t(apply(rpkmmergeA[, -c(1,2)], 2, performLinearRegression))
colnames(lmResults) <- c("Df", "SumSq", "MeanSq", "Fvalue", "Pvalue")

## calculate log10P
log10P <- -log10(lmResults[,"Pvalue"])
print("Calculated p-values")

## Read B. napus to Arabidopsis matches 
codesFile = read.csv("Marker to At_AT2018.csv", header=TRUE, na.strings = ".")

## Get everything together for writing to a file. Made an explicit
## rownames column so don't have to use row.names in next regress
## plotter step. Not a problem because I use row.names=FALSE in
## write.table. Note needed to use as.data.frame otherwise it all came
## out as a list of character vectors i.e. in quotes which obviously
## can't be order()ed.
getGeneModels = function(results) {
  unigenes = unlist(strsplit(sub("_", "~#~", results$unigene), "~#~"))
  ArabidopsisHits = codesFile[ codesFile$unigene %in% unigenes, ]
  results$tmp = unigenes
  resultsWithAGI = merge(results, ArabidopsisHits, by.x="tmp", by.y="unigene")  ## likely in the future can just merge by="unigene" as identical(as.character(results$unigene), results$tmp) == TRUE
  return(resultsWithAGI)
}

BnapusUnigene = rownames(lmResults)
traitAndUnigene = as.data.frame(cbind(trait = traitLabel, unigene = BnapusUnigene))

print("Producing final results")

finalResults <- cbind(traitAndUnigene, log10P, lmResults)
finalResults = getGeneModels(finalResults)
finalResults <- finalResults[order(finalResults$log10P, decreasing = TRUE),]

## Write these columns for now. In future if the apply() is quick
## enough probably combine this script and grapher one together.
finalResults = finalResults[,c("trait", "unigene",  "A.Chr", "C.Chr", "AGI", "log10P", "Df", "SumSq", "MeanSq", "Fvalue", "Pvalue")]
resultsFile = paste("V6All_reg-", paste(outputName, collapse="-"), ".txt", sep="")
write.table(finalResults, resultsFile, quote=FALSE, sep="," ,row.names=FALSE, col.names=TRUE)
## ? Example output in console
##     trait      unigene      A.Chr C.Chr    AGI     log10P   Df  SumSq
##  minusefvern A_JCVI_14001    A3    C2 AT5G59050.1 14.28305  1 6632.848
##  minusefvern A_JCVI_40108    A2    C2 AT1G65480.1 13.99106  1 6532.559
##    MeanSq   Fvalue     Pvalue
##  6632.848 88.62775 5.211345e-15
##  6532.559 85.99291 1.020799e-14


print ("Phew! Complete...")







