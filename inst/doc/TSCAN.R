## -----------------------------------------------------------------------------
library(TSCAN)
data(lpsdata)
procdata <- preprocess(lpsdata)

## -----------------------------------------------------------------------------
lpsmclust <- exprmclust(procdata)

## -----------------------------------------------------------------------------
plotmclust(lpsmclust)

## -----------------------------------------------------------------------------
lpsorder <- TSCANorder(lpsmclust)
lpsorder

## -----------------------------------------------------------------------------
diffval <- difftest(procdata,lpsorder)
#Selected differentially expressed genes under qvlue cutoff of 0.05
head(row.names(diffval)[diffval$qval < 0.05])

## -----------------------------------------------------------------------------
STAT2expr <- log2(lpsdata["STAT2",]+1)
singlegeneplot(STAT2expr, TSCANorder(lpsmclust,flip=TRUE,orderonly=FALSE))

## -----------------------------------------------------------------------------
subpopulation <- data.frame(cell = colnames(procdata), sub = ifelse(grepl("Unstimulated",colnames(procdata)),0,1), stringsAsFactors = FALSE)
#Comparing ordering with or without marker gene information
order1 <- TSCANorder(lpsmclust)
order2 <- TSCANorder(lpsmclust, c(1,2,3))
orders <- list(order1,order2)

## -----------------------------------------------------------------------------
orderscore(subpopulation, orders)

