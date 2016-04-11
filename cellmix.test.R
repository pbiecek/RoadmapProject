options(error = browser())
# install biocLite if not already thereg
if( !require(BiocInstaller) ){
  # enable Bioconductor repositories
  # -> add Bioc-software
  setRepositories() 
  
  install.packages('BiocInstaller')
  library(BiocInstaller)
}

if( !require(CellMix) ){
  # or alternatively do: 
  # source('http://www.bioconductor.org/biocLite.R')
  
  # install (NB: this might ask you to update some of your packages)
  biocLite('CellMix', siteRepos = 'http://web.cbio.uct.ac.za/~renaud/CRAN', type='both')
  library(CellMix)
}

if ( !require(GEOquery) ){
  
  biocLite('GEOquery')
  
  library(GEOquery)
  
}


gse <- ExpressionMix("GSE19830", verbose = 2)


# extract data for mixed samples only
mix <- mixedSamples(gse)

# extract stored known signatures (= average pure profiles)
sig <- basis(mix)



# estimate proportions from pure signatures (default: no normalization)
res <- ged(mix, sig, verbose = TRUE, log = FALSE)
## Mapping signature ids onto target ids (method: auto) ... SKIP [signatures have no annotations]
## Limit/reorder to common set of features ... OK [31099 features x 3 cell types]
## Checking data dimension compatibility ... OK [31099 features x 3 cell types]
## Using cell type signatures: ✬Brain✬, ✬Liver✬, ✬Lung✬ [3 total]
## Checking log-scale ... data:YES - signatures:YES
## Reverting log-transform on signatures (base 2) ... OK
## Reverting log-transform on data (base 2) ... OK
## Normalizing signatures and target together (method: none) ... SKIP
## Using ged algorithm: "lsfit"
## Estimating cell proportions from cell-specific signatures [lsfit: ls]
## Timing:
## user system elapsed
## 4.200 0.304 4.515
## GED final wrap up ... OK
# plot against known proportions


profplot(mix, res)






# extract data for pure samples only
pure <- pureSamples(gse)
# compute p-values for all probes
ml <- extractMarkers(pure, pure$Type, method = "Abbas")
# all probes get attributed a cell-type and a p-value
summary(ml)






# Filtering 1: show p-values histogram
hist(ml, breaks = 20)
summary(ml <= 10^-8)




# refit proportions using only the subset of markers with p-value <= 10^-8
res2 <- ged(mix, basis(gse), subset = ml <= 10^-8, log = FALSE)
# plot against known proportions
profplot(mix, res2)




# Filtering 2: select limited number of markers based on the signature
# matrix✬s condition number as proposed by Abbas et al. (2009)
sel <- screeplot(ml, basis(gse), range = 1:500)
summary(sel)
## <object of class MarkerList>
## Types: 3 [✬Brain✬, ✬Liver✬, ✬Lung✬]
## Mode: numeric
## Markers: 14
## IDtype: .Affymetrix [✬1387772_at✬, ✬1391092_at✬, ..., ✬1388033_at✬]
## Values: [5.84863829573202e-11, 5.63037162881734e-10, ..., 2.13771098034808e-09]
## Source: rat2302.db
## Breakdown:
## Brain Liver Lung
## 5 5 4
# refit proportions using the optimised set of markers
res3 <- ged(mix, basis(gse), subset = sel, log = FALSE)
# plot against known proportions
profplot(mix, res3)
