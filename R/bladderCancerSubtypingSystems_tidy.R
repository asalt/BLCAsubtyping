library(tidyverse)
library(ComplexHeatmap)

MDA.predict <- function(Exp, Gpl = NULL, Symbol = "Symbol", nmin = 2000) {
  # Create or format Gpl if null
  if (is.null(Gpl)) {
    Gpl <- tibble(Probe.ID = rownames(Exp), !!Symbol := rownames(Exp))
  }
  
  probes <- mda.training$probes
  probes. <- intersect(probes, rownames(Exp))
  d2 <- cit.quantileNormalize(Exp) # Assuming this is a defined function
  
  if (length(probes.) > nmin) {
    d1 <- mda.training$exp[probes., ]
    d2 <- Exp[probes., ]
  } else {
    G <- intersect(
      mda.training$gpl[probes, "Symbol"], 
      Gpl[rownames(Gpl) %in% rownames(Exp), Symbol]
    )
    
    probes1 <- probes[mda.training$gpl[probes, "Symbol"] %in% G]
    sd1 <- mda.training$exp[probes1, ] %>% 
      apply(1, sd)
    
    probes1 <- split(sd1, mda.training$gpl[probes1, "Symbol"]) %>%
      map(which.max) %>%
      unlist(use.names = TRUE)
    
    d1 <- mda.training$exp[probes1, ]
    rownames(d1) <- names(probes1)
    
    probes2 <- Gpl %>% 
      filter(rowname(.) %in% rownames(Exp) & .[[Symbol]] %in% G) %>%
      pull(rowname(.))
    
    sd2 <- Exp[probes2, ] %>% 
      apply(1, sd, na.rm = TRUE)
    
    probes2 <- split(sd2, Gpl[probes2, Symbol]) %>%
      map(which.max) %>%
      unlist(use.names = TRUE)
    
    d2 <- Exp[probes2, ]
    rownames(d2) <- names(probes2)
  }
  
  G <- intersect(rownames(d1), rownames(d2))
  d1 <- d1[G, ] - apply(d1[G, ], 1, median, na.rm = TRUE)
  d2 <- d2[G, ] - apply(d2[G, ], 1, median, na.rm = TRUE)
  d <- cbind(d1, d2)
  
  tmp <- as.matrix(dist(t(d)))[colnames(d1), colnames(d2)]
  cl <- mda.training$clin[apply(tmp, 2, which.min), "Cluster"]
  
  MDA.subtype <- setNames(cl, colnames(d2))
  
  # Optionally create a heatmap
  if (interactive()) {
    Heatmap(matrix = as.matrix(d), name = "Expression")
  }
  
  return(MDA.subtype)
}

