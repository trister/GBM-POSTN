# code loaded by the GSEA crawler
## written by Inn Sock Jang, PhD

preRankedTest <- function(reference,geneset,np=1000,w=1){
  
  ## output should contains  
  # es    : Enrichment Score
  # nes   : Normalized Enrichment Score
  # pv    : Nonimal Pvalue
  # ledge : Leading Edges
  # es_all
  # isInGeneset
  ## Input with 2 gene sets
  # reference   : ranked gene list and their score
  # geneset     : positive gene set (e.g. from KEGG, MSigDB, etc)
  # np          : number of permutation
  # w           : weight
  
  # check input format
  if (length(reference) < length(geneset)) print("Reference List should include given Gene Set")
  
  nes <- NaN
  pv  <- NaN
  
  # combine ranked list and score
  
  pn      <- matrix(1,ncol = length(reference));
  
  sortRef <- sort(reference, decreasing = TRUE)
  
  isInGeneset <- sign(match(names(sortRef), geneset, nomatch=0))
  
  # compute ES
  score_hit   <- cumsum((abs(sortRef * isInGeneset)) ^w)
  
  if (sum(score_hit) == 0){
    print("no significant genes in geneset. Returning ES = 0")    
    #     return(list(es_all = NaN, es = NaN,nes = NaN,p.value = NaN, hit = NaN))
    return(list(es_all = NaN, es = NaN,nes = NaN,p.value = NaN,hit = NaN,reference = NaN, geneset = NaN))
  }
  
  score_hit   <- score_hit/score_hit[length(score_hit)]
  score_miss  <- cumsum( 1 - isInGeneset)
  score_miss  <- score_miss/score_miss[length(score_miss)]
  es_all      <- score_hit - score_miss
  #   es          <- max(abs(es_all)) + min(abs(es_all))
  
  max.ES <- max(es_all)
  min.ES <- min(es_all)
  
  if (max.ES > - min.ES) {    
    es <- max.ES
  } else {
    es <-min.ES
  }
  
  #plot(1:length(es_all), es_all, type="l")
  
  # identify leading edge
  #isen <- mat.or.vec(length(es_all),1)
  #if (es < 0){
  #  ixpk <- which(es_all==min(es_all))
  #  isen[ixpk:length(isen)] <- 1
  #}
  #else{
  #  ixpk <- which(es_all==max(es_all))
  #  isen[1:ixpk] <- 1
  #}
  
  #ledge <- sortRef$x[(isen==1)&(isInGeneset==1)];
  
  # compute norminal p-value
  if (np>0) { 
    bg.es <- mat.or.vec(np,1)
    
    for (i in 1:np){
      bg.isInGeneset  = isInGeneset[sample(length(isInGeneset))]      
      bg.hit   = cumsum((abs(sortRef * bg.isInGeneset))^w);
      if(sum(bg.hit)==0){
        bg.es[i]<- NA
      }
      else{
        bg.hit   = bg.hit/bg.hit[length(bg.hit)]
        bg.miss  = cumsum(1 - bg.isInGeneset);
        bg.miss  = bg.miss/bg.miss[length(bg.miss)];
        bg.all   = bg.hit - bg.miss;
        
        max.bES <- max(bg.all)
        min.bES <- min(bg.all)
        
        if (max.bES > - min.bES) {    
          bg.es[i] <- max.bES
        } else {
          bg.es[i] <- min.bES
        }
      }
    }
    if(es<0){
      pv  <- mean(bg.es[bg.es<0] < es,na.rm = TRUE)
      nes <- es/abs(mean(bg.es[bg.es<0],na.rm = TRUE))
    }
    else {
      pv  <- mean(bg.es[bg.es>0] > es,na.rm = TRUE)
      nes <- es/abs(mean(bg.es[bg.es>0],na.rm = TRUE))
    }
  }
  
  return(list(es_all = es_all, es = es,nes = nes,p.value = pv,hit = isInGeneset,reference = sortRef, geneset = geneset))
  
}