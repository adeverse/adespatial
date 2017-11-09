# La fonction crée les MEM pour la matrice W demandée. Si plusieurs exposants (y) sont 
# fournis, toutes les matrices W résultantes sont retournées. Dans l'ancienne fonction
# test.W, seul le meilleur modèle était rendu et celui-ci était sélectionné sur base d'un
# critère de AICc.
# Cette méthode ne doit plus être utilisée dorénavant. Et le seuil de significativité
# du spatial doit être corrigé pour le nombre de matrices W comparées, y compris les
# différents exposants au sein d'une matrice de pondération. Ceci est fait par la fonction
# MEM.modsel mais pas par test.W.R2 qui ne fait que créer des matrices W.

# **************************************************************************************
test.W.R2 <- function(Y, nb, xy, style = "B", MEM.autocor = c("positive", "negative",
                                                              "all"), f = NULL, ...) 
{
  mycall <- pairlist(...)   
  res <- list()  
  MEM.autocor <- match.arg(MEM.autocor)
  control <- FALSE
  
  if (!(is.null(f))) {
    nbdist <- nbdists(nb, as.matrix(xy))
    if (!(is.null(mycall))) {   
      param <- expand.grid(as.list(mycall))
      m1 <- match(names(param), names(formals(f)))
      for (i in 1:nrow(param)) {
        formals(f)[m1] <- unclass(param[i, ])
        res[[i]] <- scores.listw(nb2listw(nb, style = style, 
                                          glist = lapply(nbdist, f), zero.policy = TRUE),
                                 MEM.autocor = MEM.autocor)
        if (i > 1) control <- TRUE
      }
    }
    else {   
      res[[1]] <- scores.listw(nb2listw(nb, style = style, glist = lapply(nbdist, f)), 
                               MEM.autocor = MEM.autocor)
    }
  }
  else {   
    res[[1]] <- scores.listw(nb2listw(nb, style = style), MEM.autocor = MEM.autocor)
  }
  
  if (control == TRUE) list(param = param, MEM = res)
  else list(MEM = res[[1]])
}                                                        # End of the function test.W.R2
# **************************************************************************************