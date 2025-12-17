suppressPackageStartupMessages(library(gtools))

built.all.models <- function(env.cov.names, env.cov.int, env.cov.2){
  all.mods <- permutations(n = 2,r = length(env.cov.names) + length(env.cov.int) +
                             length(env.cov.2), v = c(0, 1), repeats.allowed = T)
  colnames(all.mods) <- c(env.cov.names,
                          unlist(lapply(env.cov.int, function(x){paste(x,collapse=":")})),
                          unlist(lapply(env.cov.2, function(x){paste(c("I(", x, "^2)"),
                                                                     collapse="")})))
  all.problems <- c(env.cov.int, lapply(env.cov.2, function(x){c(x,x)}))
  if (length(all.problems) > 0){
    for(i in 1:length(all.problems)){
      to.del <- all.mods[ , length(env.cov.names) + i] == 1 &
        (all.mods[ ,colnames(all.mods) == (all.problems[[i]])[1]] == 0 |
           all.mods[ ,colnames(all.mods) == (all.problems[[i]])[2]] == 0
        )
      all.mods <- all.mods[to.del == F,     ]
    }
  }
  all.mods <- cbind("(Intercept)" = rep(1, nrow(all.mods)), all.mods)
  return(all.mods)
}
