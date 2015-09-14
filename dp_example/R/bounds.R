bound.error <- function(r, P, ra, Pa, gamma) {
    epsilon <- max(gamma/(1 - gamma) * inf.norm(P - Pa), 
                   inf.norm(r-ra)/inf.norm(r))
    if (epsilon < 1) 2*epsilon / (1 - epsilon)
    else Inf
   }
                        
  

compare.error.bound <- function(filename, df = 0.9, v.va = 11, v = 9, p.pa = 3, 
                  r.ra = 6, r = 4, ra = 5) {
   d <- read.table(filename, head=TRUE)
   
   rel.error <- d[,v.va] / d[,v]
   rel.bound <- array(Inf, nrow(d))
   epsilon <- pmax(df / (1 - df) * d[,p.pa], d[,r.ra] / d[,r])
   rel.bound <- 2*epsilon / (1-epsilon)
   
   ea <- d[,v.va]
   ba <- 1/(1-df) * (d[,r.ra] + d[,p.pa] * d[,r])
   ba2 <- 1/(1-df) * (d[,r] + d[,ra])
#    ob <- 1/(1-df) * (d[,r.ra] + 2*d[,r])
#    oy <- 1/(1-df) * (d[,ra] + d[,r])
    
   list(re = rel.error, rb = rel.bound, ea = ea, ba = ba, ba2 = ba2)
   }
   
compare.bounds <- function(delta.P, delta.r, max.r, gamma) {
      ep <- max(gamma / (1-gamma) * delta.P, delta.r/ max.r)
      bound.golub <- Inf
      if (ep < 1) bound.golub <- ((2 * ep) / (1 - ep))
      bound.mine <- 1/(1-gamma) * (delta.r + gamma * delta.P * max.r)
   #   print(paste(bound.golub, bound.mine))
      if (bound.golub < bound.mine) {
     #    print(paste(delta.P, delta.r, max.r, gamma))
         TRUE
         }
      else FALSE
      }
     
   
compare.bounds.batch <- function(delta.P, delta.r, max.r, gamma) {
   g <- 0
   m <- 0
   t <- 0
   for (dp in delta.P) {
      for (dr in delta.r) {
         for (mr in max.r) {
            for (g in gamma) {
               t <- t + 1
               if (compare.bounds(dp, dr, mr, g)) g <- g + 1
               else m <- m + 1
               }
            }
         }
      }
   m / t
   }
               
   
   
compare.rel.error.bound <- function(filename, cols=c(11,16,21), dfs =
                        c(0.1,0.3,0.5), v = cols-2, p.pa = 3, 
                        r.ra = 6, r = 4) {
   
   D <- read.table(filename, head=TRUE)
   
   for (d in 1:length(dfs)) {
      df <- dfs[d]
      # relative bound
      ep <- pmax(df / (1-df) * D[,p.pa], D[,r.ra] / D[,r])
      ep[ep >= 1] <- Inf
      rel.bound <- 2 * ep / (1-ep)
      
      # absolute bound
      coef <- 1/(1-df)
      bound <- coef * (D[,r.ra] + coef * df * D[,p.pa] * D[,r])
      
      ylim = c(0, max(bound, rel.bound))
      
      rel.error <- D[,cols[d]] / D[,v[d]]
      plot(rel.error, rel.bound, ylim = ylim)
      points(rel.error, bound, col=2)
      x <- c(min(rel.error),max(rel.error))
      lines(x,x)
      }
   }
