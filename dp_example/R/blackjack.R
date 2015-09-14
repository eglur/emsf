dyn.load(paste("bj_evaluate_policy", .Platform$dynlib.ext, sep = "")) #melhorar



bj.draw.card <- function(hand)
{
  # hand = list(sum.cards, has.ace)
  
   c <- sample(1:13, 1) # 1 -> ace; >= 11  -> facecards
   if (c > 10) c <- 10
   
   hand$sum.cards <- hand$sum.cards + c
   
   if (!hand$has.ace && c ==1)
   {
      hand$has.ace <- TRUE
      hand$sum.cards <- hand$sum.cards + 10
   }
   
   if (hand$has.ace && hand$sum.cards > 21)
   {
     hand$sum.cards <- hand$sum.cards - 10
     hand$has.ace <- FALSE
   }
   
   hand
}


bj.transition <- function(s, a)
{
  # a = 0 -> stick; a = 1 -> hit
  r <- 0
  terminal <- FALSE
  
  if (a == 2)
  {
    s$ph <- bj.draw.card(s$ph)
    if (s$ph$sum.cards > 21)
    {
      r <- -1
      s$idx <- 201 # agent loses
    }
  }
  else
  {
    while (s$dh$sum.cards < 17) s$dh <- bj.draw.card(s$dh)

    if (s$dh$sum.cards > 21) 
    {
      r <- 1
      s$idx <- 203 # agent wins
    }
    else
    {
      p.diff <- 21 - s$ph$sum.cards
      d.diff <- 21 - s$dh$sum.cards
      
      if (p.diff < d.diff) 
      {
	r <- 1
	s$idx <- 203 # agent wins
      }
      else if (d.diff < p.diff) 
      {
	r <- -1
	s$idx <- 201 # agent loses
      }
      else s$idx <- 202 # draw
    }
    
  }
  list(s = s, r = r)
}  

bj.state.index <- function(s)
{
  if (s$idx != 0) s$idx
  else (s$ph$has.ace * 100) + ((s$ph$sum.cards - 12) * 10) + (s$dh$sum.cards - 2) + 1
}

bj.evaluate.policy <- function(pi, num.episodes = 100)
{
  
  R <- array(0, num.episodes)
  for (i in 1:num.episodes)
  {
    ph <- bj.draw.card(list(sum.cards = 0, has.ace = FALSE))
    dh <- bj.draw.card(list(sum.cards = 0, has.ace = FALSE))
    
    while (ph$sum.cards < 12) ph <- bj.draw.card(ph) # nothing to do otherwise
      
    s <- list(ph = ph, dh = dh, idx = 0)
    
    while (s$idx == 0)
    {
     
      a <- pi[bj.state.index(s)]
      
      t <- bj.transition(s, a)
     
      s <- t$s
      if (s$idx != 0) R[i] <- t$r
    }
  }
  
  mean(R)
}	

bj.evaluate.policy.fast <- function(pi, num.episodes = 100)
{
   ret <- 0
   U <- .C("bj_evaluate_policy", as.integer(pi), as.integer(num.episodes), ret)
   U[[3]]
}    
  


