    for (u in 1:na) 
    {
      for (i in 1:n)
      {
        ss <- sum(Dh[i,,u])
        if (ss > 0) Dh[i,,u] <- Dh[i,,u] / ss
      }
    }

    for (i in 1:m)
    {
      ss <- sum(Kh[i,])
      if (ss > 0) Kh[i,] <- Kh[i,] / ss
    }
    
