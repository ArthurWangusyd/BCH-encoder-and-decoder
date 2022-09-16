function num = Cnk(n,k)
    kk = 1;
    for i = 1:k
        kk = kk*i;
    end

    n_k = 1;
    for i = (n-k+1):n
        n_k = n_k*i;
    end
    
    num = n_k/kk;
   
end