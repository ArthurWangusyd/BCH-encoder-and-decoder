function code = bch_encode(msg,n,k,g)
      code = zeros(1,n);
      code(1:k) = msg;
      A = code;
      B = g;
      g_len = length(g);
      for i = 1:n-g_len+1
          if A(i) == 1
              for j = 1:g_len
                  A(i+j-1) = xor(A(i+j-1),B(j));  
              end
          end
          
      end 
      rr = A(k+1:n);
      code(k+1:n) = code(k+1:n) + rr;
end