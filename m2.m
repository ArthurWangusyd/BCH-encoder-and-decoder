clc
clear 
close
tic

pe = 0:0.001:0.01;
pe1 = 0.02:0.01:0.1;
pe = [pe pe1];
pe_len = length(pe);

ebn0 = zeros(1,pe_len);

for i = 2:pe_len   
    ebn0(i) = -log(pe(i)); 
end
    
BER = zeros(1,pe_len);
BER1 = zeros(1,pe_len);
t = 2; 
n = 127;
m = 7;
for ii = 1:pe_len
    p = pe(ii);
    sum = 0;
    for i = 1:(t+1)
       sum = sum + (Cnk(n,(i-1)))*(p^(i-1))*((1-p)^(n-(i-1)));
    end
    BER(ii) = 1 - sum;
end

k = 113;
px = [1 0 0 0 1 0 0 1];  %x^7 + x^3 +1 
g = [1 0 0 0 0 1 1 0 1 1 1 0 1 1 1]; %x^14 + x^9 + x^8 + x^6 + x^5 + x^4 + x^2 + x + 1

pxnum = 0; %域生成多项式的10进制值
for i = 1:length(px)
    pxnum = pxnum + px(i)*2^(length(px)-i);
end
alpha = gf(2,m,pxnum);
rate = k/n;

codenum = 10000;

diff = 10*log10(2);
esn0 = ebn0;
esn0(2:pe_len) = esn0(2:pe_len) + diff;


hard = zeros(1,2*n);
hard1 = zeros(1,n);
hard2 = zeros(1,n);

len = length(ebn0);

err = zeros(1,len);
err1 = zeros(1,len);

for i = 2:length(ebn0)
    %i
    bit = 0;
    en = 10^(esn0(i)/10);   
	sigma = 1/sqrt(2*rate*en);
    
    for j = 1:1:codenum/2  
        
        msg = randi([0 1],1,k);
        msg1 = randi([0 1],1,k);
       %% 编码
        encode = bch_encode(msg,n,k,g);
        encode1 = bch_encode(msg1,n,k,g);
        
        encode_2 = [encode encode1]; %由于码长为127 不是2的倍数 这里一次生成2个码字进行QPSK传输
        
        encode_2_QPSK_I = zeros(1,n);
        encode_2_QPSK_Q = zeros(1,n);
      
       %% QPSK星座映射
       for ii = 1:1:n
           if encode_2(2*ii-1) == 0 && encode_2(2*ii) == 0
               encode_2_QPSK_I(ii) = sqrt(2)/2;
               encode_2_QPSK_Q(ii) = sqrt(2)/2;
           end
            if encode_2(2*ii-1) == 0 && encode_2(2*ii) == 1
               encode_2_QPSK_I(ii) = -sqrt(2)/2;
               encode_2_QPSK_Q(ii) = sqrt(2)/2;
            end
            if encode_2(2*ii-1) == 1 && encode_2(2*ii) == 1
               encode_2_QPSK_I(ii) = -sqrt(2)/2;
               encode_2_QPSK_Q(ii) = -sqrt(2)/2;
            end
            if encode_2(2*ii-1) == 1 && encode_2(2*ii) == 0
               encode_2_QPSK_I(ii) = sqrt(2)/2;
               encode_2_QPSK_Q(ii) = -sqrt(2)/2;
            end   
       end
       %% 加噪声
        rev_I = encode_2_QPSK_I + sigma*randn(1,n);
        rev_Q = encode_2_QPSK_Q + sigma*randn(1,n);
       %% QPSK星座解映射-硬判
        for jj = 1:n
            if rev_I(jj) > 0 && rev_Q(jj) > 0
                hard(2*jj-1) = 0;
                hard(2*jj) = 0;
            end
             if rev_I(jj) < 0 && rev_Q(jj) > 0
                hard(2*jj-1) = 0;
                hard(2*jj) = 1;
             end
             if rev_I(jj) < 0 && rev_Q(jj) < 0
                hard(2*jj-1) = 1;
                hard(2*jj) = 1;
             end
             if rev_I(jj) > 0 && rev_Q(jj) < 0
                hard(2*jj-1) = 1;
                hard(2*jj) = 0;
             end  
        end
        %% 译码
         hard1 = hard(1:n);
         hard2 = hard(n+1:2*n);
        
         decode = bch_decode_m(hard1,pxnum,t,n,k);
         bit = bit + length(find(decode ~= msg));
        
         decode1 = bch_decode_m(hard2,pxnum,t,n,k);
         bit = bit + length(find(decode1 ~= msg1));

      
 
        
    end
    err(i) = bit/k/codenum;
end

figure
semilogy(pe,BER,'-o');hold on;
semilogy(pe,err,'-o');hold on;
xlabel('Pe');
ylabel('BER');
legend('Theoretical value','real BER');
grid on;

toc