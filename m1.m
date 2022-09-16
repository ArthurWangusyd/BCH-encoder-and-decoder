tic
clc
clear 
close

n = 127;  
k = 113;
m = 7;
px = [1 0 0 0 1 0 0 1];  %x^7 + x^3 +1 
g = [1 0 0 0 0 1 1 0 1 1 1 0 1 1 1]; %x^14 + x^9 + x^8 + x^6 + x^5 + x^4 + x^2 + x + 1
t = 2;
pxnum = 0; %域生成多项式的10进制值
for i = 1:length(px)
    pxnum = pxnum + px(i)*2^(length(px)-i);
end
alpha = gf(2,m,pxnum);

rate = k/n;

codenum = 10000;
ebn0 = 1:0.5:10;

diff = 10*log10(2);
esn0 = ebn0 + diff;


hard = zeros(1,2*n);
hard1 = zeros(1,n);
hard2 = zeros(1,n);
hard3 = zeros(1,2*k);
hard4 = zeros(1,k);
hard5 = zeros(1,k); 

len = length(ebn0);

err = zeros(1,len);
err1 = zeros(1,len);
err2 = zeros(1,len);
err3 = zeros(1,len);

for i = 1:length(ebn0)
    %i
    bit = 0;
    bit1 = 0;
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
        % now find the BER of uncoded
        uncode_2 = [msg msg1];
        uncode_2_QPSK_I = zeros(1,k);
        uncode_2_QPSK_Q = zeros(1,k); % I pahse and Q phase
        %% QPSK modulation
       for ii = 1:1:k
           if uncode_2(2*ii-1) == 0 && uncode_2(2*ii) == 0
               uncode_2_QPSK_I(ii) = sqrt(2)/2;
               uncode_2_QPSK_Q(ii) = sqrt(2)/2;
           end
            if uncode_2(2*ii-1) == 0 && uncode_2(2*ii) == 1
               uncode_2_QPSK_I(ii) = -sqrt(2)/2;
               uncode_2_QPSK_Q(ii) = sqrt(2)/2;
            end
            if uncode_2(2*ii-1) == 1 && uncode_2(2*ii) == 1
               uncode_2_QPSK_I(ii) = -sqrt(2)/2;
               uncode_2_QPSK_Q(ii) = -sqrt(2)/2;
            end
            if uncode_2(2*ii-1) == 1 && uncode_2(2*ii) == 0
               uncode_2_QPSK_I(ii) = sqrt(2)/2;
               uncode_2_QPSK_Q(ii) = -sqrt(2)/2;
            end   
       end
     
      
        unrev_I = uncode_2_QPSK_I + sigma*randn(1,k);
        unrev_Q = uncode_2_QPSK_Q + sigma*randn(1,k); %pass AWGN channel
       
        %% QPSK demodulation
        for jj = 1:k
            if unrev_I(jj) > 0 && unrev_Q(jj) > 0
                hard3(2*jj-1) = 0;
                hard3(2*jj) = 0;
            end
             if unrev_I(jj) < 0 && unrev_Q(jj) > 0
                hard3(2*jj-1) = 0;
                hard3(2*jj) = 1;
             end
             if unrev_I(jj) < 0 && unrev_Q(jj) < 0
                hard3(2*jj-1) = 1;
                hard3(2*jj) = 1;
             end
             if unrev_I(jj) > 0 && unrev_Q(jj) < 0
                hard3(2*jj-1) = 1;
                hard3(2*jj) = 0;  %Hard decision. 
             end  
        end
        hard4 = hard3(1:k);
        hard5 = hard3(k+1:2*k);
        bit1 = bit1 + length(find(hard4 ~= msg));
        bit1 = bit1 + length(find(hard5 ~= msg1));
 
        
    end
    err(i) = bit/k/codenum;
    err2(i) = bit1/k/codenum; %calculate the average error
end

figure
semilogy(ebn0,err,ebn0,err2,'-o');hold on; %plot the figure 
xlabel('SNR');
ylabel('BER');
legend('code','uncode');

grid on;
















toc