clc
clear 
close
tic

pe = 0:0.001:0.01;
pe1 = 0.02:0.01:0.1;
pe = [pe pe1];
pe_len = length(pe);


BER = zeros(1,pe_len);
BER1 = zeros(1,pe_len);
t = 2;
n = 127;
e = 4;


for ii = 1:pe_len
    p = pe(ii);
    sum = 0;
    for i = 1:(e+1)
       sum = sum + (Cnk(n,(i-1)))*(p^(i-1))*((1-p)^(n-(i-1)));
    end
    BER(ii) = 1 - sum;
end





figure
semilogy(pe,BER,'-o');hold on;
xlabel('Pe');
ylabel('undetected');
grid on;

toc