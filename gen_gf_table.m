function [b] = gen_gf_table(px)

gf_m = length(px)-1;     %���m
gf_length = 2^gf_m-1;      %���Ԫ�ظ���
gf_table = zeros(1,gf_length);    %������Ԫ��ָ����10�������Ķ�Ӧ�� �±��ʾָ������Ϊ��ֵ 
                                  %��gf_table��3�� 3-2=1 ��ʾa^1 = gf_table��3��= 2 10������2
%% ������� gf_table(1)ָ�� �±��ȥ1 Ϊָ��ֵ
gf_table(1) = 1;
for i = 2:gf_m
    gf_table(i) = gf_table(i-1)*2;
end
ta = px(2:end); %��ʾ��Ķ�Ӧ��ϵ ��px = x^3 + x + 1 ,��ôa^3 = a + 1 ��ʾ011
ta_10 = 0;
for i = 1:length(ta)
    ta_10 = ta_10 + ta(i)*2^(gf_m-i);
end
gf_table(1+gf_m) = ta_10;

for i = gf_m+2:gf_length
    tt = dec2bin(gf_table(i-1),gf_m);
    tt1 = boolean(tt-'0');
    tt2 = double(tt1);
     
    ta_table = tt2;
    
    for j = 1:gf_m
        if ta_table(j) == 1
            gf_table(i) =  bitxor(gf_table(i),gf_table(gf_m+2-j));
        end
    end
end

gf_index_table = zeros(1,gf_length);
gf_index_table(1) = gf_length;

for i = 2:gf_length
    for j = 1:gf_length
        if gf_table(j) == i
            gf_index_table(i) = j-1;
        end
    end
    
end

b = gf_index_table;

end