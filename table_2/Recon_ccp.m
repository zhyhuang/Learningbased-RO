function [x_Recon] = Recon_ccp(dataset_1,dataset_2,data_rank_1,data_rank_2,c,b)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
d=length(c);
%% phase 1
[x_0] = RO_ccp(dataset_1,dataset_1,data_rank_1,c,b);

%% phase 2

t_value=   dataset_2* x_0;
[sorted_t,sort_index]=sort(t_value);
s_value=sorted_t(data_rank_2);

%% solve


    cvx_begin

        variable x_Recon(d);
        variable p(1);
        minimize(c'*x_Recon)
        subject to
            p*s_value<=b;
            p*x_0==x_Recon;
            p>=0;


    cvx_end

end

