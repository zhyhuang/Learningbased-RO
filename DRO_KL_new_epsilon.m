function [new_epsilon] = DRO_KL_new_epsilon(dataset,epsilon)
%UNTITLED7 此处显示有关此函数的摘要
%   此处显示详细说明
[N_data_total,d]=size(dataset);
N1_f_dro=ceil(N_data_total/2);
% N2_f_dro=N_data_total-N1_f_dro;
kernel_bw_f0=3*N_data_total^(-1/(d+4));
K_knn=1;
 gg=@(x,d) (exp(-d)*x.^(1-epsilon)-1)./(x-1);
gg_gradient=@(x,d) (1-epsilon*exp(-d)*x.^(1-epsilon)-exp(-d)*(1-epsilon)*x.^-epsilon)./(x-1).^2;


    data_set_1=dataset(1:N1_f_dro,:);
    data_set_2=dataset(N1_f_dro+1:end,:);
    likelihood_mean=KNN_density(data_set_2,data_set_2,K_knn)./Gaus_KDE(data_set_1,data_set_2,kernel_bw_f0);
    estimated_kl_mean=mean(log(likelihood_mean));
    
    KL_est=estimated_kl_mean;
    low_point=0;
    high_point=1;
    while high_point-low_point>10^-10
       
        mid_point=(high_point+low_point)/2;
        mid_value=gg_gradient(mid_point,KL_est);
        if mid_value>0
            high_point=mid_point;
%             high_value=gg_gradient(high_point,KL_est);
        else
            low_point=mid_point;
%             low_value=gg_gradient(low_point,KL_est);
        end
    
    end
    
    new_epsilon=1-gg(mid_point,KL_est);

    
end

function [density] = KNN_density(data,pts,k)

[N_data,n_dim]=size(data);

sorted_euclidean_dist=sort(dist(pts,data'),2);


if isequal(data,pts)
knn_dist=sorted_euclidean_dist(:,k+1);
density=(k/(N_data-1)*1/(pi^(n_dim/2)/gamma(n_dim/2+1)))./knn_dist.^n_dim;   
else
knn_dist=sorted_euclidean_dist(:,k);
density=(k/(N_data)*1/(pi^(n_dim/2)/gamma(n_dim/2+1)))./knn_dist.^n_dim;   
    
end



end

function [density] = Gaus_KDE(data,pts,bw)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[N_data,n_dim]=size(data);

celled_pts=num2cell(pts, 2);

point_density=@(x) 1/N_data/(bw^n_dim)*sum(mvnpdf(x/(bw),data/(bw),eye(n_dim)));

density=cell2mat(cellfun(point_density,celled_pts,'UniformOutput',false));

end

