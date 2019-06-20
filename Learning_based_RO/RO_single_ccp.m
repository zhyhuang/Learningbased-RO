function [x_RO] = RO_ccp(dataset_1,dataset_2,data_rank,c,b)
% [x_RO] = RO_ccp(dataset_1,dataset_2,data_rank,c,b)
% Learning-based RO with ellipsoidal set.

[n1,d]=size(dataset_1);

%% phase 1
sample_mean=mean(dataset_1);

if n1>d
    sample_cov=((cov(dataset_1)));
else
    sample_cov=diag(diag(cov(dataset_1)));
end

%% phase 2
n2=size(dataset_2,1);
sample_mean_mat=repmat(sample_mean,n2,1);  
T_value=(dataset_2-sample_mean_mat)*sample_cov^-1*(dataset_2-sample_mean_mat)';  % calculate only for diag
t_value=sqrt(diag(T_value));
[sorted_t,sort_index]=sort(t_value);
s_value=sorted_t(data_rank);

%% solve

root_sample_cov=sqrtm(sample_cov);
cvx_begin
    variable x_RO(d);
    minimize(c'*x_RO)
    subject to       
        sample_mean*x_RO+(s_value)*norm(root_sample_cov*x_RO)-b <= 0;
cvx_end

end

