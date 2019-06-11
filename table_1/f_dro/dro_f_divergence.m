% f divergence

clear
% parameter setting
addpath ..
n=11;                 % dimension
m=1;                 % number of constraints
N_data=120;          % sample size
N_test_data=10000;   % test size
n_outer=1000;        % outer test size

N1_f_dro=N_data/2;
N2_f_dro=N_data-N1_f_dro;
N_SAA_data=50000;

rng(123)

alpha=0.05;
beta=0.05;
% LP setting
load('c_sigma_for_11.mat') % d11 in paper
% load('c_sigma_for_25.mat')
% load('c_sigma_for_50.mat')
% load('c_sigma_for_100.mat')
% load('c_sigma_50_set2.mat')  % d50  in paper
% load('c_sigma_100_set2.mat')  % d100 in paper
% c=-20*(ones(n,1)+randn(n,1));


A1=-c'; %1
[A1_r A1_c]=size(A1);
b1=  [1200];
 c_t=c';
  
  A1_tsps=A1';
  miu_0=A1_tsps(:)';   %%% NOTE here!!!!!!!!!!!
  
% sigma=wishrnd(eye(A1_r*A1_c),A1_r*A1_c);
  
  obj_fun=@(x) c'*x;

    
%   
%   
fv_f_dro=zeros(n_outer,1);
vio_f_dro=zeros(n_outer,1);
kernel_bw_est=N1_f_dro^(-1/(n+1));
kernel_bw_f0=N1_f_dro^(-1/(n+4));
gg=@(x,d) (exp(-d)*x.^0.95-1)./(x-1);
gg_gradient=@(x,d) (1-0.05*exp(-d)*x.^0.95-exp(-d)*0.95*x.^-0.05)./(x-1).^2;
kl_record=zeros(n_outer,1);
new_alpha_record=zeros(n_outer,1);

for i=1:n_outer
    
    dataset_points_1=mvnrnd(miu_0,sigma,N_data);
    data_set_1=dataset_points_1(1:N1_f_dro,:);
    data_set_2=dataset_points_1(N1_f_dro+1:end,:);
    % bootstrap
    B=1000;
    estimated_kl_divergence=zeros(B,1);
    
%     for i_B=1:B
%         resample_ind=randi(N2_f_dro,N2_f_dro,1);
%         boostrap_data=data_set_2(N2_f_dro,:);
% 
%         likelihood=Gaus_KDE(data_set_2,data_set_2,kernel_bw_est)./Gaus_KDE(data_set_1,data_set_2,kernel_bw_f0);
%         % likelihood=KNN_density(boostrap_data,boostrap_data,K_knn)./mvksdensity(dataset_points_1,boostrap_data,'Bandwidth',2);
% 
%         estimated_kl_divergence(i_B)=mean(log(likelihood));
%     end
%     
%     KL_est=quantile(estimated_kl_divergence,0.95);
%     kl_record(i)=KL_est;
%     
%     
%     low_point=0;
%     high_point=1;
%     while high_point-low_point>10^-10
%        
%         mid_point=(high_point+low_point)/2;
%         mid_value=gg_gradient(mid_point,KL_est);
%         if mid_value>0
%             high_point=mid_point;
%             high_value=gg_gradient(high_point,KL_est);
%         else
%             low_point=mid_point;
%             low_value=gg_gradient(low_point,KL_est);
%         end
%     
%     end
%     
%     new_alpha_record(i)=1-gg(mid_point,KL_est);
    
    Samples_for_SAA=sample_KDE(data_set_1,kernel_bw_f0,N_SAA_data);
    
    A_f_dro_1=reshape(Samples_for_SAA',A1_c,N_SAA_data*A1_r)';
    b_f_dro_1=repmat(b1,N_SAA_data,1);
    
    A_f_dro=[A_f_dro_1];
    b_f_dro=[b_f_dro_1];
    
    cvx_begin
                variable x_f_dro(A1_c)
                minimize c_t*x_f_dro
                subject to
                A_f_dro*x_f_dro <= b_f_dro
    cvx_end
    
    fv_f_dro(i)=obj_fun(x_f_dro);
    vio_f_dro(i)=1-normcdf((b1-A1*x_f_dro)/norm(sqrtm(sigma)*x_f_dro));
    disp(['Iteration NO.' num2str(i)])
end


val_obj=mean(fv_f_dro)
epsilon_hat=mean(vio_f_dro)
delta_hat=mean(vio_f_dro>0.05)

