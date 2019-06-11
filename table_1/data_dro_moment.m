% test case
clear
% parameter setting

n=100;                 % dimension
m=1;                 % number of constraints
N_data=120;          % sample size
N_test_data=10000;   % test size
B_2=60;              % phase II budget
B_1=N_data-B_2;      % phase I budget
n_outer=1000;        % outer test size

N1_fast=61;
N2_fast=59;

% N1_fast=318;
% N2_fast=N_data-N1_fast;

% N1_fast=2326;
% N2_fast=N_data-N1_fast;

rng(123)

x_fast_0=zeros(n,1);

% rng(123) first run
rng(321) % second run

alpha=0.05;
beta=0.05;
rank_of_data=binoinv(1-beta,B_2,1-alpha); % estimated quantile
rank_of_data_p1=binoinv(1-beta,B_1,1-alpha);
% LP setting
% load('c_sigma_for_11.mat') % d11 in paper
load('c_sigma_100_set2.mat')  % d100 in paper


A1=-c'; %1
   c_t=c';
 
    [A1_r A1_c]=size(A1);

b1=  [1200];
 
  
  A1_tsps=A1';
  miu_0=A1_tsps(:)';   %%% NOTE here!!!!!!!!!!!
  
% sigma=wishrnd(eye(A1_r*A1_c),A1_r*A1_c);
  
  obj_fun=@(x) c'*x;
% 
  triu_ind=triu(true(n));
%   
 col_num=repmat(1:n,n,1);
 row_num=col_num';
%  
 row_ind_use=row_num(triu_ind);
 col_ind_use=col_num(triu_ind);
%  
% 
n_phi_input=n+n*(n+1)/2;
% x_var=sym('x',[n_phi_input,1]);
% % syms a b
% multi_mat=x_var(1:n)*x_var(1:n).';
% phi_origin=[x_var];
% phi_origin(1:n)=sqrt(alpha/(1-alpha))*phi_origin(1:n);
% phi_origin(n+1:end)=phi_origin(n+1:end)-multi_mat(triu_ind);
% 
% jacob_phi_origin1=jacobian(phi_origin, x_var);
% phi_fun_temp(x_var)=jacob_phi_origin1;
% phi_fun=matlabFunction(phi_fun_temp,'Vars',{x_var.'});


fv_ddm=zeros(n_outer,1);
ddm_vio=zeros(n_outer,1);

for i=1:n_outer

    dataset_points_1=mvnrnd(miu_0,sigma,N_data);
    tic
    dataset_for_theta=[dataset_points_1, dataset_points_1(:,row_ind_use).*dataset_points_1(:,col_ind_use)];

    theta_estimate=mean(dataset_for_theta);
    new_sigma=cov(dataset_for_theta,1);

   value_jacob= jacobian_phi_origin(theta_estimate);

    V_est=value_jacob*new_sigma*value_jacob';
    V_est_inv=inv(V_est);

    mu_hat=mean(dataset_points_1);
    Sigma_hat=cov(dataset_points_1,1);
    sqrt_V=sqrtm(V_est_inv);
    R_inv=inv(sqrt_V);
    rho=sqrt(chi2inv(1-beta,n_phi_input)/N_data);
    tilde_c=-b1;

    svec_operator=sqrt(2)*ones(n)-sqrt(2)*eye(n)+eye(n);
    svec_multiplier=svec_operator(triu_ind);
    A=[eye(n),zeros(n,n_phi_input-n);
        zeros(n_phi_input-n,n),diag(svec_multiplier)];

        cvx_begin
                    variable x_data_dro(n)
                    variable W(n,n)
                    variable vec_q(n_phi_input)
                    variable eta
                    minimize (c'*x_data_dro)
                    subject to
                    sqrt(alpha/(1-alpha))*mu_hat*x_data_dro + trace (Sigma_hat *W) + rho *norm((A*R_inv)'*(vec_q)) +sqrt(alpha/(1-alpha))* tilde_c+eta/4  <= 0
                    vec_q==[x_data_dro;W(triu_ind).*svec_multiplier];
                    [W,x_data_dro;x_data_dro',eta] == semidefinite(n+1)

        cvx_end
        fv_ddm(i)=obj_fun(x_data_dro);
    
    
        ddm_vio(i)=1-normcdf((b1-A1*x_data_dro)/norm(sqrtm(sigma)*x_data_dro));

       time_dro(i)=toc;
    
    disp(i)
    
    
    
    
end
    
rest=[mean(fv_ddm);
    mean(ddm_vio);
    sum(ddm_vio>0.05)/1000];