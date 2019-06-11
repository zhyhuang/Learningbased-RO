function [outputArg1,outputArg2] = moment_DRO_ccp(dataset_points_1,inputArg2)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
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


end

