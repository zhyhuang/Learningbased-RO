function [x_mo_DRO] = moment_DRO_ccp(dataset,c,b,epsilon,delta)
%[x_mo_DRO] = moment_DRO_ccp(dataset,c,b,epsilon,delta)
%   moment_based DRO

[N_data,d]=size(dataset);

 triu_ind=triu(true(d));
%   
 col_num=repmat(1:d,d,1);
 row_num=col_num';
%  
 row_ind_use=row_num(triu_ind);
 col_ind_use=col_num(triu_ind);
%  
% 
d_phi_input=d+d*(d+1)/2;

dataset_for_theta=[dataset, dataset(:,row_ind_use).*dataset(:,col_ind_use)];

    theta_estimate=mean(dataset_for_theta);
    sigma_estimate=cov(dataset_for_theta,1);

   value_jacob= jacobian_phi_origin(theta_estimate,delta);

    V_est=value_jacob*sigma_estimate*value_jacob';
    V_est_inv=inv(V_est);

    mu_hat=mean(dataset);
    Sigma_hat=cov(dataset,1);
    sqrt_V=sqrtm(V_est_inv);
    R_inv=inv(sqrt_V);
    rho=sqrt(chi2inv(1-delta,d_phi_input)/N_data);
    tilde_c=-b;

    svec_operator=sqrt(2)*ones(d)-sqrt(2)*eye(d)+eye(d);
    svec_multiplier=svec_operator(triu_ind);
    A=[eye(d),zeros(d,d_phi_input-d);
        zeros(d_phi_input-d,d),diag(svec_multiplier)];

        cvx_begin
                    variable x_mo_DRO(d)
                    variable W(d,d)
                    variable vec_q(d_phi_input)
                    variable eta
                    minimize (c'*x_mo_DRO)
                    subject to
                    sqrt(epsilon/(1-epsilon))*mu_hat*x_mo_DRO + trace (Sigma_hat *W) + rho *norm((A*R_inv)'*(vec_q)) +sqrt(epsilon/(1-epsilon))* tilde_c+eta/4  <= 0
                    vec_q==[x_mo_DRO;W(triu_ind).*svec_multiplier];
                    [W,x_mo_DRO;x_mo_DRO',eta] == semidefinite(d+1)

        cvx_end


end

function [jac_mat] = jacobian_phi_origin(x,delta)
%only for alpha=0.05
n_phi_input=length(x);
% n_phi_input=n_variable+n_variable*(n_variable+1)/2;
roos_v=roots([1/2 3/2 -n_phi_input]);
n_variable=round(roos_v(roos_v>0));
jac_mat=eye(n_phi_input);
jac_mat(1:n_variable,1:n_variable)=sqrt(delta/(1-delta))*eye(n_variable);
row_num=n_variable;
for i=1:n_variable
    for j=1:i
        row_num=row_num+1;
        if j==i
            jac_mat(row_num,j)=-2*x(i);
        else
            jac_mat(row_num,j)=-x(i);
            jac_mat(row_num,i)=-x(j);
        end
    end
    
end
end



