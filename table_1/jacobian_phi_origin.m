function [jac_mat] = jacobian_phi_origin(x)
%only for alpha=0.05
n_phi_input=length(x);
% n_phi_input=n_variable+n_variable*(n_variable+1)/2;
roos_v=roots([1/2 3/2 -n_phi_input]);
n_variable=round(roos_v(roos_v>0));
jac_mat=eye(n_phi_input);
jac_mat(1:n_variable,1:n_variable)=sqrt(0.05/(1-0.05))*eye(n_variable);
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

