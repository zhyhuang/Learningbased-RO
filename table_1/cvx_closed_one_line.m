function [ cvx_status cvx_optval x ] = cvx_closed_one_line(  c,phi_quantile,rt_sigma,miu_0,b1)
%UNTITLED2 Summary of this function goes here
% [ cvx_status cvx_optval ] = cvx_elli_LP( s1_1,M_spl,A_m,b1,A2,b2 )
d=length(c);
cvx_begin
    variable x(d);
    minimize(c'*x)
    subject to
   
    phi_quantile*norm(rt_sigma*x)-b1+miu_0*x <= 0;
 
cvx_end

end

