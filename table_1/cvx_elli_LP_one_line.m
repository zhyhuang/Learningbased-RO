function [ cvx_status cvx_optval x ] = cvx_elli_LP_one_line( c,s1_1,M_spl,A_m,b1,n)
%UNTITLED2 Summary of this function goes here
% [ cvx_status cvx_optval ] = cvx_elli_LP( s1_1,M_spl,A_m,b1,A2,b2 )

cvx_begin
    variable x(n);
    minimize(c'*x)
    subject to
   
    
    A_m(1,:)*x+sqrt(s1_1)*norm(M_spl(:,:,1)*x)-b1(1) <= 0;
 
cvx_end

end

