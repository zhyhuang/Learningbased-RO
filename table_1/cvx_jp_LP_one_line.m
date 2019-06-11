function [ cvx_status cvx_optval x ] = cvx_jp_LP_one_line( f,t_jp,x0_for_cut,b1,n)
%UNTITLED2 Summary of this function goes here
% [ cvx_status cvx_optval ] = cvx_elli_LP( s1_1,M_spl,A_m,b1,A2,b2 )


cvx_begin

    variable x(n);
    variable p(1);
    minimize(f'*x)
    subject to
        p*t_jp<=b1;
        p*x0_for_cut==x;
        p>=0;

     
cvx_end
end
