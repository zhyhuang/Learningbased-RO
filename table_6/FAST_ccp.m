function [x_FAST] = FAST_ccp(dataset_1,dataset_2,c,b,x_0)
% [x_FAST] = FAST_ccp(dataset_1,dataset_2,c,b,x_0)
% Using FAST algorithm.
% Provide solution for minimize c'x s.t. a'x <= b with two datasets for a
% and an initial solution x_0.
    c_t=c';
    d=length(c);
    [N1_fast,n_var]=size(dataset_1);
    N2_fast=size(dataset_2,1);
    A_r=n_var/d;
    A_c=d;
    
    A_fast_1=reshape(dataset_1',A_c,N1_fast*A_r)';
    b_fast_1=repmat(b,N1_fast,1);
        
    cvx_begin
    	variable x_fast_1(d)
    	minimize c_t*x_fast_1
        subject to
        	A_fast_1*x_fast_1 <= b_fast_1
        	x_fast_1 >= 0
	cvx_end
          
	A_fast_2=reshape(dataset_2',A_c,N2_fast*A_r)';
	b_fast_2=repmat(b,N2_fast,1);
        
          
	cvx_begin
        variable alpha_fast
        minimize c_t*( (1-alpha_fast)*x_fast_1  + alpha_fast*x_0   )
        subject to
        	A_fast_2*( (1-alpha_fast)*x_fast_1  + alpha_fast*x_0   ) <= b_fast_2
        	( (1-alpha_fast)*x_fast_1  + alpha_fast*x_0   ) >= 0
            alpha_fast >= 0
            alpha_fast <= 1
	cvx_end
          
	x_FAST=(1-alpha_fast)*x_fast_1  + alpha_fast*x_0 ;
end