function [x_FAST] = FAST_ccp(dataset_1,dataset_2,c,b,x_0)
% [x_FAST] = FAST_ccp(dataset_1,dataset_2,c,b,x_0)
% Using FAST algorithm.
% Provide solution for minimize c'x s.t. a'x <= b with two datasets for a
% and an initial solution x_0.
    c_t=c';
    m=length(c);
    
    cvx_begin
        variable x_1(m)
        minimize c_t*x_1
        subject to
            dataset_1*x_1 <= b
    cvx_end

    if strcmp(cvx_status,'Solved')
        cvx_begin
            variable alpha_x
            minimize c_t*( (1-alpha_x)*x_1  + alpha_x*x_0   )
            subject to
                dataset_2*( (1-alpha_x)*x_1  + alpha_x*x_0   ) <= b
                alpha_x >= 0
                alpha_x <= 1
        cvx_end
        x_FAST=(1-alpha_x)*x_1  + alpha_x*x_0;
    else
        x_FAST=nan(m,1);
    end
end

