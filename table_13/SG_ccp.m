function [x_SG]=SG_ccp(A_gen,c,b_gen)


    c_t=c';
    d=length(c);
	cvx_begin
        variable x_SG(d)
        minimize c_t*x_SG
        subject to
            A_gen*x_SG <= b_gen
            x_SG>=0
	cvx_end
    
end

