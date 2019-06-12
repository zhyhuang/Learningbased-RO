function [x_SG] = SG_ccp(dataset,c,b)
% [x_SG] = SG_ccp(dataset,c,b)
% SG algorithm.

    c_t=c';
    d=length(c);
    cvx_begin
        variable x_SG(d)
        minimize c_t*x_SG
            subject to
            dataset*x_SG <= b
	cvx_end
    
    if ~strcmp(cvx_status,'Solved')
        x_SG=nan(d,1);
    end
end

