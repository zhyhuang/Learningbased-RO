function [x_SCA] = SCA_ccp(c,b,mu_0,sigma,epsilon)
%[x_SCA] = SCA_ccp(c,b,mu,sigma,epsilon)
% safe convex approximation
rt_sigma=sqrtm(sigma);  
phi_quantile=sqrt(2*log(1/epsilon));
d=length(c);
cvx_begin
    variable x_SCA(d);
    minimize(c'*x_SCA)
    subject to
   
    phi_quantile*norm(rt_sigma*x_SCA)-b+mu_0*x_SCA <= 0;
 
cvx_end
end

