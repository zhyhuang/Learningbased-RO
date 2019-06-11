function [density] = Gaus_KDE(data,pts,bw)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[N_data,n_dim]=size(data);

celled_pts=num2cell(pts, 2);

point_density=@(x) 1/N_data/(bw^n_dim)*sum(mvnpdf(x/(bw),data/(bw),eye(n_dim)));

%test_density=@(x) 1/N_data*sum(mvnpdf(x,data,eye(n_dim).*(bw)^2));

density=cell2mat(cellfun(point_density,celled_pts,'UniformOutput',false));

%test_dddensity=cell2mat(cellfun(test_density,celled_pts,'UniformOutput',false));
end

