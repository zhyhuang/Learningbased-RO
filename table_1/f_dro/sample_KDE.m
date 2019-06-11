function [samples] = sample_KDE(data,bw,N_SAA_data)

[N_data,n_dim]=size(data);

centers_ind=randi(N_data,N_SAA_data,1);

centers=data(centers_ind,:);


samples=mvnrnd(centers,eye(n_dim)*bw^2);

end