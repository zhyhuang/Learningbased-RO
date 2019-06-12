    d=100;                 % dimension
    m=1;                 % number of constraints
    N_data=2331;          % sample size
    n_outer=1000;        % outer test size
    delta=0.05;
    epsilon=0.05;

    rng(123)
    % LP setting
    load('c_sigma_100.mat') % d11 in paper
    A=-c'; %
    [A_r A_c]=size(A);
    b=  [1200];

    % parameters for data 
    miu_0=A;

    % setting for RO and Recon
    B_2=1013;              % phase II budget
    B_1=N_data-B_2;      % phase I budget
    rank_of_data=binoinv(1-delta,B_2,1-epsilon); % estimated quantile
    rank_of_data_p1=binoinv(1-delta,B_1,1-epsilon); % estimated quantile for recon phase 1

    % setting for FAST
    N1_fast=2326;
    N2_fast=N_data-N1_fast;
    x_fast_0=zeros(d,1);

    % result record
    fv_fast=zeros(n_outer,1);
    fv_sg=zeros(n_outer,1);
    fv_ro=zeros(n_outer,1);
    fv_recon=zeros(n_outer,1);
    fv_mo_dro=zeros(n_outer,1);

    violation_fast=zeros(n_outer,1);
    violation_sg=zeros(n_outer,1);
    violation_ro=zeros(n_outer,1);
    violation_recon=zeros(n_outer,1);
    violation_mo_dro=zeros(n_outer,1);
    
    time_fast=zeros(n_outer,1);
    time_sg=zeros(n_outer,1);
    time_ro=zeros(n_outer,1);
    time_recon=zeros(n_outer,1);
    time_mo_dro=zeros(n_outer,1);

    for i=1:n_outer
       dataset=mvnrnd(miu_0,sigma,N_data);
        
       %% FAST
       tic 
       dataset_fast_1=dataset(1:N1_fast,:);
       dataset_fast_2=dataset(N1_fast+1:end,:);
       [x_FAST] = FAST_ccp(dataset_fast_1,dataset_fast_2,c,b,x_fast_0);
       time_fast(i)=toc;
       fv_fast(i)=c_t*x_FAST;
       violation_fast(i)=1-normcdf((b-A*x_FAST)/norm(sqrtm(sigma)*x_FAST));
       %% SG
       tic
       A_gen=reshape(dataset',A_c,N_data*A_r)';
       [x_SG]=SG_ccp(A_gen,c,b);
       time_sg(i)=toc;
       fv_sg(i)=c_t*x_SG;
       violation_sg(i)=1-normcdf((b-A*x_SG)/norm(sqrtm(sigma)*x_SG));
       %% RO
       tic
       dataset_ro_1=dataset(1:B_1,:);
       dataset_ro_2=dataset(B_1+1:end,:);
       [x_RO] = RO_ccp(dataset_ro_1,dataset_ro_2,rank_of_data+1,c,b);
       time_ro(i)=toc;
       fv_ro(i)=c_t*x_RO;
       violation_ro(i)=1-normcdf((b-A*x_RO)/norm(sqrtm(sigma)*x_RO));
       %% Reconstructed RO
       tic
       dataset_recon_1=dataset(1:B_1,:);
       dataset_recon_2=dataset(B_1+1:end,:);
       [x_Recon] = Recon_ccp(dataset_recon_1,dataset_recon_2,rank_of_data_p1+1,rank_of_data+1,c,b);
       time_recon(i)=toc;
       fv_recon(i)=c_t*x_Recon;
       violation_recon(i)=1-normcdf((b-A*x_Recon)/norm(sqrtm(sigma)*x_Recon));
       %% Moment-based DRO
       tic
        [x_mo_DRO] = moment_DRO_ccp(dataset,c,b,epsilon,delta)
        time_mo_dro(i)=toc;
        fv_mo_dro(i)=c_t*x_mo_DRO;
        violation_mo_dro(i)=1-normcdf((b-A*x_mo_DRO)/norm(sqrtm(sigma)*x_mo_DRO));
    end
    % Safe Convex Approximation
    [x_SCA] = SCA_ccp(c,b,mu,sigma,epsilon)
    fv_sca=c'*x_SCA;
    violation_sca=1-normcdf((b-A*x_SCA)/norm(sqrtm(sigma)*x_SCA));
    
    % ture solution
    phi_quantile=norminv(1-epsilon,0,1);
    rt_sigma=sqrtm(sigma);     
    [x_true] =cvx_closed_one_line(c,phi_quantile,rt_sigma,miu_0,b);
    fv_true=c'*x_true;
    
    result_table=cell(7,7);
    result_table(1,:)={'','RO','Recon','SG','FAST','DRO Mo','SCA'};
    result_table(2,:)={'n',N_data,N_data,N_data,N_data,N_data,'-'};
    result_table(3,:)={'n1',B_1,B_1,'-',N1_fast,'-','-'};
    result_table(4,:)={'n2',B_2,B_2,'-',N2_fast,'-','-'};
    result_table(5,:)={'ov',mean(fv_ro),mean(fv_recon),mean(fv_sg),mean(fv_fast),mean(fv_mo_dro),fv_sca};
    result_table(6,:)={'eps',mean(violation_ro),mean(violation_recon),mean(violation_sg),mean(violation_fast),mean(violation_mo_dro),violation_sca};
    result_table(7,:)={'delta',sum(violation_ro>delta)/n_outer,sum(violation_recon>delta)/n_outer,sum(violation_sg>delta)/n_outer,sum(violation_fast>delta)/n_outer,sum(violation_mo_dro>delta)/n_outer,0};
    disp('Results')
    disp(result_table)
    
    disp('True solution')
    disp(fv_true)
    
    computation_time=[mean(time_ro),mean(time_recon),mean(time_sg),mean(time_fast),mean(time_mo_dro)];
    disp('Average Computation Time')
    disp(computation_time)
