clear  
    d=11;                 % dimension
    m=1;                 % number of constraints
    N_data=336;          % sample size
    n_outer=1000;        % outer test size
    N_test_data=10000;   % test size
    delta=0.05;
    epsilon=0.05;

    rng(123)
    % LP setting
    load('c_sigma_for_11.mat') % d11 in paper
    A=-c'; %
    [A_r A_c]=size(A);
    b=  [1200];

    % parameters for data 
    mu_0=A;
    [Cor_r,std_sigma]=corrcov(sigma+eye(d)*100); 
    nu_t=5;
    cov_mat=nu_t/(nu_t-2)*Cor_r;
    real_mean=mu_0;

    % setting for RO and Recon
    B_2=124;              % phase II budget
    B_1=N_data-B_2;      % phase I budget
    rank_of_data=binoinv(1-delta,B_2,1-epsilon); % estimated quantile
    rank_of_data_p1=binoinv(1-delta,B_1,1-epsilon); % estimated quantile for recon phase 1

    % setting for FAST
    N1_fast=318;
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
       
        dataset=mvtrnd(Cor_r,nu_t,N_data)+repmat(mu_0,N_data,1);
        
       %% FAST
       tic 
       dataset_fast_1=dataset(1:N1_fast,:);
       dataset_fast_2=dataset(N1_fast+1:end,:);
       [x_FAST] = FAST_ccp(dataset_fast_1,dataset_fast_2,c,b,x_fast_0);
       time_fast(i)=toc;
       fv_fast(i)=c'*x_FAST;
       
       %% SG
       tic
       A_gen=reshape(dataset',A_c,N_data*A_r)';
       [x_SG]=SG_ccp(A_gen,c,b);
       time_sg(i)=toc;
       fv_sg(i)=c'*x_SG;
       
       %% RO
       tic
       dataset_ro_1=dataset(1:B_1,:);
       dataset_ro_2=dataset(B_1+1:end,:);
       [x_RO] = RO_ccp(dataset_ro_1,dataset_ro_2,rank_of_data+1,c,b);
       time_ro(i)=toc;
       fv_ro(i)=c'*x_RO;
       
       %% Reconstructed RO
       tic
       dataset_recon_1=dataset(1:B_1,:);
       dataset_recon_2=dataset(B_1+1:end,:);
       [x_Recon] = Recon_ccp(dataset_recon_1,dataset_recon_2,rank_of_data_p1+1,rank_of_data+1,c,b);
       time_recon(i)=toc;
       fv_recon(i)=c'*x_Recon;
       
       %% Moment-based DRO
       tic
        [x_mo_DRO] = moment_DRO_ccp(dataset,c,b,epsilon,delta);
        time_mo_dro(i)=toc;
        fv_mo_dro(i)=c'*x_mo_DRO;
        
        
        %% violation test
        test_data=mvtrnd(Cor_r,nu_t,N_test_data)+repmat(mu_0,N_test_data,1);
        
        violate_num_sg=0;
        violate_num_recon=0;
        violate_num_ro=0;
        violate_num_mo_dro=0;
        violate_num_fast=0;
       
        for j=1:N_test_data
           A_test=reshape(test_data(j,:),A_c,A_r)';
           violate_num_sg=violate_num_sg+(sum(A_test*x_SG-b >= 0)>0);
           violate_num_recon=violate_num_recon+(sum(A_test*x_Recon-b >= 0)>0);
           violate_num_ro=violate_num_ro+(sum(A_test*x_RO-b >= 0)>0);
           violate_num_mo_dro=violate_num_mo_dro+(sum(A_test*x_mo_DRO-b >= 0)>0);
           violate_num_fast=violate_num_fast+(sum(A_test*x_FAST-b >= 0)>0);
        end
       
        violation_fast(i)=violate_num_fast/N_test_data;
        violation_sg(i)=violate_num_sg/N_test_data;
        violation_ro(i)=violate_num_ro/N_test_data;
        violation_recon(i)=violate_num_recon/N_test_data;
        violation_mo_dro(i)=violate_num_mo_dro/N_test_data;
    end
   
  
    result_table=cell(7,6);
    result_table(1,:)={'','RO','Recon','SG','FAST','DRO Mo'};
    result_table(2,:)={'n',N_data,N_data,N_data,N_data,N_data};
    result_table(3,:)={'n1',B_1,B_1,'-',N1_fast,'-'};
    result_table(4,:)={'n2',B_2,B_2,'-',N2_fast,'-'};
    result_table(5,:)={'ov',mean(fv_ro),mean(fv_recon),mean(fv_sg),mean(fv_fast),mean(fv_mo_dro)};
    result_table(6,:)={'eps',mean(violation_ro),mean(violation_recon),mean(violation_sg),mean(violation_fast),mean(violation_mo_dro)};
    result_table(7,:)={'delta',sum(violation_ro>delta)/n_outer,sum(violation_recon>delta)/n_outer,sum(violation_sg>delta)/n_outer,sum(violation_fast>delta)/n_outer,sum(violation_mo_dro>delta)/n_outer};
    disp('Results')
    disp(result_table)
        
    computation_time=[mean(time_ro),mean(time_recon),mean(time_sg),mean(time_fast),mean(time_mo_dro)];
    disp('Average Computation Time for RO,Recon,SG,FAST,DRO_Mo:')
    disp(computation_time)

