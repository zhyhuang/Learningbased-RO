clear
    d=11;                 % dimension
    m=15;                 % number of constraints
    N_data=120;          % sample size
    n_outer=1000;        % outer test size
    delta=0.05;
    epsilon=0.05;
    N_test_data=10000;   % test size

    rng(123)
    % LP setting
    c=[-139 -88 -133 -137 -165    -80     -120    -150      -100     -150     -110    ]';

    A=[ 14    18    15   14   14   10     15    18      26     15    34  ; %1
        98    24    14   14   14   0    0    30      50      0   25   ;
        40    4    14   14   14   0    40    0     30      0     30  ;
        28    6   13   14    15  55   0    20     40      40     0 ;                                       
        22   4    13    15   12   0    55    0   44      40     30;
        22   4   13    13   12   0     0     55   0      0     0;
        22   6    13   13   12   0      10     0    30    0   0;
        26   58    15    15   12  0    0     10      0     30  0;
        6   40    13    0      0    30    0     0    30      0     45;
        6   12     13   13    26  0   0       0   40       0   0   ;
        6    18   12    12   12   0   20       0   40       40   50   ;
        6    18    15   14   14   0   0       0   10       0   0   ;
        16      0     0     12    35   50    50      0   50       0   20   ;
        20    0       10      36   50   0     0    50    50     0   0   ;
        16     0      0       12   35  0     10     0     0     50    30];

    [A_r A_c]=size(A);

    b=  [1200;     
         2150;
         1550;
         1450;
         1600;
         950;
         1800;
         1950;
         1700;
         2300;
         1100;
         1200;
         1200;
         1200;
         3600];
     
     
    % parameters for data 
    A_tsps=A';
    mu_0=A_tsps(:)';
    L=165;
    A_l=[eye(d*m)*2+rand(d*m)*0.1; rand(L-d*m,d*m)]*100;
    beta_dist_a=10;
    beta_dist_b=10;
    
    % setting for RO and Recon
    B_2=60;              % phase II budget
    B_1=N_data-B_2;      % phase I budget
    rank_of_data=binoinv(1-delta,B_2,1-epsilon); % estimated quantile
    rank_of_data_p1=binoinv(1-delta,B_1,1-epsilon); % estimated quantile for recon phase 1

    % setting for FAST
    N1_fast=61;
    N2_fast=59;
    x_fast_0=zeros(d,1);
    
    rng(123)
     
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
        beta_rnd=betarnd(beta_dist_a,beta_dist_b,N_data,L)*2-1;

        al_mat=zeros(N_data,d*m,L);
            for l_i=1:L

                al_mat(:,:,l_i)=beta_rnd(:,l_i)*A_l(l_i,:);

            end

        al_sum=sum(al_mat,3);
        A1_mat=repmat(mu_0,N_data,1);
       dataset=A1_mat+al_sum;
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
       b_gen=repmat(b,N_data,1);
       [x_SG]=SG_ccp(A_gen,c,b_gen);
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
%        tic
%        [x_mo_DRO] = moment_DRO_ccp(dataset,c,b,epsilon,delta);
%        time_mo_dro(i)=toc;
%        fv_mo_dro(i)=c'*x_mo_DRO;
        
        
       %% violation test
       beta_rnd_test=betarnd(beta_dist_a,beta_dist_b,N_test_data,L)*2-1;

        al_mat_test=zeros(N_test_data,d*m,L);
            for l_i_test=1:L

                al_mat_test(:,:,l_i_test)=beta_rnd_test(:,l_i_test)*A_l(l_i_test,:);

            end

        al_sum_test=sum(al_mat_test,3);
        A1_mat_test=repmat(mu_0,N_test_data,1);
       test_data=A1_mat_test+al_sum_test;
       
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
%           violate_num_mo_dro=violate_num_mo_dro+(sum(A_test*x_mo_DRO-b >= 0)>0);
          violate_num_fast=violate_num_fast+(sum(A_test*x_FAST-b >= 0)>0);
       end
       
       violation_fast(i)=violate_num_fast/N_test_data;
       violation_sg(i)=violate_num_sg/N_test_data;
       violation_ro(i)=violate_num_ro/N_test_data;
       violation_recon(i)=violate_num_recon/N_test_data;
       violation_mo_dro(i)=violate_num_mo_dro/N_test_data;
       
       
       
       
    end
    
    
     
     
     
     
