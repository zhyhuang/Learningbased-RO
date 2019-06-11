clear
% parameter setting

n=100;                 % dimension
m=1;                 % number of constraints
N_data=2331;          % sample size
N_test_data=10000;   % test size
B_2=1013;              % phase II budget
B_1=N_data-B_2;      % phase I budget
n_outer=1000;        % outer test size

% N1_fast=61;
% N2_fast=59;



% N1_fast=318;
% N2_fast=N_data-N1_fast;

N1_fast=2326;
N2_fast=N_data-N1_fast;

x_fast_0=zeros(n,1);

rng(123)

alpha=0.05;
beta=0.05;
rank_of_data=binoinv(1-beta,B_2,1-alpha); % estimated quantile
rank_of_data_p1=binoinv(1-beta,B_1,1-alpha);
% LP setting
% load('c_sigma_for_11.mat') % d11 in paper
load('c_sigma_100_set2.mat')

A1=-c'; %1
c_t=c';
[A1_r A1_c]=size(A1);
b1=  [1200];
A1_tsps=A1';
miu_0=A1_tsps(:)';   %%% NOTE here!!!!!!!!!!!
  
  
obj_fun=@(x) c'*x;

    
fv_el=zeros(n_outer,1);
fv_gen=zeros(n_outer,1);
gen_vio=zeros(n_outer,1);
el_vio=zeros(n_outer,1);
fv_jp=zeros(n_outer,1);
jp_vio=zeros(n_outer,1);
scb_fv=zeros(n_outer,1);
scb_vio=zeros(n_outer,1);
dro_est_fv=zeros(n_outer,1);
dro_est_vio=zeros(n_outer,1);
fv_FAST=zeros(n_outer,1);
vio_FAST=zeros(n_outer,1);

time_fast=zeros(n_outer,1);
time_sg=zeros(n_outer,1);
time_ro=zeros(n_outer,1);
time_recon=zeros(n_outer,1);

for i=1:n_outer
  
  
dataset_points_1=mvnrnd(miu_0,sigma,N_data);


 % seperate data for phase I and II

        
        
        
        %% FAST x_fast_0
       tic 
         dataset_fast_1=dataset_points_1(1:N1_fast,:);
        dataset_fast_2=dataset_points_1(N1_fast+1:end,:);
        
          cvx_begin
                variable x_fast_1(A1_c)
                minimize c_t*x_fast_1
                subject to
                    dataset_fast_1*x_fast_1 <= b1
          cvx_end
          
          if strcmp(cvx_status,'Solved')
                cvx_begin
                    variable alpha_fast
                    minimize c_t*( (1-alpha_fast)*x_fast_1  + alpha_fast*x_fast_0   )
                    subject to
                        dataset_fast_2*( (1-alpha_fast)*x_fast_1  + alpha_fast*x_fast_0   ) <= b1
                        alpha_fast >= 0
                        alpha_fast <= 1
              cvx_end

              x_FAST=(1-alpha_fast)*x_fast_1  + alpha_fast*x_fast_0 ;
              fv_FAST(i)=c_t*x_FAST;
              vio_FAST(i)=1-normcdf((b1-A1*x_FAST)/norm(sqrtm(sigma)*x_FAST));
          else
                fv_FAST(i)=nan;
          end
          time_fast(i)=toc;
%           [x_FAST_2] = FAST_ccp(dataset_fast_1,dataset_fast_2,c,b1,x_fast_0);
          
          
        
           %% ellipsoidal set for violated data
           tic
           
        dataset_p1_1=dataset_points_1(1:B_1,:);
        dataset_p2_1=dataset_points_1(B_1+1:end,:);
        % decide size
        sample_mean_1=mean(dataset_p1_1);
        
%         sample_cov_1=sigma;
        sample_cov_1=diag(diag(cov(dataset_p1_1)));
%         sample_cov_1=((cov(dataset_p1_1)));
         mean_mat_b2_1=repmat(sample_mean_1,B_2,1);  
              T_value_b2_elli=(dataset_p2_1-mean_mat_b2_1)*sample_cov_1^-1*(dataset_p2_1-mean_mat_b2_1)';  % calculate only for diag
        t_value_b2_elli=sqrt(diag(T_value_b2_elli));
      [sorted_t_old,sort_indexold]=sort(t_value_b2_elli);
        s_old=sorted_t_old(rank_of_data+1);
         A_m=reshape(sample_mean_1,A1_c,A1_r)';
       
       M_total=sqrtm(sample_cov_1);
       [M_r,M_c]=size(M_total);
       M_spl=zeros(M_r,A1_c,A1_r);
       for j=1:A1_r
         M_spl(:,:,j)=M_total(:,(j-1)*A1_c+1:j*A1_c);
       end
       
      [ cvx_status1,f_val_old,x_elli_old] = cvx_elli_LP_one_line( c,(s_old)^2,M_spl,A_m,b1,n);  
      time_ro(i)=toc;
%        sample_cov_1=eye(m*n); 

%         sample_cov_2=eye(m*n); 
%         sample_cov_2=diag(diag(cov(dataset_p1_1)));
%% recon
           tic
           
        dataset_p1_1=dataset_points_1(1:B_1,:);
        dataset_p2_1=dataset_points_1(B_1+1:end,:);
        sample_mean_1=mean(dataset_p1_1);
        sample_cov_1=((cov(dataset_p1_1)));

        mean_mat_b1_for_scale=repmat(sample_mean_1,B_1,1);
       
        T_value_b1_for_cut=(dataset_p1_1-mean_mat_b1_for_scale)*sample_cov_1^-1*(dataset_p1_1-mean_mat_b1_for_scale)';  % calculate only for diag
        t_value_b1_for_cut=diag(T_value_b1_for_cut);

        [sorted_t_for_cut,sort_index_t_for_cut]=sort(t_value_b1_for_cut);

        s1_1_for_solve_x0=sorted_t_for_cut(rank_of_data_p1+1);
   

       % constraints for ellipsoid

       A_m=reshape(sample_mean_1,A1_c,A1_r)';
       
       M_total=sqrtm(sample_cov_1);
       [M_r,M_c]=size(M_total);
       M_spl=zeros(M_r,A1_c,A1_r);
       for j=1:A1_r
         M_spl(:,:,j)=M_total(:,(j-1)*A1_c+1:j*A1_c);
       end
       
                          
       [ cvx_status1 f_val_cut x0_for_cut] = cvx_elli_LP_one_line( c,s1_1_for_solve_x0,M_spl,A_m,b1,n);          
       
       

%         s1_1_for_cut=sorted_t_for_cut(binoinv(1-beta,B_1,1-alpha));
        
%         mean_mat_b2_1=repmat(sample_mean_1,B_2,1);  
         T_j=   dataset_p2_1* x0_for_cut;
       [t_j_plane_sort,j_ind]=     sort(T_j);
        t_jp=t_j_plane_sort(rank_of_data+1);
        
       [ cvx_status1 f_val_jp x_jp] = cvx_jp_LP_one_line( c,t_jp,x0_for_cut,b1,n); 
       
       time_recon(i)=toc;
%        A_cutline=x0_for_cut'/(x0_for_cut' *in_elli_cut_data(cut_point,:)'-x0_for_cut'*sample_mean_1');
       
%        t_cut_line=(dataset_p2_1-mean_mat_b2_1)*A_cutline';
       
       
       
%        ell_built_fun=@(x,y) ([x y]-sample_mean_1)*sample_cov_1^-1*([x y]-sample_mean_1)'/s1_1_for_cut-1;
%        line_built_fun=@(x,y) A_cutline*([x y]-sample_mean_1)'-1;
%        
%        
 
       
 
%       
%         
%         T_value_b2_scb=(dataset_p2_1-mean_mat_b2_1)*sample_cov_2^-1*(dataset_p2_1-mean_mat_b2_1)'/s1_1_for_cut;  % calculate only for diag
%         t_value_b2_scb=sqrt(diag(T_value_b2_scb));
       
%        scale_r=1;
%        t_combine=max(t_cut_line*scale_r,t_value_b2_elli);
       
%         [sorted_t_new,sort_indexnew]=sort(t_combine);
%        s_combine=sorted_t_new(rank_of_data+1);
  
    
        
        %%  point ball set
%         t_value_mat_pbs=zeros(B_2,B_1);
%     for ii=1:B_1
%         mean_c1=dataset_p1_1(ii,:);
%         mean_mat_b2_c1=repmat(mean_c1,B_2,1);
%         T_v_c1=(dataset_p2_1-mean_mat_b2_c1)*(dataset_p2_1-mean_mat_b2_c1)';
%         t_value_mat_pbs(:,ii)=diag(T_v_c1);
%         
%         
%     end
%    
%     t_value_pbs=min(t_value_mat_pbs');  % calculate only for diag
% 
%     [sorted_t_pbs,sort_index]=sort(t_value_pbs);
%     
%     
%     s_pbs=sorted_t_pbs(rank_of_data+1);
        
        
        
        
        
  %% j halfplane
       
        
        %%  compare fast
        
%         T_j=   dataset_fast_2* x_fast_1;
%        [t_j_plane_sort,j_ind]=     sort(T_j);
%         t_jp=t_j_plane_sort(rank_of_data_p1+1);
%         
%        [ cvx_status1fast f_val_jpfast x_jpjpast] = cvx_jp_LP_one_line( c,t_jp,x_fast_1,b1,n); 
%         
%         fv_compare_fast(i)=f_val_jpfast;
%         vio_compare_fast(i)=1-normcdf((b1-A1*x_jpjpast)/norm(sqrtm(sigma)*x_jpjpast));
        
       %% solve
       
       
       
%          A_m=reshape(sample_mean_1,A1_c,A1_r)';
%        
%        M_total=sqrtm(sample_cov_1);
%        [M_r,M_c]=size(M_total);
%        M_spl=zeros(M_r,A1_c,A1_r);
%        for j=1:A1_r
%          M_spl(:,:,j)=M_total(:,(j-1)*A1_c+1:j*A1_c);
%        end

% for one line actually     M_spl=sqrtm(sample_cov_1)  A_m=sample_mean_1
% n=2 here
       
%          [ cvx_status1 scb_fv(i) x_elli_scb] = cvx_elli_LP_one_line( c,(s_scb)^2*s1_1_for_cut,M_spl_scb,A_m,b1,n);  
        
    


%         [ cvx_status1 f_val_cut x_cut ] = cvx_cut_LP_one_line( c,s_combine,sample_cov_1,s1_1_for_cut,A_cutline,sample_mean_1,b1,n);
        
         % generating
         tic
         A_gen_1=reshape(dataset_points_1',A1_c,N_data*A1_r)';
         b_gen_1=repmat(b1,N_data,1);
         
         A_gen=[A_gen_1];
         b_gen=[b_gen_1];
         
         
        
          
          cvx_begin
                variable x_gen(A1_c)
                minimize c_t*x_gen
                subject to
                A_gen*x_gen <= b1
          cvx_end
          
          
          fv_el(i)=f_val_old;
        if strcmp(cvx_status,'Solved')
        fv_gen(i)=obj_fun(x_gen);
        else
            fv_gen(i)=nan;
        end
        time_sg(i)=toc;
         %% dro
%     estimated_mean=mean(dataset_points_1);
%     estimated_covariance=cov(dataset_points_1);
%     c_t=c';
%         cvx_begin
%                 variable x_dro_el(A1_c)
%                 minimize c_t*x_dro_el
%                 subject to
%                 sqrt((1-alpha)/alpha) *norm(sqrtm(estimated_covariance)*x_dro_el)+ estimated_mean* x_dro_el -b1 <=0;
%         cvx_end
%     dro_est_fv(i)=obj_fun(x_dro_el);
%     dro_est_vio(i)=1-normcdf((b1-A1*x_dro_el)/norm(sqrtm(sigma)*x_dro_el));
   %%  
       
   
   
   
    
         
    gen_vio(i)=1-normcdf((b1-A1*x_gen)/norm(sqrtm(sigma)*x_gen));
    el_vio(i)=1-normcdf((b1-A1*x_elli_old)/norm(sqrtm(sigma)*x_elli_old));
%     
%     fv_cut(i)=f_val_cut;
%     cut_vio(i)=1-normcdf((b1-A1*x_cut)/norm(sqrtm(sigma)*x_cut));
%     
    fv_jp(i)=f_val_jp;
    jp_vio(i)=1-normcdf((b1-A1*x_jp)/norm(sqrtm(sigma)*x_jp));
    
%     [ cvx_status pointball_fv(i) x_pbs ] = cvx_elli_LP_bps_lines( c,s_pbs,dataset_p1_1,b1,n,B_1);
    
    
%     pointball_vio(i)=1-normcdf((b1-A1*x_pbs)/norm(sqrtm(sigma)*x_pbs));
%     scb_vio(i)=1-normcdf((b1-A1*x_elli_scb)/norm(sqrtm(sigma)*x_elli_scb));
%          fv_el=zeros(n_outer,1);
% fv_gen=zeros(n_outer,1);
% gen_vio=zeros(n_outer,1);
% el_vio=zeros(n_outer,1);
% % 
% fv_cut=zeros(n_outer,1);
% fv_jp=zeros(n_outer,1);
% cut_vio=zeros(n_outer,1);
% jp_vio=zeros(n_outer,1);
%          
%         
%         f_val_old
%         f_val_cut
%         f_val_jp
end    
        
 
ov_eli=mean(fv_el)
ov_eli_std=std(fv_el)
vio_eil=mean(el_vio)
vio_eil_std=std(el_vio)
elli_beta_est=sum(el_vio>0.05)/1000
elli_beta_ci=[elli_beta_est-1.98*sqrt(elli_beta_est*(1-elli_beta_est)/n_outer) elli_beta_est+1.98*sqrt(elli_beta_est*(1-elli_beta_est)/n_outer)];






% sg for bounded
% ov_gen=mean(fv_gen)
% ov_gen_std=std(fv_gen)
% vio_gen=mean(gen_vio)
% vio_gen_std=std(gen_vio)
% gen_beta_est=sum(gen_vio>0.05)/1000      
% gen_beta_ci=[gen_beta_est-1.98*sqrt(gen_beta_est*(1-gen_beta_est)/n_outer) gen_beta_est+1.98*sqrt(gen_beta_est*(1-gen_beta_est)/n_outer)];

% scenerio without unbounded
ov_gen=mean(fv_gen(fv_gen>-2e5))
ov_gen_std=std(fv_gen(fv_gen>-2e5))
vio_gen=mean(gen_vio(fv_gen>-2e5))
vio_gen_std=std(gen_vio(fv_gen>-2e5))
gen_beta_est=sum(gen_vio(fv_gen>-2e5)>0.05)/sum((fv_gen>-2e5))      
gen_beta_ci=[gen_beta_est-1.98*sqrt(gen_beta_est*(1-gen_beta_est)/sum((fv_gen>-2e5))  ) gen_beta_est+1.98*sqrt(gen_beta_est*(1-gen_beta_est)/sum((fv_gen>-2e5))  )];
  



% 
% ov_cut=mean(fv_cut)
% ov_cut_std=std(fv_cut)
% vio_cut=mean(cut_vio)
% vio_cut_std=std(cut_vio)
% cut_beta_est=sum(cut_vio>0.05)/1000
% cut_beta_ci=[cut_beta_est-1.98*sqrt(cut_beta_est*(1-cut_beta_est)/n_outer) cut_beta_est+1.98*sqrt(cut_beta_est*(1-cut_beta_est)/n_outer)];

ov_jp=mean(fv_jp)
ov_jp_std=std(fv_jp)
vio_jp=mean(jp_vio)
vio_jp_std=std(jp_vio)
jp_beta_est=sum(jp_vio>0.05)/1000  
jp_beta_ci=[jp_beta_est-1.98*sqrt(jp_beta_est*(1-jp_beta_est)/n_outer) jp_beta_est+1.98*sqrt(jp_beta_est*(1-jp_beta_est)/n_outer)]


% ov_pbs=mean(pointball_fv)
% ov_pbs_std=std(pointball_fv)
% vio_pbs=mean(pointball_vio)
% vio_pbs_std=std(pointball_vio)
% pbs_beta_est=sum(pointball_vio>0.05)/1000  
% pbs_beta_ci=[pbs_beta_est-1.98*sqrt(pbs_beta_est*(1-pbs_beta_est)/n_outer) pbs_beta_est+1.98*sqrt(pbs_beta_est*(1-pbs_beta_est)/n_outer)]

ov_scb=mean(scb_fv)
ov_scb_std=std(scb_fv)
vio_scb=mean(scb_vio)
vio_scb_std=std(scb_vio)
scb_beta_est=sum(scb_vio>0.05)/1000  
scb_beta_ci=[scb_beta_est-1.98*sqrt(scb_beta_est*(1-scb_beta_est)/n_outer) scb_beta_est+1.98*sqrt(scb_beta_est*(1-scb_beta_est)/n_outer)]

ov_FAST=mean(fv_FAST)
ov_FAST_std=std(fv_FAST)
vio_FASTtable=mean(vio_FAST)
vio_FASTtable_std=std(vio_FAST)
FAST_beta_est=sum(vio_FAST>0.05)/1000  
FAST_beta_ci=[FAST_beta_est-1.98*sqrt(FAST_beta_est*(1-FAST_beta_est)/n_outer) FAST_beta_est+1.98*sqrt(FAST_beta_est*(1-FAST_beta_est)/n_outer)]


% ov_compare_fast=mean(fv_compare_fast)
% ov_compare_fast_std=std(fv_compare_fast)
% vio_compare_fast=mean(vio_compare_fast)
% vio_compare_fast_std=std(vio_compare_fast)
% compare_fast_beta_est=sum(vio_compare_fast>0.05)/1000  
% compare_fast_beta_ci=[compare_fast_beta_est-1.98*sqrt(compare_fast_beta_est*(1-compare_fast_beta_est)/n_outer) compare_fast_beta_est+1.98*sqrt(compare_fast_beta_est*(1-compare_fast_beta_est)/n_outer)]


ov_dro_est=mean(dro_est_fv)
ov_dro_est_std=std(dro_est_fv)
vio_dro_est=mean(dro_est_vio)
vio_dro_est_std=std(dro_est_vio)
dro_est_beta_est=sum(dro_est_vio>0.05)/1000  
dro_est_beta_ci=[dro_est_beta_est-1.98*sqrt(dro_est_beta_est*(1-dro_est_beta_est)/n_outer) dro_est_beta_est+1.98*sqrt(dro_est_beta_est*(1-dro_est_beta_est)/n_outer)]




result_table=[
    ov_gen,ov_gen_std,vio_gen,vio_gen_std,gen_beta_est, gen_beta_ci ; 
    ov_eli,ov_eli_std,vio_eil,vio_eil_std,elli_beta_est,elli_beta_ci;
%     ov_cut,ov_cut_std,vio_cut,vio_cut_std,cut_beta_est,cut_beta_ci;
    ov_jp,ov_jp_std,vio_jp,vio_jp_std,jp_beta_est, jp_beta_ci;
%     ov_pbs,ov_pbs_std,vio_pbs,vio_pbs_std,pbs_beta_est,pbs_beta_ci;
    ov_scb,ov_scb_std,vio_scb,vio_scb_std,scb_beta_est,scb_beta_ci;
    ov_dro_est,ov_dro_est_std,vio_dro_est,vio_dro_est_std,dro_est_beta_est,dro_est_beta_ci;
    ov_FAST,ov_FAST_std,vio_FASTtable,vio_FASTtable_std,FAST_beta_est,FAST_beta_ci]';
%     ov_compare_fast,ov_compare_fast_std,vio_compare_fast,vio_compare_fast_std,compare_fast_beta_est,compare_fast_beta_ci]';


sum((fv_gen>-2e5)) 


phi_quantile=norminv(1-alpha,0,1);
rt_sigma=sqrtm(sigma);
     
[ cvx_status cvx_optval x_closed ] = cvx_closed_one_line( c,phi_quantile,rt_sigma,miu_0,b1,n);   
cvx_optval


% approximation 1
phi_quantile_2=sqrt(2*log(1/alpha));
[ cvx_status cvx_optval_2 x_closed_2 ] = cvx_closed_one_line( c,phi_quantile_2,((rt_sigma)),miu_0,b1,n);   
cvx_optval_2
1-normcdf((b1-A1*x_closed_2)/norm(sqrtm(sigma)*x_closed_2))

% % approximation 2
% phi_quantile_3=sqrt(2*log(1/alpha));
% [ cvx_status cvx_optval_3 x_closed_3 ] = cvx_closed_one_line( c,phi_quantile_3,dataset_points_1(1:120,:),miu_0,b1,n);   
% cvx_optval_3
% 1-normcdf((b1-A1*x_closed_3)/norm(sqrtm(sigma)*x_closed_3))



% scenerio without unbounded
% ov_gen=mean(fv_gen(fv_gen>-2e5))
% ov_gen_std=std(fv_gen(fv_gen>-2e5))
% vio_gen=mean(gen_vio(fv_gen>-2e5))
% vio_gen_std=std(gen_vio(fv_gen>-2e5))
% gen_beta_est=sum(gen_vio(fv_gen>-2e5)>0.05)/sum((fv_gen>-2e5))      
% gen_beta_ci=[gen_beta_est-1.98*sqrt(gen_beta_est*(1-gen_beta_est)/sum((fv_gen>-2e5))  ) gen_beta_est+1.98*sqrt(gen_beta_est*(1-gen_beta_est)/sum((fv_gen>-2e5))  )];
% sum((fv_gen>-2e5))   
% 


        cvx_begin
                variable x_dro_el(A1_c)
                minimize c_t*x_dro_el
                subject to
                sqrt((1-alpha)/alpha) *norm(sqrtm(sigma)*x_dro_el)+ miu_0* x_dro_el -b1 <=0;
        cvx_end
    dro_el_fv=obj_fun(x_dro_el);
    dro_el_vio=1-normcdf((b1-A1*x_dro_el)/norm(sqrtm(sigma)*x_dro_el));








