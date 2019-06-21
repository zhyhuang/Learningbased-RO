function [x_Recon] = Recon_ccp(dataset_1,dataset_2,data_rank_p1,data_rank,c,b);


    [n1,num_rv]=size(dataset_1);
    n2=size(dataset_2,1);
    d=length(c);
    A_c=d;
    A_r=num_rv/d;

    
    
    % phase 1
    sample_mean=mean(dataset_1);

%     dataset_row_split=zeros(n1+n2,A_c,A_r);
    dataset_p1_row_split=zeros(n1,A_c,A_r);
    dataset_p2_row_split=zeros(n2,A_c,A_r);
    cov_row_split=zeros(A_c,A_c,A_r);
    mean_row_split=zeros(A_c,A_r);

%     s_boffer_cov=zeros(A_r,1);

    t_value_cov=zeros(n1,A_r);

    cov_cov=[];
    for i_row=1:A_r
%             dataset_row_split(:,:,i_row)=dataset_1(:,((i_row-1)*A_c+1):(i_row*A_c));
            dataset_p1_row_split(:,:,i_row)=dataset_1(:,((i_row-1)*A_c+1):(i_row*A_c));
            dataset_p2_row_split(:,:,i_row)=dataset_2(:,((i_row-1)*A_c+1):(i_row*A_c));

            mean_row_split(:,i_row)=mean(dataset_p1_row_split(:,:,i_row))';
            cov_row_split(:,:,i_row)=((cov(dataset_p1_row_split(:,:,i_row))));


            mean_mat=repmat(mean_row_split(:,i_row)',n1,1);

            T_value_mat_Split_cov=(dataset_p1_row_split(:,:,i_row)-mean_mat)*cov_row_split(:,:,i_row)^-1*(dataset_p1_row_split(:,:,i_row)-mean_mat)';  % calculate only for diag
         
            t_value_cov(:,i_row)=diag(T_value_mat_Split_cov);%%

%             [t_value_boffer_cov index_temp]=sort(t_value_cov(:,i_row));
% 
%             s_boffer_cov(i_row)=t_value_boffer_cov(data_rank_p1);

            cov_cov=blkdiag(cov_cov,cov_row_split(:,:,i_row));

     end

        [t_value_max_cov,index_temp]=sort(max(t_value_cov,[],2));
        s_max_cov=t_value_max_cov(data_rank_p1);

         A_m=reshape(sample_mean,A_c,A_r)';
       

       
       M_total_cov=sqrtm((cov_cov));
       [M_r,M_c]=size(M_total_cov);
       M_spl_cov=zeros(M_r,A_c,A_r);
       for j=1:A_r
         M_spl_cov(:,:,j)=M_total_cov(:,(j-1)*A_c+1:j*A_c);
       end
       
               cvx_begin
            variable x_0(d);
            minimize(c'*x_0)
            subject to

            for j_cons=1:15
                A_m(j_cons,:)*x_0+sqrt(s_max_cov)*norm(M_spl_cov(:,:,j_cons)*x_0)-b(j_cons) <= 0;
            end

            x_0 >= 0;
        cvx_end
                
        T_mat_sc2=zeros(n2,A_r);
        D_sc2=zeros(A_c,A_r);  %%%%%%%%% NOTICE: D fixed here!!!
        D_sc2_seperate_x0=zeros(A_r,1);
        d_sc2=zeros(A_r,1);
        for k=1:A_r
            
            if (b(k)-sample_mean((k-1)*d+1:k*d)*x_0)>0
                D_sc2(:,k)=x_0/(b(k)-sample_mean((k-1)*d+1:k*d)*x_0);
                D_sc2_seperate_x0(k)=1/(b(k)-sample_mean((k-1)*d+1:k*d)*x_0);%%%%%%%%% NOTICE: D fixed here!!!
                d_sc2(k)=sample_mean((k-1)*d+1:k*d)*x_0/(b(k)-sample_mean((k-1)*d+1:k*d)*x_0);
            else
                D_sc2(:,k)=-x_0/(b(k)-sample_mean((k-1)*d+1:k*d)*x_0);
                D_sc2_seperate_x0(k)=-1/(b(k)-sample_mean((k-1)*d+1:k*d)*x_0);%%%%%%%%% NOTICE: D fixed here!!!
                d_sc2(k)=-sample_mean((k-1)*d+1:k*d)*x_0/(b(k)-sample_mean((k-1)*d+1:k*d)*x_0);
            end
            T_mat_sc2(:,k)=dataset_2(:,(k-1)*d+1:k*d)*D_sc2(:,k)-d_sc2(k);
            
        end
        t_sc2=max(T_mat_sc2,[],2);
       
        
        [sorted_t_sc2,sort_index_sc2]=sort(t_sc2);
        s_sc2=sorted_t_sc2(data_rank);
        
%         [ cvx_status fv_sc1 x_sc1 ] = cvx_hp_LP_sc1( c,s_sc1,x_0,b1,n,m)
%         [ cvx_status fv_sc2 x_sc2 ] = cvx_hp_LP_sc2( c,s_sc2,x_0,D_sc2_seperate_x0,d_sc2,b1,n,m);%%%%%%%%% NOTICE: D fixed here!!!
        cvx_begin
            variable p;
            minimize(p*c'*x_0)
            subject to


                p*(s_sc2+d_sc2)./D_sc2_seperate_x0<=b;


                p>=0;

            p*x_0 >= 0;
        cvx_end
        x_Recon=p*x_0;
end

