function [x_RO] = RO_ccp(dataset_1,dataset_2,data_rank,c,b)

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

    t_value_cov=zeros(n2,A_r);

    cov_cov=[];
    for i_row=1:A_r
%             dataset_row_split(:,:,i_row)=dataset_1(:,((i_row-1)*A_c+1):(i_row*A_c));
            dataset_p1_row_split(:,:,i_row)=dataset_1(:,((i_row-1)*A_c+1):(i_row*A_c));
            dataset_p2_row_split(:,:,i_row)=dataset_2(:,((i_row-1)*A_c+1):(i_row*A_c));
            
            mean_row_split(:,i_row)=mean(dataset_p1_row_split(:,:,i_row))';
            cov_row_split(:,:,i_row)=((cov(dataset_p1_row_split(:,:,i_row))));


            mean_mat=repmat(mean_row_split(:,i_row)',n2,1);

            T_value_mat_Split_cov=(dataset_p2_row_split(:,:,i_row)-mean_mat)*cov_row_split(:,:,i_row)^-1*(dataset_p2_row_split(:,:,i_row)-mean_mat)';  % calculate only for diag
         
            t_value_cov(:,i_row)=diag(T_value_mat_Split_cov);%%

%             [t_value_boffer_cov index_temp]=sort(t_value_cov(:,i_row));
% 
%             s_boffer_cov(i_row)=t_value_boffer_cov(data_rank);

            cov_cov=blkdiag(cov_cov,cov_row_split(:,:,i_row));

        end

        [t_value_max_cov,index_temp]=sort(max(t_value_cov,[],2));
        s_max_cov=t_value_max_cov(data_rank);

         A_m=reshape(sample_mean,A_c,A_r)';
       

       
       M_total_cov=sqrtm((cov_cov));
       [M_r,M_c]=size(M_total_cov);
       M_spl_cov=zeros(M_r,A_c,A_r);
       for j=1:A_r
         M_spl_cov(:,:,j)=M_total_cov(:,(j-1)*A_c+1:j*A_c);
       end
  
%         
        cvx_begin
            variable x_RO(d);
            minimize(c'*x_RO)
            subject to

            for j_cons=1:15
                A_m(j_cons,:)*x_RO+sqrt(s_max_cov)*norm(M_spl_cov(:,:,j_cons)*x_RO)-b(j_cons) <= 0;
            end

            x_RO >= 0;
        cvx_end

end

