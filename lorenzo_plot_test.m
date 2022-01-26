% clc
% %clear all
% close all
% %load('lorenzo_dyn_9_26_opt_V_01.mat')
% %load('lorenzo_dyn_9_26_opt_finer.mat')
% %  Num_E=10; %Number of eNBs per GW
% %  Num_G=2; %Number of GWs per OP
% %  Num_O=2; %Number of operators
% %  M=100e3;
% %  %alpha_vec=[1e-2,5e-2,1e-1,5e-1,1,5,10];
% %   alpha_vec=[2e-1,5e-1,8e-1];
% %   iter_vec=[0:M/40:M/2]; % I am using as input peak distance for now.
% %  %V_vec=[1e1,5e1,1e2,5e2,1e3,5e3,1e4];
% %  V_vec=[1,1e1,1e2,1e3];
% %  % tot_params=length(V_vec)*length(iter_vec)*length(alpha_vec);
% % % input_params=cell(tot_params,1);
% % % for aa=1:length(alpha_vec)
% % %     for ii=1:length(iter_vec)
% % %         for vv=1:length(V_vec)
% % %             input_params{vv+(ii-1)*length(V_vec)+(aa-1)*length(V_vec)*length(iter_vec)}=[V_vec(vv),iter_vec(ii),alpha_vec(aa)];
% % %         end
% % %     end
% % % end
% % %load('lorenzo_plot_ready_28_9')
% demandsOp1=cell(length(iter_vec),1);
% demandsOp2=cell(length(iter_vec),1);
% load('lorenzo_28_9_waiting_scheme_finer_count_50.mat')
% for tt=1:length(iter_vec)
% %         mu1=4;
% %         mu2=8;
% %         demandsOp1 = poissrnd(mu1,Num_E,Num_G)/Num_E;
% %         demandsOp2 = poissrnd(mu1,Num_E,Num_G)/Num_E;
% demandsOp1{tt}=zeros(M/100,1);
% demandsOp2{tt}=zeros(M/100,1);
% dist=iter_vec(tt);
%         for mm = 1:100:M        
%             if mm>(M/10) && mm < (2*M/10)+1
%                 demandsOp1{tt}(floor((mm-1)/100)+1)= 16;
%             else
%                 demandsOp1{tt}(floor((mm-1)/100)+1)= 8;
%             end
%             if mm>(M/10+dist) && mm < (2*M/10)+dist+1
%                 demandsOp2{tt}(floor((mm-1)/100)+1)= 16;
%             else
%                 demandsOp2{tt}(floor((mm-1)/100)+1)= 8;
%             end
%         end
% end
%load('lorenzo_28_9_waiting_scheme_finer_count_50.mat')
%load('lorenzo_28_9_waiting_scheme_noSDN_count_50.mat')
close all
%load('nurullah_10_4_updated.mat')
%load('lorenzo_10_7_pre_scheme.mat')
%load('lorenzo_10_7_pre_scheme_test')
%load('lorenzo_10_7_retuned')
%%
cardV=length(V_vec);
cardD=length(iter_vec);
cardA=length(alpha_vec);
aa=1;
% filelist=[{'allocation_1_V_1_d_0'},{'allocation_2_V_1_d_0'},{'allocation_1_V_1e2_d_0'},{'allocation_2_V_1e2_d_0'};
%           {'allocation_1_V_1_d_7'},{'allocation_2_V_1_d_7'},{'allocation_1_V_1e2_d_7'},{'allocation_2_V_1e2_d_7'};
%           {'allocation_1_V_1_d_13'},{'allocation_2_V_1_d_13'},{'allocation_1_V_1e2_d_13'},{'allocation_2_V_1e2_d_13'}];
% filelist=[{'allocation_1_optimal_d_0'},{'allocation_2_optimal_d_0'};
%           {'allocation_1_optimal_d_7'},{'allocation_2_optimal_d_7'};
%           {'allocation_1_optimal_d_13'},{'allocation_2_optimal_d_13'}];
 filelist3=[{'inst_queues_1_noSDN_d_0'},{'inst_queues_2_noSDN_d_0'};
           {'inst_queues_1_noSDN_d_7'},{'inst_queues_2_noSDN_d_7'};
           {'inst_queues_1_noSDN_d_13'},{'inst_queues_2_noSDN_d_13'}];      
%filelist3=[{'inst_queues_1_V_1_d_0'},{'inst_queues_2_V_1_d_0'},{'inst_queues_1_V_1e2_d_0'},{'inst_queues_2_V_1e2_d_0'};
%          {'inst_queues_1_V_1_d_7'},{'inst_queues_2_V_1_d_7'},{'inst_queues_1_V_1e2_d_7'},{'inst_queues_2_V_1e2_d_7'};
%          {'inst_queues_1_V_1_d_13'},{'inst_queues_2_V_1_d_13'},{'inst_queues_1_V_1e2_d_13'},{'inst_queues_2_V_1e2_d_13'}];
dist_idx_vec=[1,7,13];
%dist_idx_vec=13;
 %for aa=1:cardA;
   % figure(aa)
%figure(1)
    %idx_to_test=sub2ind([cardV,cardD,cardA],[1:cardV],dist_idx*ones(1,4),aa*ones(1,4));
    for dd=1:length(dist_idx_vec)
    dist_idx=dist_idx_vec(dd);
  %  idx_to_test=sub2ind([cardV,cardD,cardA],[1,3],dist_idx*ones(1,2),aa*ones(1,2));
    idx_to_test=sub2ind([cardV,cardD,cardA],1,dist_idx,aa);
    T=1e1;
  %  figure(dd)
    %plot([1:100:M],movmean(demandsOp1{dist_idx},T),'--')
    %plot([1:100:M],reshape(sum(sum(demandsOp1{idx_to_test(1)},1),2),[1,M/100]))
%    hold on 
    %plot([100:100:M],reshape(sum(sum(demandsOp2{idx_to_test(1)},1),2),[1,M/100]))
    %plot([1:100:M],movmean(demandsOp2{dist_idx},T),'--')
        for jj=1:length(idx_to_test)
    %    plot([10:10:M],movmean(allocations_database{idx_to_test(jj)}(:,1),T),'LineWidth',2.0);
     %    pippo=movmean(allocations_database{idx_to_test(jj)}(:,1),T);
      %   saveDat([[10:500:M]',pippo(1:50:end)],char(filelist(dd,2*(jj-1)+1)))
      %   plot([10:500:M],pippo(1:50:end))
     %    hold on
     %    saveDat([[10:100:M]',pippo(1:10:end)],char(filelist(dd,2*(jj-1)+1)))
        pippo=sum(qG1_Opt_database{idx_to_test(jj)},2);
        saveDat([[10:100:M]',pippo(1:10:end)],char(filelist3(dd,2*(jj-1)+1)));
%        plot([10:10:M],movmean(allocations_database{idx_to_test(jj)}(:,2),T),'LineWidth',2.0);
 %        pippo=movmean(allocations_database{idx_to_test(jj)}(:,2),T);
         %plot([10:500:M],pippo(1:50:end))
 %        saveDat([[10:500:M]',pippo(1:50:end)],char(filelist(dd,2*jj)))
 %        saveDat([[10:100:M]',pippo(1:10:end)],char(filelist(dd,2*jj)))
        pippo=sum(qG2_Opt_database{idx_to_test(jj)},2);
        saveDat([[10:100:M]',pippo(1:10:end)],char(filelist3(dd,2*jj)));
        end
    end
%     idx_to_test_no=sub2ind([1,cardD,cardA],1,dist_idx,aa);
%     plot([10:10:M],movmean(allocations_database_NOsdn{idx_to_test_no}(:,1),T),'LineWidth',2.0);
%     plot([10:10:M],movmean(allocations_database_NOsdn{idx_to_test_no}(:,2),T),'LineWidth',2.0);
%end 
%%
% pippo=zeros(length(idx_to_test),1);
% for jj=1:length(idx_to_test)
%     pippo(jj)=mean(utility_database{jj});
% end
%%
aa=1;
avg_over_d1_Opt=zeros(21,cardV);
avg_over_d2_Opt=zeros(21,cardV);
%avg_over_d1_Opt=zeros(6,cardV);
%avg_over_d2_Opt=zeros(6,cardV);
%avg_over_d1_OptNo=zeros(21,1);
%avg_over_d2_OptNo=zeros(21,1);
%for vv=1:4
    for ii=1:21
       avg_over_d1_Opt(ii)=mean(mean(qG1_Opt_database{vv+(ii-1)*cardV+(aa-1)*cardV*length(iter_vec)}));
       avg_over_d2_Opt(ii)=mean(mean(qG2_Opt_database{vv+(ii-1)*cardV+(aa-1)*cardV*length(iter_vec)}));
 %       if(vv==1)
 %           avg_over_d1_OptNo(ii)=mean(mean(qG1_Opt_database_NOsdn{ii+(aa-1)*length(iter_vec)}));
 %           avg_over_d2_OptNo(ii)=mean(mean(qG2_Opt_database_NOsdn{ii+(aa-1)*length(iter_vec)}));
 %       end
    end
%end
% filelist2=[{'queues_1_V_1'},{'queues_2_V_1'};
%           {'queues_1_V_1e1'},{'queues_2_V_1e1'};
%           {'queues_1_V_1e2'},{'queues_2_V_1e2'};
%           {'queues_1_V_1e3'},{'queues_2_V_1e3'}];
% filelist2=[{'queues_1_optimal'},{'queues_2_optimal'}];
% saveDat([[0:1:20]',avg_over_d1_Opt(:)],char(filelist2(1)))
% saveDat([[0:1:20]',avg_over_d2_Opt(:)],char(filelist2(2)))
%          
% figure(3)
%  plot(avg_over_d1_OptNo,'--')
%  hold on
%  plot(avg_over_d2_OptNo,'--')
% hold on
figure(4)
%for vv=1:4
% plot(avg_over_d1_Opt(1:2:end,vv))
% hold on
% saveDat([[0:1:20]',avg_over_d1_Opt(:,vv)],char(filelist2(vv,1)))
% saveDat([[0:1:20]',avg_over_d2_Opt(:,vv)],char(filelist2(vv,2)))
% plot(avg_over_d2_Opt(1:2:end,vv))
%end
% 
% %%
% %idx_to_test=[13*(dist_idx-1)+1,13*(dist_idx-1)+4,13*(dist_idx-1)+7,13*(dist_idx-1)+10,13*(dist_idx-1)+11,13*dist_idx];
% %idx_to_test=[13*(dist_idx-1)+1,13*(dist_idx-1)+7,13*(dist_idx-1)+10,13*dist_idx];
% %idx_to_test=[6:11:50];
% %idx_to_test=[13*(dist_idx-1)+10,13*dist_idx];
% % arrivals1=movmean(reshape(sum(sum(demandsOp1,1),2),[1,M]),T);
% % arrivals2=movmean(reshape(sum(sum(demandsOp2,1),2),[1,M]),T);  
% % figure(1)
% % plot([1:M],arrivals1,'b--','LineWidth',2.0);
% % hold on
% % plot([1:M],arrivals2,'r--','LineWidth',2.0);
% % figure(1)
% % hold on
% % for jj=1:length(idx_to_test)
% % plot([10:10:M],movmean(allocations_database{idx_to_test(jj)}(:,1),T),'LineWidth',2.0);
% % plot([10:10:M],movmean(allocations_database{idx_to_test(jj)}(:,2),T),'LineWidth',2.0);
% % end
% %%
%%
figure(2)
for jj=1:length(idx_to_test)
    %plot([10:10:M],movmean(sum(qG1_Opt_database{idx_to_test(jj)},2),T))
    pippo=sum(qG1_Opt_database{idx_to_test(jj)},2);
    plot([10:100:M],pippo(1:10:end));
    hold on
    pippo=sum(qG2_Opt_database{idx_to_test(jj)},2);
    %plot([10:10:M],movmean(sum(qG2_Opt_database{idx_to_test(jj)},2),T))
    plot([10:100:M],pippo(1:10:end))
    %plot([1:M],movmean(sum(qG1_OptNo_database{idx_to_test(jj)},2),T))
    %plot([1:M],movmean(sum(qG2_OptNo_database{idx_to_test(jj)},2),T))    
end
%
% % V_idx_test=[1,7,10,13];
% % avg_over_d1_Opt=zeros(11,length(V_idx_test));
% % avg_over_d2_Opt=zeros(11,length(V_idx_test));
% % avg_over_d1_OptNo=zeros(11,1);
% % avg_over_d2_OptNo=zeros(11,1);
% % for jj=1:length(V_idx_test)
% % for ii=1:11
% % avg_over_d1_Opt(ii,jj)=mean(mean(qG1_Opt_database{V_idx_test(jj)+13*(ii-1)}));
% % avg_over_d2_Opt(ii,jj)=mean(mean(qG2_Opt_database{V_idx_test(jj)+13*(ii-1)}));
% %     if(jj==1)    
% %     avg_over_d1_OptNo(ii)=mean(mean(qG1_OptNo_database{ii}));
% %     avg_over_d2_OptNo(ii)=mean(mean(qG2_OptNo_database{ii}));
% %     end
% % end
% % end
% %%
% figure(4)
% plot(avg_over_d1_OptNo)
% hold on
% plot(avg_over_d2_OptNo)
% % for jj=1:length(V_idx_test)
% % plot(avg_over_d1_Opt(:,jj))
% % plot(avg_over_d2_Opt(:,jj))
% % end
% % legend('static 1','static 2','Lay V 1','Lay V 1/2','Lay V 7','Lay V 7/2','Lay V 10','Lay V 10/2','Lay V 13','Lay V 13/2')
% %%
% % figure(5)
% % legend('static','Lay V 1','Lay V 7','Lay V 10','Lay V 13')
% % %%
% % averagesLB=zeros(1,6);
% % averagesNoLB=zeros(1,6);
% % averagesLB2=zeros(1,6);
% % averagesNoLB2=zeros(1,6);
% % for i=1:6
% %     meanqG1={i};
% %     meanqG2=qG2_Opt_database{i};
% %     meanqG1No=qG1_OptNo_database{i};
% %     meanqG2No=qG2_OptNo_database{i};
% %     sum_qG1=movmean(sum(meanqG1,2),k);
% %     sum_qG2=movmean(sum(meanqG2,2),k);
% %     sum_qG1No=movmean(sum(meanqG1No,2),k);
% %     sum_qG2No=movmean(sum(meanqG2No,2),k);
% %     averagesLB(i)=mean((sum_qG1));
% %     averagesLB2(i)=mean((sum_qG2));
% %     averagesNoLB(i)=mean(sum_qG1No);
% %     averagesNoLB2(i)=mean(sum_qG2No);    
% % end    
% % figure;
% % plot(iter_vec,averagesLB,iter_vec,averagesLB2,iter_vec,averagesNoLB,iter_vec,averagesNoLB2)
% % legend('OP1-With LB','OP2-With LB','OP1-Without LB','OP2-Without LB')
% % xlabel('Input Peak Difference Time');
% % ylabel('Average Queue Length')
% % % sum_qG1=movmean(sum(qG1_Opt_database{idx_to_test(jj)},2),T);
% % % sum_qG2=movmean(sum(qG1_Opt_database{idx_to_test(jj)},2),T);
% % % sum_qG1No=movmean(sum(qG1_OptNo_database{idx_to_test(jj)},2),T);
% % % sum_qG2No=movmean(sum(qG2_OptNo_database{idx_to_test(jj)},2),T);
% % % load('lorenzo_dyn_9_26_opt_V_1e4_1e6.mat')
% % % idx_to_test=[11,12];
% % % for jj=1:length(idx_to_test)
% % % figure(jj+3)
% % % plot([1:T:M],allocations_database{idx_to_test(jj)}(1:T:M,1));
% % % hold on
% % % plot([1:T:M],allocations_database{idx_to_test(jj)}(1:T:M,2));
% % % end
