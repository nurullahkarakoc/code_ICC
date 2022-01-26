figure;
%load('nurullah_dyn_9_25_opt_V_07.mat')
%load('nurullah_dyn_9_25_opt_V_06_NoSDN.mat')
load('lorenzo_dyn_9_26_opt_V_01.mat')
%load('lorenzo_dyn_9_26_opt_NoSDN.mat')
k=100;
arrivals1=movmean(reshape(sum(sum(demandsOp1,1),2),[1,M]),k);
arrivals2=movmean(reshape(sum(sum(demandsOp2,1),2),[1,M]),k);
% meanqG1=qG1_Opt_database{6};
% meanqG2=qG2_Opt_database{6};
% meanqG1No=qG1_OptNo_database{6};
% meanqG2No=qG2_OptNo_database{6};
% meanAlloc=allocations_database{6};
% meanAllocNo=allocationsNo_database{6};
%%
sum_qG1=movmean(sum(meanqG1,2),k);
sum_qG2=movmean(sum(meanqG2,2),k);
sum_qG1No=movmean(sum(meanqG1No,2),k);
sum_qG2No=movmean(sum(meanqG2No,2),k);
plot((1:M),arrivals1,(1:M),arrivals2,(1:M),movmean(meanAlloc(:,1),k),(1:M),movmean(meanAlloc(:,2),k),(1:M),movmean(meanAllocNo(:,1),k),(1:M),movmean(meanAllocNo(:,2),k))
xlabel('Time slots');
ylabel('Total Demand Rate per Operator (packets/slot)');

yyaxis right
plot((1:M),sum_qG1,(1:M),sum_qG2,(1:M),sum_qG1No,(1:M),sum_qG2No)
ylabel('Queue Length (packets)');
legend('input - OP1','input - OP2','OP1-total allocations with LB','OP2-total allocations with LB','OP1-total allocations w/o LB','OP2-total allocations w/o LB','OP1-queue with LB','OP2-queue with LB','OP1-queue w/o LB','OP2-queue w/o LB')

averagesLB=zeros(1,6);
averagesNoLB=zeros(1,6);
averagesLB2=zeros(1,6);
averagesNoLB2=zeros(1,6);
for i=1:6
    meanqG1=qG1_Opt_database{i};
    meanqG2=qG2_Opt_database{i};
    meanqG1No=qG1_OptNo_database{i};
    meanqG2No=qG2_OptNo_database{i};
    sum_qG1=movmean(sum(meanqG1,2),k);
    sum_qG2=movmean(sum(meanqG2,2),k);
    sum_qG1No=movmean(sum(meanqG1No,2),k);
    sum_qG2No=movmean(sum(meanqG2No,2),k);
    averagesLB(i)=mean((sum_qG1));
    averagesLB2(i)=mean((sum_qG2));
    averagesNoLB(i)=mean(sum_qG1No);
    averagesNoLB2(i)=mean(sum_qG2No);    
end
    
    
figure;
plot(iter_vec,averagesLB,iter_vec,averagesLB2,iter_vec,averagesNoLB,iter_vec,averagesNoLB2)
legend('OP1-With LB','OP2-With LB','OP1-Without LB','OP2-Without LB')
xlabel('Input Peak Difference Time');
ylabel('Average Queue Length')