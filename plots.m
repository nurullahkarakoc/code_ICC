figure;
inputlambda=[6*ones(1,5) 16*ones(1,5) 6*ones(1,15)];
inputlambda2=[6*ones(1,25)];
load('lorenzo_test_9_13_no_opt.mat')
load('lorenzo_test_9_13_opt.mat')
plot((1:25),inputlambda,(1:25),inputlambda2,(1:25),sum(qG1_Opt,2),(1:25),sum(qG2_Opt,2),(1:25),sum(qG1,2),(1:25),sum(qG2,2))
legend('input - OP1','input - OP2','OP1-queue with LB','OP2-queue with LB','OP1-queue w/o LB','OP2-queue w/o LB')