clc
clear all
close all
load('input_function_test')
demands(:,:,1)=demandsOp1(:,:,1);
demands(:,:,2)=demandsOp2(:,:,1);
avg_demands_new=0;
avg_w_s_new=0;
gamma=.4*ones(Num_E,Num_G,Num_O);
for mm=1:M     
            if(mod(mm,1e1)==0)                
                demands(:,:,1)=max(demands(:,:,1)-gamma(:,:,1)+sum(demandsOp1(:,:,mm+1:mm+11),3)/10,0);
                demands(:,:,2)=max(demands(:,:,2)-gamma(:,:,2)+sum(demandsOp1(:,:,mm+1:mm+11),3)/10,0);
            end
  sum(demands,1)
end
%%
T=10;
arrivals1=movmean(reshape(sum(sum(demandsOp1,1),2),[1,M]),T);
arrivals2=movmean(reshape(sum(sum(demandsOp2,1),2),[1,M]),T);
plot((arrivals1+arrivals2)/2)
% hold on 
% plot(arrivals2)
