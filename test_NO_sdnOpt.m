%% initialization
clc
clear all
close all
linear_utility=true;
V=1; %Drift plus penalty constant
Num_E=10; %Number of eNBs per GW
Num_G=2; %Number of GWs
Num_O=2; %Number of operators
Num_iter=10;
Num_iter2=10;
Num_iter3=10;
Num_iter4=10;
alpha=0.1;
alpha_2=0.1;
alpha_3=0.1;
alpha_4=0.1;
Kbar=.4*Num_E*Num_G*Num_O;
Ko=.4*(Num_E*Num_G)*ones(1,Num_O);
Kbar_o =Ko;
% G=Num_user*0.5;
% lambda=rand;
gammas = .4*ones(Num_E,Num_G,Num_O);
Gs=.4*Num_E*ones(Num_G,Num_O);
lambda_gamma = rand(Num_G,Num_O);
lambda_g=rand(Num_O,1);
lambda_k=rand();
Q=zeros(Num_O,1);
w_s = rand(Num_E,Num_G,Num_O);
%utility = (min(gammas,demands))'*w_s-(lambda*ones(1,Num_user))*gammas;
%% input configuration
M=25;
%Demand arrival rates
mu1=3;
mu2=8;
demandsOp1 = poissrnd(mu1,Num_E,Num_G)/10;
demandsOp2 = poissrnd(mu1,Num_E,Num_G)/10;
for nn = 2:M
    if nn>(5) && nn < (10)+1
        demandsOp1(:,:,nn)= poissrnd(mu2,Num_E,Num_G)/10;
    else
        demandsOp1(:,:,nn)= poissrnd(mu1,Num_E,Num_G)/10;
    end
    demandsOp2(:,:,nn)= poissrnd(mu1,Num_E,Num_G)/10;
end
allocations=zeros(M,Num_O);
%% optimization
utilities=zeros(Num_G,Num_iter);
q1=zeros(Num_E,Num_G);
q2=zeros(Num_E,Num_G);
qG1=zeros(M,2);
qG2=zeros(M,2);
for mm=1:M
    w_s = rand(Num_E,Num_G,Num_O);   %generate new state variable every m
    disp(mm)
    demands(:,:,1)=demandsOp1(:,:,mm)+q1;
    demands(:,:,2)=demandsOp2(:,:,mm)+q2;
  %  for t=1:Num_iter
        for oo=1:Num_O
            %for t2=1:Num_iter2
                for t3=1:Num_iter3  % Updating nu values 
                    for ss=1:Num_G % For each eNB's
                        lambda=lambda_gamma(ss,oo);
                        %         G=Gs(l);
                        for t4=1:Num_iter4 % Updating G's
                            if(linear_utility)
                                %case with linear utility
                                [~,idx_decreasing_state]=sort(w_s(:,ss,oo),'descend');
                                idx_stop_assign=find(cumsum(demands(idx_decreasing_state,ss,oo))>=Gs(ss,oo),1);
                                if(isempty(idx_stop_assign))
                                lambda_gamma(ss,oo)=0;
                                %gammas(:,ss,oo)=demands(:,ss,oo);
                                else
                                lambda_gamma(ss,oo)=w_s(find(cumsum(demands(idx_decreasing_state,ss,oo))>=Gs(ss,oo),1),ss,oo);
                                end
                            else 
                                %case with log utility
                                error('Not implemented yet')
                            end                            
                            Gs(ss,oo)=max(Gs(ss,oo)+alpha_4*(lambda_gamma(ss,oo)-lambda_g(oo)),0);
                        end
                    end
                    lambda_g(oo)=max(lambda_g(oo)-alpha_3*(Ko(oo)-sum(Gs(:,oo))),0);   
                end
              %  Ko(oo)=max(Ko(oo)+alpha_2*(lambda_g(oo)-1/V*Q(oo)-lambda_k),0);
           % end
       % end        
      %  lambda_k=max(lambda_k-alpha*(Kbar-sum(Ko)),0);
    end
    %% Projection for intermediate assignment
   
    for oo=1:Num_O
    Gs_proj(:,oo)=[Ko(oo)/2+(Gs(1,oo)-Gs(2,oo))/2;Ko(oo)/2-(Gs(1,oo)-Gs(2,oo))/2];
    end
                      for oo=1:Num_O
                        for ss=1:Num_G
                            if(linear_utility)
                                %case with linear utility
                                [~,idx_decreasing_state]=sort(w_s(:,ss,oo),'descend');
                                idx_stop_assign=find(cumsum(demands(idx_decreasing_state,ss,oo))>=Gs_proj(ss,oo),1);
                                if(isempty(idx_stop_assign))
                                gammas(:,ss,oo)=demands(:,ss,oo);
                                else
                                    if(idx_stop_assign==1)
                                    gammas(idx_decreasing_state(1),ss,oo)=Gs_proj(ss,oo);
                                    gammas(idx_decreasing_state(2:end),ss,oo)=0;
                                    else
                                    gammas(idx_decreasing_state(1:idx_stop_assign-1),ss,oo)=demands(idx_decreasing_state(1:idx_stop_assign-1),ss,oo);
                                    gammas(idx_decreasing_state(idx_stop_assign),ss,oo)=Gs_proj(ss,oo)-sum(demands(idx_decreasing_state(1:idx_stop_assign-1),ss,oo));
                                    gammas(idx_decreasing_state(idx_stop_assign+1:end),ss,oo)=0;
                                    end   
                                end
                            else 
                                %case with log utility
                                error('Not implemented yet')
                            end 
                          end
                        end
    q=demands-gammas;
    %disp(mm)
    qG1(mm,:)=sum(q(:,:,1));
    qG2(mm,:)=sum(q(:,:,2));
    allocations(mm,:)=reshape(sum(sum(gammas,1),2),[1,Num_O]);  
end
save('lorenzo_test_9_13_no_opt.mat')
%save('firstresults_sdn6_9_11.mat')
% plot(utilities)
% figure;
% plot(lambdas)
% figure;
% plot(Gss(1,:))
% sum(demands)
% sum(gammas)
% Gs
% sum(Gs)
% Ko
% Kbar_o
% sum(Ko)
% Kbar
