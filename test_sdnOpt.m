%% initialization
clc
close all
linear_utility=true;
Num_E=10; %Number of eNBs per GW
Num_G=2; %Number of GWs
Num_O=2; %Number of operators
alpha_vec=[.1,1,10];
iter_vec=[10,20,30];
V_vec=[.1,1,10];
tot_params=length(V_vec)*length(iter_vec)*length(alpha_vec);
input_params=cell(tot_params,1);
for aa=1:length(alpha_vec)
    for ii=1:length(iter_vec)
        for vv=1:length(V_vec)
            input_params{vv+(ii-1)*length(V_vec)+(aa-1)*length(V_vec)*length(iter_vec)}=[V_vec(vv),iter_vec(ii),alpha_vec(aa)];
        end
    end
end
allocations_database=cell(tot_params,1);
qG1_Opt_database=cell(tot_params,1);
qG2_Opt_database=cell(tot_params,1);
utility_database=cell(tot_params,1);
M=25;
parfor tt=1:tot_params
    V=input_params{tt}(1);
    Num_iter=input_params{tt}(2);
    alpha=input_params{tt}(3);
    Num_iter2=Num_iter;
    Num_iter3=Num_iter;
    Num_iter4=Num_iter;
    alpha_2=alpha;
    alpha_3=alpha;
    alpha_4=alpha;
    utility=zeros(M,1);
    demands=zeros(Num_E,Num_G,Num_O);
    Gs_proj=zeros(Num_G,Num_O);
    Kbar=.4*Num_E*Num_G*Num_O;
    Ko=.4*(Num_E*Num_G)*ones(1,Num_O);
    Kbar_o =Ko;
    % G=Num_user*0.5;
    % lambda=rand;
    gamma = .4*ones(Num_E,Num_G,Num_O);
    Gs=.4*Num_E*ones(Num_G,Num_O);
    lambda_gamma = rand(Num_G,Num_O);
    lambda_g=rand(Num_O,1);
    lambda_k=rand();
    Q=zeros(Num_O,1);
    %utility = (min(gamma,demands))'*w_s-(lambda*ones(1,Num_user))*gamma;
    %% input configuration
    %Demand arrival rates
    mu1=3;
    mu2=8;
    demandsOp1 = poissrnd(mu1,Num_E,Num_G)/10;
    demandsOp2 = poissrnd(mu1,Num_E,Num_G)/10;
    for mm = 2:M
        if mm>(5) && mm < (10)+1
            demandsOp1(:,:,mm)= poissrnd(mu2,Num_E,Num_G)/10;
        else
            demandsOp1(:,:,mm)= poissrnd(mu1,Num_E,Num_G)/10;
        end
        demandsOp2(:,:,mm)= poissrnd(mu1,Num_E,Num_G)/10;
    end
    allocations=zeros(M,Num_O);
    %% optimization
    utilities=zeros(Num_G,Num_iter);
    q1=zeros(Num_E,Num_G);
    q2=zeros(Num_E,Num_G);
    qG1_Opt=zeros(M,2);
    qG2_Opt=zeros(M,2);
    for mm=1:M
        w_s = rand(Num_E,Num_G,Num_O);   %generate new state variable every m
        %disp(mm)
        demands(:,:,1)=demandsOp1(:,:,mm)+q1;
        demands(:,:,2)=demandsOp2(:,:,mm)+q2;
        for t=1:Num_iter
            for oo=1:Num_O
                for t2=1:Num_iter2
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
                                        %gamma(:,ss,oo)=demands(:,ss,oo);
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
                    Ko(oo)=max(Ko(oo)+alpha_2*(lambda_g(oo)-1/V*Q(oo)-lambda_k),0);
                end
            end
            lambda_k=max(lambda_k-alpha*(Kbar-sum(Ko)),0);
        end
        %% Projection for intermediate assignment
        Ko_proj=[Kbar/2+(Ko(1)-Ko(2))/2,Kbar/2-(Ko(1)-Ko(2))/2];
        if(~isempty(find(Ko_proj<0)))
            Ko_proj(find(Ko_proj<0))=0;
            Ko_proj(find(Ko_proj>0))=Kbar;
        else
        end
        for oo=1:Num_O
            Gs_proj(:,oo)=[Ko_proj(oo)/2+(Gs(1,oo)-Gs(2,oo))/2;Ko_proj(oo)/2-(Gs(1,oo)-Gs(2,oo))/2];
        end
        for oo=1:Num_O
            for ss=1:Num_G
                if(linear_utility)
                    %case with linear utility
                    [~,idx_decreasing_state]=sort(w_s(:,ss,oo),'descend');
                    idx_stop_assign=find(cumsum(demands(idx_decreasing_state,ss,oo))>=Gs_proj(ss,oo),1);
                    if(isempty(idx_stop_assign))
                        gamma(:,ss,oo)=demands(:,ss,oo);
                    else
                        if(idx_stop_assign==1)
                            gamma(idx_decreasing_state(1),ss,oo)=Gs_proj(ss,oo);
                            gamma(idx_decreasing_state(2:end),ss,oo)=0;
                        else
                            gamma(idx_decreasing_state(1:idx_stop_assign-1),ss,oo)=demands(idx_decreasing_state(1:idx_stop_assign-1),ss,oo);
                            gamma(idx_decreasing_state(idx_stop_assign),ss,oo)=Gs_proj(ss,oo)-sum(demands(idx_decreasing_state(1:idx_stop_assign-1),ss,oo));
                            gamma(idx_decreasing_state(idx_stop_assign+1:end),ss,oo)=0;
                        end
                    end
                else
                    %case with log utility
                    error('Not implemented yet')
                end
            end
        end
        q=demands-gamma;
        %%disp(mm)
        qG1_Opt(mm,:)=sum(q(:,:,1));
        qG2_Opt(mm,:)=sum(q(:,:,2));
        allocations(mm,:)=reshape(sum(sum(gamma,1),2),[1,Num_O]);
        utility(mm)=sum(sum(sum(w_s.*gamma)));
    end
    allocations_database{tt}=allocations;
    qG1_Opt_database{tt}=qG1_Opt;
    qG2_Opt_database{tt}=qG2_Opt;
    utility_database{tt}=utility;
end
%save('lorenzo_test_9_13_opt_V_01.mat')
% plot(utilities)
% figure;
% plot(lambdas)
% figure;
% plot(Gss(1,:))
% sum(demands)
% sum(gamma)
% Gs
% sum(Gs)
% Ko
% Kbar_o
% sum(Ko)
% Kbar
