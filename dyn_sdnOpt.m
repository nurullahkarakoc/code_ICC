%% initialization
clc
close all
linear_utility=false;  % if false logarithmic utility for fairness
Num_E=10; %Number of eNBs per GW
Num_G=2; %Number of GWs per OP
Num_O=2; %Number of operators
M=100000;
% alpha_vec=[.1,1,10];
alpha_vec=0.2;
% iter_vec=[10,20,30];
iter_vec=[0:M/40:M/4]; % I am using as input peak distance for now.
% V_vec=[.1,1,10];
%V_vec=[1e-2,1,1e2];
V_vec=[1,2,5,10,20,50,100,200,500,1e3,2e3,5e3,1e4];
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
% Simulation length
for tt=1:tot_params % parfor for parallel computation
    meanAlloc=zeros(M,Num_O);
    meanqG1=zeros(M,Num_O);
    meanqG2=zeros(M,Num_O);
    meanUti=zeros(M,1);
    for count=1:5        
        Num_iter=4;
        V=input_params{tt}(1);  % Lyapunov drift plus penalty parameter
        dist=input_params{tt}(2); % Number of iterations of algorithms
        alpha=input_params{tt}(3); % Step sizes of algorithms
        Num_iter2=Num_iter;
        Num_iter3=Num_iter;
        Num_iter4=Num_iter;
        alpha_2=alpha;
        alpha_3=alpha;
        alpha_4=alpha;
        utility=zeros(M,1);
        demands=zeros(Num_E,Num_G,Num_O);
        q=zeros(Num_E,Num_G,Num_O);
        Gs_proj=zeros(Num_G,Num_O);
        Kbar=.4*Num_E*Num_G*Num_O;
        Ko=.4*(Num_E*Num_G)*ones(1,Num_O);
        Ko_proj=.4*(Num_E*Num_G)*ones(1,Num_O);
        Kbar_o=Ko;
        gamma = .4*ones(Num_E,Num_G,Num_O);
        Gs=.4*Num_E*ones(Num_G,Num_O);
        lambda_gamma = rand(Num_G,Num_O);
        lambda_g=rand(Num_O,1);
        lambda_k=rand();
        QQ=zeros(Num_O,1);      
        %% input configuration
        %Demand arrival rates
        mu1=3;
        mu2=7;
        demandsOp1 = poissrnd(mu1,Num_E,Num_G)/10;
        demandsOp2 = poissrnd(mu1,Num_E,Num_G)/10;
        for mm = 2:M
            if mm>(M/20) && mm < (2*M/20)+1
                demandsOp1(:,:,mm)= poissrnd(mu2,Num_E,Num_G)/10;
            else
                demandsOp1(:,:,mm)= poissrnd(mu1,Num_E,Num_G)/10;
            end
            if mm>(M/20+dist) && mm < (2*M/20)+dist+1
                demandsOp2(:,:,mm)= poissrnd(mu2,Num_E,Num_G)/10;
            else
                demandsOp2(:,:,mm)= poissrnd(mu1,Num_E,Num_G)/10;
            end
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
            demands(:,:,1)=demandsOp1(:,:,mm)+q(:,:,1);
            demands(:,:,2)=demandsOp2(:,:,mm)+q(:,:,2);
            if(mod(mm,Num_iter2*Num_iter3*Num_iter4)==0)
                lambda_k=max(lambda_k-alpha*(Kbar-sum(Ko)),0);
            end
            for oo=1:Num_O % For each OP's
                if(mod(mm,Num_iter3*Num_iter4)==0) % It occurs every Num_iter3*Num_iter4*mm
                    Ko(oo)=max(Ko(oo)+alpha_2*(lambda_g(oo)-(1/V)*QQ(oo)-lambda_k),0);
                    QQ(oo)=max(0,QQ(oo)+Ko(oo)-Kbar_o(oo));
                end
                if(mod(mm,Num_iter4)==0) % It occurs every Num_iter4*mm
                    lambda_g(oo)=max(lambda_g(oo)-alpha_3*(Ko(oo)-sum(Gs(:,oo))),0);
                end
                for ss=1:Num_G % For each Gw's
                    Gs(ss,oo)=max(Gs(ss,oo)+alpha_4*(lambda_gamma(ss,oo)-lambda_g(oo)),0);
                    [~,idx_decreasing_state]=sort(w_s(:,ss,oo),'descend');
                    idx_stop_assign=find(cumsum(demands(idx_decreasing_state,ss,oo))>=Gs(ss,oo),1);
                    if(isempty(idx_stop_assign))
                        lambda_gamma(ss,oo)=0;
                    else
                        if(linear_utility)
                            %case with linear utility
                            lambda_gamma(ss,oo)=w_s(find(cumsum(demands(idx_decreasing_state,ss,oo))>=Gs(ss,oo),1),ss,oo);
                        else
                            %case with log utility
                            %run bisection
                            margin=1e-6;
                            lambda_a=1e-6;
                            lambda_b=1-1e-6;
                            err=1;
                            while(err>margin)
                                lambda_th=mean([lambda_a,lambda_b]);
                                if(sum(max(min(w_s(:,ss,oo)./lambda_th-1,demands(:,ss,oo)),0))>Gs(ss,oo))
                                    lambda_a=lambda_th;
                                else
                                    lambda_b=lambda_th;
                                end
                                err=abs(lambda_a-lambda_b);
                            end
                            lambda_gamma(ss,oo)=mean([lambda_a,lambda_b]);
                        end
                    end
                end
                
            end
            %% Projection for intermediate assignment
            if(mod(mm,Num_iter2*Num_iter3*Num_iter4)==0) % It occurs every Num_iter3*Num_iter4*mm
                Ko_proj=[Kbar/2+(Ko(1)-Ko(2))/2,Kbar/2-(Ko(1)-Ko(2))/2];
                if(~isempty(find(Ko_proj<0)))
                    Ko_proj(find(Ko_proj<0))=0;
                    Ko_proj(find(Ko_proj>0))=Kbar;
                end
            end
            if(mod(mm,Num_iter3*Num_iter4)==0)
            for oo=1:Num_O
                Gs_proj(:,oo)=[Ko_proj(oo)/2+(Gs(1,oo)-Gs(2,oo))/2;Ko_proj(oo)/2-(Gs(1,oo)-Gs(2,oo))/2];
            end
            end
            if(mod(mm,Num_iter4)==0)
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
                          %run bisection (run again with the projected value)
                            margin=1e-6;
                            lambda_a=1e-6;
                            lambda_b=1-1e-6;
                            err=1;
                            while(err>margin)
                                lambda_th=mean([lambda_a,lambda_b]);
                                if(sum(max(min(w_s(:,ss,oo)./lambda_th-1,demands(:,ss,oo)),0))>Gs_proj(ss,oo))
                                    lambda_a=lambda_th;
                                else
                                    lambda_b=lambda_th;
                                end
                                err=abs(lambda_a-lambda_b);
                            end
                            gamma(:,ss,oo)=max(min(w_s(:,ss,oo)./mean([lambda_a,lambda_b])-1,demands(:,ss,oo)),0);
                    end
                end
            end
            end
            q=demands-gamma;
            %%disp(mm)
            qG1_Opt(mm,:)=sum(q(:,:,1));
            qG2_Opt(mm,:)=sum(q(:,:,2));
            sum(q,1)
            allocations(mm,:)=reshape(sum(sum(gamma,1),2),[1,Num_O]);
            if(linear_utility)
                utility(mm)=sum(sum(sum(w_s.*gamma)));
            else
                utility(mm)=sum(sum(sum(w_s.*log(1+gamma))));
            end
        end
        meanAlloc= (meanAlloc*(count-1)+allocations)/count;
        meanqG1= (meanqG1*(count-1)+qG1_Opt)/count;
        meanqG2= (meanqG2*(count-1)+qG2_Opt)/count;
        meanUti= (meanUti*(count-1)+utility)/count;
    end
    allocations_database{tt}=meanAlloc;
    qG1_Opt_database{tt}=meanqG1;
    qG2_Opt_database{tt}=meanqG2;
    utility_database{tt}=meanUti;
    %     disp(tt)
end
save('lorenzo_dyn_9_26_opt_finer.mat')
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
