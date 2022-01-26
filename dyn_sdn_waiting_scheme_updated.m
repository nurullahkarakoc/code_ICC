clc
clear all
close all
linear_utility=true;  % if false logarithmic utility for fairness
Num_E=10; %Number of eNBs per GW
Num_G=2; %Number of GWs per OP
Num_O=2; %Number of operators
M=100e3;
num_avg_each_case=1;
%alpha_vec=[1e-2,5e-2,1e-1,5e-1,1,5,10];
%alpha_vec=[2e-1,5e-1,8e-1];
alpha_vec=5e-1;
%alpha_vec=[];
%alpha_vec=[1e-2,1e-1,5e-1,1];
%alpha_vec=5e-1;
%iter_vec=[0:M/20:M/2]; % I am using as input peak distance for now.
%iter_vec=[0:M/10:M/2];
iter_vec=[M/4];
V_vec=[1e1,1e3];
%iter_vec=[0:M/10:M/2];
%V_vec=[1e1,5e1,1e2,5e2,1e3,5e3,1e4];
%V_vec=[1,1e1,1e2,1e3];
%V_vec=[1e-2,1e-1,1];
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
for tt=1:tot_params % parfor for parallel computation\
    disp(tt)
    meanAlloc=zeros(M/10,Num_O);
    meanqG1=zeros(M/10,Num_O);
    meanqG2=zeros(M/10,Num_O);
    meanUti=zeros(M/10,1);
    for count=1:num_avg_each_case
        count_lambda_Z=0;
        count_x_0=0;
        count_lambda_x_o=0;
        count_y_g=0;
        V=input_params{tt}(1);  % Lyapunov drift plus penalty parameter
        dist=input_params{tt}(2); % Number of iterations of algorithms
        alpha=input_params{tt}(3); % Step sizes of algorithms
        alpha_2=alpha;
        alpha_3=alpha; 
        alpha_4=alpha;
        demands=zeros(Num_E,Num_G,Num_O);
        q=zeros(Num_E,Num_G,Num_O);
        Kbar=.5*Num_E*Num_G*Num_O;
        Ko=.5*(Num_E*Num_G)*ones(Num_O,1);
        Ko_proj=Ko;
        Kbar_o=Ko;
        gamma =.5*ones(Num_E,Num_G,Num_O);
        Gs=.5*Num_E*ones(Num_G,Num_O);
        Gs_proj=Gs;
        lambda_gamma=zeros(Num_G,Num_O);
        lambda_g=zeros(Num_O,1);
        lambda_k=0;
        QQ=zeros(Num_O,1);
        %% input configuration
        %Demand arrival rates
%        mu1=4;
%        mu2=8;
%         demandsOp1 = poissrnd(mu1,Num_E,Num_G)/Num_E;
%         demandsOp2 = poissrnd(mu1,Num_E,Num_G)/Num_E;
%         for mm = 2:M
%             if mm>(M/10) && mm < (2*M/10)+1
%                 demandsOp1(:,:,mm)= poissrnd(mu2,Num_E,Num_G)/Num_E;
%             else
%                 demandsOp1(:,:,mm)= poissrnd(mu1,Num_E,Num_G)/Num_E;
%             end
%             if mm>(M/10+dist) && mm < (2*M/10)+dist+1
%                 demandsOp2(:,:,mm)= poissrnd(mu2,Num_E,Num_G)/Num_E;
%             else
%                 demandsOp2(:,:,mm)= poissrnd(mu1,Num_E,Num_G)/Num_E;
%             end
%         end
load('demandsOp_input_database')
        %% optimization
        demands_received=gamma;
        w_s_received=zeros(Num_E,Num_G,Num_O);
        avg_w_s_new=zeros(Num_E,Num_G,Num_O);
        utility=zeros(M/10,1);
        allocations=zeros(M/10,Num_O);
        q1=zeros(Num_E,Num_G);
        q2=zeros(Num_E,Num_G);
        qG1_Opt=zeros(M/10,2);
        qG2_Opt=zeros(M/10,2);
        for mm=1:M
            if(mm>M/10)
                pippo=1;  %for debug purpose
            end
            demands(:,:,1)=q(:,:,1)+demandsOp1(:,:,mm);
            demands(:,:,2)=q(:,:,2)+demandsOp2(:,:,mm);
            %w_s = rand(Num_E,Num_G,Num_O);   %generate new state variable every m
            w_s=demands;
            %Upper level
          %  if(mod(floor(mm/1e3),2)==0)
                %Update K
                %new lambda_K in
                if(mod(mm,2e3)==0)
                    lambda_k=max(lambda_k-alpha*(Kbar-sum(Ko)),0);
                    count_lambda_Z=count_lambda_Z+1;
                    %new intermediate assignment of resources K and update of queues
                    if(sum(Ko)~=0)
                    Ko_proj=Kbar*Ko/sum(Ko)
                    else
                    Ko_proj(:)=Kbar/length(Ko)
                    end
                    %can be distributed at the operators
                    QQ=max(0,QQ+Ko_proj-Kbar_o)
                    %reinitialize into the projected value to push back?
                    %Ko=Ko_proj;
                end
                if(mod(mm,2e2)==0)
                    for oo=1:Num_O % For each OP's
                        Ko(oo)=max(Ko(oo)+alpha_2*(lambda_g(oo)-(1/V)*QQ(oo)-lambda_k),0);
                        %QQ(oo)=max(0,QQ(oo)+Ko(oo)-Kbar_o(oo));  %not sure if queues should be updated with the K_proj values
                    end
                    count_x_0=count_x_0+1;
                end
         %   else
                %Keep K constant
        %    end
            %Lower level (shifted by 50 ms)
            %   if(mod(floor((mm-50)/1e2),2)==0)   %check when mm negative
                if(mod(mm-50,2e2)==0)
                    %new lambda_G in
                    for oo=1:Num_O
                        lambda_g(oo)=max(lambda_g(oo)-alpha_3*(Ko_proj(oo)-sum(Gs(:,oo))),0);
                    end
                    count_lambda_x_o=count_lambda_x_o+1;
                    %new intermediate assignment of resources G
                    for oo=1:Num_O
                        if(sum(Gs(:,oo))~=0)
                        Gs_proj(:,oo)=Ko_proj(oo)*Gs(:,oo)/sum(Gs(:,oo));
                        else
                        Gs_proj(:,oo)=Ko_proj(oo)/length(Gs(:,oo));
                        end
                    end
                    %reinitialize into the projected value to push back?
                    %Gs=Gs_proj;
                end
                %Update G
                if(mod(mm,1e1)==0)
                    for oo=1:Num_O
                        for ss=1:Num_G
                            Gs(ss,oo)=max(Gs(ss,oo)+alpha_4*(lambda_gamma(ss,oo)-lambda_g(oo)),0);
                        end
                    end
                    count_y_g=count_y_g+1;
                end
        %    else
                %Keep G constant
        %    end
            %Bottom level
            if(mod(mm,1e1)==0)
                %REAL ALLOCATION PROBLEM
                for oo=1:Num_O
                    for ss=1:Num_G
                        if(Gs_proj(ss,oo)==0)
                        
                        else
                            
                        end
                        HD_weights=min(demands_received(:,ss,oo),Gs_proj(ss,oo));
                        HD_weights=2*demands_received(:,ss,oo).*HD_weights-HD_weights.^2;   %BUG FIXED
                        HD_weights=HD_weights./sum(HD_weights);
                        [~,idx_decreasing_state]=sort(HD_weights,'descend');
                        idx_stop_assign=find(cumsum(demands_received(idx_decreasing_state,ss,oo))>=Gs_proj(ss,oo),1);
                        if(isempty(idx_stop_assign))
                            %lambda_gamma(ss,oo)=0;
                            gamma(:,ss,oo)=demands_received(:,ss,oo);
                        else
                            if(linear_utility)
                                %case with linear utility
                                %lambda_gamma(ss,oo)=w_s_received(find(cumsum(demands_received(idx_decreasing_state,ss,oo))>=Gs_proj(ss,oo),1),ss,oo);
                                %lambda_gamma(ss,oo)=HD_weights(idx_decreasing_state(find(cumsum(demands_received(idx_decreasing_state,ss,oo))>=Gs_proj(ss,oo),1)));
                                if(idx_stop_assign==1)
                                    gamma(idx_decreasing_state(1),ss,oo)=Gs_proj(ss,oo);
                                    gamma(idx_decreasing_state(2:end),ss,oo)=0;
                                else
                                    gamma(idx_decreasing_state(1:idx_stop_assign-1),ss,oo)=demands_received(idx_decreasing_state(1:idx_stop_assign-1),ss,oo);
                                    gamma(idx_decreasing_state(idx_stop_assign),ss,oo)=Gs_proj(ss,oo)-sum(demands_received(idx_decreasing_state(1:idx_stop_assign-1),ss,oo));
                                    gamma(idx_decreasing_state(idx_stop_assign+1:end),ss,oo)=0;
                                end
                            else
                                %case with log utility
                                %run bisection
                                margin=1e-6;
                                lambda_a=1e-6;
                                lambda_b=1-1e-6;
                                err=1;
                                while(err>margin)
                                    lambda_th=mean([lambda_a,lambda_b]);
                                    if(sum(max(min(w_s_received(:,ss,oo)./lambda_th-1,demands_received(:,ss,oo)),0))>Gs_proj(ss,oo))
                                        lambda_a=lambda_th;
                                    else
                                        lambda_b=lambda_th;
                                    end
                                    err=abs(lambda_a-lambda_b);
                                end
                                %lambda_gamma(ss,oo)=mean([lambda_a,lambda_b]);
                                gamma(:,ss,oo)=max(min(w_s_received(:,ss,oo)./mean([lambda_a,lambda_b])-1,demands_received(:,ss,oo)),0);
                            end
                        end
                    end
                end
                % PSEUDO-ALLOCATION TO DETERMINE THE LAGRANGIAN
                for oo=1:Num_O
                    for ss=1:Num_G
                        HD_weights=min(demands_received(:,ss,oo),Gs(ss,oo));
                        HD_weights=2*demands_received(:,ss,oo).*HD_weights-HD_weights.^2;   %BUG FIXED
                        HD_weights=HD_weights./sum(HD_weights)*5; %check the *5 to look at stability
                        [~,idx_decreasing_state]=sort(HD_weights,'descend');
                        idx_stop_assign=find(cumsum(demands_received(idx_decreasing_state,ss,oo))>=Gs(ss,oo),1);
                        if(isempty(idx_stop_assign))
                            lambda_gamma(ss,oo)=0;
                            %gamma(:,ss,oo)=demands_received(:,ss,oo);
                        else
                            if(linear_utility)
                                %case with linear utility
                                lambda_gamma(ss,oo)=HD_weights(idx_decreasing_state(find(cumsum(demands_received(idx_decreasing_state,ss,oo))>=Gs(ss,oo),1)));
                            else
                                %case with log utility
                                %run bisection
                                margin=1e-6;
                                lambda_a=1e-6;
                                lambda_b=1-1e-6;
                                err=1;
                                while(err>margin)
                                    lambda_th=mean([lambda_a,lambda_b]);
                                    if(sum(max(min(w_s_received(:,ss,oo)./lambda_th-1,demands_received(:,ss,oo)),0))>Gs(ss,oo))
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
            end
            avg_w_s_new=avg_w_s_new+w_s;
            if(mod(mm-5,1e1)==0)
                if(mm==5)
                    avg_w_s_new=avg_w_s_new/5;
                    demands_received(:,:,1)=max(demands_received(:,:,1)-gamma(:,:,1)+sum(demandsOp1(:,:,mm-5+1:mm),3)/5,0);
                    demands_received(:,:,2)=max(demands_received(:,:,2)-gamma(:,:,2)+sum(demandsOp2(:,:,mm-5+1:mm),3)/5,0);
                else
                    avg_w_s_new=avg_w_s_new/10;
                    demands_received(:,:,1)=max(demands_received(:,:,1)-gamma(:,:,1)+sum(demandsOp1(:,:,mm-10+1:mm),3)/10,0);
                    demands_received(:,:,2)=max(demands_received(:,:,2)-gamma(:,:,2)+sum(demandsOp2(:,:,mm-10+1:mm),3)/10,0);
                end
%              for oo=1:Num_O
%                     for ss=1:Num_G
%                w_s_received(:,ss,oo)=avg_w_s_new(:,ss,oo)/sum(avg_w_s_new(:,ss,oo));
%                     end   
%              end
                w_s_received=avg_w_s_new;
                %sum(demands_received,1);
                avg_w_s_new=zeros(Num_E,Num_G,Num_O);
            end
            q=max(demands-gamma,0);
            %%disp(mm)
            if(mod(mm,1e1)==0)
                qG1_Opt(floor(mm/10),:)=sum(q(:,:,1));
                qG2_Opt(floor(mm/10),:)=sum(q(:,:,2));
                allocations(floor(mm/10),:)=reshape(sum(sum(gamma,1),2),[1,Num_O]);
                if(linear_utility)
                    utility(floor(mm/10))=sum(sum(sum(w_s.*gamma)));
                else
                    utility(floor(mm/10))=sum(sum(sum(w_s.*log(1+gamma))));
                end
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
%save('nurullah_10_4_updated_v2')
save('lorenzo_10_7_new')