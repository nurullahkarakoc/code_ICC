clc
clear all
close all
linear_utility=true;  % if false logarithmic utility for fairness
Num_E=10; %Number of eNBs per GW
Num_G=2; %Number of GWs per OP
Num_O=2; %Number of operators
M=100e3;
%alpha_vec=[1e-2,5e-2,1e-1,5e-1,1,5,10];
%alpha_vec=[4e-1,5e-1,6e-1];
%alpha_vec=5e-1;
alpha_vec=4e-1;
%V_vec=20;
%iter_vec=M/4;
iter_vec=[0:M/40:M/2]; % I am using as input peak distance for now.
%iter_vec=M/4;
%iter_vec=linspace(0,M/2,100);
%V_vec=[1e1,5e1,1e2,5e2,1e3,5e3,1e4];
%iter_vec=M/4;
V_vec=1;
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
%utility_database=cell(tot_params,1);
% Simulation length
% parfor tt=1:tot_params % parfor for parallel computation
parfor tt=1:tot_params
    w_s_received=zeros(Num_E,Num_G,Num_O);
    disp(tt)
    meanAlloc=zeros(M/10,Num_O);
    meanqG1=zeros(M/10,Num_O);
    meanqG2=zeros(M/10,Num_O);
  %  meanUti=zeros(M/10,1);
  dist=input_params{tt}(2);
    for count=1:50      
        demands=zeros(Num_E,Num_G,Num_O);
        q=zeros(Num_E,Num_G,Num_O);
        Kbar=.5*Num_E*Num_G*Num_O;
        gamma =.5*ones(Num_E,Num_G,Num_O);
        %% input configuration
        %Demand arrival rates
        mu1=4;
        mu2=8;
        demandsOp1 = poissrnd(mu1,Num_E,Num_G)/Num_E;
        demandsOp2 = poissrnd(mu1,Num_E,Num_G)/Num_E;
        for mm = 2:M
            if mm>(M/10) && mm < (2*M/10)+1
                demandsOp1(:,:,mm)= poissrnd(mu2,Num_E,Num_G)/Num_E;
            else
                demandsOp1(:,:,mm)= poissrnd(mu1,Num_E,Num_G)/Num_E;
            end
            if mm>(M/10+dist) && mm < (2*M/10)+dist+1
                demandsOp2(:,:,mm)= poissrnd(mu2,Num_E,Num_G)/Num_E;
            else
                demandsOp2(:,:,mm)= poissrnd(mu1,Num_E,Num_G)/Num_E;
            end
        end
        %% optimization
        demands_received=gamma;
        avg_w_s_new=zeros(Num_E,Num_G,Num_O);
        utility=zeros(M/10,1);
        allocations=zeros(M/10,Num_O);
        q1=zeros(Num_E,Num_G);
        q2=zeros(Num_E,Num_G);
        qG1_Opt=zeros(M/10,2);
        qG2_Opt=zeros(M/10,2);
        for mm=1:M
            if(mm>M/10)
                pippo=1;
            end
            if(mm>M/10+dist)
            pippo=1;
            end
            %w_s = rand(Num_E,Num_G,Num_O);   %generate new state variable every m
            demands(:,:,1)=q(:,:,1)+demandsOp1(:,:,mm);
            demands(:,:,2)=q(:,:,2)+demandsOp2(:,:,mm);
            w_s=demands;   %1ST LINE ADDED
            %Bottom level
            if(mod(mm,1e1)==0)
                %send allocations
                        global_w_s_received=reshape(w_s_received,[Num_E*Num_G*Num_O,1]);
                        global_demands_received=reshape(demands_received,[Num_E*Num_G*Num_O,1]);
                        [~,idx_decreasing_state]=sort(global_w_s_received,'descend');
                        %idx_stop_assign=find(cumsum(global_demands_received(idx_decreasing_state)>=Kbar),1);
                        idx_stop_assign=find(cumsum(global_demands_received(idx_decreasing_state))>=Kbar,1);
                        if(isempty(idx_stop_assign))
                       %     lambda_gamma(ss,oo)=0;
                            global_gamma=global_demands_received;
                        else
                                %case with linear utility
                       %         lambda_gamma(ss,oo)=w_s_received(idx_decreasing_state(find(cumsum(demands_received(idx_decreasing_state,ss,oo))>=Gs_proj(ss,oo),1)),ss,oo);
                        %        lambda_gamma(ss,oo)=lambda_gamma(ss,oo)/(sum(w_s_received(:,ss,oo)))*5;    %2ND LINE ADDED
                                if(idx_stop_assign==1)
                                    global_gamma(idx_decreasing_state(1))=Kbar;
                                    global_gamma(idx_decreasing_state(2:end))=0;
                                else
                                    global_gamma(idx_decreasing_state(1:idx_stop_assign-1))=global_demands_received(idx_decreasing_state(1:idx_stop_assign-1));
                                    global_gamma(idx_decreasing_state(idx_stop_assign))=Kbar-sum(global_demands_received(idx_decreasing_state(1:idx_stop_assign-1)));
                                    global_gamma(idx_decreasing_state(idx_stop_assign+1:end))=0;
                                end
                        end
                        gamma=reshape(global_gamma,[Num_E,Num_G,Num_O]);
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
                 %   utility(floor(mm/10))=sum(sum(sum(w_s.*gamma)));
                else
                 %   utility(floor(mm/10))=sum(sum(sum(w_s.*log(1+gamma))));
                end
            end
        end
        meanAlloc= (meanAlloc*(count-1)+allocations)/count;
        meanqG1= (meanqG1*(count-1)+qG1_Opt)/count;
        meanqG2= (meanqG2*(count-1)+qG2_Opt)/count;
      %  meanUti= (meanUti*(count-1)+utility)/count;
    end
    allocations_database{tt}=meanAlloc;
    qG1_Opt_database{tt}=meanqG1;
    qG2_Opt_database{tt}=meanqG2;
  %  utility_database{tt}=meanUti;
    %     disp(tt)
end
save('lorenzo_10_11_optimal_LCQ')
