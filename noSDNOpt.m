%% initialization
clear;
clear clc;
close all;
Num_E=10;
Num_G=2;
Num_O=2;
Num_iter=100;
Num_iter2=10;
Num_iter3=10;
Num_iter4=10;
Num_iter5=10;
Num_iter6=15;
alpha=0.01;
alpha2=0.1;
alpha3=0.1;
alpha4=0.1;
alpha5=0.1;
alpha6=0.1;
K=Num_E*Num_G;
Ko=[K/2 K/2];
Kind =[K K];
% G=Num_user*0.5;
% lamda=rand;
demands = 2*rand(Num_E,Num_G);  % input demands
demands(:,:,2) = 2*rand(Num_E,Num_G);


gammas = rand(Num_E,Num_G);
gammas(:,:,2) = rand(Num_E,Num_G);
Gs=rand(Num_G,Num_O)*Num_E;
lamdaa = rand(Num_G,Num_O);
nuu=rand(Num_O,1);
queue=zeros(Num_O,1);
mu=rand();

w_s = ones(Num_E,Num_G);
w_s(:,:,2)= ones(Num_E,Num_G);

%utility = (min(gammas,demands))'*w_s-(lamda*ones(1,Num_user))*gammas;
%% input configuration
interval=25;
lambda1=3;
lambda2=8;
demandsOp1 = poissrnd(lambda1,Num_E,Num_G)/10;
demandsOp2 = poissrnd(lambda1,Num_E,Num_G)/10;
for i = 2:interval
    if i>(5) && i < (10)+1
        demandsOp1(:,:,i)= poissrnd(lambda2,Num_E,Num_G)/10;
    else
        demandsOp1(:,:,i)= poissrnd(lambda1,Num_E,Num_G)/10;
    end
    demandsOp2(:,:,i)= poissrnd(lambda1,Num_E,Num_G)/10;
end


%% optimization
utilities=zeros(Num_G,Num_iter);
lamdas=zeros(1,Num_iter);
Gss=zeros(Num_G,Num_iter4);
q1=zeros(Num_E,Num_G);
q2=zeros(Num_E,Num_G);
qG1=zeros(interval,2);
qG2=zeros(interval,2);
for s=1:interval
    demands(:,:,1)=demandsOp1(:,:,s)+q1;
    demands(:,:,2)=demandsOp2(:,:,s)+q2;
    for r=1:Num_iter6
        for p=1:Num_O
            nu=nuu(p);
                 for n=1:Num_iter4  % Updating nu values
                    for l=1:Num_G % For each eNB's
                        lamda=lamdaa(l,p);
                        %         G=Gs(l);
                        for m=1:Num_iter3 % Updating G's
                            for j=1:Num_iter % Updating lamdas
                                %                                 gammass=zeros(Num_E,Num_iter2);
                                lamdas(j)=lamda;
                                for i=1:Num_E  % For each user end
                                    for k=1:Num_iter2  % Updating gamma
                                        gammas(i,l,p)=gammas(i,l,p)+alpha2*(w_s(i,l,p)-lamda);
                                        gammas(i,l,p)=max(min(demands(i,l,p),gammas(i,l,p)),0);
                                        %                         gammas(i)=max(gammas(i),0);
                                        %                         gammass(i,k)=gammas(i);
                                    end
                                end
                                
                                %                 utility2 = (min(gammas,demands))'*w_s-(lamda*ones(1,Num_user))*gammas;
                                
                                
                                lamda=max(lamda-alpha*(Gs(l,p)-sum(gammas(:,l,p))),0);
                                
                                
                                utility3 = (min(gammas(:,l,p),demands(:,l,p)))'*w_s(:,l,p)-(lamda*ones(1,Num_E))*gammas(:,l,p);
                                utilities(l,j)=utility3;
                            end
                            Gs(l,p)=max(Gs(l,p)+alpha3*(lamda-nu),0);
                        end
                    end
                    nu=max(nu-alpha4*(Ko(p)-sum(Gs(:,p))),0);
                end
        end
%         mu=max(mu-alpha6*(K-sum(Ko)),0);
        disp(r)
    end
    A=demands-gammas;
    q1= A(:,:,1);
    q2= A(:,:,2);
    disp(s)
    qG1(s,:)=sum(q1);
    qG2(s,:)=sum(q2);
end
save('firstresults_nosdn4.mat')
% plot(utilities)
% figure;
% plot(lamdas)
% figure;
% plot(Gss(1,:))
sum(demands)
sum(gammas)
Gs
sum(Gs)
Ko
Kind
sum(Ko)
K
