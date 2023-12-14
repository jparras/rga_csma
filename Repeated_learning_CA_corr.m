%% REPEATED GAME LEARNING USING CA-CORRELATED
% Juan Parras, GAPS-UPM, May 2017
clear all; clc; close all;
%% Initial parameters

n_iter=2e3; %Number of iterations for RM
n_avg=50; %Number of averages per equilibrium

%Save parameter
sa=1;

%Plot parameter
pl=0;

% Load prestored data
load('Data_network_cost_simulations');
N=5; %Number of players
% Cost parameters
ks=1;
kc=1;
kd=0.1;

N_com=100;
delta=0.99;

%% Main loop
actions_avg_cell=cell(4,3); %n2 x case (RM, lambda=0.5, lambda=1)
u_avg_cell=cell(4,3);
idcell=1;
region_out=[];
for n2=[1 2 3 4]
    display(['Case n2 = ' num2str(n2)]);
    n1=N-n2; 

    npl=1+n2; %Number of players
    
    u=obtain_u(npl,S_1,S_2,n1,ks,kc,kd);
    
    actions_avg=zeros(n_avg,2^npl,3);
    u_avg=zeros(n_avg,npl,3);

    for i=1:n_avg
        display(['Iteration ' num2str(i) ' of ' num2str(n_avg)])
        %RM
        [u_out,a]=regret_min_n(npl,u,n_iter);
        a_aux=squeeze(mean(a));
        ac_aux=a_aux(1,:); 
        % Obtain phi vector
        phi=zeros(2^npl,1);
        for k=1:2^npl %Encode actions as binary: each player only has 2!!
            action_string=dec2bin(k-1,npl);
            action=action_string-'0';
            action=action+1;
            action_aux=num2cell(action); %For indexing purposes
            pr=zeros(1,npl);
            for j=1:npl %Update probabilities
                if action(j)==1
                    pr(j)=ac_aux(j);
                else
                    pr(j)=1-ac_aux(j);
                end
            end
            phi(k)=prod(pr);
        end
        actions_rm=phi;
        payoff_rm=obtain_payoff_corr(npl,u,actions_rm);
        u_avg(i,:,1)=payoff_rm;
        actions_avg(i,:,1)=actions_rm;
        % CA case 1
        lambda=0.5;
        if n2==1 && i==1
            [actions,v,region_out]=obtain_CA_corr(N_com,delta,u,npl,actions_rm,payoff_rm,lambda);
        else
            [actions,v]=obtain_CA_corr(N_com,delta,u,npl,actions_rm,payoff_rm,lambda);
        end            
        u_avg(i,:,2)=v;
        actions_avg(i,:,2)=actions;
        % CA case 2
        lambda=1;
        [actions,v]=obtain_CA_corr(N_com,delta,u,npl,actions_rm,payoff_rm,lambda);
        u_avg(i,:,3)=v;
        actions_avg(i,:,3)=actions;
    end
    
    for cas=[1 2 3]
        actions_avg_cell{idcell,cas}=actions_avg(:,:,cas);
        u_avg_cell{idcell,cas}=u_avg(:,:,cas);
    end
    idcell=idcell+1;
end
%% Plot 
if pl==1
    dif_u_s=zeros(4,3,2); %n2 x (mean, minimum, maximum) x lambda
    dif_u_c=zeros(4,3,2);
    for i=1:4 %n2
        % Case: lambda = 0.5
        aux=mean(u_avg_cell{i,2}-u_avg_cell{i,1});
        dif_u_s(i,1,1)=aux(1);
        dif_u_c(i,1,1)=mean(aux(2:end));
        aux=min(u_avg_cell{i,2}-u_avg_cell{i,1});
        dif_u_s(i,2,1)=aux(1);
        dif_u_c(i,2,1)=min(aux(2:end));
        aux=max(u_avg_cell{i,2}-u_avg_cell{i,1});
        dif_u_s(i,3,1)=aux(1);
        dif_u_c(i,3,1)=max(aux(2:end));
        % Case: lambda = 1
        aux=mean(u_avg_cell{i,3}-u_avg_cell{i,1});
        dif_u_s(i,1,2)=aux(1);
        dif_u_c(i,1,2)=mean(aux(2:end));
        aux=min(u_avg_cell{i,3}-u_avg_cell{i,1});
        dif_u_s(i,2,2)=aux(1);
        dif_u_c(i,2,2)=min(aux(2:end));
        aux=max(u_avg_cell{i,3}-u_avg_cell{i,1});
        dif_u_s(i,3,2)=aux(1);
        dif_u_c(i,3,2)=max(aux(2:end));
    end
%     figure();
%     plot(1:4,dif_u_s(:,1,1),'b-o',1:4,dif_u_s(:,2,1),'b-^',1:4,dif_u_s(:,3,1),'b-v');
%     grid on; hold on;
%     plot(1:4,dif_u_s(:,1,2),'r--o',1:4,dif_u_s(:,2,2),'r--^',1:4,dif_u_s(:,3,2),'r--v');
%     figure();
%     plot(1:4,dif_u_c(:,1,1),'b-o',1:4,dif_u_c(:,2,1),'b-^',1:4,dif_u_c(:,3,1),'b-v');
%     grid on; hold on;
%     plot(1:4,dif_u_c(:,1,2),'r--o',1:4,dif_u_c(:,2,2),'r--^',1:4,dif_u_c(:,3,2),'r--v');
    figure();
    plot(region_out{4}(:,1),region_out{4}(:,2),'c.');
    hold on;
    grid on;
    plot(region_out{1}(:,1),region_out{1}(:,2),'bo');
    if isempty(region_out{2})==0
        plot(region_out{2}(:,1),region_out{2}(:,2),'kx');
    end
    plot(region_out{3}(:,1),region_out{3}(:,2),'rs');
%     figure();
%     errorbar(1:4,dif_u_s(:,1,1),dif_u_s(:,2,1)-dif_u_s(:,1,1),dif_u_s(:,1,1)-dif_u_s(:,3,1),'bo-');
%     hold on; grid on;
%     errorbar(1:4,dif_u_s(:,1,2),dif_u_s(:,2,2)-dif_u_s(:,1,2),dif_u_s(:,1,2)-dif_u_s(:,3,2),'rs--');
%     figure();
%     errorbar(1:4,dif_u_c(:,1,1),dif_u_c(:,2,1)-dif_u_c(:,1,1),dif_u_c(:,1,1)-dif_u_c(:,3,1),'bo-');
%     hold on; grid on;
%     errorbar(1:4,dif_u_c(:,1,2),dif_u_c(:,2,2)-dif_u_c(:,1,2),dif_u_c(:,1,2)-dif_u_c(:,3,2),'rs--');
%     % Payoff without difference
    payoff_u_s_ca=zeros(4,3,2); %n2 x (mean, minimum, maximum) x lambda
    payoff_u_c_ca=zeros(4,3,2);
    payoff_u_s_rm=zeros(4,3,1); %n2 x (mean, minimum, maximum) x 1
    payoff_u_c_rm=zeros(4,3,1);
    for i=1:4 %n2
        % Case: RM
        aux=mean(u_avg_cell{i,1});
        payoff_u_s_rm(i,1,1)=aux(1);
        payoff_u_c_rm(i,1,1)=mean(aux(2:end));
        aux=min(u_avg_cell{i,1});
        payoff_u_s_rm(i,2,1)=aux(1);
        payoff_u_c_rm(i,2,1)=min(aux(2:end));
        aux=max(u_avg_cell{i,1});
        payoff_u_s_rm(i,3,1)=aux(1);
        payoff_u_c_rm(i,3,1)=max(aux(2:end));
        % Case: lambda = 0.5
        aux=mean(u_avg_cell{i,2});
        payoff_u_s_ca(i,1,1)=aux(1);
        payoff_u_c_ca(i,1,1)=mean(aux(2:end));
        aux=min(u_avg_cell{i,2});
        payoff_u_s_ca(i,2,1)=aux(1);
        payoff_u_c_ca(i,2,1)=min(aux(2:end));
        aux=max(u_avg_cell{i,2});
        payoff_u_s_ca(i,3,1)=aux(1);
        payoff_u_c_ca(i,3,1)=max(aux(2:end));
        % Case: lambda = 1
        aux=mean(u_avg_cell{i,3});
        payoff_u_s_ca(i,1,2)=aux(1);
        payoff_u_c_ca(i,1,2)=mean(aux(2:end));
        aux=min(u_avg_cell{i,3});
        payoff_u_s_ca(i,2,2)=aux(1);
        payoff_u_c_ca(i,2,2)=min(aux(2:end));
        aux=max(u_avg_cell{i,3});
        payoff_u_s_ca(i,3,2)=aux(1);
        payoff_u_c_ca(i,3,2)=max(aux(2:end));
    end
    print_error_area(payoff_u_s_rm, payoff_u_s_ca, ['r', 'k', 'b'], 'corr_s_tikz')
    print_error_area(payoff_u_c_rm, payoff_u_c_ca, ['r', 'k', 'b'], 'corr_c_tikz')
end

%% Save
if sa==1
    save('Values_CA_corr_learning_paper');
end