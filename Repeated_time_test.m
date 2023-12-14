%% REPEATED GAME LEARNING USING CA: TIMES ANALYSIS
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
l_grid=30; 

%% Main loops
nc=5; %Number of cases: 1 RM, 2 nash 0.5, 3 nash 1, 4 corr 0.5, 5 corr 1
time=zeros(4,n_avg,nc);
for n2=[1 2 3 4]
    n1=N-n2; 
    npl=1+n2; %Number of players
    u=obtain_u(npl,S_1,S_2,n1,ks,kc,kd);
    for i=1:n_avg
        display(['Iteration ' num2str(i) ' of ' num2str(n_avg)]);
        % Case 1
        tic;
        [u_out,a]=regret_min_n(npl,u,n_iter);
        time(n2,i,1)=toc;
        a_aux=squeeze(mean(a));
        actions_rm=a_aux(1,:);
        ac_aux=actions_rm;
        payoff_rm=obtain_payoff(npl,u,actions_rm);
        a_aux=squeeze(mean(a));
        % Case 2
        lambda=0.5; 
        tic;
        [actions,v]=obtain_CA_nash(N_com,delta,u,npl,actions_rm,payoff_rm,lambda,l_grid);
        time(n2,i,2)=toc;
        % Case 3
        lambda=1; 
        tic;
        [actions,v]=obtain_CA_nash(N_com,delta,u,npl,actions_rm,payoff_rm,lambda,l_grid);
        time(n2,i,3)=toc;
        % Prepare for corr data obtention
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
        % Case 4
        lambda=0.5; 
        tic;
        [actions,v]=obtain_CA_corr(N_com,delta,u,npl,actions_rm,payoff_rm,lambda);
        time(n2,i,4)=toc;
        % Case 5
        lambda=1; 
        tic;
        [actions,v]=obtain_CA_corr(N_com,delta,u,npl,actions_rm,payoff_rm,lambda);
        time(n2,i,5)=toc;
    end
end

%% Plot
if pl==1
    fig=figure;
    for i=1:5
        time_mean=mean(time(:,:,i),2);
        time_std=std(time(:,:,i),0,2);
        errorbar(1:4,time_mean,1.96*time_std);
        hold on;
    end
    ax = get(fig,'CurrentAxes');
    set(ax,'YScale','log');
    grid on;
    legend('RM','CA Nash \lambda=0.5','CA Nash \lambda=1','CA Corr \lambda=0.5','CA Corr \lambda=1');
end

%% Save

if sa
    save('Data_time_CA');
end
