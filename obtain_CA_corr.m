function[a,v_out,varargout]=obtain_CA_corr(N_com,delta,u,npl,threat_phi,threat_v,lambda)

% Declare variables

valid_grid=[];
vc=[];
v_out=obtain_payoff_corr(npl,u,threat_phi); %Initialize to worst case
global param; % Declared as global
param=struct('lambda',lambda,'valid_grid',valid_grid,'vc',vc,'threat_v',threat_v,'delta',delta,'npl',npl,'player_optimizing',0);
global param_u;
param_u=struct('u',u);
% Use SOO to sample
settings.dim = 2^npl-1; % The hypercube will be mapped to the simplex!
settings.type = 'det';
for i=1:npl %Optimize as many times as players, one per player!
    param.player_optimizing=i;
    [jitter]=oo(@function_cost_soo_corr,N_com,settings);
end

% Update global variables
par=getGlobalx;
vc=par.vc;
valid_grid=par.valid_grid;

if nargout==3
    region_out=cell(3,1);
    region_out{1}=vc;
    region_out{2}=[];
    region_out{3}=v_out; %Threat v!

    % Replace valid grid by valid payoffs!!
    for i=1:size(valid_grid,1);
        region_out{2}(i,:)=obtain_payoff_corr(npl,u,valid_grid(i,:));
    end
    % Only works for npl==2
    nsample=30;
    actions_sample=linspace(0,1,nsample);
    region_out{4}=[];
    for i=1:nsample
        for j=1:nsample
            for k=1:nsample
                phi_t=[actions_sample(i) actions_sample(j) actions_sample(k)];
                if sum(phi_t)<=1
                    region_out{4}(end+1,:)=obtain_payoff_corr(npl,u,[phi_t 1-sum(phi_t)]);
                end
            end
        end
    end
    varargout{1}=region_out;
end


if ~isempty(valid_grid)
    player_grid=cell(npl,1);
    pl=length(threat_phi);
    % Each player sorts its payoffs
    for i=1:npl
        player_grid{i}=valid_grid;
        for j=1:size(player_grid{i},1)
            aux_phi=player_grid{i}(j,:);
            payoff=obtain_payoff_corr(npl,u,aux_phi);
            player_grid{i}(j,pl+1)=payoff(i);
        end
        player_grid{i}=sortrows(player_grid{i},-(pl+1)); %Sort by payoff value
    end
    % Find a pareto value
    pareto_search=1;
    a=threat_phi; %Security intialization: update if possible

    while pareto_search
        %Jointly controlled lottery (here, just random selection for implementation issues)
        w=rand(1);
        
        id=round(1+(size(valid_grid,1)-1)*w);
        a=valid_grid(id,:); % Starting point for search
        v_out=obtain_payoff_corr(npl,u,a);
        % Each player erases dominated strategies
        for i=1:npl
            id=find(player_grid{i}(:,pl+1)==v_out(i),1);
            player_grid{i}=player_grid{i}(1:id,:); %Update grid
        end
        % Intersect and reorder
        for i=1:npl
            valid_grid=intersect(valid_grid(:,1:pl),player_grid{i}(:,1:pl),'rows');
        end
        if isempty(valid_grid) || size(valid_grid,1)==1 % Dominating strategy found
            pareto_search=0;
        else
            % Each player sorts its payoffs
            for i=1:npl
                player_grid{i}=valid_grid;
                for j=1:size(player_grid{i},1)
                    aux_phi=player_grid{i}(j,:);
                    payoff=obtain_payoff_corr(npl,u,aux_phi);
                    player_grid{i}(j,pl+1)=payoff(i);
                end
                player_grid{i}=sortrows(player_grid{i},-(pl+1)); %Sort by payoff value
            end
        end

    end
else %No valid equilibrium for both players
    a=threat_phi;
    v_out=threat_v;
end