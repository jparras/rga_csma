function[a,v_out,varargout]=obtain_CA_nash(N_com,delta,u,npl,threat_str,threat_v,lambda,l_grid)

% Declare variables

valid_grid=[];
ac_dev=linspace(0,1,l_grid);
vc=[];
v_out=obtain_payoff(npl,u,threat_str); %Initialize to worst case
global param; % Declared as global
param=struct('lambda',lambda,'valid_grid',valid_grid,'ac_dev',ac_dev,'vc',vc,'threat_v',threat_v,'delta',delta,'npl',npl,'player_optimizing',0);
global param_u;
param_u=struct('u',u);
% Use SOO to sample
settings.dim = npl; %As many dimensions as players!
settings.type = 'det';
for i=1:npl %Optimize as many times as players, one per player!
    param.player_optimizing=i;
    [jitter]=oo(@function_cost_soo,N_com,settings);
end

% Update global variables
par=getGlobalx;
vc=par.vc;
valid_grid=par.valid_grid;

if nargout==3
    region_out=cell(3,1);
    region_out{1}=vc;
    region_out{2}=valid_grid;
    region_out{3}=v_out; %Threat v!

    % Replace valid grid by valid payoffs!!
    for i=1:size(valid_grid,1);
        region_out{2}(i,:)=obtain_payoff(npl,u,valid_grid(i,:));
    end
    % Only works for npl==2
    nsample=100;
    actions_sample=linspace(0,1,nsample);
    region_out{4}=zeros(nsample^2,2);
    idx=1;
    for i=1:nsample
        for j=2:nsample
            aux_ac=[actions_sample(i) actions_sample(j)];
            region_out{4}(idx,:)=obtain_payoff(npl,u,aux_ac);
            idx=idx+1;
        end
    end
    varargout{1}=region_out;
end


if ~isempty(valid_grid)
    player_grid=cell(npl,1);
    % Each player sorts its payoffs
    for i=1:npl
        player_grid{i}=valid_grid;
        for j=1:size(player_grid{i},1)
            aux_act=player_grid{i}(j,:);
            payoff=obtain_payoff(npl,u,aux_act);
            player_grid{i}(j,3)=payoff(i);
        end
        player_grid{i}=sortrows(player_grid{i},-3); %Sort by payoff value
    end
    % Find a pareto value
    pareto_search=1;
    a=threat_str; %Security intialization: update if possible

    while pareto_search
        %Jointly controlled lottery (here, just random selection for implementation issues)
        w=rand(1);
        
        id=round(1+(size(valid_grid,1)-1)*w);
        a=valid_grid(id,:); % Starting point for search
        v_out=obtain_payoff(npl,u,a);
        % Each player erases dominated strategies
        for i=1:npl
            id=find(player_grid{i}(:,3)==v_out(i),1);
            player_grid{i}=player_grid{i}(1:id,:); %Update grid
        end
        % Intersect and reorder
        for i=1:npl
            valid_grid=intersect(valid_grid(:,1:2),player_grid{i}(:,1:2),'rows');
        end
        if isempty(valid_grid) || size(valid_grid,1)==1 % Dominating strategy found
            pareto_search=0;
        else
            % Each player sorts its payoffs
            for i=1:npl
                player_grid{i}=valid_grid;
                for j=1:size(player_grid{i},1)
                    aux_act=player_grid{i}(j,:);
                    payoff=obtain_payoff(npl,u,aux_act);
                    player_grid{i}(j,3)=payoff(i);
                end
                player_grid{i}=sortrows(player_grid{i},-3); %Sort by payoff value
            end
        end

    end
else %No valid equilibrium for both players
    a=threat_str;
    v_out=obtain_payoff(npl,u,a);
end