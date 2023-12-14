function [value_of_function]=function_cost_soo(s)

param=getGlobalx;
lambda=param.lambda;
valid_grid=param.valid_grid;
ac_dev=param.ac_dev;
vc=param.vc;
delta=param.delta;
threat_v=param.threat_v;
i=param.player_optimizing;
npl=param.npl;
param_u=getGlobalu;
u=struct2cell(param_u);

gamma_i=[];
gamma_rest=0;
valid=1;
% Check action for player i
ac_test=s;
payoff=obtain_payoff(npl,u,ac_test);
ui=payoff(i);
% Add to region tested
vc(end+1,:)=payoff;
% Check if satisfies pi restrictions
ac_dev(ac_dev==ac_test(i))=[]; %All actions but current
dev_v=zeros(1,length(ac_dev));
for j=1:length(dev_v) %Obtain all possible deviations
    ac_aux=ac_test;
    ac_aux(i)=ac_dev(j); %Player i deviates and plays j!!
    payoff=obtain_payoff(npl,u,ac_aux);
    dev_v(j)=payoff(i);
end
dev=max(dev_v);
idx=(ui-(1-delta)*dev-delta*threat_v(i)<0);
if sum(idx)==0 && ui>=threat_v(i) %No profitable deviation
    gamma_i=norm(ui-threat_v(i));
else
    gamma_i=-norm(ui-threat_v(i));
    valid=0;
end

% Check action for the rest of the players
rest_of_players=1:npl;
rest_of_players(i)=[]; %All players but player i
for k=1:length(rest_of_players)
    idp=rest_of_players(k); %Player to analyze
    payoff=obtain_payoff(npl,u,ac_test);
    uidp=payoff(idp);
    % Check if satisfies player restrictions
    ac_dev(ac_dev==ac_test(idp))=[]; %All actions but current
    dev_v=zeros(1,length(ac_dev));
    for j=1:length(dev_v) %Obtain all possible deviations
        ac_aux=ac_test;
        ac_aux(idp)=ac_dev(j); %Player i deviates and plays j!!
        payoff=obtain_payoff(npl,u,ac_aux);
        dev_v(j)=payoff(idp);
    end
    dev=max(dev_v);
    idx=(uidp-(1-delta)*dev-delta*threat_v(idp)<0);
    if sum(idx)==0 && uidp>=threat_v(idp) %No profitable deviation
        gamma_rest=gamma_rest+norm(uidp-threat_v(idp));
    else
        gamma_rest=gamma_rest-norm(uidp-threat_v(idp));
        valid=0;
    end
end
if valid==1
    valid_grid(end+1,:)=ac_test;
end
%Update global
setGlobalgrid(valid_grid,vc)

%output value

value_of_function=lambda*gamma_i+(1-lambda)*gamma_rest;

return;