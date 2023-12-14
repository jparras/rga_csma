function [value_of_function]=function_cost_soo_corr(s)

param=getGlobalx;
lambda=param.lambda;
valid_grid=param.valid_grid;
vc=param.vc;
delta=param.delta;
threat_v=param.threat_v;
i=param.player_optimizing;
npl=param.npl;
param_u=getGlobalu;
u=struct2cell(param_u);

% Map s to the unit simplex
sm=sum(s);
m=max(s);
s=s*m/sm;

if sum(s)<=1 % valid simplex!

    gamma_i=[];
    gamma_rest=0;
    valid=1;
    % Obtain correlated actions indexes
    action=zeros(2^npl,npl);
    for k=1:2^npl %Encode actions as binary: each player only has 2!!
        action_string=dec2bin(k-1,npl);
        action_s=action_string-'0';
        action(k,:)=action_s+1;
    end
    % Check action for player i
    phi_test=[s 1-sum(s)];
    payoff=obtain_payoff_corr(npl,u,phi_test);
    ui=payoff(i);
    % Add to region tested
    vc(end+1,:)=payoff;
    % Check if satisfies pi restrictions
    r1=0; % One accumulator per action (2 per player!)
    r2=0;
    for j=1:size(action,1)
        pure_action_dev=action(j,:); %Action from which deviate
        pa_ndev_index=num2cell(pure_action_dev);
        if pure_action_dev(i)==1 %Player i plays api=1
            pure_action_dev(i)=2;
            pa_index=num2cell(pure_action_dev); %For indexing purposes
            udev=(1-delta)*u{i}(pa_index{:})+delta*threat_v(i);
            undev=(1-delta)*u{i}(pa_ndev_index{:})+delta*ui;
            aux=phi_test(j)*(undev-udev);
            r1=r1+aux;
        else
            pure_action_dev(i)=1;
            pa_index=num2cell(pure_action_dev); %For indexing purposes
            udev=(1-delta)*u{i}(pa_index{:})+delta*threat_v(i);
            undev=(1-delta)*u{i}(pa_ndev_index{:})+delta*ui;
            aux=phi_test(j)*(undev-udev);
            r2=r2+aux;
        end
    end
    
    if r1>=0 && r2>=0 && ui>=threat_v(i) %No profitable deviation
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
        payoff=obtain_payoff_corr(npl,u,phi_test);
        uidp=payoff(idp);
        % Add to region tested
        vc(end+1,:)=payoff;
        % Check if satisfies other players restrictions
        r1=0; % One accumulator per action (2 per player!)
        r2=0;
        for j=1:size(action,1)
            pure_action_dev=action(j,:); %Action from which deviate
            pa_ndev_index=num2cell(pure_action_dev);
            if pure_action_dev(i)==1 %Player i plays api=1
                pure_action_dev(i)=2;
                pa_index=num2cell(pure_action_dev); %For indexing purposes
                udev=(1-delta)*u{i}(pa_index{:})+delta*threat_v(i);
                undev=(1-delta)*u{i}(pa_ndev_index{:})+delta*ui;
                aux=phi_test(j)*(undev-udev);
                r1=r1+aux;
            else
                pure_action_dev(i)=1;
                pa_index=num2cell(pure_action_dev); %For indexing purposes
                udev=(1-delta)*u{i}(pa_index{:})+delta*threat_v(i);
                undev=(1-delta)*u{i}(pa_ndev_index{:})+delta*ui;
                aux=phi_test(j)*(undev-udev);
                r2=r2+aux;
            end
        end
        if r1>=0 && r2>=0 && uidp>=threat_v(idp) %No profitable deviation
            gamma_rest=gamma_rest+norm(uidp-threat_v(idp));
        else
            gamma_rest=gamma_rest-norm(uidp-threat_v(idp));
            valid=0;
        end
    end
    if valid==1
        valid_grid(end+1,:)=phi_test;
    end
    %Update global
    setGlobalgrid(valid_grid,vc)

    %output value

    value_of_function=lambda*gamma_i+(1-lambda)*gamma_rest;
else
    display('Fail to find simplex');
    value_of_function=-2*abs(sum(s)-1)*max(max(max(max(max([u{:}]))))); %High value!!
end

return;