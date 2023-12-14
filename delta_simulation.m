%% REPEATED GAME delta comparison in the CSMA/CA game
% Juan Parras, GAPS-UPM, November 2019
clear all; clc; close all;
%% Main parameters
R1 = [-0.3668 0; 0.2668 -0.1];
R2 = [0.3608 0; -0.1617 0];
u1=@(y,z) [y, 1-y]*R1*[z, 1-z]';
u2=@(y,z) [y, 1-y]*R2*[z, 1-z]';
v1n=-R1(1,1)*R1(2,2)/sum(sum(abs(R1)));
v2n=0;
yn = -R2(2,1) / (R2(1,1) - R2(2,1));
zn = -R1(2,2) /sum(sum(abs(R1)));

n_delta = 100;
n_ac = 100;
delta_v = linspace(0.5, 1, n_delta);
y_v = linspace(0, 1, n_ac);
z_v = linspace(0, 1, n_ac);
regions = zeros(n_delta, 2, n_ac * n_ac); % Delta x (y,z) x payoffs

for id=1:n_delta
    delta = delta_v(id);
    display(['Case delta = ' num2str(delta)]);
    p1 = v1n * ones(length(y_v) * length(z_v), 1);
    p2 = v2n * ones(length(y_v) * length(z_v), 1);
    idp = 1;
    for iy=1:n_ac
        y = y_v(iy);
        for iz=1:n_ac
            z = z_v(iz);
            u = [u1(y,z), u2(y,z)];
            [um1, um2] = max_payoffs(R1, R2, y, z, yn, zn);
            cond1 = u(1) >= (1 - delta) * um1 + delta * v1n;
            cond2 = u(2) >= (1 - delta) * um2 + delta * v2n;
            if cond1 && cond2
                p1(idp) = u(1);
                p2(idp) = u(2);
            end
            idp = idp + 1;
        end
    end
    regions(id, 1, :) = p1;
    regions(id, 2, :) = p2;
end
% Plot the maximum payoff obtained for each player
vals = [max(regions(:,1, :),[], 3), max(regions(:,2, :),[], 3)];
plot(delta_v, v1n * ones(n_delta), 'b--') % Player 1, static
hold on
plot(delta_v, v2n * ones(n_delta), 'r--') % Player 2, static
plot(delta_v, vals(:, 1), 'b') % Player 1, repeated
plot(delta_v, vals(:, 2), 'r') % Player 2, repeated
axis([delta_v(1) delta_v(n_delta) -0.1 0.1]);
xlabel('\delta');
ylabel('Payoff');
hold off
matlab2tikz('delta_sim.tikz');
save('Values_delta_paper');

% %% Plot the regions obtained
% cm = colormap(jet(n_delta)); 
% for id=n_delta:-1:1
%     delta = delta_v(id);
%     plot(squeeze(regions(id,1,:)), squeeze(regions(id,2,:)), 'o', 'Color', cm(id, :));
%     hold on
% end
% 
% %% Plot the regions obtained (convex hull)
% cm = colormap(jet(n_delta)); 
% for id=1:n_delta
%     delta = delta_v(id);
%     p = [squeeze(regions(id,1,:)), squeeze(regions(id,2,:))];
%     cond = length(unique(p)) == 2;
%     display(['Case id = ' num2str(id) ', delta = ' num2str(delta) ', cond = ' num2str(cond)]);
%     if cond
%         plot(p(1,1), p(1,2), 'o', 'Color', cm(id, :));
%     else
%         k = convhull(p);
%         plot(p(k,1), p(k,2), 'Color', cm(id, :));
%     end
%     hold on
% end