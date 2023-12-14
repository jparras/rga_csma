function [] = print_error_area(datas, datac, color, name)
% Data RM
figure();
x_vector = [(1:4), fliplr((1:4))];
patch = fill(x_vector, [datas(:,2,1)',fliplr(datas(:,3,1)')], color(1));
set(patch, 'edgecolor', 'none');
set(patch, 'FaceAlpha', 0.3);
hold on;
plot((1:4), datas(:,1,1)', 'color', color(1),'LineWidth', 2);
% Data CA
patch = fill(x_vector, [datac(:,2,1)',fliplr(datac(:,3,1)')], color(2));
set(patch, 'edgecolor', 'none');
set(patch, 'FaceAlpha', 0.3);
plot((1:4), datac(:,1,1)', 'color', color(2),'LineWidth', 2);
patch = fill(x_vector, [datac(:,2,2)',fliplr(datac(:,3,2)')], color(3));
set(patch, 'edgecolor', 'none');
set(patch, 'FaceAlpha', 0.3);
plot((1:4), datac(:,1,2)', 'color', color(3),'LineWidth', 2);
hold off
matlab2tikz(name)