%% Plot table in a competable way to the rest of the simulations:

t=readtable('MNIST_results.csv');
t.num_steps

figure(410);clf
pos=find(t.num_steps<501); % we ignore the 1000 conditions that were run, because this simulation had a techincal problem
mclr='k';
plot(t.num_steps(pos),t.FID(pos),'k-','MarkerFaceColor',mclr,'MarkerSize',10, 'LineWidth',2,'Color','k');hold on;
for I=1:length(pos)
    II=pos(I);
    mclr=0.9*mod([II*123213,(II+123)*12312,II*II+II+123+ II*1232131],256)/256+ [0.1,0.1,0.1];
    plot(t.num_steps(II),t.FID(II),'o','MarkerFaceColor',mclr,'MarkerSize',10, 'LineWidth',2,'Color',mclr);hold on;

end
set(gca,'FontSize',14);
ylabel('FID')
xlabel('Number of steps');
title('Performance (FID)')