% plot_psr_dist.m

clear, clf

load 'psr_dist_XYZ.mat'

for idx = 1:33;
plot3([0,x(idx)], [0,y(idx)], [0,z(idx)], '-*')
text(x(idx)+0.05,y(idx)+0.05,z(idx)+0.05, num2str(idx), 'FontSize', 14)
psrl(idx,:) = strcat(num2str(idx, '%02.f'),' : ', psr(idx,:));
hold on

end

plot3([0,0], [0,0], [0,0], 'ko', 'MarkerSize',20)
plot3([0,0], [0,0], [0,0], 'y*', 'MarkerSize',18)

grid on



legend(psrl, 'FontSize', 14, 'Location', 'northeastoutside')
xlabel('SSB Distance in kpc', 'FontSize', 16)
ylabel('SSB Distance in kpc', 'FontSize', 16)
zlabel('SSB Distance in kpc', 'FontSize', 16)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',18 , 'FontWeight','bold')


for az = 37.5:0.5:220;
    view(az, 30)
    pause(0.1)
end

for elv = 30:0.5:90;
    view(az, elv)
    pause(0.1)
end

for elv = elv:-0.5:37.5;
    view(az, elv)
    pause(0.1)
end

%view(37.5, 30)

% plt.savefig("movie%d.png" % ii)
