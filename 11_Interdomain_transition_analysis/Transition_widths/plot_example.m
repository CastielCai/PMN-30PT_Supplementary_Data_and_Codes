% Fig. S15, representative D3-D6. Display only.
e=example_profile([3 6]); colors=[.06 .53 .53;.47 .30 .53];
figure('Color','w'); tiledlayout(2,3);
for j=1:2
    ax=nexttile(j); hold(ax,'on'); x=e.qr_A/10;
    plot(x,e.observed(:,j),'-o','Color',colors(j,:),'MarkerFaceColor',colors(j,:));
    plot(x,e.baseline(:,j),'--','Color',[.61 .34 .10]);
    plot(x,e.half_level(:,j),':k'); xline(e.bounds_A(j,:)/10,':k');
    if j==1, ylabel('POP'); else, ylabel('Polarization magnitude (µC/cm²)'); end
    title(sprintf('w = %.2f nm',diff(e.bounds_A(j,:))/10)); box on
    nexttile(j+3); plot(x,e.reduction(:,j),'-','Color',colors(j,:)); hold on
    yline(50,':k');xline(e.bounds_A(j,:)/10,':k');ylabel('Reduction (%)');xlabel('Distance (nm)');box on
end
nexttile(3);hold on; colors=[.84 .48 .12;.16 .48 .66];
for j=1:2
    scatter(e.qr_A/10,e.observed(:,j+2),12,colors(j,:),'filled');
    plot(e.grid_A/10,e.angle_fit(:,j),'Color',colors(j,:));
end
xline(e.bounds_A(3,:)/10,':k');ylabel('Angle (deg)');title(sprintf('w = %.2f nm',diff(e.bounds_A(3,:))/10));box on
nexttile(6);plot(e.grid_A/10,e.cumulative,'k-'); hold on; xline(e.bounds_A(3,:)/10,':k');
yline([12.5 87.5],':k');ylabel('Cumulative change (%)');xlabel('Distance (nm)');box on
