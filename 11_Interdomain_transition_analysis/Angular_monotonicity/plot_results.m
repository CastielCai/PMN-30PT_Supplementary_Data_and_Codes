% Fig. 4F. Run Angular_monotonicity/main.m to regenerate the statistics.
here=fileparts(mfilename('fullpath'));
data=load(fullfile(here,'output','angular_monotonicity.mat'));
keep=data.in_figure; p=data.trend_parameters; theta=linspace(0,180,361);
figure('Color','w'); hold on
scatter(data.delta_theta_deg(keep),data.M(keep),16,[.13 .52 .79],'filled');
plot(theta,1-(1-p(1))*exp(-theta/p(2)),'k--','LineWidth',1.7);
xlim([0 180]); ylim([0 1.1]); axis square; box on
xlabel('\Delta\theta (°)'); ylabel('Monotonicity');
