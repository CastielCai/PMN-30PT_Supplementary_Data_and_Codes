% Equation 2 uses the two fits already integrated in Transition_widths.
calculateM(fileparts(mfilename('fullpath')));

function calculateM(here)
root=fileparts(fileparts(here));
a=load(fullfile(root,'10_Polar_order_and_spatial_correlations','Domain_adjacency','output','adjacent_pairs.mat'),'pairs');
f=load(fullfile(root,'08_polar&domain&POP','output','polar_field.mat'),'core_direction');
w=load(fullfile(fileparts(here),'Transition_widths','output','transition_widths.mat'),'signed_net_deg','total_variation_deg');
denominator=sum(w.total_variation_deg,2);
M=sum(abs(w.signed_net_deg),2)./denominator;
M(denominator<=1e-6)=NaN;
directions=f.core_direction;
delta_theta_deg=acosd(max(-1,min(1,sum(directions(a.pairs(:,1),:).*directions(a.pairs(:,2),:),2))));
% D11-D14 was omitted in the manuscript's fixed Fig. 4F population.
in_figure=all(isfinite([delta_theta_deg M]),2) & ~ismember(a.pairs,[11 14],'rows');
spearman_rho=corr(delta_theta_deg(in_figure),M(in_figure),'Type','Spearman');
% Unweighted nonlinear least squares on the same 233 pairs as Fig. 4F.
% trend_parameters = [M at zero misorientation, angular scale in degrees].
trend=@(b,theta) 1-(1-b(1))*exp(-theta/b(2));
options=optimoptions('lsqcurvefit','Display','off','FunctionTolerance',1e-13,'OptimalityTolerance',1e-11);
trend_parameters=lsqcurvefit(trend,[.05 60],delta_theta_deg(in_figure),M(in_figure),[0 eps],[1 Inf],options);
out=fullfile(here,'output'); if ~isfolder(out), mkdir(out); end
save(fullfile(out,'angular_monotonicity.mat'),'M','delta_theta_deg','in_figure','spearman_rho','trend_parameters');
fprintf('Fig. 4F: %d pairs, Spearman rho = %.7f.\n',nnz(in_figure),spearman_rho);
end
