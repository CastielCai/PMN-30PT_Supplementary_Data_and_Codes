% Distribution of measurable POP half-depth widths (Fig. 4C).
here=fileparts(mfilename('fullpath'));
data=load(fullfile(here,'output','transition_widths.mat'),'width_A');
width=data.width_A(:,1); width=width(isfinite(width));
figure('Color','w');
histogram(width,0:4:72,'FaceColor',[.22 .62 .64],'EdgeColor','k');
axis square; box on; xlabel('w_{POP} (Å)'); ylabel('Counts');
