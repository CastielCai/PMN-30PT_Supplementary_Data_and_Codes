function window_A=valley_window(x_A,angles,count,L_A,cell_A)
% Reproduce the Fig. 4C valley-localization window from the observed angles.
% This auxiliary detector selects a bin; angular widths and M use B-splines.
x=x_A/cell_A; L=L_A/cell_A; weight=count/sum(count);
opts=optimset('Display','off','MaxIter',500,'TolX',1e-8);
best=Inf;
for fraction=[.25 .5 .75]
    [model,score]=fminsearch(@objective,[fraction*L,log(L/8)],opts);
    if score<best, best=score; chosen=model; end
end
window_A=cell_A*(chosen(1)+[-1 1]*sqrt(2)*erfinv(.8)*exp(chosen(2)));
    function loss=objective(z)
        c=z(1); sigma=exp(max(-15,min(15,z(2))));
        t=.5*(1+erf((x-c)/(sqrt(2)*sigma)));
        H=[ones(size(x)),t];
        if rcond(H'*(H.*weight))<1e-12, loss=1e12+sum(z.^2); return; end
        B=(H.*sqrt(weight))\(angles.*sqrt(weight));
        loss=sum(weight.*mean((angles-H*B).^2,2))...
            +1e4*(max(0,-c)^2+max(0,c-L)^2+max(0,sigma-2*L)^2+max(0,.02-sigma)^2)...
            +100*(max(0,-B(2,1))^2+max(0,B(2,2))^2);
    end
end
