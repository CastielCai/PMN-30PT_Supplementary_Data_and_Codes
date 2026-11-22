function [width,net,tv,edges,cuts,mass] = angular_change(curve,lo,hi)
% Total absolute change on exactly monotone subintervals of both cubics.
cuts=[lo;hi];
for j=1:2
    [br,co,np]=unmkpp(curve{j});
    for i=1:np
        r=roots(polyder(co(i,:)));
        t=real(r(abs(imag(r))<1e-9));
        t=t(t>0 & t<br(i+1)-br(i));
        cuts=[cuts;t(:)+br(i)]; %#ok<AGROW>
    end
end
cuts=unique(cuts(cuts>=lo & cuts<=hi));
vals=[ppval(curve{1},cuts),ppval(curve{2},cuts)];
net=vals(end,:)-vals(1,:); tv=sum(abs(diff(vals,1,1)),1);
mass=[0;cumsum(sum(abs(diff(vals,1,1)),2))];
width=NaN; edges=[NaN NaN];
if mass(end)<=1e-6, return; end
edges=zeros(1,2); fractions=[.125 .875];
for j=1:2
    target=fractions(j)*mass(end);
    i=find(mass>=target,1);
    if mass(i)==target, edges(j)=cuts(i); continue; end
    a=i-1;
    advance=@(s) sum(abs([ppval(curve{1},s) ppval(curve{2},s)]-vals(a,:)));
    edges(j)=fzero(@(s) mass(a)+advance(s)-target,[cuts(a) cuts(i)]);
end
width=diff(edges);
end
