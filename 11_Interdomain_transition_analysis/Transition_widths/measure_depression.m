function [bounds,method] = measure_depression(x,y,n,window,hasVoid,binWidth)
% Fig. 4C rule: observe the depression and interpolate crossings first.
% A count-weighted single-valley fit may supply otherwise unresolved widths.
bounds=[NaN NaN]; method=uint8(0);
good=n>0 & isfinite(y); ids=find(good);
if numel(ids)<3, return; end
base=@(s)y(ids(1))+(y(ids(end))-y(ids(1)))*(s-x(ids(1)))/(x(ids(end))-x(ids(1)));
if min(base(x(ids)))<=0, return; end
depth=1-y./base(x); candidate=find(good & x>=window(1) & x<=window(2));
if ~isempty(candidate)
    [height,j]=max(depth(candidate)); peak=candidate(j); half=height/2;
    lo=peak; hi=peak;
    while lo>1 && good(lo-1) && depth(lo-1)>half, lo=lo-1; end
    while hi<numel(x) && good(hi+1) && depth(hi+1)>half, hi=hi+1; end
    outside=good; outside(max(1,lo-1):min(numel(x),hi+1))=false;
    isolated=~any(depth(outside)>=.8*height);
    if height>=.10 && isolated && lo>1 && hi<numel(x) && good(lo-1) && good(hi+1)
        left=crossing(x(lo-1:lo),depth(lo-1:lo),half);
        right=crossing(x(hi:hi+1),depth(hi:hi+1),half);
        if right-left>=binWidth
            bounds=[left right];method=uint8(1);return
        end
    end
end
x=x(ids); v=y(ids); w=n(ids);
% Freeze endpoints before calling the interpolant on the reduced observations.
reference=@(s)v(1)+(v(end)-v(1))*(s-x(1))/(x(end)-x(1));
if hasVoid || max(diff(x))>2*binWidth+1e-7, return; end
[z,loss]=singleValley(v,w); amp=max(v)-min(v);
if sqrt(loss/sum(w))/max(amp,1e-12)>.10, return; end
valley=find(abs(z-min(z))<=1e-10+1e-9*abs(min(z))); arms={1:valley(1),valley(end):numel(x)};
for j=1:2
    arm=arms{j}; endpoint=z(1+(j-1)*(numel(z)-1)); contrast=endpoint-min(z);
    if numel(arm)<3 || contrast<=1e-10 || contrast/max(abs(endpoint),1e-8)<.10, return; end
    if sqrt(sum(w(arm).*(z(arm)-v(arm)).^2)/sum(w(arm)))/contrast>.20, return; end
end
fun=pchip(x,z); reduction=@(s)1-ppval(fun,s)./reference(s);
dense=linspace(x(1),x(end),10001); dd=reduction(dense); [~,k]=max(dd);
[height,at]=max(reduction(x)); peak=x(at);
if k>1 && k<numel(dense)
    [s,yneg]=fminbnd(@(t)-reduction(t),dense(k-1),dense(k+1),optimset('TolX',1e-12));
    if -yneg>height, height=-yneg; peak=s; end
end
if height<.10, return; end
half=height/2; mid=find(dense>=peak,1);
a=mid;b=mid;
while a>1 && dd(a)>half, a=a-1; end
while b<numel(dd) && dd(b+1)>half, b=b+1; end
if a==1 || b==numel(dd), return; end
outside=true(size(dd)); outside(max(1,a-1):min(numel(dd),b+1))=false;
if any(dd(outside)>=.8*height), return; end
left=fzero(@(s)reduction(s)-half,dense([a a+1]));
right=fzero(@(s)reduction(s)-half,dense([b b+1]));
if right-left>=binWidth, bounds=[left right]; method=uint8(2); end
end

function value=crossing(x,y,h)
value=x(1)+(h-y(1))*(x(2)-x(1))/(y(2)-y(1));
end

function [best,loss]=singleValley(y,w)
loss=Inf; n=numel(y);
for k=1:n
    % Pool each branch in outward order, then pool both shared valley blocks.
    branches={pava((k-1):-1:1,y,w),pava(k+1:n,y,w)};
    z=y; blocks=k; total=y(k)*w(k); weight=w(k);
    for j=1:2
        arm=branches{j};
        for m=1:numel(arm), z(arm(m).ids)=arm(m).value; end
    end
    while true
        near=[Inf Inf];
        for j=1:2, if ~isempty(branches{j}), near(j)=branches{j}(1).value; end, end
        [v,j]=min(near); if v>=total/weight-1e-13, break; end
        entry=branches{j}(1); branches{j}(1)=[];
        total=total+entry.value*entry.weight; weight=weight+entry.weight;
        blocks=[blocks entry.ids]; %#ok<AGROW>
    end
    z(blocks)=total/weight;
    score=sum(w.*(y-z).^2);
    if score<loss, loss=score; best=z; end
end
end

function b=pava(order,y,w)
b=struct('ids',{},'value',{},'weight',{});
for i=order
    b(end+1)=struct('ids',i,'value',y(i),'weight',w(i)); %#ok<AGROW>
    while numel(b)>1 && b(end-1).value>=b(end).value
        a=b(end-1); c=b(end); ww=a.weight+c.weight;
        b(end-1)=struct('ids',[a.ids c.ids],'value',(a.value*a.weight+c.value*c.weight)/ww,'weight',ww);
        b(end)=[];
    end
end
end
