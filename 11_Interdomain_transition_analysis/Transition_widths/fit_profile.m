function splines = fit_profile(x,y,count,bounded)
% Penalized cubic B-spline. Select a common lambda across columns by
% count-weighted leave-one-bin-out CV, refitting active coefficient bounds.
n=numel(x); w=count(:)/mean(count);
nb=min([16,max(4,floor(n/2)+2),n-1]);
if nb<4 || any(diff(x)<=0) || any(w<=0)
    error('Fit needs at least five distinct, occupied bins.');
end
p=linspace(0,1,nb-2);
inner=interp1(1:n,x,1+(n-1)*p(2:end-1),'linear');
knots=[repmat(x(1),1,4),inner(:)',repmat(x(end),1,4)];
basis=spmak(knots,eye(nb));
B=fnval(basis,x')';
breaks=unique(knots); half=diff(breaks)/2; mid=breaks(1:end-1)+half;
nodes=[-sqrt(3/5);0;sqrt(3/5)]; weights=[5/9;8/9;5/9];
zx=mid+nodes*half; zw=weights*half;
D=sqrt(zw(:)).*fnval(fnder(basis,2),zx(:)')'; D=D/norm(D,2);
G=B'*(w.*B); K=D'*D; rhs=B'*(w.*y);
bestScore=Inf; cap=min([12,.6*n,n-2]);
for lambda=logspace(-4,6,41)
    Q=G+lambda*K;
    if trace(Q\G)>cap+1e-8, continue; end
    C=Q\rhs; U=Q\(B'.*w');
    h=sum(B.*U',2);
    loo=(y-B*C)./max(1-h,1e-10);
    if bounded
        for i=1:n
            held=C-U(:,i)*loo(i,:);
            outside=find(any(held< -1e-8 | held>180+1e-8,1));
            for j=outside
                A=Q-w(i)*(B(i,:)'*B(i,:));
                b=rhs(:,j)-w(i)*B(i,:)'*y(i,j);
                coef=boxSolve(A,b,held(:,j));
                loo(i,j)=y(i,j)-B(i,:)*coef;
            end
        end
    end
    score=sum(w.*mean(loo.^2,2))/sum(w);
    if score<bestScore
        bestScore=score; bestQ=Q;
    end
end
if ~isfinite(bestScore), error('No eligible smoothing parameter.'); end
coef=bestQ\rhs; splines=cell(1,size(y,2));
for j=1:size(y,2)
    c=coef(:,j);
    if bounded, c=boxSolve(bestQ,rhs(:,j),c); end
    splines{j}=fn2fm(spmak(knots,c'),'pp');
end
end

function c=boxSolve(Q,b,unbounded)
% Active-set least squares on [0,180]; hessian is positive definite.
if all(unbounded>=0 & unbounded<=180), c=unbounded; return; end
c=max(0,min(180,unbounded)); active=zeros(size(c));
active(unbounded<=0)=-1; active(unbounded>=180)=1;
for step=1:200
    free=active==0; v=c;
    if any(free)
        v(free)=Q(free,free)\(b(free)-Q(free,~free)*c(~free));
    end
    low=free & v< -1e-10; high=free & v>180+1e-10;
    if any(low|high)
        ratio=inf(size(c));
        ratio(low)=-c(low)./(v(low)-c(low)); ratio(high)=(180-c(high))./(v(high)-c(high));
        [alpha,i]=min(ratio); c=c+max(0,min(1,alpha))*(v-c);
        if low(i), c(i)=0; active(i)=-1; else, c(i)=180; active(i)=1; end
        continue
    end
    c=max(0,min(180,v)); grad=Q*c-b;
    violation=max(0,grad.*active); [largest,i]=max(violation);
    if largest<=1e-9, return; end
    active(i)=0;
end
error('Bounded least-squares active set did not converge.');
end
