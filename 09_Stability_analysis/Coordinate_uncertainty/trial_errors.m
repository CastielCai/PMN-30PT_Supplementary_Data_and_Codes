function [lead,bsite] = trial_errors(stream,k,nLead,nB)
% Replay trial k from the archived PCG64 starting state. Arithmetic below
% is native MATLAB; the normal layers and stream constants are numeric data.
% Each lead atom is sampled once, including when shared by several B sites.
nNormal=3*(nLead+nB); block=32768;
nBlocks=ceil((1.08*nNormal+nB)/block);
state=double(stream.initial(k,:)); increment=double(stream.increment(k,:));
[words,state]=blocks(state,increment,double(stream.multiplier),block,nBlocks);
while true
    [lead,bsite,enough]=decode(words,stream,nLead,nB);
    if enough, return; end
    [extra,state]=blocks(state,increment,double(stream.multiplier),block,1);
    words=[words;extra]; %#ok<AGROW>
end
end

function [raw,state]=blocks(state,inc,mult,block,nBlocks)
% Advance the 128-bit LCG modulo 2^128, using eight exact base-65536 limbs.
% The powers and geometric sums for all offsets follow by doubling.
persistent powers sums cacheMult
if isempty(powers) || ~isequal(cacheMult,mult) || size(powers,1)~=block
    powers=mult; sums=[1 zeros(1,7)]; cacheMult=mult;
    while size(powers,1)<block
        endP=powers(end,:); endS=sums(end,:);
        nextP=carry(powers*productMatrix(endP));
        nextS=carry(powers*productMatrix(endS)+sums);
        powers=[powers;nextP]; sums=[sums;nextS]; %#ok<AGROW>
    end
end
incrementPart=sums*productMatrix(inc);
raw=zeros(block*nBlocks,1,'uint64');
for b=1:nBlocks
    values=carry(powers*productMatrix(state)+incrementPart);
    lo=pack64(values(:,1:4)); hi=pack64(values(:,5:8)); rot=double(bitshift(hi,-58));
    u=bitxor(hi,lo);
    raw((b-1)*block+(1:block))=bitor(bitshift(u,-rot),bitshift(u,64-rot));
    state=values(end,:);
end
end

function T=productMatrix(a)
T=triu(toeplitz([a(1);zeros(7,1)],a));
end
function a=carry(a)
for j=1:7
    c=floor(a(:,j)/65536); a(:,j)=a(:,j)-65536*c; a(:,j+1)=a(:,j+1)+c;
end
a(:,8)=mod(a(:,8),65536);
end
function z=pack64(v)
lo=uint64(v(:,1)+65536*v(:,2)); hi=uint64(v(:,3)+65536*v(:,4));
z=bitor(lo,bitshift(hi,32));
end

function [a,b,ok]=decode(words,s,nA,nB)
u=double(bitshift(words,-11))*2^-53;
layer=double(bitand(words,uint64(255)))+1;
absBits=bitand(bitshift(words,-9),uint64(2^52-1));
v=double(absBits).*s.wi(layer);
negative=bitget(words,9)~=0; v(negative)=-v(negative);
fast=absBits<s.ki(layer);
slow=find(~fast); slow(end+1)=numel(words)+1;
[a,index,next,ok]=normals(3*nA,1,1);
if ~ok || index+nB-1>numel(words), b=[];ok=false;return; end
large=u(index:index+nB-1)>=.30; index=index+nB;
[b,~,~,ok]=normals(3*nB,index,next);
if ok
    a=10.4*reshape(a,3,nA)';
    b=reshape(b,3,nB)'.*(7.7+10.4*large);
end
    function [out,j,is,success]=normals(n,j,is)
        out=zeros(n,1); written=0; success=false; R=3.6541528853610088;
        while written<n
            while is<numel(slow) && slow(is)<j, is=is+1; end
            amount=min(n-written,slow(is)-j);
            if amount>0
                out(written+(1:amount))=v(j:j+amount-1);
                written=written+amount;j=j+amount;
            end
            if written==n, success=true; return; end
            if j>numel(v)-2, return; end
            current=j; j=j+1; is=is+1; l=layer(current); x=v(current);
            if l==1
                while true
                    if j+1>numel(u), return; end
                    tail=-log1p(-u(j))/R; h=-log1p(-u(j+1)); j=j+2;
                    if 2*h>tail*tail, break; end
                end
                x=R+tail; if bitget(absBits(current),9), x=-x; end
            else
                accept=(s.fi(l-1)-s.fi(l))*u(j)+s.fi(l)<exp(-0.5*x*x);
                j=j+1; if ~accept, continue; end
            end
            written=written+1; out(written)=x;
        end
        success=true;
    end
end
