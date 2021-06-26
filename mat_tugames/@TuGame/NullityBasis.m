function NB=NullityBasis(clv,tol)
% NULLITYBASIS determines a basis of the game space for which n is
% a null player. 
%
% Source: L. Hernández-Lamoneda et al. (2007), Dissection of solutions in cooperative game theory
% using representation techniques, IJGT 35:395–426.
% DOI: 10.1007/s00182-006-0036-3
%
% Usage: NB=clv.NullityBasis(tol)
%
%  output:
%  Define variables of the structure element:
%  bVj       -- Basis of the nullity game space N inside of C + T.
%  bMub      -- Some Basis of space of the orthogonal complement of N inside of space of C + T. 
%  bXi       -- Remaining basis vector of the orthogonal complement of N inside of space of C + T.
%  spN       -- Subspace N of C+T.
%  spNC      -- Orthogonal complement of N.
%  bvj       -- Alternative basis of the nullity game space N inside of C + T.
%  bMub      -- Alternative basis of space of the orthogonal complement of N inside of space of C + T. 
%  bxi       -- Alternative basis vector of the orthogonal complement of N inside of space of C + T.
%  spN2      -- Alternative basis of subspace N of C+T.
%  spNC2     -- Alternative orthogonal complement of N.
%  vn        -- Projection of game v on N inside of space of C + T.
%  vb        -- Projection of game v on the orthogonal complement of N inside of space of C + T.
%  vcn       -- Projection of game v on the orthogonal complement of N..
%  egQ       -- Returns true (1), if v=vb+vn+vcn is satisfied, otherwise false (0).
%  v         -- Original game v.
%
%  input:
%  clv      -- TuGame class object.
%  tol      -- Tolerance value. Its default value is set to 10^5*eps.
%


%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   11/20/2015        0.8             hme
%

N=clv.tusize;
n=clv.tuplayers;
v=clv.tuvalues;


if nargin < 2
   tol=10^5*eps; 
end    

vj=zeros(n-1,N);
S=1:N;
ov=ones(n,1);
it=0:-1:1-n;
matC=rem(floor(S(:)*pow2(it)),2);
sc=matC*ov;
Cj=zeros(n,N);

exg=cell(1,n);
for ii=1:n
    Cj(ii,:)=sc==ii; 
    for jj=1:n
        exom=zeros(1,N);
        sK=S(sc==jj);
        sKn=sK(bitget(sK,n)==1);
        sKwn=sK(bitget(sK,n)==0);
        exom(sKwn)=jj;
        exom(sKn)=jj-n;
        exg{ii}(jj,:)=exom;
    end
end
Xi=zeros(1,N);
k=1:n;
ki=2.^(k-1);
Xi(ki(n))=1;
vxi=Xi*(v*Xi')/(Xi*Xi');
Mu=zeros(n-1,N);
vj=zeros(n-1,N);
vm=zeros(1,N);
vn=zeros(1,N);
for jj=1:n-1
    Vj(jj,:)= (n-jj)*Cj(jj,:)/n + (jj+1)*Cj(jj+1,:)/n + (1/n)*(exg{n}(jj,:)-exg{n}(jj+1,:));
    sK2=S(sc==jj);
    sK3=S(sc==jj+1);
    sKn3=sK3(bitget(sK3,n)==1);
    sKwn2=sK2(bitget(sK2,n)==0);
    Mu(jj,sKwn2)=1;
    Mu(jj,sKn3)=-1;
    vj(jj,sKwn2)=1;
    vj(jj,sKn3)=1;
    vn=Vj(jj,:)*(v*Vj(jj,:)')/(Vj(jj,:)*Vj(jj,:)')+vn;
    vm=Mu(jj,:)*(v*Mu(jj,:)')/(Mu(jj,:)*Mu(jj,:)')+vm;
    Mub(jj,:)= (n-jj)*Cj(jj,:)/n - (jj+1)*Cj(jj+1,:)/n + (1/n)*(exg{n}(jj,:)+exg{n}(jj+1,:));
end
vb=vxi+vm;
Xin=(Cj(1,:)-exg{n}(1,:))*1/n;
N=[Vj;Mub;Xin];
cN=null(N);
vcn=(pinv(cN')*cN'*v')';
N2=[vj;Mu;Xi];
cN2=null(N2);
v2=vb+vn+vcn;
egQ=all(abs(v2-v)<tol);
%% Formatting output.
NB=struct('bVj',Vj,'bMub',Mub,'bXi',Xi,'spN',N,'spCN',cN,'bvj',vj,'bMu',Mu,'bxi',Xin,'spN2',N2,'spCN2',cN2,'vn',vn,'vb',vb,'vcn',vcn,'egQ',egQ,'v',v);
    
