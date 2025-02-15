function x=nucl_formula(v,tol)
% NUCL_FORMULA computes the nucleolus of a three-person super-additive game v from a formula. 
% If the game is not zero-normalized, it will be done internally. 
%    
% Resource: Leng and Parlar (2010) 
%    
% Usage: x=nucl_formula(v,tol)
% Define variables:
%  output:
%  x        -- Nucleolus of a three-person super-additive game and zero-normalized game.
%
%  input:
%  v        -- A three person super-additive TU-game of length 2^n-1.
%
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/15/2020        1.9             hme
%                
narginchk(1,2); % check for legal number of input arguments.

if nargin<2
 tol=10^8*eps;
end
N=length(v);
[~, n]=log2(N);
if n>3
     error('NuclThree:Exit','Game has not the correct size! Only applicable for three-person!');
    return
end

%% Zero-Normalization of the Game 
w=v;
v=zero_normalization(v);
n=3;
k=1:n;
N=2^n-1;
Nk=N-2.^(k-1);
x=-inf(1,3);
v1=v;
iS=2.^(k-1);
S=1:N;
if super_additiveQ(v)==0
     error('NuclThree:Exit','Game is not super-additive!'); 
    return
end

%% Checking conditions 
if CddCoreQ(v,tol)==0
   cdQ="Case: Empty Core"
    for ii=1:3
        vnjk=v(Nk(bitget(Nk,ii)==0));
        vik=v(Nk(bitget(Nk,ii)==1));    
        x(ii)=(v(N)+sum(vik)-2*vnjk)./3;
    end	
else
cd1=3*v(Nk)<=v(N);    
cd2=~cd1; 

 if any(cd2)
   cij=Nk(cd2);
   lc=length(cij);
   if lc==1
      ck=k(bitget(cij,k)==1);
      pk=k(bitget(cij,k)==0);
      cmk=Nk(bitget(Nk,pk)==1);
      vmk=v(cmk);        
      cd2_grQ=all(v(N)>=v(cij)+2.*vmk);
      if cd2_grQ==0
          ck=k(bitget(cij,k)==1);
          pk=k(bitget(cij,k)==0);
          cmk=Nk(bitget(Nk,pk)==1);
          vmk=v(cmk);
          cd3_grQ=v(N)>=v(cij)+2.*vmk;
          cd4_grQ=v(N)+v(cij)>=2*sum(vmk);
      else    
      end
   else
      ll=1:lc; 
      cd2_grQ0=false(1,lc);
      ck0=cell(lc,1);
      pk0=zeros(1,lc);
      for ii=1:lc
          kpl=k(bitget(cij(ii),k)==1);
          nkpl=k(bitget(cij(ii),k)==0);
          cmk=Nk(bitget(Nk,nkpl)==1);
          vmk=v(cmk);
          cd2_grQ0(ii)=all(v(N)>=v(cij(ii))+2.*vmk);
          ck0{ii}=kpl;
          pk0(ii)=nkpl;
      end
      ps=ll(cd2_grQ0==1);
      if isempty(ps)
          cd2_grQ=false;
          [mx,idx]=max(v(cij));
          ck=k(bitget(cij(idx),k)==1);
          pk=k(bitget(cij(idx),k)==0);
          cmk=Nk(bitget(Nk,pk)==1);
          vmk=v(cmk);
          % kij=cij(idx)
          v(cij(idx))+2.*vmk;
          cd3_grQ=v(N)>=v(cij(idx))+2.*vmk;
          cd4_grQ=v(N)+v(cij(idx))>=2*sum(vmk);
          cij=cij(idx);
      elseif length(ps)==1    
         cij=cij(ps);
         ck=ck0{ps};
         pk=pk0(ps);
         cd2_grQ=cd2_grQ0(ps);
      elseif length(ps)==2   
      end   
   end    
 end


%%% Assigning the payoffs to the players. 
    if all(cd1)
        cdQ="Case 1"
        x=ones(1,n)*v(N)/3;
    elseif any(cd2) & cd2_grQ
        cdQ="Case 2"
        x(ck)=(v(N)+v(cij))/4;
        x(pk)=(v(N)-v(cij))/2;
    elseif cd2_grQ == 0 & any(cd3_grQ)
        cdQ="Case 3"
        slc=cmk(cd3_grQ==0);
        ip=k(bitget(bitand(cij,slc),k)==1);
        cip=ck(ismember(ck,ip)==0);
        x(cip)=(v(N)-v(slc))/2;
        x(ip)=(v(cij)+v(slc))/2;
        x(pk)=(v(N)-v(cij))/2;
    elseif cd2_grQ == 0 & all(cd3_grQ == 0) & cd4_grQ
        cdQ="Case 4"
        x(ck(1))=(v(N)+v(cij)-2*(vmk(1)-vmk(2)))/4;
        x(ck(2))=(v(N)+v(cij)+2*(vmk(1)-vmk(2)))/4;
        x(pk)=(v(N)-v(cij))/2;
    elseif cd4_grQ==0
        cdQ="Case 5"
        x(ck(1))=(v(N)+v(cij)+vmk(2)-2*vmk(1))/3;
        x(ck(2))=(v(N)+v(cij)+vmk(1)-2*vmk(2))/3;
        x(pk)=(v(N)+sum(vmk)-2*v(cij))/3;
    else    
    end
end    
%%% Reformatting solution vector x w.r.t. the original game setting. 
x=x+w(iS);

end

