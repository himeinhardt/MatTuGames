function [v,wpn,wpk]=game_Wsys(pS,str,tol)
% GAME_WSYS creates a set of games from an asymmetric weight system.
%  
%  Usage: v=game_Wsys(pS,str,tol)
%
% Define variables:
% output:
%  v        -- A set of Tu game of length 2^n-1. Output is a cell containing the games.
%  wpn      -- The set of weighted nucleoli (cell output).
%  wpk      -- The set of weighted Pre-Kernels (cell output).
%
% input: 
%  pS       -- An asymmetric weight sytem of length 2^n-1.
%  str      -- Permissible strings are 
%              'type1'   game type given by Kleppe et al (2013).
%              'type2'   game type given by Solymosi (2015) of case 1.
%              'type3'   game type given by Solymosi (2015) of case 1 and 2.
%                        the last game is of case 2.
%               Default is 'type3'.
%
%  tol      -- Tolerance value. Its default value is set to 10^6*eps.
%
%  Example:
%          Generates a set of five person games.
%
%          pS=[1 1 1 1 1 1 1 1 1 1 7 1 7 7 1 1 1 1 7 1 7 7 1 1 1 1 1 1 1 1 1];
%          [v,wpn,wpk]=game_Wsys(pS,'type1')
%
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   06/26/2015        0.7             hme
%

if nargin < 1
   error('At least an asymmetric weight system pS must be given!');
elseif nargin < 2
   game_type='type3';
   tol=10^6*eps;
elseif nargin < 3
   tol=10^6*eps;
   game_type=str;
else
   game_type='type3';
end

if symmetricQ(pS)
  error('An asymmetric weight system pS must be given!');
else
  AsC=FindAsymmetricCoalition(pS,tol);
end
nA=numel(AsC);
N=length(pS);
[~, n]=log2(N);


Tc=cell(1,nA);
ltc=zeros(1,nA);
ll=1:nA;

for k=1:nA
    if isempty(AsC{k})==0
       Tc{k}=unique(bitand(AsC{k},AsC{k}(end)));
       ltc(k)=length(Tc{k});
   end
end
idx=ll(ltc>0);
ml=max(ltc);
li=length(idx);
v=cell(li,ml);
wpn=cell(li,ml);
wpk=cell(li,ml);

switch game_type
     case 'type1'
         for k=1:li
             if isempty(Tc{idx(k)})==0
                nt=length(Tc{idx(k)});
                for ii=1:nt
                   v{k,ii}=game_wghs1(Tc{idx(k)}(1,ii),pS,n);
                   wpn{k,ii}=cplex_weightedPreNucl(v{k,ii},pS);
                   wpk{k,ii}=weightedPreKernel(v{k,ii},pS);
                end
             end
         end
     case 'type2' % {'type2','type3'}
         for k=1:li
             if isempty(Tc{idx(k)})==0
                nt=length(Tc{idx(k)});
                T=Tc{idx(k)};
                for ii=1:nt
                   v{k,ii}=game_wghs2(T(ii),pS,n);
                   wpn{k,ii}=cplex_weightedPreNucl(v{k,ii},pS);
                   wpk{k,ii}=weightedPreKernel(v{k,ii},pS);
                end
             end
         end
     case 'type3'
         for k=1:li
             if isempty(Tc{idx(k)})==0
                nt=length(Tc{idx(k)})+1;
                T=[Tc{idx(k)},0];
                for ii=1:nt
                   v{k,ii}=game_wghs3(T(ii),pS,n);
                   wpn{k,ii}=cplex_weightedPreNucl(v{k,ii},pS);
                   wpk{k,ii}=weightedPreKernel(v{k,ii},pS);
                end
             end
         end
end
