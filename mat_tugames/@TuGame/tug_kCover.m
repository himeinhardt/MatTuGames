function SOL=tug_kCover(clv,int_k)
% TUG_KCOVER computes the k-cover of game v with the Mathematica Package TuGames.
%
% Usage: SOL=tug_kCover(v,int_k)
% Define variables:
%  output:
%  SOL        -- The k-cover of game v.
%                Field variable gives result in Matlab and Mathematica format.
%  input:
%  clv        -- TuGame class object.
%  int_k      -- An integer samller than or equal to n.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   11/02/2012        0.3             hme
%   07/02/2021        1.9             hme
%


n=clv.tuplayers;

if nargin<2
  error('An integer not larger than n must be given!')
else
   if isnumeric(int_k)
      if length(int_k)>1
        error('Second argument must be an integer not larger than n!')
      end
   end
end

math('quit')
pause(1)
math('$Version')
try 
    math('{Needs["TUG`"] }'); 
catch 
    math('{Needs["coop`CooperativeGames`"],Needs["VertexEnum`"],Needs["TuGames`"],Needs["TuGamesAux`"] }'); 
end
disp('Passing Game to Mathematica ...')
w=gameToMama(clv);
math('matlab2math','mg1',w);
math('matlab2math','n1',n);
math('matlab2math','intk',int_k);
math('k=Rationalize[Flatten[intk,1][[1]]]');
math('bds=Flatten[n1][[1]]');
math('T=Flatten[Range[n1]]');
math('expg=Flatten[mg1,1]');
math('{T,mg=PrependTo[expg,0]}');
math('ExpGame:=(DefineGame[T,mg];);');
madw=math('adw=kCover[ExpGame,k]');
adw_v=math('math2matlab','adw');
adw_v(:,1)=[];
adw_g=gameToMatlab(adw_v);
SOL=struct('kCover',adw_g,'MkCover',madw);
math('quit')
