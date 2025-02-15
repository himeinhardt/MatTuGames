function [bzv theta]=p_BanzhafColeman(clv);
% P_BANZHAFCOLEMAN computes the Banzhaf/Coleman index of a simple game sv 
% using MATLAB's PCT (normalized Banzhaf value by 2^n-1).
%
% Usage: [bzv theta]=clv.BanzhafColeman()
%
% Define variables:
%  output:
%  bzv       -- The Banzhaf/Coleman index.
%  theta     -- A vector of swings.
%
%  input:
%  clv      -- TuGame class object.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/13/2021        1.9.1           hme
%                

sv=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
gt=clv.tutype;
if strcmp(gt,'sv')
else
  error('Wrong game type!. Game must be a simple game!')
end

S=1:N;
bzv=zeros(1,n);
theta=zeros(1,n);

parfor k=1:n
  plk=bitget(S,k)==1;
  nplk=plk==0;
  Clk=S(plk);
  Clnk=S(nplk);
  w_clk=sv(Clk);
  nw_clk=sv(Clnk);
  sum_wk=w_clk*ones(length(w_clk),1);
  sum_nwk=nw_clk*ones(length(nw_clk),1);
  theta(k)=sum_wk-sum_nwk;
end
bzv=theta/N;
