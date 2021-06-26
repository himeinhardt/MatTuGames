function ct=contentment(clv,x);
% Computes the contentment vector of game v w.r.t. x. It is a
% negative excess vector.
%
% Usage: ct=clv.contentment(x)
%
% Define variables:
% output:
% ex         -- contentment vector of game v w.r.t. x.
%
% input: 
%  clv        -- TuGame class object.
%  x        -- payoff vector of size(1,n).
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/11/2017        0.9             hme
%                

% Computing the excess vector w.r.t. x.
v=clv.tuvalues;
n=clv.tuplayers;
% Borrowed from J. Derks
Xm(1)=x(1);
for ii=2:n
    Xm=[Xm x(ii) Xm+x(ii)];
end
ct=Xm-v;
