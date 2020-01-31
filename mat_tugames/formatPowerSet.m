function pwm=formatPowerSet(A)
% FORMATPOWERSET formats the Matlab cell output that contains the 
% representation of coalitions into matrix form.
%
%
% Usage: pwm=formatPowerSet(pws)
%
% Define variables:
%  output:
%  pwm      -- all subsets of  
%
%  input:
%  A        -- a set A like A=[2 3 4], or the power set pws.
% 
% Example:
% A=[2 3 4];
% pws = formatPowerSet(A)
% >> pwm=formatPowerSet([2 3 4])
%
% pwm =
%
%     0     0     0     4
%     0     0     3     0
%     0     0     3     4
%     0     2     0     0
%     0     2     0     4
%     0     2     3     0
%     0     2     3     4
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/01/2018        0.9             hme

if iscell(A)==0
   pws=PowerSet(A);
else 
   pws=A;
end

np=numel(pws);
n=length(pws{np});
pwm=zeros(np,n);
for k=1:np
    pwm(k,pws{k})=pws{k};
end


