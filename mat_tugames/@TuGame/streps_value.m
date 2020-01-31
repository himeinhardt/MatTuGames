function v_eps=streps_value(clv,t);
% STREPS_VALUE computes the game v-t.
% 
% Usage: v_eps=streps_value(v,t)
%

% Define variables:
%  output:
%  v_eps    -- Returns the characteristic function of v-t.
%
%  input:
%  clv      -- TuGame class object.
%  t        -- epsilon value (cost) to form a proper subcoalition. 


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   06/19/2013        0.4             hme
%                

v=clv.tuvalues;

v_eps=v-t;
v_eps(end)=v(end);


