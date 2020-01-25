function PlayCoreMovie(fm,mp,fps)
% PLAYCOREMOVIE plays a movie from a collection of frames stored in a
% structure element fm.
%
%   Usage: PlayCoreMovie(fm,mp,fps) 
%
%  output     
%             -- A movie sequence of the strong epsilon cores
%                to illuminate an important (pre-)kernel property.
%  input:
%  fm         -- A frame or a sequence of frames like
%      
%                MCrMov (structure element)
%
%                which contains the frames of a core movie.
%  mp         -- Plays the movie mp times, default is 1.
%  fps        -- Plays the movie at fps frames per second, 
%                default is 12 (Matlab).
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/05/2014        0.6             hme
%

if nargin<2
   mp=1;
   fps=25;
else
   fps=25;
end

fig1=figure(1);
winsize = get(fig1,'Position');
winsize(1:2) = [0 0];
movie(gcf,fm,mp,fps,winsize);
