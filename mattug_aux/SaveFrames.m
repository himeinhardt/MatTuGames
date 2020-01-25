function SaveFrames(varargin)
% SAVEFRAMES saves the frames of a movie to different file formats
%  
%   Usage: SaveFrames(fm,'png')
%
%  output:
%             -- Saves sequence of frames to different file formats  
%
%  input:     At most two input arguments are admissible. Any order of
%             the input arguments below is allowed. Nevertheless, at
%             least a frame or a sequence of frames must be given.
%
%  fm         -- A frame or a sequence of frames like
%      
%                MCrMov (structure element)
%
%                which contains the frames of a core movie.
%
%  format     -- A file format to save the frames. 
%                Permissible formats are:
%                'jpg' to save the sequence of images in jpg-format. 
%                   Useful to produce your own QuickTime movie. 
%                'eps' to save the sequence of images in eps-format.
%                'png' to save the sequence of images in png-format.
%                   Useful to produce your own movie with Blender.
%                'pdf' to save the sequence of images in pdf-format.
%                'tiff' to save the sequence of images in tiff-format.
%
%
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   04/19/2014        0.5             hme
%

narginchk(1,2);

fm='';
format=''; 

for i=1:nargin
   
    if isstruct(varargin{i})
       fm=varargin{i};
    elseif strcmp(varargin{i},'jpg')
       format='jpg';
    elseif strcmp(varargin{i},'eps')
       format='eps';
    elseif strcmp(varargin{i},'png')
       format='png';
    elseif strcmp(varargin{i},'tiff')
       format='tiff';
    elseif strcmp(varargin{i},'pdf')
       format='pdf';
    else
       error(['Input argument at position ' int2str(i) ' not recognized']);
    end

end

if isempty(fm)
  error(['A frame or a sequence of frames must be given!']);;
end
if isempty(format)
  format='png';
end

if strcmp(format,'pdf')
   msg01=['You are going to produce ', num2str(size(fm,2)),' PDF images!  This can take some while to complete!'];
   warning('Mov:Images',msg01);
   msg02='Waiting 15 seconds to interrupt...';
   warning('Mov:pause',msg02);
   pause(15);
elseif strcmp(format,'eps')
   msg01=['You are going to produce ', num2str(size(fm,2)),' EPS images! This can take some while to complete!'];
   warning('Mov:Images',msg01);
   msg02='Waiting 15 seconds to interrupt...';
   warning('Mov:pause',msg02);
   pause(15);
end

if strcmp(format,'jpg')
    for k=1:size(fm,2),
          fname=sprintf('core3Dim_%.04d.jpg',k); 
          imwrite(fm(k).cdata,fname,'jpg');
    end
elseif strcmp(format,'png') 
    for k=1:size(fm,2), 
           fname=sprintf('core3Dim_%.04d.png',k); 
           imwrite(fm(k).cdata,fname,'png');
    end
elseif strcmp(format,'tiff')
    for k=1:size(fm,2), 
           fname=sprintf('core3Dim_%.04d.tiff',k); 
           imwrite(fm(k).cdata,fname,'tiff');
    end
elseif strcmp(format,'eps')
    for k=1:size(fm,2),
       [im,map] = frame2im(fm(k));
       image(im);
       axis off;
       Fig=myaa('publish');
       fname=sprintf('core3Dim_%.04d.eps',k);
       saveas(Fig,fname, 'psc2')
    end
elseif strcmp(format,'pdf')
    for k=1:size(fm,2),
       [im,map] = frame2im(fm(k));
       image(im);
       axis off;
       Fig=myaa('publish');
       fname=sprintf('core3Dim_%.04d.pdf',k);
       saveas(Fig,fname,'pdf')
    end
else
    for k=1:size(fm,2),
       fname=sprintf('core3Dim_%.04d.png',k);
       saveas(fm(k),fname,'png')
    end
end
