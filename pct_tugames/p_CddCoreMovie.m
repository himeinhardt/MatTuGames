function p_CddCoreMovie(varargin)
% P_CDDCOREMOVIE creates a movie w.r.t. the strong epsilon-cores. 
% The Multi-Parametric Toolbox 3 is needed.
% http://control.ee.ethz.ch/~mpt/3/Main/HomePage
%
% Usage: p_CddCoreMovie(v,'all','3','2.sec')
% Define variables:
%  output:
%             -- A movie of the strong epsilon cores of game v. 
%
%  input:     At most six input arguments are admissible. Any order of
%             the input arguments below is allowed. Nevertheless, at
%             least a game v must be specified. 
%
%
%  v          -- A Tu-Game v of length 2^n-1.
%  crit_val   -- A number as a string to specify the largest epsilon core.
%                Use the functions:
%                critical_value1,critical_value2, or critical_value_star
%                to find an useful number.
%                For instance:
%                crv=critical_value1(v)
%
%                 crv =
%                       8
%                and invoke
%                CddCoreMovie(v,'8')                
% 
%  add_sol    -- A string to invoke additional solutions into the plot.
%                Permissible solutions are:
%                'none', this is the default value.
%                'prk', a pre-kernel element will be incorporated.
%                'prn', the pre-nucleolus will be incorporated.
%                'shap', the Shapley value will be incorporated.
%                'all', all three solutions above will be incorporated.
%  zoom       -- A string argument to zoom into the movie. Permissible 
%                arguments are: 
%                'zoom'
%                '' 
%                The empty set is default, which means no zooming.
%  format     -- A file format to save images. 
%                Permissible formats are:
%                'avi' to save the sequence of pictures as a movie.
%                'jpg' to save the sequence of images in jpg-format. 
%                 Useful to produce your own QuickTime movie. 
%                'eps' to save the sequence of images in eps-format.
%  dur        -- Specify the duration of the movie by a string argument. A number
%                as a string and the keyword 'sec' must be supplied. The default value
%                is one second, that means 25 images will be produced. 
%                The number of seconds determines the multiple of it.
%                Permissible strings are for instance:
%                '3.sec'
%                '3 sec'
%                '3sec'
%                to produce 75 images.  
%  tol        -- A positive tolerance value. Its default value is set to 10^9*eps.
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

narginchk(1,7);

v='';
add_sol='';
tol='';
format=''; % saves bad avi-file!
crit_val='';
dur=1;  % length of the movie, default is one sec, that is, 25 images.
zm='';

for i=1:nargin

  if ischar(varargin{i})
    if strcmp(varargin{i},'none')
       add_sol='none';
    elseif strcmp(varargin{i},'prk')
       add_sol='prk';
    elseif strcmp(varargin{i},'prn')
       add_sol='prn';
    elseif strcmp(varargin{i},'shap')
       add_sol='shap';
    elseif strcmp(varargin{i},'all')
       add_sol='all';
    elseif strcmp(varargin{i},'avi')
       format='avi';
    elseif strcmp(varargin{i},'jpg')
       format='jpg';
    elseif strcmp(varargin{i},'eps')
       format='eps';
    elseif strcmp(varargin{i},'png')
       format='png';
    elseif strcmp(varargin{i},'tiff')
       format='svg';
    elseif strcmp(varargin{i},'pdf')
       format='pdf';
    elseif strcmp(varargin{i},'zoom')
       zm='zoom';
    elseif isempty(regexp(varargin{i},'sec'))==0
       str=varargin{i};
       p=regexp(varargin{i},'sec')-1;
       dur=str2num(str(1:p));
    else
       crit_val=str2num(varargin{i});
       if isempty(crit_val)==0
          crit_val=str2num(varargin{i});
       end
    end
  else
    if varargin{i}>0 & varargin{i}<1
      tol=varargin{i};
    elseif length(varargin{i})==7
      v=varargin{i};
    elseif length(varargin{i})==15
      v=varargin{i};
    else
     N=15;
     if varargin{i}>1
      if length(varargin{i})<N
        error('Game has not the correct size!');
      elseif length(varargin{i})>N
        error('Game has dimension larger than four!');
      else
        error(['Input argument at position ' int2str(i) ' not recognized']);
      end
     else
       error(['Input argument at position ' int2str(i) ' not recognized']);
     end
    end
  end
end


if isempty(v)
   error('At least the game must be given!');
end
if isempty(add_sol)
  add_sol='none';
end
if isempty(tol)
  tol=10^9*eps;
end
if isempty(format)
  format='jpg';
end

N=length(v);
[~, n]=log2(N);

clf;
fmin=CddLeastCore(v);
if isempty(crit_val)
   ctv1=critical_value1(v);
   ctv2=critical_value2(v);
   ctv3=critical_value_star(v);
   vc=[ctv1,ctv2,ctv3];
   ctv=max(vc);
   if ctv<=0
      ctv=2;
   end
else
   ctv=crit_val;
end
div=min(1000,dur*25);
if div>=1000
   msg01='You are going to produce more than 1000 images! Interrupt if you do not want to continue!';
   warning('Mov:Images',msg01);
   msg02='Waiting 5 seconds to interrupt...';
   warning('Mov:pause',msg02);
   pause(5);
end
sz=abs(fmin-ctv)/div;
t=ctv:-sz:fmin;
t(end+1)=fmin;
[crv_vert,~,~,Pv]=CddCoreVertices(v);
y=range(crv_vert);
[~, idx]=min(y);
[imp_vert,~,~,Pip]=CddImputationVertices(v,idx);
v_ctv=streps_value(v,ctv);
[cr_ctv,crst_ctv,vol_ctv,Pctv]=CddCoreVertices(v_ctv,idx);
v_prk=PreKernel(v);
v_prn=CddPrenucl(v);
v_sh=ShapleyValue(v);
ms1=min(Pip.V);
ml1=max(Pip.V);
if n==4
   v_prk(:,idx)=[];
   v_prn(:,idx)=[];
   v_sh(:,idx)=[];
else
   [X1,X2]=ToSimplex(v_prk);
   v_prk=[X1,X2];
   [X1,X2]=ToSimplex(v_prn);
   v_prn=[X1,X2];
   [X1,X2]=ToSimplex(v_sh);
   v_sh=[X1,X2];
end
if n==4
%   for k=1:length(t)
  parfor k=1:8
    v_eps=streps_value(v,t(k));
    [cr_eps,crst_eps,vol_eps,P]=CddCoreVertices(v_eps,idx);
    if strcmp(zm,'zoom')
       ms2=min(P.V);
       ml2=max(P.V);
       sm=floor(min(ms1,ms2))-2;
       lr=ceil(max(ml1,ml2))+2;
    else
       sm=floor(min(Pctv.V))-2;
       lr=ceil(max(Pctv.V))+2;
    end
%    figure(k);
%    figure('Renderer','opengl')
%    set(gcf,'Renderer','opengl');
%    get(gcf)
    hp=plot3(v_prk(1),v_prk(2),v_prk(3));
    grid on;
    set(hp,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
    hold on
    hn=plot3(v_prn(1),v_prn(2),v_prn(3));
    hs=plot3(v_sh(1),v_sh(2),v_sh(3));
    set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
    set(hs,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
    h=P.plot('linewidth', 1.2);
    set(h,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0.7 0 0]);
    Pv.plot('alpha',0);
    Pip.plot('alpha',0,'linewidth', 0.8);
    title(['Strong ',num2str(t(k)),'-Core']);
    view(120, 25);
    if strcmp(zm,'zoom')
%       set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
       axis tight;
       axis off;
    else
       set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
    end
%    get(gcf)
    hold off;
    camlight('headlight');
    lighting phong;
    material shiny;
    MCrMov(k)=getframe(gcf(k));
%    close(k);
  end
else
    for k=1:length(t)
    v_eps(k,:)=streps_value(v,t(k));
    [cr_eps{k},crst_eps{k},vol_eps(k),P(k)]=CddCoreVertices(v_eps(k,:));
    if strcmp(zm,'zoom')
       sm=floor(min(Pip.V))-ctv;
       lr=ceil(max(Pip.V))+ctv;
    else
       sm=floor(min(P(1).V))-0.5;
       lr=ceil(max(P(1).V))+0.5;
    end
    hp=plot(v_prk(1),v_prk(2));
    grid on;
    set(hp,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
    hold on
    hn=plot(v_prn(1),v_prn(2));
    hs=plot(v_sh(1),v_sh(2));
    set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
    set(hs,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
    P(k).plot('alpha',0.5);
    Pv.plot('alpha',0);
    Pip.plot('alpha',0);
    title(['Strong ',num2str(t(k)),'-Core']);
    if strcmp(zm,'zoom')
%       set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
       axis tight;
       axis off;
    else
       set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
    end
    hold off;
    if strcmp(format,'jpg');
       fname=sprintf('core2Dim_%.04d.jpg',k);
       print('-djpeg100',fname);
    elseif strcmp(format,'eps')
       fname=sprintf('core2Dim_%.04d.eps',k);
       print('-depsc2','-r0',fname);
    elseif strcmp(format,'png')
       fname=sprintf('core3Dim_%.04d.png',k);
       print('-dpng','-r0',fname);
    elseif strcmp(format,'svg')
       fname=sprintf('core3Dim_%.04d.svg',k);
       print('-dsvg','-r0',fname);
    elseif strcmp(format,'pdf')
       fname=sprintf('core3Dim_%.04d.pdf',k);
       print('-dpdf','-r0',fname);
    end
  end 
end

if strcmp(format,'jpg')
    parfor k=1:size(MCrMov,2),
          fname=sprintf('core3Dim_%.04d.jpg',k); 
          imwrite(MCrMov(k).cdata,fname,'jpg');
    end
elseif strcmp(format,'png') 
    parfor k=1:size(MCrMov,2), 
           fname=sprintf('core3Dim_%.04d.png',k); 
           imwrite(MCrMov(k).cdata,fname,'png');
    end
elseif strcmp(format,'tiff')
    parfor k=1:size(MCrMov,2), 
           fname=sprintf('core3Dim_%.04d.tiff',k); 
           imwrite(MCrMov(k).cdata,fname,'tiff');
    end
end
