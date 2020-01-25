function CddCoreMovie(varargin)
% CDDCOREMOVIE creates a movie w.r.t. the strong epsilon-cores. 
% The Multi-Parametric Toolbox 3 is needed.
% http://control.ee.ethz.ch/~mpt/3/Main/HomePage
%
% Usage:       CddCoreMovie(v) 
%                    which is the default  usage. 
%              CddCoreMovie(v,'all','3','2.sec','zoom','jpg')
%                    which is a full usage.
% Define variables:
%  output:
%             -- A movie of the strong epsilon cores of game v. 
%
%  input:     At most eight input arguments are admissible. Any order of
%             the input arguments below is allowed. Nevertheless, at
%             least a game v must be supplied. 
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
%                'none', no solution below will be plotted.
%                'prk', a pre-kernel element will be incorporated.
%                'prn', the pre-nucleolus will be incorporated.
%                'shap', the Shapley value will be incorporated.
%                'all', all three solutions above will be incorporated.
%                       This is the default value.
%  vw_pt       -- A string command to determine the view point.
%                The default view point is
%                [120, 25]
%
%                To specify a different view point type, for instance
%                vw_pt='view(130,35)'
%                
%                Then invoke
%
%                CddCoreMovie(v,'all',vw_pt)
%
%  zoom       -- A string argument to zoom into the movie. Permissible 
%                arguments are: 
%                'zoom' 
%                '' 
%                The empty set is default, which means no zooming. 
%  format     -- A file format to save images. 
%                Permissible formats are:
%                'avi' to save the sequence of pictures as a movie.
%                'jpg' to save the sequence of images in jpg-format. 
%                   Useful to produce your own QuickTime movie. 
%                'eps' to save the sequence of images in eps-format.
%                'png' to save the sequence of images in png-format.
%                   Useful to produce your own movie with Blender.
%                'pdf' to save the sequence of images in pdf-format.
%                'tiff' to save the sequence of images in tiff-format.
%                'mat' to save the frames of the pictures.
%                   To replay the movie load the mat-file by
%
%                   load CoreMovie.mat
%
%                   and type in the command line
%                                       
%                   movie(gcf,MCrMov,10);  to play the movie ten times, for instance.
% 
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

narginchk(1,8);

v='';
add_sol='';
tol='';
format=''; % saves bad avi-file!
crit_val='';
dur=1;  % length of the movie, default is one sec, that is, 25 images.
zm='';
vw_pt='';

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
       format='tiff';
    elseif strcmp(varargin{i},'pdf')
       format='pdf';
    elseif strcmp(varargin{i},'mat')
       format='mat';
    elseif strcmp(varargin{i},'zoom')
       zm='zoom';
    elseif isempty(regexp(varargin{i},'sec'))==0
       str=varargin{i};
       p=regexp(varargin{i},'sec')-1;
       dur=str2num(str(1:p));
    elseif isempty(regexp(varargin{i},'view'))==0
       vw_pt=varargin{i};
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
  add_sol='all';
end
if isempty(tol)
  tol=10^9*eps;
end
if isempty(format)
  format='play';
end
if isempty(vw_pt)
  vw_pt='view(120,25)';
end

N=length(v);
[~, n]=log2(N);

clf;
fmin=CddLeastCore(v);
crQ=CddCoreQ(v);
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
vfmax=streps_value(v,ctv);
[crv_vfm,~,~,~]=CddCoreVertices(vfmax);
y=range(crv_vfm);
[~, idx]=min(y);
if crQ==1
  [crv_vert,~,~,Pv]=CddCoreVertices(v,idx);
else
  vfmin=streps_value(v,fmin);
  [crv_vert,~,~,Pv]=CddCoreVertices(vfmin,idx);
end
[imp_vert,~,~,Pip]=CddImputationVertices(v,idx);
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
winsize1=zeros(1,n);
if n==4
  if strcmp(add_sol,'all')
%   for k=1:8
    for k=1:length(t)
    v_eps(k,:)=streps_value(v,t(k));
    [cr_eps{k},crst_eps{k},vol_eps(k),P(k)]=CddCoreVertices(v_eps(k,:),idx);
    ms2=min(P(1).V);
    ml2=max(P(1).V);
    mrg=min(max(ml2)/4,2);
    sm=floor(min(P(1).V))-mrg;
    lr=ceil(max(P(1).V))+mrg;
    hp=plot3(v_prk(1),v_prk(2),v_prk(3));
    grid on;
    set(hp,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
    hold on
    hn=plot3(v_prn(1),v_prn(2),v_prn(3));
    hs=plot3(v_sh(1),v_sh(2),v_sh(3));
    set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
    set(hs,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
    h=P(k).plot('linewidth', 1.2);
    set(h,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0.7 0 0]);
    Pv.plot('alpha',0);
    Pip.plot('alpha',0,'linewidth', 0.8);
    title(['Strong ',num2str(t(k)),'-Core']);
    if strcmp(zm,'zoom')
       axis tight;
       axis off;
    else
       set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
    end
    eval(vw_pt); % setting view point.
    hold off;
    camlight('headlight');
    lighting phong;
    material shiny;
    winsize = get(gcf,'Position');
    winsize1(3:4)=winsize(3:4);
    MCrMov(k)= getframe(gcf,winsize1);
%    s = size(MCrMov(k).cdata);
%    fprintf('%d %d\n', s(2), s(1));
    if strcmp(format,'eps')
       fname=sprintf('core3Dim_%.04d.eps',k);
       print('-depsc2','-r0',fname);
    elseif strcmp(format,'pdf')
       fname=sprintf('core3Dim_%.04d.pdf',k);
       print('-dpdf','-r0',fname);
    end
   end
% Pre_Kernel
  elseif strcmp(add_sol,'prk')
   for k=1:length(t)
    v_eps(k,:)=streps_value(v,t(k));
    [cr_eps{k},crst_eps{k},vol_eps(k),P(k)]=CddCoreVertices(v_eps(k,:),idx);
    ms2=min(P(1).V);
    ml2=max(P(1).V);
    mrg=min(max(ml2)/4,2);
    sm=floor(min(P(1).V))-mrg;
    lr=ceil(max(P(1).V))+mrg;
    hp=plot3(v_prk(1),v_prk(2),v_prk(3));
    grid on;
    set(hp,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
    hold on
    h=P(k).plot('linewidth', 1.2);
    set(h,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0.7 0 0]);
    Pv.plot('alpha',0);
    Pip.plot('alpha',0,'linewidth', 0.8);
    title(['Strong ',num2str(t(k)),'-Core']);
    eval(vw_pt);
    if strcmp(zm,'zoom')
       axis tight;
       axis off;
    else
       set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
    end
    hold off;
    camlight('headlight');
    lighting phong;
    material shiny;
    winsize = get(gcf,'Position');
    winsize1(3:4)=winsize(3:4);
    MCrMov(k)= getframe(gcf,winsize1);
    if strcmp(format,'eps')
       fname=sprintf('core3Dim_%.04d.eps',k);
       print('-depsc2','-r0',fname);
    elseif strcmp(format,'pdf')
       fname=sprintf('core3Dim_%.04d.pdf',k);
       print('-dpdf','-r0',fname);
    end
   end
% Pre-Nucleolus
  elseif strcmp(add_sol,'prn')
   for k=1:length(t)
    v_eps(k,:)=streps_value(v,t(k));
    [cr_eps{k},crst_eps{k},vol_eps(k),P(k)]=CddCoreVertices(v_eps(k,:),idx);
    ms2=min(P(1).V);
    ml2=max(P(1).V);
    mrg=min(max(ml2)/4,2);
    sm=floor(min(P(1).V))-mrg;
    lr=ceil(max(P(1).V))+mrg;
    hn=plot3(v_prn(1),v_prn(2),v_prn(3));
    grid on;
    set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
    hold on
    h=P(k).plot('linewidth', 1.2);
    set(h,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0.7 0 0]);
    Pv.plot('alpha',0);
    Pip.plot('alpha',0,'linewidth', 0.8);
    title(['Strong ',num2str(t(k)),'-Core']);
    eval(vw_pt);
    if strcmp(zm,'zoom')
       axis tight;
       axis off;
    else
       set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
    end
    hold off;
    camlight('headlight');
    lighting phong;
    material shiny;
    winsize = get(gcf,'Position');
    winsize1(3:4)=winsize(3:4);
    MCrMov(k)= getframe(gcf,winsize1);
    if strcmp(format,'eps')
       fname=sprintf('core3Dim_%.04d.eps',k);
       print('-depsc2','-r0',fname);
    elseif strcmp(format,'pdf')
       fname=sprintf('core3Dim_%.04d.pdf',k);
       print('-dpdf','-r0',fname);
    end
   end
% Shapley value
  elseif strcmp(add_sol,'shap')
   for k=1:length(t)
    v_eps(k,:)=streps_value(v,t(k));
    [cr_eps{k},crst_eps{k},vol_eps(k),P(k)]=CddCoreVertices(v_eps(k,:),idx);
    ms2=min(P(1).V);
    ml2=max(P(1).V);
    mrg=min(max(ml2)/4,2);
    sm=floor(min(P(1).V))-mrg;
    lr=ceil(max(P(1).V))+mrg;
    hs=plot3(v_sh(1),v_sh(2),v_sh(3));
    grid on;
    set(hs,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
    hold on
    h=P(k).plot('linewidth', 1.2);
    set(h,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0.7 0 0]);
    Pv.plot('alpha',0);
    Pip.plot('alpha',0,'linewidth', 0.8);
    title(['Strong ',num2str(t(k)),'-Core']);
    eval(vw_pt);
    if strcmp(zm,'zoom')
       axis tight;
       axis off;
    else
       set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
    end
    hold off;
    camlight('headlight');
    lighting phong;
    material shiny;
    winsize = get(gcf,'Position');
    winsize1(3:4)=winsize(3:4);
    MCrMov(k)= getframe(gcf,winsize1);
    if strcmp(format,'eps')
       fname=sprintf('core3Dim_%.04d.eps',k);
       print('-depsc2','-r0',fname);
    elseif strcmp(format,'pdf')
       fname=sprintf('core3Dim_%.04d.pdf',k);
       print('-dpdf','-r0',fname);
    end
   end
% No solution
  elseif strcmp(add_sol,'none')
   for k=1:length(t)
    v_eps(k,:)=streps_value(v,t(k));
    [cr_eps{k},crst_eps{k},vol_eps(k),P(k)]=CddCoreVertices(v_eps(k,:),idx);
    ms2=min(P(1).V);
    ml2=max(P(1).V);
    mrg=min(max(ml2)/4,2);
    sm=floor(min(P(1).V))-mrg;
    lr=ceil(max(P(1).V))+mrg;
    h=P(k).plot('linewidth', 1.2);
    grid on;
    set(h,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0.7 0 0]);
    hold on
    Pv.plot('alpha',0);
    Pip.plot('alpha',0,'linewidth', 0.8);
    title(['Strong ',num2str(t(k)),'-Core']);
    eval(vw_pt);
    if strcmp(zm,'zoom')
       axis tight;
       axis off;
    else
       set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
    end
    hold off;
    camlight('headlight');
    lighting phong;
    material shiny;
    winsize = get(gcf,'Position');
    winsize1(3:4)=winsize(3:4);
    MCrMov(k)= getframe(gcf,winsize1);
    if strcmp(format,'eps')
       fname=sprintf('core3Dim_%.04d.eps',k);
       print('-depsc2','-r0',fname);
    elseif strcmp(format,'pdf')
       fname=sprintf('core3Dim_%.04d.pdf',k);
       print('-dpdf','-r0',fname);
    end
   end
% Default all solutions
  else
   for k=1:length(t)
    v_eps(k,:)=streps_value(v,t(k));
    [cr_eps{k},crst_eps{k},vol_eps(k),P(k)]=CddCoreVertices(v_eps(k,:),idx);
    ms2=min(P(1).V);
    ml2=max(P(1).V);
    mrg=min(max(ml2)/4,2);
    sm=floor(min(P(1).V))-mrg;
    lr=ceil(max(P(1).V))+mrg;
    hp=plot3(v_prk(1),v_prk(2),v_prk(3));
    grid on;
    set(hp,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
    hold on
    hn=plot3(v_prn(1),v_prn(2),v_prn(3));
    hs=plot3(v_sh(1),v_sh(2),v_sh(3));
    set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
    set(hs,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
    h=P(k).plot('linewidth', 1.2);
    set(h,'FaceLighting','phong','FaceAlpha',0.5,'FaceColor',[0.7 0 0]);
    Pv.plot('alpha',0);
    Pip.plot('alpha',0,'linewidth', 0.8);
    title(['Strong ',num2str(t(k)),'-Core']);
    eval(vw_pt);
    if strcmp(zm,'zoom')
       axis tight;
       axis off;
    else
       set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)],'ZLim',[sm(3) lr(3)]);
    end
    hold off;
    camlight('headlight');
    lighting phong;
    material shiny;
    winsize = get(gcf,'Position');
    winsize1(3:4)=winsize(3:4);
    MCrMov(k)= getframe(gcf,winsize1);
    if strcmp(format,'eps')
       fname=sprintf('core3Dim_%.04d.eps',k);
       print('-depsc2','-r0',fname);
    elseif strcmp(format,'pdf')
       fname=sprintf('core3Dim_%.04d.pdf',k);
       print('-dpdf','-r0',fname);
    end
   end  
  end
% 3-Persons 
% Plot All Solutions
else
  if strcmp(add_sol,'all') 
   for k=1:length(t)
    v_eps(k,:)=streps_value(v,t(k));
    [cr_eps{k},crst_eps{k},vol_eps(k),P(k)]=CddCoreVertices(v_eps(k,:));
    ms2=min(P(1).V);
    ml2=max(P(1).V);
    mrg=min(max(ml2)/4,2);
    sm=floor(min(P(1).V))-mrg;
    lr=ceil(max(P(1).V))+mrg;
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
       axis tight;
       axis off;
    else
       set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
    end
    hold off;
    winsize = get(gcf,'Position');
    winsize1(3:4)=winsize(3:4);
    MCrMov(k)= getframe(gcf,winsize1);
    if strcmp(format,'eps')
       fname=sprintf('core2Dim_%.04d.eps',k);
       print('-depsc2','-r0',fname);
    elseif strcmp(format,'pdf')
       fname=sprintf('core3Dim_%.04d.pdf',k);
       print('-dpdf','-r0',fname);
    end
   end
% Pre-Kernel
  elseif strcmp(add_sol,'prk')
   for k=1:length(t)
    v_eps(k,:)=streps_value(v,t(k));
    [cr_eps{k},crst_eps{k},vol_eps(k),P(k)]=CddCoreVertices(v_eps(k,:));
    ms2=min(P(1).V);
    ml2=max(P(1).V);
    mrg=min(max(ml2)/4,2);
    sm=floor(min(P(1).V))-mrg;
    lr=ceil(max(P(1).V))+mrg;
    hp=plot(v_prk(1),v_prk(2));
    grid on;
    set(hp,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
    hold on
    P(k).plot('alpha',0.5);
    Pv.plot('alpha',0);
    Pip.plot('alpha',0);
    title(['Strong ',num2str(t(k)),'-Core']);
    if strcmp(zm,'zoom')
       axis tight;
       axis off;
    else
       set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
    end
    hold off;
    winsize = get(gcf,'Position');
    winsize1(3:4)=winsize(3:4);
    MCrMov(k)= getframe(gcf,winsize1);
    if strcmp(format,'eps')
       fname=sprintf('core2Dim_%.04d.eps',k);
       print('-depsc2','-r0',fname);
    elseif strcmp(format,'pdf')
       fname=sprintf('core3Dim_%.04d.pdf',k);
       print('-dpdf','-r0',fname);
    end
   end
% Pre-Nucleolus
  elseif strcmp(add_sol,'prn')
   for k=1:length(t)
    v_eps(k,:)=streps_value(v,t(k));
    [cr_eps{k},crst_eps{k},vol_eps(k),P(k)]=CddCoreVertices(v_eps(k,:));
    ms2=min(P(1).V);
    ml2=max(P(1).V);
    mrg=min(max(ml2)/4,2);
    sm=floor(min(P(1).V))-mrg;
    lr=ceil(max(P(1).V))+mrg;
    hn=plot(v_prn(1),v_prn(2));
    grid on;
    set(hn,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
    hold on
    P(k).plot('alpha',0.5);
    Pv.plot('alpha',0);
    Pip.plot('alpha',0);
    title(['Strong ',num2str(t(k)),'-Core']);
    if strcmp(zm,'zoom')
       axis tight;
       axis off;
    else
       set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
    end
    hold off;
    winsize = get(gcf,'Position');
    winsize1(3:4)=winsize(3:4);
    MCrMov(k)= getframe(gcf,winsize1);
    if strcmp(format,'eps')
       fname=sprintf('core2Dim_%.04d.eps',k);
       print('-depsc2','-r0',fname);
    elseif strcmp(format,'pdf')
       fname=sprintf('core3Dim_%.04d.pdf',k);
       print('-dpdf','-r0',fname);
    end
   end
% Shapley Value
  elseif strcmp(add_sol,'shap')
   for k=1:length(t)
    v_eps(k,:)=streps_value(v,t(k));
    [cr_eps{k},crst_eps{k},vol_eps(k),P(k)]=CddCoreVertices(v_eps(k,:));
    ms2=min(P(1).V);
    ml2=max(P(1).V);
    mrg=min(max(ml2)/4,2);
    sm=floor(min(P(1).V))-mrg;
    lr=ceil(max(P(1).V))+mrg;
    hs=plot(v_sh(1),v_sh(2));
    grid on;
    set(hs,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
    hold on
    P(k).plot('alpha',0.5);
    Pv.plot('alpha',0);
    Pip.plot('alpha',0);
    title(['Strong ',num2str(t(k)),'-Core']);
    if strcmp(zm,'zoom')
       axis tight;
       axis off;
    else
       set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
    end
    hold off;
    winsize = get(gcf,'Position');
    winsize1(3:4)=winsize(3:4);
    MCrMov(k)= getframe(gcf,winsize1);
    if strcmp(format,'eps')
       fname=sprintf('core2Dim_%.04d.eps',k);
       print('-depsc2','-r0',fname);
    elseif strcmp(format,'pdf')
       fname=sprintf('core3Dim_%.04d.pdf',k);
       print('-dpdf','-r0',fname);
    end
   end
% No solution
  elseif strcmp(add_sol,'none')
   for k=1:length(t)
    v_eps(k,:)=streps_value(v,t(k));
    [cr_eps{k},crst_eps{k},vol_eps(k),P(k)]=CddCoreVertices(v_eps(k,:));
    ms2=min(P(1).V);
    ml2=max(P(1).V);
    mrg=min(max(ml2)/4,2);
    sm=floor(min(P(1).V))-mrg;
    lr=ceil(max(P(1).V))+mrg;
    P(k).plot('alpha',0.5);
    grid on;
    if strcmp(zm,'zoom')
       axis tight;
       axis off;
    else
       set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
    end
    hold on
    Pv.plot('alpha',0);
    Pip.plot('alpha',0);
    title(['Strong ',num2str(t(k)),'-Core']);
    hold off;
    winsize = get(gcf,'Position');
    winsize1(3:4)=winsize(3:4);
    MCrMov(k)= getframe(gcf,winsize1);
    if strcmp(format,'eps')
       fname=sprintf('core2Dim_%.04d.eps',k);
       print('-depsc2','-r0',fname);
    elseif strcmp(format,'pdf')
       fname=sprintf('core3Dim_%.04d.pdf',k);
       print('-dpdf','-r0',fname);
    end
   end
% Default, all solutions
  else
   for k=1:length(t)
    v_eps(k,:)=streps_value(v,t(k));
    [cr_eps{k},crst_eps{k},vol_eps(k),P(k)]=CddCoreVertices(v_eps(k,:));
    ms2=min(P(1).V);
    ml2=max(P(1).V);
    mrg=min(max(ml2)/4,2);
    sm=floor(min(P(1).V))-mrg;
    lr=ceil(max(P(1).V))+mrg;
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
       axis tight;
       axis off;
    else
       set(gca,'XLim',[sm(1) lr(1)],'YLim',[sm(2) lr(2)]);
    end
    hold off;
    winsize = get(gcf,'Position');
    winsize1(3:4)=winsize(3:4);
    MCrMov(k)= getframe(gcf,winsize1);
    if strcmp(format,'eps')
       fname=sprintf('core2Dim_%.04d.eps',k);
       print('-depsc2','-r0',fname);
    elseif strcmp(format,'pdf')
       fname=sprintf('core3Dim_%.04d.pdf',k);
       print('-dpdf','-r0',fname);
    end
   end
  end
end
% Create AVI file.
if strcmp(format,'avi')
   movie2avi(MCrMov, 'CoreMovie.avi', 'compression', 'None');
elseif strcmp(format,'jpg')
    for k=1:size(MCrMov,2),
          fname=sprintf('core3Dim_%.04d.jpg',k); 
          imwrite(MCrMov(k).cdata,fname,'jpg');
    end
    save CoreMovie02.mat MCrMov v;
elseif strcmp(format,'png') 
    for k=1:size(MCrMov,2), 
           fname=sprintf('core3Dim_%.04d.png',k); 
           imwrite(MCrMov(k).cdata,fname,'png');
    end
    save CoreMovie02.mat MCrMov v;
elseif strcmp(format,'tiff')
    for k=1:size(MCrMov,2), 
           fname=sprintf('core3Dim_%.04d.tiff',k); 
           imwrite(MCrMov(k).cdata,fname,'tiff');
    end
    save CoreMovie02.mat MCrMov v;
elseif strcmp(format,'mat')
   save CoreMovie.mat MCrMov v;
end
% Play the movie three times
% movie(MCrMov,3);
movie(gcf,MCrMov,1,25,winsize1);
