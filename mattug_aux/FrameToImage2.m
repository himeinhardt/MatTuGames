function FrameToImage2(fm)
%
%
%
for k=1:size(fm,2),
    [im,map] = frame2im(fm(k));
    image(im);
    axis off;
    Fig=myaa('publish');
    fname=sprintf('strcore3Dim_%.04d.eps',k);
    saveas(Fig,fname,'eps')
end
