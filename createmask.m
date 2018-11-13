function [mask,mask2,xstr,xstp,ystr,ystp]=createmask(dsm,xm,ym,wi)
warning('off')
[ni,mi]=size(dsm);
% wi=10;
dsm2=dsm;
for i=1:length(xm)
    %     disp(i)
    dsm(ym(i),xm(i))=1e5;
    xstr=xm(i)-wi; xstr=max(1,xstr);
    xstp=xm(i)+wi; xstp=min(mi,xstp);
    ystr=ym(i)-wi; ystr=max(1,ystr);
    ystp=ym(i)+wi; ystp=min(ni,ystp);
    mask=dsm(ystr:ystp,xstr:xstp);
    mask2=dsm2(ystr:ystp,xstr:xstp);
    
end
% hold on; plot(xm2, ym2,'g*','MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor','r');
