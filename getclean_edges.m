function [edges_clean, edges] = getclean_edges (inputIm, cannySigma, removeEdgesVarianceThr)

edges = edge(inputIm, 'Canny', [], cannySigma);

C = bwconncomp(edges);
S = regionprops(C,'Image','BoundingBox');

edges_clean = zeros(size(inputIm));
numberOfEdgeSegments = size(S,1);
counter = 0;

wi = 1;
[ns, ms] = size(inputIm);
for i = 1 : numberOfEdgeSegments
    try
    XY = ceil(S(i).BoundingBox);
    xstr = XY(1); xstp = xstr + XY(3);
    ystr = XY(2); ystp = ystr + XY(4);
    
    xstr0 = xstr; ystr0 = ystr;
    xstp0 = xstp; ystp0 = ystp;
    
    xstr = max(xstr-wi, 1); xstp = min(xstp+wi, ms);
    ystr = max(ystr-wi, 1); ystp = min(ystp+wi, ns);
    
    maskDsm = inputIm(ystr:ystp, xstr:xstp);    
   
    dY = ystp0 - ystr0;
    dX = xstp0 - xstr0;
    
    if xstr0 > wi; xstr1 = wi + 1; else xstr1 = xstr0; end
    if ystr0 > wi; ystr1 = wi + 1; else ystr1 = ystr0; end
    
    xstp1 = xstr1 + dX - 1;
    ystp1 = ystr1 + dY - 1;
    
    edgeSegment = zeros(ystp-ystr+1, xstp-xstr+1);
    edgeSegment(ystr1:ystp1, xstr1:xstp1) = S(i).Image;
    
    thickenEdgeSegment = bwmorph(edgeSegment, 'bridge');
    thickenEdgeSegment = bwmorph(thickenEdgeSegment, 'diag');
    thickenEdgeSegment = bwmorph(thickenEdgeSegment, 'thicken',1);
        
    if var( maskDsm(thickenEdgeSegment) ) > removeEdgesVarianceThr
        edges_clean(ystr:ystp, xstr:xstp) = edges_clean(ystr:ystp, xstr:xstp) | edgeSegment;
        counter = counter + 1;
    end
    catch
       disp('Remove edges, sth went wrong!') 
    end
end







