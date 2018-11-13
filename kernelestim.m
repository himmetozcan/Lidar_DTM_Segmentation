function kernelMap=kernelestim(inputIm,voteLocations,kernelSigma)

[n, m] = size(inputIm);

kernelMap = zeros(size(inputIm));

Ni = size(voteLocations, 1);

for i = 1 : Ni
    try
    sigmaIn = kernelSigma;
    win = round(8 * sigmaIn);
    gaussWin = fspecial('gaussian', win + 1, sigmaIn);
    gaussWin = gaussWin / max(gaussWin(:));
    
    xstr = max(voteLocations(i,1) - win/2, 1);
    xend = min(voteLocations(i,1) + win/2, m);
    ystr = max(voteLocations(i,2) - win/2, 1);
    yend = min(voteLocations(i,2) + win/2, n);
    
    if xstr > 1 && xend < m
        x1 = 1;
        x2 = win + 1;
    elseif xstr == 1 && xend < m
        x1 = win + 1 - xend + 1;
        x2 = win + 1;
    elseif xstr > 1 && xend == m
        x1 = 1;
        x2 = m - xstr + 1;
    end
    
    if ystr > 1 && yend < n
        y1 = 1;
        y2 = win + 1;
    elseif ystr == 1 && yend < n
        y1 = win + 1 - yend + 1;
        y2 = win + 1;
    elseif ystr > 1 && yend == n
        y1 = 1;
        y2 = n - ystr + 1;
    end

    gaussWinSub = gaussWin(y1 : y2, x1 : x2);
    kernelMap(ystr : yend, xstr : xend) = kernelMap(ystr : yend, xstr : xend) + gaussWinSub;  
    catch
              disp('Kernel est.,  sth went wrong!')  
    end
end






