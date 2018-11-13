function locs2 = spatialvoting(inputIm, locMax, voteMaxWinSize, voteHeightDifference)

[n, m] = size(inputIm);

counter = 0;

locs2 = [];

numOfLocMaxs = size(locMax,1);

for i = 1 : numOfLocMaxs
    try
    x = locMax(i, 1); x = round(x); x = max(x, 1); x = min(x, m);
    y = locMax(i, 2); y = round(y); y = max(y, 1); y = min(y, n);
    
    strx2 = max(1, x - voteMaxWinSize); stpx2 = min(x + voteMaxWinSize, m);
    stry2 = max(1, y - voteMaxWinSize); stpy2 = min(y + voteMaxWinSize, n);
    
    dsm_sub = inputIm(stry2 : stpy2, strx2 : stpx2);
    
    maskKose = ones(size(dsm_sub));
    tmp_dsm = dsm_sub .* maskKose;
    
    maxHeight = max(tmp_dsm(:));
    heightDif = abs(min(dsm_sub(:)) - maxHeight);
    
    if heightDif >= voteHeightDifference
        
        [y2, x2] = find(dsm_sub == maxHeight);
        
        y2 = y2(1);
        x2 = x2(1);
        
        counter = counter + 1;
        
        locs2(counter, :) = [x2 + strx2 - 1, y2 + stry2 - 1];
        
    end
    catch
               disp('Spatial voting, sth went wrong!') 
    end
end




