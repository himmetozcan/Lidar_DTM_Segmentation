function buis=segmentation_process(dsm,dsm2,xm,ym,segWinSize,hDif,Ni,...
    thickStep, tUpper,tLower,muFinal, compareOuters)

buis = zeros(size(dsm));
numberOfSegments = numel(xm);
buiSegments=cell(numberOfSegments,1);

parfor i = 1 : numberOfSegments
    try
        
        [maskDsm, maskDsm2, xstr, xstp, ystr, ystp] = ...
            createmask(dsm2, xm(i), ym(i), segWinSize);
        
        cHeight = dsm2(ym(i),xm(i));
        J = ( maskDsm2 > (cHeight-hDif) ) & ( maskDsm2 < (cHeight+hDif*1.5) );
        
        J2 = imopen(J, strel('disk', 1));
        
        [amax, bmax] = find(maskDsm == max(maskDsm(:)) );
        J3 = bwselect(J2, bmax, amax);
        
        if ~ any( J3(:) )
            continue;
        end
        
        counter = 1;
        
        %- calculate the mean elevation value
        muX = mean(maskDsm2(J3));
        while( counter < Ni + 1 )
            
            Jprev = J3;
            for kk = 1 : Ni
                
                %- Thicken the current segment
                J = bwmorph(Jprev, 'thicken', thickStep);
                
                %- Extract the boundary of thickened segment
                ind = find((J-Jprev) > 0);
                
                %- calculate the elevation diffrence for each boundary pixel
                hdifs = muX - maskDsm2(ind);
                
                %                 indPos = [indPos(:); ind(hdifs > 0)];
                %                 indNeg = [indNeg(:); ind(hdifs < 0)];
                
                %- If any boundary pixel’s elevation value is too low
                J( ind (hdifs > tUpper ) ) = 0;
                
                %- If any boundary pixel’s elevation value is too high
                J( ind (hdifs < tLower ) ) = 0;
                
                J = bwselect(J, bmax, amax);
                
                tmpDif=1;
                if kk>1
                    tmpDif=J-Jprev;
                end
                if ~ any( J(:) ) || ~any(tmpDif(:))
                    counter = Ni + 1;
                    continue;
                end
                
                Jprev = J;
                muX = mean(maskDsm2(Jprev));
                counter = counter + 1;
                
            end
        end
        
        indInner = bwperim( bwmorph (J, 'thin', 1 ) );
        indOuter = bwperim( bwmorph( J, 'thicken', compareOuters ) );
        
        muX=mean( maskDsm2(indInner) ); % mean elevation value of current segment
        muB=mean( maskDsm2(indOuter) ); % mean elevation value of outer boundary
        deltaMu = muX - muB; % elevation height difference of segment and outer boundary
        
        flagLabel=deltaMu > muFinal;
        
        if flagLabel
            tmp=zeros(size(dsm));
            tmp(ystr:ystp,xstr:xstp)=J;
            buiSegments{i}.seg = find(tmp);
        end
    catch
        disp('Segmentation, sth went wrong!')
    end
    
end

%% Combine all segmentation results
for i = 1 : numel(buiSegments)
    if ~ isempty( buiSegments{i} )
        Jlast = buiSegments{i}.seg;
        buis(Jlast)=1;
    end
end

