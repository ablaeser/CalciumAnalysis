function fiber = MakeFibers(expt, fiber, ROI, varargin)
Nfiber = numel(fiber);
if nargin > 3
    show = varargin{1};
else
    show = false;
end
if show, figure('Units','normalized', 'OuterPosition', [0,0,1,1], 'Color','w', 'PaperOrientation','landscape'); end
for f = 1:Nfiber
    fiber(f).Nroi = numel( fiber(f).ROI );
    fiber(f).Nvox = sum([ROI(fiber(f).ROI).Nvox]);
    fiber(f).ind = vertcat( ROI(fiber(f).ROI).ind );
    fiber(f).neuropil = vertcat( ROI(fiber(f).ROI).neuropil );
    tempCrop = vertcat( ROI(fiber(f).ROI).crop );
    fiber(f).crop = [min(tempCrop(:,1)), max(tempCrop(:,2)), min(tempCrop(:,3)), max(tempCrop(:,4))];     
    fiber(f).width = fiber(f).crop(2) - fiber(f).crop(1); 
    fiber(f).height = fiber(f).crop(4) - fiber(f).crop(3);
    fiber(f).zUse = unique( horzcat(ROI(fiber(f).ROI).zUse) );
    fiber(f).cent = vertcat(ROI(fiber(f).ROI).cent);
    fiber(f).apSep = expt.umPerPixel*squareform(pdist(fiber(f).cent(:,1)));
    fiber(f).mlSep = expt.umPerPixel*squareform(pdist(fiber(f).cent(:,2)));
    fiber(f).footSep = expt.umPerPixel*squareform(pdist(fiber(f).cent(:,1:2)));
    fiber(f).length = sum( diag(fiber(f).footSep,1) );
    if expt.Nplane > 1
        fiber(f).volSep = expt.umPerPixel*squareform(pdist(fiber(f).cent(:,1:3)));
    else
        fiber(f).volSep = [];
    end
    %fiber(f).similarity = stillCosSim(fiber(f).ROI, fiber(f).ROI); %nan(fiber(f).Nroi); %
    %fiber(f).correlation = stillCorr(fiber(f).ROI, fiber(f).ROI); %nan(fiber(f).Nroi); %
    fiber(f).labelFoot = zeros(expt.Nrow, expt.Ncol);
    fiber(f).labelVol = zeros(expt.Nrow, expt.Ncol, expt.Nplane);
    for r = fiber(f).ROI
        %axonLabelFoot(ROI(r).footprintInd) = a;
        %axonLabelVol(ROI(r).ind) = a;
        fiber(f).labelFoot(ROI(r).footprintInd) = r;
        fiber(f).labelVol(ROI(r).ind) = r;
    end
    if fiber(f).Nroi > 1
        tempHull = bwconvhull(fiber(f).labelFoot );
        if show
            subplot(2,1,1); imshow(fiber(f).labelFoot); title(sprintf('fiber %i', f)); % sprintf('[x,f] = [%i, %i]', x,f)
            subplot(2,1,2); imshow(tempHull ); title('Convex Hull'); %pause(0.1);
            impixelinfo;
            pause
        end
        tempHullProps = regionprops(tempHull, 'Orientation', 'Eccentricity', 'MajoraxisLength', 'MinorAxisLength'); % 
        fiber(f).orientation = deg2rad( tempHullProps.Orientation );
        fiber(f).eccentricity = tempHullProps.Eccentricity;
        fiber(f).MajorAxisLength = tempHullProps.MajorAxisLength;
        fiber(f).MinorAxisLength = tempHullProps.MinorAxisLength;
        %fiber(f).PrinAxLength = [tempHullProps.MajorAxisLength, tempHullProps.MinorAxisLength];
    else
        fiber(f).orientation = deg2rad( ROI(fiber(f).ROI).orientation ); %median( [] );
        fiber(f).eccentricity = ROI(fiber(f).ROI).eccentricity;
        fiber(f).MajorAxisLength = ROI(fiber(f).ROI).PrinAxLength(1);
        fiber(f).MinorAxisLength = ROI(fiber(f).ROI).PrinAxLength(2);
        %fiber(f).PrinAxLength = ROI(fiber(f).ROI).PrinAxLength;%(1:2)
    end
end
%fluor = GetAxonFluor(expt, catInfo, fiber, fluor, deform, 'window',find(Tscan{1}<=32,1,'last'), 'overwrite',true);
end