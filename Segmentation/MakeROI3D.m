function [ROI, preROI] = MakeROI3D(expt, varargin) % , run exptDir, exptName %ROIlabel, 
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addParameter( IP, 'corrPrct', 90, @isnumeric)  % 68
addParameter( IP, 'minVol', 50, @isnumeric) 
addParameter( IP, 'minFoot', 50, @isnumeric) 
addParameter( IP, 'neuropilRad', [8,21], @isnumeric ) % inner and outer radii (in xy plane) to draw neuropil ROI{z}
addParameter( IP, 'overlap', 0.5, @isnumeric )
addParameter( IP, 'corr', false, @islogical)
addParameter( IP, 'show', false, @islogical)
addParameter( IP, 'overwrite', false, @islogical)
parse( IP, expt, varargin{:} );  %  exptDir, exptName fluor, deform, T, 
params.neuropilRad = IP.Results.neuropilRad;
params.minCorrPrct = IP.Results.corrPrct; 
params.minVol = IP.Results.minVol; 
params.minFoot = IP.Results.minFoot; 
params.overlap = IP.Results.overlap;
removeCorr = ~IP.Results.corr;
overwrite = IP.Results.overwrite;
show = IP.Results.show;
roiDir = strcat(expt.dir, 'ROI\'); mkdir(roiDir);
savePath = sprintf('%s%s_ROI.mat', expt.dir, expt.name ); 
opt = {[0.04, 0.03], [0.06, 0.05], [0.03,0.05]};  % {[vert, horz], [bottom, top], [left, right] }
tic;
if ~exist(savePath,'file') || overwrite
    nanVol = nan( expt.Nrow, expt.Ncol, expt.Nplane );  falseVol = false( expt.Nrow, expt.Ncol, expt.Nplane );
    % Get and process fluorescence signals
    [~,seg3Dpath] = FileFinder(expt.dir, 'contains','final');
    %seg3Dpath = sprintf('%s%s_seg_data_final.mat', expt.dir, expt.name);  % tracePath = sprintf('%s%s_seg_traces.mat', expt.dir, expt.name);   
    Nz = numel(seg3Dpath); % # of projections used for segmentation
    ROI = cell(1,Nz); preROI = cell(1,Nz); Nroi = nan(1,Nz);
    for z = 1:Nz
        fprintf('\nLoading %s... ', seg3Dpath{z}); 
        load(seg3Dpath{z}, 'seg_data'); 
        preROI{z} = seg_data([seg_data.keep]); %seg_data_keep;
        if expt.Nplane == 1       
            % {
            % Get intermediate 2D masks by breaking preROI{z} (PCA/ICA results) up into blobs and thresholding on max correlation
            fprintf('\nBreaking up preliminary ROIs')
            c = 0; intROI = struct();
            for p = 1:numel(preROI{z})
                tempConnComp = bwconncomp( preROI{z}(p).mask_2d, 8 );
                tempSubProps = regionprops( tempConnComp, 'Area','PixelIdxList'); % , preROI{z}(p).correlation
                goodSub = find([tempSubProps.Area] > params.minFoot );
                preROI{z}(p).Nsub = numel(goodSub);
                tempBinary = preROI{z}(p).mask_2d; %false( expt.Nrow, expt.Ncol ); tempBinary(ROI{z}(r).footprintInd) = true; 
                tempBinary = imdilate(tempBinary, strel('rectangle',[3,3]) ); % imshow( tempBinary )
                tempEdge = edge(tempBinary);
                [preROI{z}(p).footprintEdge(:,1), preROI{z}(p).footprintEdge(:,2)] = ind2sub( [expt.Nrow, expt.Ncol], find( tempEdge ) );
                for s = goodSub
                    c = c + 1;
                    % Retain some info from preliminary ROI{z}
                    intROI(c).mask_2d = falseVol; %preROI{z}(p).mask_2d;
                    intROI(c).mask_2d(tempSubProps(s).PixelIdxList) = true;
                    intROI(c).ind = tempSubProps(s).PixelIdxList;
                    intROI(c).Nvox = tempSubProps(s).Area;
                    intROI(c).correlation = [];
                    intROI(c).preROI{z} = p; %keep track of which mask this ROI{z} is derived from
                    %{
                    sp(1) = subplot(1,2,1);
                    imshow( preROI{z}(p).mask_2d, [] ); axis image;
                    title(sprintf('Min area = %i', params.minFoot));
                    colorbar; %caxis
                    sp(2) = subplot(1,2,2);
                    imshow( intROI(c).mask_2d, [] ); axis image;
                    title(sprintf('c = %i', c));
                    impixelinfo;
                    linkaxes(sp,'xy')
                    pause;
                    %}
                end
            end
            Nint = numel(intROI);
            %VisualizeSegmentation( expt, intROI, 'overwrite',true );
    
            % Exclude pixels that overlap multiple intermediate ROI{z}, or merge intermediate ROI{z} with high degree of overlap
            % Determine pairwise overlap
            overlapInd = cell(Nint,Nint); overlapFrac = nan(Nint,Nint); Noverlap = nan(Nint,Nint);
            for r = 1:Nint-1
                for R = r+1:Nint
                    overlapInd{r,R} = intersect(intROI(r).ind, intROI(R).ind);
                    overlapInd{R,r} = overlapInd{r,R};
                    Noverlap(r,R) = numel(overlapInd{r,R});
                    Noverlap(R,r) = numel(overlapInd{r,R});
                    overlapFrac(r,R) = Noverlap(r,R)/intROI(r).Nvox;
                    overlapFrac(R,r) = Noverlap(r,R)/intROI(R).Nvox;
                end
            end
            %figure; imagesc(overlapFrac); impixelinfo; axis square;
            %figure; imagesc(Noverlap); impixelinfo; axis square;
            % Exclude overlapping voxels, or merge ROI{z}
            %figure('Units','normalized','OuterPosition',[0,0,1,1]);
            for r = 1:Nint
                for R = find(Noverlap(r,:) > 0)
                    fprintf('\n[r,R] = [%i, %i]: overlaps fractions = [%2.1f, %2.1f] pct (%i pixels)... ', r, R, 100*overlapFrac(r,R), 100*overlapFrac(R,r), Noverlap(r,R) );
                    %{
                    sp(1) = subplot(2,2,1); imshow( intROI(r).mask_2d, [] ); axis image;
                    title(sprintf('r = %i: %2.2f percent overlap', r, 100*overlapFrac(r,R)));
                    sp(2) = subplot(2,2,2); imshow(intROI(R).mask_2d, [] ); axis image;
                    title(sprintf('r = %i: %2.2f percent overlap', R, 100*overlapFrac(R,r)) )
                    %}   
                    if overlapFrac(r,R) < params.overlap && overlapFrac(R,r) < params.overlap % exclude overlapping pixels from both ROI{z}
                        fprintf('Excluding overlap from both ROI{z}');
                        intROI(r).ind = setdiff(intROI(r).ind, overlapInd{r,R});
                        intROI(r).mask_2d = falseVol; intROI(r).mask_2d(intROI(r).ind) = true;
                        intROI(R).ind = setdiff(intROI(R).ind, overlapInd{r,R});
                        intROI(R).mask_2d = falseVol; intROI(R).mask_2d(intROI(R).ind) = true;
                        %{
                        sp(3) = subplot(2,2,3); imshow( intROI(r).mask_2d, [] ); axis image;
                        title('Post-separation')
                        sp(4) = subplot(2,2,4); imshow(intROI(R).mask_2d, [] ); axis image;
                        title('Post-separation')
                        impixelinfo;
                        linkaxes(sp,'xy')
                        pause;
                        %}
                    else 
                        % Absorb the ROI{z} with higher overlapFrac into the other one
                        rPair = [r, R];
                        [~,absorbInd] = min([overlapFrac(r,R), overlapFrac(R,r)]);
                        rAbsorb = rPair(absorbInd);
                        rDelete = setdiff(rPair, rAbsorb);
                        fprintf('Absorbing R = %i into r = %i', rDelete, rAbsorb);
             
                        %intROI(rAbsorb).ind = unique(vertcat(intROI(rAbsorb).ind, intROI(rDelete).ind));
                        distinctInd = setdiff(intROI(rDelete).ind, intROI(rAbsorb).ind );
                        intROI(rAbsorb).ind = vertcat(intROI(rAbsorb).ind, distinctInd); %unique(vertcat(intROI(rAbsorb).ind, intROI(rDelete).ind));
                        intROI(rAbsorb).mask_2d(intROI(rAbsorb).ind) = true;
                        intROI(rDelete).ind = [];
                        intROI(rDelete).mask_2d = falseVol;
                        %{
                        sp(3) = subplot(2,2,3); imshow( intROI(r).mask_2d, [] ); axis image;
                        title('Post-absorption')
                        sp(4) = subplot(2,2,4); imshow(intROI(R).mask_2d, [] ); axis image;
                        title('Post-absorption')
                        impixelinfo;
                        linkaxes(sp,'xy')
                        pause;
                        %}
                    end
                    intROI(r).Nvox = numel(intROI(r).ind);
                    intROI(R).Nvox = numel(intROI(R).ind);
                    %pause;
                end
            end
            intROI([intROI.Nvox] < params.minFoot) = []; % remove intROI that have fallen below minimum area
            Nint = numel(intROI);
            
            r = 0; ROI{z} = struct();
            if show, figure('Units','normalized','OuterPosition',[0,0,1,1]); end
            for p = 1:Nint %numel(preROI{z})
                % Suppress overlapping pixels
                if show
                    subtightplot(1,3,1,opt{:}); imshow(preROI{z}(intROI(p).preROI{z}).mask_2d); title(sprintf('Initial ROI{z} %i', intROI(p).preROI{z}));   
                    subtightplot(1,3,2,opt{:}); imshow(intROI(p).mask_2d); title(sprintf('p = %i: Intermediate ROI{z}', p)); 
                end
                tempConnComp = bwconncomp( intROI(p).mask_2d, 8 );
                tempProps = regionprops( tempConnComp, 'Area', 'Centroid', 'BoundingBox','Eccentricity','Orientation', 'PixelList','PixelIdxList', 'MajoraxisLength', 'MinorAxisLength'); % projVol, ,'MeanIntensity', 'MaxIntensity'
                for s = find(cellfun(@numel, tempConnComp.PixelIdxList) > params.minFoot) %1:tempConnComp.NumObjects
                    r = r+1;
                    ROI{z}(r).ind = tempConnComp.PixelIdxList{s};
                    ROI{z}(r).mask_2d = falseVol; ROI{z}(r).mask_2d(ROI{z}(r).ind) = true; 
                    if show, subtightplot(1,3,3,opt{:}); imshow(ROI{z}(r).mask_2d); title(sprintf('r = %i: Final ROI{z}', r)); impixelinfo; pause; end
                    ROI{z}(r).mask_3d = []; % 
                    ROI{z}(r).correlation = [];
                    ROI{z}(r).Nvox = tempProps(s).Area;
                    ROI{z}(r).box_x = [tempProps(s).BoundingBox(1), tempProps(s).BoundingBox(1)+tempProps(s).BoundingBox(3)];%tempSubRegion.BoundingBox(:)
                    ROI{z}(r).box_y = [tempProps(s).BoundingBox(2), tempProps(s).BoundingBox(2)+tempProps(s).BoundingBox(4)];
                    ROI{z}(r).box_z = [1,1];
                    ROI{z}(r).cent = tempProps(s).Centroid;
                    ROI{z}(r).solidity = NaN;
                    ROI{z}(r).orientation = tempProps(s).Orientation;
                    ROI{z}(r).eccentricity = tempProps(s).Eccentricity;
                    ROI{z}(r).corrMax = NaN;
                    ROI{z}(r).corrMean = NaN;
                    ROI{z}(r).PrinAxLength = [tempProps(s).MajorAxisLength, tempProps(s).MinorAxisLength, NaN]; %NaN;
                    ROI{z}(r).footprintInd = ROI{z}(r).ind;
                    ROI{z}(r).footprint = ROI{z}(r).Nvox; 
                    [footRow, footCol] = ind2sub( [expt.Nrow,expt.Ncol], ROI{z}(r).footprintInd );
                    ROI{z}(r).footprintXY = [footRow, footCol];
                    tempBinary = ROI{z}(r).mask_2d; %false( expt.Nrow, expt.Ncol ); tempBinary(ROI{z}(r).footprintInd) = true; 
                    tempBinary = imdilate(tempBinary, strel('rectangle',[3,3]) ); % imshow( tempBinary )
                    tempEdge = edge(tempBinary);
                    [ROI{z}(r).footprintEdge(:,1), ROI{z}(r).footprintEdge(:,2)] = ind2sub( [expt.Nrow, expt.Ncol], find( tempEdge ) );
                    % How many voxels in each plane?
                    ROI{z}(r).zProfile = ROI{z}(r).Nvox;
                    ROI{z}(r).zUse = 1; %find( squeeze(sum(sum(ROI{z}(r).mask_3d, 1), 2)) ); % which planes were used for this ROI{z}?
                    ROI{z}(r).preROI{z} = intROI(p).preROI{z}; %p;
                    ROI{z}(r).axon = 0;
                    % Get cropped summary projection images (may be overwritten usign WriteROIproj later)
                    ROI{z}(r).crop = [floor(ROI{z}(r).box_x(1)), ceil(ROI{z}(r).box_x(2)), floor(ROI{z}(r).box_y(1)), ceil(ROI{z}(r).box_y(2))];
                    %ROI{z}(k).crop = [floor(ROI{z}(k).box_x(1)), ceil(ROI{z}(k).box_x(2)), floor(ROI{z}(k).box_y(1)), ceil(ROI{z}(k).box_y(2))];
                    ROI{z}(r).projPath = ''; %projPath;
                    ROI{z}(r).cropMean = expt.meanProj(ROI{z}(r).crop(3)+1:expt.Nrow-ROI{z}(r).crop(4), ROI{z}(r).crop(1)+1:expt.Ncol-ROI{z}(r).crop(2)); % imshow( ROI{z}(r).cropMean, [] )  
                    ROI{z}(r).cropMax = []; 
                    ROI{z}(r).cropStd = []; 
                    %{
                    close all;
                    imshow( imadjust(ROI{z}(r).zProj), [] ); hold on; % 
                    plot( ROI{z}(r).footprintEdge(:,2), ROI{z}(r).footprintEdge(:,1), '.b' ); hold on;
                    %}
                end
            end
            Nroi(z) = numel( ROI{z} );
            %VisualizeSegmentation( expt, ROI{z}, 'overwrite',true );
            %}
            %{
            preROIclaim = cat(3, preROI{z}.mask_2d); % which preROIs claim a given pixel?
            preROIclaimTot = sum(preROIclaim,3); % how many preROIs claim a given pixel?
            opt = {[0.04, 0.03], [0.06, 0.05], [0.03,0.05]};  % {[vert, horz], [bottom, top], [left, right] }
            subtightplot(1,3,1,opt{:}); imshow(preROIclaimTot,[]);
            overlapPix = find( preROIclaimTot > 1); %unique(preROIoccTot);
            allClaimInd = find(preROIclaimTot); %find( allMask(:) ); 
    
            r = 0;
            for p = 1:numel(preROI{z})
                % Suppress overlapping pixels
                tempMask = preROI{z}(p).mask_2d; 
                subtightplot(1,3,2,opt{:}); imshow(preROI{z}(p).mask_2d)
                
                tempMask(intersect(find(preROI{z}(p).mask_2d),overlapPix )) = false;      
                tempConnComp = bwconncomp( tempMask, 8 );
                tempProps = regionprops( tempConnComp, 'Area', 'Centroid', 'BoundingBox','Eccentricity','Orientation', 'PixelList','PixelIdxList', 'MajoraxisLength', 'MinorAxisLength'); % projVol, ,'MeanIntensity', 'MaxIntensity'
                %ROI{z}(r).mean_trace_raw = preROI{z}(p).mean_trace_raw;
                %ROI{z}(r).mean_trace_hp = preROI{z}(p).mean_trace_hp;
                for s = find(cellfun(@numel, tempConnComp.PixelIdxList) > params.minFoot) %1:tempConnComp.NumObjects
                    r = r+1;
                    ROI{z}(r).ind = tempConnComp.PixelIdxList{s};
                    ROI{z}(r).mask_2d = falseVol; ROI{z}(r).mask_2d(ROI{z}(r).ind) = true; 
                    subtightplot(1,3,3,opt{:}); imshow(ROI{z}(r).mask_2d); title(sprintf('[p, r] = [%i, %i]',p,r)); impixelinfo; pause;
                    ROI{z}(r).mask_3d = []; % 
                    ROI{z}(r).correlation = [];
                    ROI{z}(r).Nvox = tempProps(s).Area;
                    ROI{z}(r).box_x = [tempProps(s).BoundingBox(1), tempProps(s).BoundingBox(1)+tempProps(s).BoundingBox(3)];%tempSubRegion.BoundingBox(:)
                    ROI{z}(r).box_y = [tempProps(s).BoundingBox(2), tempProps(s).BoundingBox(2)+tempProps(s).BoundingBox(4)];
                    ROI{z}(r).box_z = [1,1];
                    ROI{z}(r).cent = tempProps(s).Centroid;
                    ROI{z}(r).solidity = NaN;
                    ROI{z}(r).orientation = tempProps(s).Orientation;
                    ROI{z}(r).eccentricity = tempProps(s).Eccentricity;
                    ROI{z}(r).corrMax = NaN;
                    ROI{z}(r).corrMean = NaN;
                    ROI{z}(r).PrinAxLength = [tempProps(s).MajorAxisLength, tempProps(s).MinorAxisLength, NaN]; %NaN;
                    ROI{z}(r).footprintInd = ROI{z}(r).ind;
                    ROI{z}(r).footprint = ROI{z}(r).Nvox; 
                    [footRow, footCol] = ind2sub( [expt.Nrow,expt.Ncol], ROI{z}(r).footprintInd );
                    ROI{z}(r).footprintXY = [footRow, footCol];
                    tempBinary = ROI{z}(r).mask_2d; %false( expt.Nrow, expt.Ncol ); tempBinary(ROI{z}(r).footprintInd) = true; 
                    tempBinary = imdilate(tempBinary, strel('rectangle',[3,3]) ); % imshow( tempBinary )
                    tempEdge = edge(tempBinary);
                    [ROI{z}(r).footprintEdge(:,1), ROI{z}(r).footprintEdge(:,2)] = ind2sub( [expt.Nrow, expt.Ncol], find( tempEdge ) );
                    % How many voxels in each plane?
                    ROI{z}(r).zProfile = ROI{z}(r).Nvox;
                    ROI{z}(r).zUse = 1; %find( squeeze(sum(sum(ROI{z}(r).mask_3d, 1), 2)) ); % which planes were used for this ROI{z}?
                    ROI{z}(r).preROI{z} = p;
                    ROI{z}(r).axon = 0;
                    % Get cropped summary projection images (may be overwritten usign WriteROIproj later)
                    ROI{z}(r).crop = [floor(ROI{z}(r).box_x(1)), ceil(ROI{z}(r).box_x(2)), floor(ROI{z}(r).box_y(1)), ceil(ROI{z}(r).box_y(2))];
                    %ROI{z}(k).crop = [floor(ROI{z}(k).box_x(1)), ceil(ROI{z}(k).box_x(2)), floor(ROI{z}(k).box_y(1)), ceil(ROI{z}(k).box_y(2))];
                    ROI{z}(r).projPath = ''; %projPath;
                    ROI{z}(r).cropMean = expt.meanProj(ROI{z}(r).crop(3)+1:expt.Nrow-ROI{z}(r).crop(4), ROI{z}(r).crop(1)+1:expt.Ncol-ROI{z}(r).crop(2)); % imshow( ROI{z}(r).cropMean, [] )  
                    ROI{z}(r).cropMax = []; 
                    ROI{z}(r).cropStd = []; 
                    %{
                    close all;
                    imshow( imadjust(ROI{z}(r).zProj), [] ); hold on; % 
                    plot( ROI{z}(r).footprintEdge(:,2), ROI{z}(r).footprintEdge(:,1), '.b' ); hold on;
                    %}
                end
            end
            %}
            % Construct neuropil regions, excluding any preROI{z}'s pixels
            preROIclaim = cat(3, preROI{z}.mask_2d); % which preROIs claim a given pixel?
            preROIclaimTot = sum(preROIclaim,3); % how many preROIs claim a given pixel?
            allPreInd = find(preROIclaimTot);
            if show
                allMask = any(cat(3, ROI{z}(:).mask_2d), 3); % Combined mask from all surviving ROI{z}
                figure;
                subtightplot(1,2,1,opt{:}); imshow(preROIclaimTot,[]);
                subtightplot(1,2,2,opt{:}); imshow(allMask,[]);
            end
            fprintf('\nGenerating neuropil ROI{z}...');
            for r = 1:Nroi(z)
                tempSubNP = imdilate(ROI{z}(r).mask_2d, strel('cube',params.neuropilRad(2))) & ~imdilate(ROI{z}(r).mask_2d,strel('cube',params.neuropilRad(1))); % imshow( max(tempSubNP,[],3), [] );
                tempNPind = find(tempSubNP);
                [~,tempNPbad, ~] = intersect(tempNPind, allPreInd );
                tempNPind(tempNPbad) = []; 
                ROI{z}(r).neuropil = tempNPind;
                % Construct labeled projection image
                ROI{z}(r).labelVol = zeros(expt.Nrow, expt.Ncol, expt.Nplane, 'uint8');
                ROI{z}(r).labelVol(ROI{z}(r).neuropil) = 1; ROI{z}(r).labelVol(ROI{z}(r).ind) = 2; % saveastiff(ROI{z}(r).labelVol, 'D:\2photon\DL102\180426_DL102\ROI1_label.tif')
            end
        else
            allPreCorr = cat(4,preROI{z}.correlation); %histogram(allPreCorr(:) )
            corrThresh = prctile(allPreCorr(:), params.minCorrPrct); 
            allPreCorrMax = max(allPreCorr,[],4);
            allPreMask3D = falseVol;  allPreMask3D(allPreCorrMax > corrThresh) = true;
            %saveastiff( uint16(allPreMask3D), 'D:\2photon\DL75\170525\allPreMask3D.tif');
            %close all; figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
            % Get intermediate 3D masks by breaking preROI{z} (PCA/ICA results) up into blobs and thresholding on max correlation
            c = 0;
            for p = 1:numel(preROI{z})
                tempConnComp = bwconncomp( preROI{z}(p).correlation > corrThresh, 18 );
                tempSubProps = regionprops3( tempConnComp, preROI{z}(p).correlation, 'Volume','MaxIntensity','VoxelIdxList'); % , preROI{z}(p).correlation
                goodSub = find(tempSubProps.Volume > params.minVol & tempSubProps.MaxIntensity > corrThresh)';
                preROI{z}(p).Nsub = numel(goodSub);
                for s = goodSub
                    c = c + 1;
                    % Retain some info from preliminary ROI{z}
                    intROI(c).mask_2d = preROI{z}(p).mask_2d;
                    intROI(c).ind = tempSubProps.VoxelIdxList{s};
                    intROI(c).Nvox = tempSubProps.Volume(s);
                    intROI(c).correlation = single(nanVol);
                    intROI(c).correlation(intROI(c).ind) = single(preROI{z}(p).correlation(intROI(c).ind));
                    intROI(c).preROI{z} = p; %keep track of which mask this ROI{z} is derived from
                    %{
                    sp(1) = subplot(1,2,1);
                    imagesc( max(preROI{z}(p).correlation,[],3) ); axis image;
                    title(sprintf('Min volume = %i,  Correlation threshold = %1.2f', params.minVol, corrThresh));
                    colorbar; %caxis
                    sp(2) = subplot(1,2,2);
                    imagesc( max(intROI(c).correlation,[],3) ); axis image;
                    title(sprintf('c = %i', c));
                    impixelinfo;
                    linkaxes(sp,'xy')
                    pause;
                    %}
                end
            end
            Nint = numel(intROI);
            %VisualizeSegmentation( expt, intROI, 'overwrite',true );
    
            % Exclude voxels that overlap multiple intermediate ROI{z}, or merge ROI{z} with high degree of overlap
            % Determine pairwise overlap
            overlapInd = cell(Nint,Nint); overlapFrac = nan(Nint,Nint); Noverlap = nan(Nint,Nint);
            for r = 1:Nint-1
                for R = r+1:Nint
                    overlapInd{r,R} = intersect(intROI(r).ind, intROI(R).ind);
                    overlapInd{R,r} = overlapInd{r,R};
                    Noverlap(r,R) = numel(overlapInd{r,R});
                    Noverlap(R,r) = numel(overlapInd{r,R});
                    overlapFrac(r,R) = Noverlap(r,R)/intROI(r).Nvox;
                    overlapFrac(R,r) = Noverlap(r,R)/intROI(R).Nvox;
                end
            end
            %imagesc(overlapFrac); impixelinfo; axis square;
            %imagesc(Noverlap); impixelinfo; axis square;
            % Exclude overlapping voxels, or merge ROI{z}
            for r = 1:Nint
                for R = find(Noverlap(r,:) > 0)
                    fprintf('\n[r,R] = [%i, %i]: overlaps fractions = [%2.1f, %2.1f] pct (%i voxels)... ', r, R, 100*overlapFrac(r,R), 100*overlapFrac(R,r), Noverlap(r,R) );
                    if overlapFrac(r,R) < params.overlap && overlapFrac(R,r) < params.overlap % exclude from both ROI{z}
                        fprintf('Excluding overlap from both ROI{z}');
                        intROI(r).ind = setdiff(intROI(r).ind, overlapInd{r,R});
                        intROI(R).ind = setdiff(intROI(R).ind, overlapInd{r,R});
                    else % overlapFrac(r,R) >= params.overlap && overlapFrac(R,r) >= params.overlap
                        % Absorb the ROI{z} with higher overlapFrac into the other one
                        rPair = [r, R];
                        [~,absorbInd] = min([overlapFrac(r,R), overlapFrac(R,r)]);
                        rAbsorb = rPair(absorbInd);
                        rDelete = setdiff(rPair, rAbsorb);
                        fprintf('Absorbing R = %i into r = %i', rDelete, rAbsorb);
                        %{
                        sp(1) = subplot(1,3,1); imagesc( max(intROI(r).correlation,[],3) ); axis image;
                        title(sprintf('r = %i', r));
                        sp(2) = subplot(1,3,2); imagesc( max(intROI(R).correlation,[],3) ); axis image;
                        title(sprintf('r = %i', R))
                        %}                  
                        %intROI(rAbsorb).ind = unique(vertcat(intROI(rAbsorb).ind, intROI(rDelete).ind));
                        distinctInd = setdiff(intROI(rDelete).ind, intROI(rAbsorb).ind );
                        intROI(rAbsorb).ind = vertcat(intROI(rAbsorb).ind, distinctInd); %unique(vertcat(intROI(rAbsorb).ind, intROI(rDelete).ind));
                        intROI(rAbsorb).correlation(distinctInd) = intROI(rDelete).correlation(distinctInd);
                        intROI(rDelete).ind = [];
                        %{
                        sp(3) = subplot(1,3,3); imagesc(max(intROI(rAbsorb).correlation,[],3)); axis image; 
                        title(sprintf('r = %i post-absorbtion', rAbsorb))
                        impixelinfo;
                        linkaxes(sp,'xy')
                        %}
                    end
                    intROI(r).Nvox = numel(intROI(r).ind);
                    intROI(R).Nvox = numel(intROI(R).ind);
                    %pause;
                end
            end
            %VisualizeSegmentation( expt, intROI, 'overwrite',true );
                  
            k = 0;
            for r = find([intROI.Nvox] >= params.minVol)%1:Nint
                tempMask3D = falseVol; tempMask3D(intROI(r).ind) = true; % imshow( max(tempMask3D,[],3), [])
                tempConnComp = bwconncomp( tempMask3D, 18 );
                tempProps3D = regionprops3( tempConnComp, intROI(r).correlation, 'Volume', 'Centroid', 'BoundingBox', 'MeanIntensity', 'MaxIntensity',...
                    'VoxelList', 'VoxelIdxList','PrincipalAxisLength', 'Solidity' ); % 
                for s = find(tempProps3D.Volume > params.minVol)'
                    k = k+1;
                    ROI{z}(k).preROI{z} = intROI(r).preROI{z}; 
                    % Update mask/indices/correlation
                    ROI{z}(k).mask_3d = falseVol;
                    ROI{z}(k).mask_3d(tempConnComp.PixelIdxList{s}) = true; % imshow( max(ROI{z}(k).mask_3d,[],3), [])
                    ROI{z}(k).ind = find(ROI{z}(k).mask_3d);
                    ROI{z}(k).Nvox = numel(ROI{z}(k).ind);
                    ROI{z}(k).correlation = intROI(r).correlation;
                    ROI{z}(k).correlation(~ROI{z}(k).mask_3d) = NaN; % imshow( max(ROI{z}(k).correlation,[],3), [])               
                    % Get additional info about the 3D ROI{z}
                    ROI{z}(k).box_x = [tempProps3D.BoundingBox(s,1), tempProps3D.BoundingBox(s,1)+tempProps3D.BoundingBox(s,4)];%tempSubRegion.BoundingBox(s,:)
                    ROI{z}(k).box_y = [tempProps3D.BoundingBox(s,2), tempProps3D.BoundingBox(s,2)+tempProps3D.BoundingBox(s,5)];
                    ROI{z}(k).box_z = [tempProps3D.BoundingBox(s,3), tempProps3D.BoundingBox(s,3)+tempProps3D.BoundingBox(s,6)];
                    ROI{z}(k).cent = tempProps3D.Centroid(s,:);
                    ROI{z}(k).solidity = tempProps3D.Solidity(s); %how solid was this volume before filling?
                    ROI{z}(k).corrMax = tempProps3D.MaxIntensity(s);
                    ROI{z}(k).corrMean = tempProps3D.MeanIntensity(s);
                    ROI{z}(k).PrinAxLength = tempProps3D.PrincipalAxisLength(s,:);
                    % Get info about the ROI{z}'s footprint
                    ROI{z}(k).mask_2d = max(ROI{z}(k).mask_3d, [], 3); % imshow( max(ROI{z}(k).mask_2d,[],3), [])
                    ROI{z}(k).footprint = sum(ROI{z}(k).mask_2d(:)); %numel(ROI{z}(k).footprintInd); 
                    ROI{z}(k).footprintInd = find( ROI{z}(k).mask_2d );
                    [footRow, footCol] = ind2sub( [expt.Nrow,expt.Ncol], ROI{z}(k).footprintInd );
                    ROI{z}(k).footprintXY = [footCol, footRow];
                    tempFootProp = regionprops(ROI{z}(k).mask_2d, 'Orientation', 'Eccentricity');
                    ROI{z}(k).orientation = tempFootProp.Orientation;
                    ROI{z}(k).eccentricity = tempFootProp.Eccentricity;
                    %close all; figure('Units','normalized', 'OuterPosition',[0,0,1,1]); 
                    %imshow(tempSubBinary, []); title( sprintf('Orientation = %2.2f,  Ecc = %2.2f', ROI{z}(k).orientation, ROI{z}(k).eccentricity) ); pause;
                    tempFootBinary = imdilate(ROI{z}(k).mask_2d, strel('rectangle',[3,3]) ); % imshow( tempSubBinary )
                    tempFootEdge = edge(tempFootBinary);
                    [ROI{z}(k).footprintEdge(:,1), ROI{z}(k).footprintEdge(:,2)] = ind2sub( [expt.Nrow, expt.Ncol], find( tempFootEdge ) );
                    ROI{z}(k).axon = 0;
                    % How many voxels in each plane?
                    ROI{z}(k).zProfile = squeeze(sum(ROI{z}(k).mask_3d, [1,2])); %zeros(1, expt.Nplane);
                    ROI{z}(k).zUse = find(ROI{z}(k).zProfile); %find( squeeze(sum(sum(ROI{z}(k).mask_3d, 1), 2)) ); % which planes were used for this ROI{z}?
    
                    % Setup cropped summary projection images (see WriteROIproj)
                    ROI{z}(k).projPath = '';
                    ROI{z}(k).crop = [floor(ROI{z}(k).box_x(1)), ceil(ROI{z}(k).box_x(2)), floor(ROI{z}(k).box_y(1)), ceil(ROI{z}(k).box_y(2))]; %[]; expt.Ncol-  expt.Nrow-
                    ROI{z}(k).cropMean = []; 
                    ROI{z}(k).cropStd = [];
                    ROI{z}(k).cropMax = []; 
                end
            end
            Nroi(z) = numel(ROI{z});
            
            % Sort ROI{z} from deep to superficial
            roiCent = vertcat(ROI{z}.cent);
            [~,rSort] = sort( roiCent(:,3), 'ascend' ); % 1 is deepest, Nplane is most superficial
            ROI{z} = ROI{z}(rSort);
            %VisualizeSegmentation( expt, intROI, 'overwrite',true );
            
            % Construct neuropil regions, excluding any intROI's voxels
            allPreMask = allPreMask3D; %any( cat(4, ROI{z}(:).mask_3d), 4 ); 
            allPreInd = find( allPreMask(:) );
            zExclude = find(squeeze(~any(allPreMask,[1,2])));  % find(squeeze(~any(mask_3d_all, [1,2,3])))'; % determine which planes were excluded from analysis
            for r = 1:Nroi(z)
                tempSubNP = imdilate(ROI{z}(r).mask_3d, strel('cuboid',[params.neuropilRad(2),params.neuropilRad(2),2])) & ...
                    ~imdilate(ROI{z}(r).mask_3d,strel('cuboid',[params.neuropilRad(1),params.neuropilRad(1),1])); % imshow( max(tempSubNP,[],3), [] );
                tempSubNP(:,:,zExclude) = false; % suppress planes that were excluded from segmentation due to poor registration
                tempNPind = find(tempSubNP);
                [~,tempNPbad, ~] = intersect(tempNPind, allPreInd );
                tempNPind(tempNPbad) = []; 
                ROI{z}(r).neuropil = tempNPind;
                % Construct labeled projection image
                ROI{z}(r).labelVol = zeros(expt.Nrow, expt.Ncol, expt.Nplane, 'uint8');
                ROI{z}(r).labelVol(ROI{z}(r).neuropil) = 1; ROI{z}(r).labelVol(ROI{z}(r).ind) = 2; % 
            end
            clearvars tempSubNP tempNPind
        end
        ROI{z}( [ROI{z}.footprint] < params.minFoot | [ROI{z}.Nvox] < params.minVol ) = [];

        if isfield(ROI{z}, 'correlation') && removeCorr
            fprintf('\nRemoving correlation field...')
            ROI{z} = rmfield(ROI{z}, 'correlation'); % correlation field is memory-intensive and not useful
        end
    end

    % Save the results
    fprintf('\nSaving %s', savePath);
    save(savePath, 'ROI','preROI','intROI','params', '-v7.3'); 
    toc
else
    fprintf('\nLoading %s...  ', savePath);
    if nargout > 1
        load(savePath, 'ROI', 'preROI');   % ,'params'
    else
        load(savePath, 'ROI');
    end
end

if iscell(ROI), ROI = [ROI{:}]; end

if isfield(ROI, 'correlation') && removeCorr
    fprintf('\nRemoving correlation field...')
    ROI = rmfield(ROI, 'correlation'); % correlation field is memory-intensive and not useful
end
toc


%{
tempCent = vertcat(ROI{z}.cent)
find( tempCent(:,1) > 480 & tempCent(:,2) > 360)
ROI{z}(find( tempCent(:,1) > 480 & tempCent(:,2) > 360)) = [];

runningOverlap = zeros(size(preROI{z}(1).mask_2d));
figure;
for p = [51,53,56]% 1:numel(preROI{z})
    subplot(2,1,1);
    imshow(preROI{z}(p).mask_2d)
    title(sprintf('p = %i',p)); 
    
    runningOverlap = runningOverlap + preROI{z}(p).mask_2d;
    subplot(2,1,2);
    imshow(runningOverlap, [])
    impixelinfo; pause;
end

%}