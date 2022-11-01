function [ROI, preROI] = MakeROI3D(expt, varargin) % , run exptDir, exptName %ROIlabel, 
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addParameter( IP, 'corrPrct', 90, @isnumeric)  % 68
addParameter( IP, 'minVol', 50, @isnumeric) 
addParameter( IP, 'minFoot', 50, @isnumeric) 
addParameter( IP, 'neuropilRad', [8,21], @isnumeric ) % inner and outer radii (in xy plane) to draw neuropil ROI
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
    seg3Dpath = sprintf('%s%s_seg_data_final.mat', expt.dir, expt.name);  % tracePath = sprintf('%s%s_seg_traces.mat', expt.dir, expt.name);   
    fprintf('\nLoading %s... ', seg3Dpath); 
    load(seg3Dpath, 'seg_data'); 
    preROI = seg_data([seg_data.keep]); %seg_data_keep;
    if expt.Nplane == 1       
        % {
        % Get intermediate 2D masks by breaking preROI (PCA/ICA results) up into blobs and thresholding on max correlation
        fprintf('\nBreaking up preliminary ROIs')
        c = 0; intROI = struct();
        for p = 1:numel(preROI)
            tempConnComp = bwconncomp( preROI(p).mask_2d, 8 );
            tempSubProps = regionprops( tempConnComp, 'Area','PixelIdxList'); % , preROI(p).correlation
            goodSub = find([tempSubProps.Area] > params.minFoot );
            preROI(p).Nsub = numel(goodSub);
            tempBinary = preROI(p).mask_2d; %false( expt.Nrow, expt.Ncol ); tempBinary(ROI(r).footprintInd) = true; 
            tempBinary = imdilate(tempBinary, strel('rectangle',[3,3]) ); % imshow( tempBinary )
            tempEdge = edge(tempBinary);
            [preROI(p).footprintEdge(:,1), preROI(p).footprintEdge(:,2)] = ind2sub( [expt.Nrow, expt.Ncol], find( tempEdge ) );
            for s = goodSub
                c = c + 1;
                % Retain some info from preliminary ROI
                intROI(c).mask_2d = falseVol; %preROI(p).mask_2d;
                intROI(c).mask_2d(tempSubProps(s).PixelIdxList) = true;
                intROI(c).ind = tempSubProps(s).PixelIdxList;
                intROI(c).Nvox = tempSubProps(s).Area;
                intROI(c).correlation = [];
                intROI(c).preROI = p; %keep track of which mask this ROI is derived from
                %{
                sp(1) = subplot(1,2,1);
                imshow( preROI(p).mask_2d, [] ); axis image;
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

        % Exclude pixels that overlap multiple intermediate ROI, or merge intermediate ROI with high degree of overlap
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
        % Exclude overlapping voxels, or merge ROI
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
                if overlapFrac(r,R) < params.overlap && overlapFrac(R,r) < params.overlap % exclude overlapping pixels from both ROI
                    fprintf('Excluding overlap from both ROI');
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
                    % Absorb the ROI with higher overlapFrac into the other one
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
        
        r = 0; ROI = struct();
        if show, figure('Units','normalized','OuterPosition',[0,0,1,1]); end
        for p = 1:Nint %numel(preROI)
            % Suppress overlapping pixels
            if show
                subtightplot(1,3,1,opt{:}); imshow(preROI(intROI(p).preROI).mask_2d); title(sprintf('Initial ROI %i', intROI(p).preROI));   
                subtightplot(1,3,2,opt{:}); imshow(intROI(p).mask_2d); title(sprintf('p = %i: Intermediate ROI', p)); 
            end
            tempConnComp = bwconncomp( intROI(p).mask_2d, 8 );
            tempProps = regionprops( tempConnComp, 'Area', 'Centroid', 'BoundingBox','Eccentricity','Orientation', 'PixelList','PixelIdxList', 'MajoraxisLength', 'MinorAxisLength'); % projVol, ,'MeanIntensity', 'MaxIntensity'
            for s = find(cellfun(@numel, tempConnComp.PixelIdxList) > params.minFoot) %1:tempConnComp.NumObjects
                r = r+1;
                ROI(r).ind = tempConnComp.PixelIdxList{s};
                ROI(r).mask_2d = falseVol; ROI(r).mask_2d(ROI(r).ind) = true; 
                if show, subtightplot(1,3,3,opt{:}); imshow(ROI(r).mask_2d); title(sprintf('r = %i: Final ROI', r)); impixelinfo; pause; end
                ROI(r).mask_3d = []; % 
                ROI(r).correlation = [];
                ROI(r).Nvox = tempProps(s).Area;
                ROI(r).box_x = [tempProps(s).BoundingBox(1), tempProps(s).BoundingBox(1)+tempProps(s).BoundingBox(3)];%tempSubRegion.BoundingBox(:)
                ROI(r).box_y = [tempProps(s).BoundingBox(2), tempProps(s).BoundingBox(2)+tempProps(s).BoundingBox(4)];
                ROI(r).box_z = [1,1];
                ROI(r).cent = tempProps(s).Centroid;
                ROI(r).solidity = NaN;
                ROI(r).orientation = tempProps(s).Orientation;
                ROI(r).eccentricity = tempProps(s).Eccentricity;
                ROI(r).corrMax = NaN;
                ROI(r).corrMean = NaN;
                ROI(r).PrinAxLength = [tempProps(s).MajorAxisLength, tempProps(s).MinorAxisLength, NaN]; %NaN;
                ROI(r).footprintInd = ROI(r).ind;
                ROI(r).footprint = ROI(r).Nvox; 
                [footRow, footCol] = ind2sub( [expt.Nrow,expt.Ncol], ROI(r).footprintInd );
                ROI(r).footprintXY = [footRow, footCol];
                tempBinary = ROI(r).mask_2d; %false( expt.Nrow, expt.Ncol ); tempBinary(ROI(r).footprintInd) = true; 
                tempBinary = imdilate(tempBinary, strel('rectangle',[3,3]) ); % imshow( tempBinary )
                tempEdge = edge(tempBinary);
                [ROI(r).footprintEdge(:,1), ROI(r).footprintEdge(:,2)] = ind2sub( [expt.Nrow, expt.Ncol], find( tempEdge ) );
                % How many voxels in each plane?
                ROI(r).zProfile = ROI(r).Nvox;
                ROI(r).zUse = 1; %find( squeeze(sum(sum(ROI(r).mask_3d, 1), 2)) ); % which planes were used for this ROI?
                ROI(r).preROI = intROI(p).preROI; %p;
                ROI(r).axon = 0;
                % Get cropped summary projection images (may be overwritten usign WriteROIproj later)
                ROI(r).crop = [floor(ROI(r).box_x(1)), ceil(ROI(r).box_x(2)), floor(ROI(r).box_y(1)), ceil(ROI(r).box_y(2))];
                %ROI(k).crop = [floor(ROI(k).box_x(1)), ceil(ROI(k).box_x(2)), floor(ROI(k).box_y(1)), ceil(ROI(k).box_y(2))];
                ROI(r).projPath = ''; %projPath;
                ROI(r).cropMean = expt.meanProj(ROI(r).crop(3)+1:expt.Nrow-ROI(r).crop(4), ROI(r).crop(1)+1:expt.Ncol-ROI(r).crop(2)); % imshow( ROI(r).cropMean, [] )  
                ROI(r).cropMax = []; 
                ROI(r).cropStd = []; 
                %{
                close all;
                imshow( imadjust(ROI(r).zProj), [] ); hold on; % 
                plot( ROI(r).footprintEdge(:,2), ROI(r).footprintEdge(:,1), '.b' ); hold on;
                %}
            end
        end
        Nroi = numel( ROI );
        %VisualizeSegmentation( expt, ROI, 'overwrite',true );
        %}
        %{
        preROIclaim = cat(3, preROI.mask_2d); % which preROIs claim a given pixel?
        preROIclaimTot = sum(preROIclaim,3); % how many preROIs claim a given pixel?
        opt = {[0.04, 0.03], [0.06, 0.05], [0.03,0.05]};  % {[vert, horz], [bottom, top], [left, right] }
        subtightplot(1,3,1,opt{:}); imshow(preROIclaimTot,[]);
        overlapPix = find( preROIclaimTot > 1); %unique(preROIoccTot);
        allClaimInd = find(preROIclaimTot); %find( allMask(:) ); 

        r = 0;
        for p = 1:numel(preROI)
            % Suppress overlapping pixels
            tempMask = preROI(p).mask_2d; 
            subtightplot(1,3,2,opt{:}); imshow(preROI(p).mask_2d)
            
            tempMask(intersect(find(preROI(p).mask_2d),overlapPix )) = false;      
            tempConnComp = bwconncomp( tempMask, 8 );
            tempProps = regionprops( tempConnComp, 'Area', 'Centroid', 'BoundingBox','Eccentricity','Orientation', 'PixelList','PixelIdxList', 'MajoraxisLength', 'MinorAxisLength'); % projVol, ,'MeanIntensity', 'MaxIntensity'
            %ROI(r).mean_trace_raw = preROI(p).mean_trace_raw;
            %ROI(r).mean_trace_hp = preROI(p).mean_trace_hp;
            for s = find(cellfun(@numel, tempConnComp.PixelIdxList) > params.minFoot) %1:tempConnComp.NumObjects
                r = r+1;
                ROI(r).ind = tempConnComp.PixelIdxList{s};
                ROI(r).mask_2d = falseVol; ROI(r).mask_2d(ROI(r).ind) = true; 
                subtightplot(1,3,3,opt{:}); imshow(ROI(r).mask_2d); title(sprintf('[p, r] = [%i, %i]',p,r)); impixelinfo; pause;
                ROI(r).mask_3d = []; % 
                ROI(r).correlation = [];
                ROI(r).Nvox = tempProps(s).Area;
                ROI(r).box_x = [tempProps(s).BoundingBox(1), tempProps(s).BoundingBox(1)+tempProps(s).BoundingBox(3)];%tempSubRegion.BoundingBox(:)
                ROI(r).box_y = [tempProps(s).BoundingBox(2), tempProps(s).BoundingBox(2)+tempProps(s).BoundingBox(4)];
                ROI(r).box_z = [1,1];
                ROI(r).cent = tempProps(s).Centroid;
                ROI(r).solidity = NaN;
                ROI(r).orientation = tempProps(s).Orientation;
                ROI(r).eccentricity = tempProps(s).Eccentricity;
                ROI(r).corrMax = NaN;
                ROI(r).corrMean = NaN;
                ROI(r).PrinAxLength = [tempProps(s).MajorAxisLength, tempProps(s).MinorAxisLength, NaN]; %NaN;
                ROI(r).footprintInd = ROI(r).ind;
                ROI(r).footprint = ROI(r).Nvox; 
                [footRow, footCol] = ind2sub( [expt.Nrow,expt.Ncol], ROI(r).footprintInd );
                ROI(r).footprintXY = [footRow, footCol];
                tempBinary = ROI(r).mask_2d; %false( expt.Nrow, expt.Ncol ); tempBinary(ROI(r).footprintInd) = true; 
                tempBinary = imdilate(tempBinary, strel('rectangle',[3,3]) ); % imshow( tempBinary )
                tempEdge = edge(tempBinary);
                [ROI(r).footprintEdge(:,1), ROI(r).footprintEdge(:,2)] = ind2sub( [expt.Nrow, expt.Ncol], find( tempEdge ) );
                % How many voxels in each plane?
                ROI(r).zProfile = ROI(r).Nvox;
                ROI(r).zUse = 1; %find( squeeze(sum(sum(ROI(r).mask_3d, 1), 2)) ); % which planes were used for this ROI?
                ROI(r).preROI = p;
                ROI(r).axon = 0;
                % Get cropped summary projection images (may be overwritten usign WriteROIproj later)
                ROI(r).crop = [floor(ROI(r).box_x(1)), ceil(ROI(r).box_x(2)), floor(ROI(r).box_y(1)), ceil(ROI(r).box_y(2))];
                %ROI(k).crop = [floor(ROI(k).box_x(1)), ceil(ROI(k).box_x(2)), floor(ROI(k).box_y(1)), ceil(ROI(k).box_y(2))];
                ROI(r).projPath = ''; %projPath;
                ROI(r).cropMean = expt.meanProj(ROI(r).crop(3)+1:expt.Nrow-ROI(r).crop(4), ROI(r).crop(1)+1:expt.Ncol-ROI(r).crop(2)); % imshow( ROI(r).cropMean, [] )  
                ROI(r).cropMax = []; 
                ROI(r).cropStd = []; 
                %{
                close all;
                imshow( imadjust(ROI(r).zProj), [] ); hold on; % 
                plot( ROI(r).footprintEdge(:,2), ROI(r).footprintEdge(:,1), '.b' ); hold on;
                %}
            end
        end
        %}
        % Construct neuropil regions, excluding any preROI's pixels
        preROIclaim = cat(3, preROI.mask_2d); % which preROIs claim a given pixel?
        preROIclaimTot = sum(preROIclaim,3); % how many preROIs claim a given pixel?
        allPreInd = find(preROIclaimTot);
        if show
            allMask = any(cat(3, ROI(:).mask_2d), 3); % Combined mask from all surviving ROI
            figure;
            subtightplot(1,2,1,opt{:}); imshow(preROIclaimTot,[]);
            subtightplot(1,2,2,opt{:}); imshow(allMask,[]);
        end
        fprintf('\nGenerating neuropil ROI...');
        for r = 1:Nroi
            tempSubNP = imdilate(ROI(r).mask_2d, strel('cube',params.neuropilRad(2))) & ~imdilate(ROI(r).mask_2d,strel('cube',params.neuropilRad(1))); % imshow( max(tempSubNP,[],3), [] );
            tempNPind = find(tempSubNP);
            [~,tempNPbad, ~] = intersect(tempNPind, allPreInd );
            tempNPind(tempNPbad) = []; 
            ROI(r).neuropil = tempNPind;
            % Construct labeled projection image
            ROI(r).labelVol = zeros(expt.Nrow, expt.Ncol, expt.Nplane, 'uint8');
            ROI(r).labelVol(ROI(r).neuropil) = 1; ROI(r).labelVol(ROI(r).ind) = 2; % saveastiff(ROI(r).labelVol, 'D:\2photon\DL102\180426_DL102\ROI1_label.tif')
        end
    else
        allPreCorr = cat(4,preROI.correlation); %histogram(allPreCorr(:) )
        corrThresh = prctile(allPreCorr(:), params.minCorrPrct); 
        allPreCorrMax = max(allPreCorr,[],4);
        allPreMask3D = falseVol;  allPreMask3D(allPreCorrMax > corrThresh) = true;
        %saveastiff( uint16(allPreMask3D), 'D:\2photon\DL75\170525\allPreMask3D.tif');
        %close all; figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
        % Get intermediate 3D masks by breaking preROI (PCA/ICA results) up into blobs and thresholding on max correlation
        c = 0;
        for p = 1:numel(preROI)
            tempConnComp = bwconncomp( preROI(p).correlation > corrThresh, 18 );
            tempSubProps = regionprops3( tempConnComp, preROI(p).correlation, 'Volume','MaxIntensity','VoxelIdxList'); % , preROI(p).correlation
            goodSub = find(tempSubProps.Volume > params.minVol & tempSubProps.MaxIntensity > corrThresh)';
            preROI(p).Nsub = numel(goodSub);
            for s = goodSub
                c = c + 1;
                % Retain some info from preliminary ROI
                intROI(c).mask_2d = preROI(p).mask_2d;
                intROI(c).ind = tempSubProps.VoxelIdxList{s};
                intROI(c).Nvox = tempSubProps.Volume(s);
                intROI(c).correlation = single(nanVol);
                intROI(c).correlation(intROI(c).ind) = single(preROI(p).correlation(intROI(c).ind));
                intROI(c).preROI = p; %keep track of which mask this ROI is derived from
                %{
                sp(1) = subplot(1,2,1);
                imagesc( max(preROI(p).correlation,[],3) ); axis image;
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

        % Exclude voxels that overlap multiple intermediate ROI, or merge ROI with high degree of overlap
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
        % Exclude overlapping voxels, or merge ROI
        for r = 1:Nint
            for R = find(Noverlap(r,:) > 0)
                fprintf('\n[r,R] = [%i, %i]: overlaps fractions = [%2.1f, %2.1f] pct (%i voxels)... ', r, R, 100*overlapFrac(r,R), 100*overlapFrac(R,r), Noverlap(r,R) );
                if overlapFrac(r,R) < params.overlap && overlapFrac(R,r) < params.overlap % exclude from both ROI
                    fprintf('Excluding overlap from both ROI');
                    intROI(r).ind = setdiff(intROI(r).ind, overlapInd{r,R});
                    intROI(R).ind = setdiff(intROI(R).ind, overlapInd{r,R});
                else % overlapFrac(r,R) >= params.overlap && overlapFrac(R,r) >= params.overlap
                    % Absorb the ROI with higher overlapFrac into the other one
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
                ROI(k).preROI = intROI(r).preROI; 
                % Update mask/indices/correlation
                ROI(k).mask_3d = falseVol;
                ROI(k).mask_3d(tempConnComp.PixelIdxList{s}) = true; % imshow( max(ROI(k).mask_3d,[],3), [])
                ROI(k).ind = find(ROI(k).mask_3d);
                ROI(k).Nvox = numel(ROI(k).ind);
                ROI(k).correlation = intROI(r).correlation;
                ROI(k).correlation(~ROI(k).mask_3d) = NaN; % imshow( max(ROI(k).correlation,[],3), [])               
                % Get additional info about the 3D ROI
                ROI(k).box_x = [tempProps3D.BoundingBox(s,1), tempProps3D.BoundingBox(s,1)+tempProps3D.BoundingBox(s,4)];%tempSubRegion.BoundingBox(s,:)
                ROI(k).box_y = [tempProps3D.BoundingBox(s,2), tempProps3D.BoundingBox(s,2)+tempProps3D.BoundingBox(s,5)];
                ROI(k).box_z = [tempProps3D.BoundingBox(s,3), tempProps3D.BoundingBox(s,3)+tempProps3D.BoundingBox(s,6)];
                ROI(k).cent = tempProps3D.Centroid(s,:);
                ROI(k).solidity = tempProps3D.Solidity(s); %how solid was this volume before filling?
                ROI(k).corrMax = tempProps3D.MaxIntensity(s);
                ROI(k).corrMean = tempProps3D.MeanIntensity(s);
                ROI(k).PrinAxLength = tempProps3D.PrincipalAxisLength(s,:);
                % Get info about the ROI's footprint
                ROI(k).mask_2d = max(ROI(k).mask_3d, [], 3); % imshow( max(ROI(k).mask_2d,[],3), [])
                ROI(k).footprint = sum(ROI(k).mask_2d(:)); %numel(ROI(k).footprintInd); 
                ROI(k).footprintInd = find( ROI(k).mask_2d );
                [footRow, footCol] = ind2sub( [expt.Nrow,expt.Ncol], ROI(k).footprintInd );
                ROI(k).footprintXY = [footCol, footRow];
                tempFootProp = regionprops(ROI(k).mask_2d, 'Orientation', 'Eccentricity');
                ROI(k).orientation = tempFootProp.Orientation;
                ROI(k).eccentricity = tempFootProp.Eccentricity;
                %close all; figure('Units','normalized', 'OuterPosition',[0,0,1,1]); 
                %imshow(tempSubBinary, []); title( sprintf('Orientation = %2.2f,  Ecc = %2.2f', ROI(k).orientation, ROI(k).eccentricity) ); pause;
                tempFootBinary = imdilate(ROI(k).mask_2d, strel('rectangle',[3,3]) ); % imshow( tempSubBinary )
                tempFootEdge = edge(tempFootBinary);
                [ROI(k).footprintEdge(:,1), ROI(k).footprintEdge(:,2)] = ind2sub( [expt.Nrow, expt.Ncol], find( tempFootEdge ) );
                ROI(k).axon = 0;
                % How many voxels in each plane?
                ROI(k).zProfile = zeros(1, expt.Nplane);
                for z = 1:expt.Nplane
                     tempPlane = ROI(k).mask_3d(:,:,z);
                     ROI(k).zProfile(z) = sum( tempPlane(:) ); 
                end
                ROI(k).zUse = find(ROI(k).zProfile); %find( squeeze(sum(sum(ROI(k).mask_3d, 1), 2)) ); % which planes were used for this ROI?

                % Setup cropped summary projection images (see WriteROIproj)
                ROI(k).projPath = '';
                ROI(k).crop = [floor(ROI(k).box_x(1)), ceil(ROI(k).box_x(2)), floor(ROI(k).box_y(1)), ceil(ROI(k).box_y(2))]; %[]; expt.Ncol-  expt.Nrow-
                ROI(k).cropMean = []; 
                ROI(k).cropStd = [];
                ROI(k).cropMax = []; 
            end
        end
        Nroi = numel(ROI);
        
        % Sort ROI from deep to superficial
        roiCent = vertcat(ROI.cent);
        [~,rSort] = sort( roiCent(:,3), 'ascend' ); % 1 is deepest, Nplane is most superficial
        ROI = ROI(rSort);
        %VisualizeSegmentation( expt, intROI, 'overwrite',true );
        
        % Construct neuropil regions, excluding any intROI's voxels
        allPreMask = allPreMask3D; %any( cat(4, ROI(:).mask_3d), 4 ); 
        allPreInd = find( allPreMask(:) );
        zExclude = find(squeeze(~any(allPreMask,[1,2])));  % find(squeeze(~any(mask_3d_all, [1,2,3])))'; % determine which planes were excluded from analysis
        for r = 1:Nroi
            tempSubNP = imdilate(ROI(r).mask_3d, strel('cuboid',[params.neuropilRad(2),params.neuropilRad(2),2])) & ...
                ~imdilate(ROI(r).mask_3d,strel('cuboid',[params.neuropilRad(1),params.neuropilRad(1),1])); % imshow( max(tempSubNP,[],3), [] );
            tempSubNP(:,:,zExclude) = false; % suppress planes that were excluded from segmentation due to poor registration
            tempNPind = find(tempSubNP);
            [~,tempNPbad, ~] = intersect(tempNPind, allPreInd );
            tempNPind(tempNPbad) = []; 
            ROI(r).neuropil = tempNPind;
            % Construct labeled projection image
            ROI(r).labelVol = zeros(expt.Nrow, expt.Ncol, expt.Nplane, 'uint8');
            ROI(r).labelVol(ROI(r).neuropil) = 1; ROI(r).labelVol(ROI(r).ind) = 2; % 
        end
        clearvars tempSubNP tempNPind
    end
    ROI( [ROI.footprint] < params.minFoot | [ROI.Nvox] < params.minVol ) = [];

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
if isfield(ROI, 'correlation') && removeCorr
    fprintf('\nRemoving correlation field...')
    ROI = rmfield(ROI, 'correlation'); % correlation field is memory-intensive and not useful
end
toc


%{
tempCent = vertcat(ROI.cent)
find( tempCent(:,1) > 480 & tempCent(:,2) > 360)
ROI(find( tempCent(:,1) > 480 & tempCent(:,2) > 360)) = [];

runningOverlap = zeros(size(preROI(1).mask_2d));
figure;
for p = [51,53,56]% 1:numel(preROI)
    subplot(2,1,1);
    imshow(preROI(p).mask_2d)
    title(sprintf('p = %i',p)); 
    
    runningOverlap = runningOverlap + preROI(p).mask_2d;
    subplot(2,1,2);
    imshow(runningOverlap, [])
    impixelinfo; pause;
end

%}