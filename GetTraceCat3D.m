function [fluor, ROI, expt] = GetTraceCat3D(expt, fluor, deform, T, varargin) % , run exptDir, exptName %ROIlabel, 
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'fluor', @isstruct )
addRequired( IP, 'deform', @isstruct )
addRequired( IP, 'T', @iscell )
addParameter( IP, 'corrPrct', 68, @isnumeric) 
addParameter( IP, 'minVol', 0, @isnumeric) 
addParameter( IP, 'minFoot', 150, @isnumeric) 
addParameter( IP, 'window', 101, @isnumeric) %rolling percentile window
addParameter( IP, 'basePrct', 10, @isnumeric) %rolling percentile value
addParameter( IP, 'overwrite', false, @islogical)
parse( IP, expt, fluor, deform, T, varargin{:} );  %  exptDir, exptName
windowSize = IP.Results.window;
if mod(windowSize,2) ==0,  windowSize = windowSize+1; end
minCorrPrct = IP.Results.corrPrct; 
minVol = IP.Results.minVol; 
minFoot = IP.Results.minFoot; 
basePrct = IP.Results.basePrct; 
overwrite = IP.Results.overwrite;
roiDir = strcat(expt.dir, 'ROI\'); mkdir(roiDir);
savePath = sprintf('%s%s_ROI.mat', expt.dir, expt.name ); 
tic;
if ~exist(savePath,'file') || overwrite
    % Get and process fluorescence signals
    seg3Dpath = sprintf('%s%s_seg_data_final.mat', expt.dir, expt.name);  % tracePath = sprintf('%s%s_seg_traces.mat', expt.dir, expt.name);   
    tifPath = [expt.dir,expt.name,'_affineProj.tif'];
    fprintf('\nLoading %s  and  %s... ', seg3Dpath,  tifPath); 
    load(seg3Dpath); %load( tracePath );
    projStack = loadtiff(tifPath);  toc % Load projection stack
    preROI = seg_data([seg_data.keep]); %seg_data_keep;
    nanVol = nan( expt.Nx, expt.Ny, expt.Nplane );
    falseVol = false( expt.Nx, expt.Ny, expt.Nplane );
    c = 0;
    for r = 1:numel(preROI)
        % Get preliminary 3D mask from thresholding correlation
        threshold = max(0, prctile( preROI(r).correlation(:), minCorrPrct, 'all') );    %histogram( preROI(r).correlation(:))
        temp3Dmask = preROI(r).correlation > threshold;
        tempConnComp = bwconncomp( temp3Dmask, 26 );
        tempSubProps = regionprops3( tempConnComp, preROI(r).correlation, 'Volume', 'Centroid', 'BoundingBox', 'MeanIntensity', 'MaxIntensity', 'VoxelList', 'VoxelIdxList','PrincipalAxisLength' ); % 
        [volSort, volSortInd] = sort( tempSubProps.Volume, 'descend' );
        preROI(r).Nsub = numel( find(volSort >= minVol) );
        if preROI(r).Nsub > 0
            for s = 1:preROI(r).Nsub
                c = c + 1;
                % Retain some info from preliminary ROI
                ROI(c).mask_2d = preROI(r).mask_2d;
                ROI(c).mean_trace_raw = preROI(r).mean_trace_raw;
                ROI(c).mean_trace_hp = preROI(r).mean_trace_hp;
                %
                ROI(c).ind = tempConnComp.PixelIdxList{volSortInd(s)};
                ROI(c).mask_3d = falseVol; % 
                ROI(c).mask_3d(ROI(c).ind) = true;
                ROI(c).correlation = nanVol;
                ROI(c).correlation(ROI(c).ind) = preROI(r).correlation(ROI(c).ind);
                %ROI(c).xproj
                %ROI(c).yproj
                %ROI(c).zproj
                                
                % Get info about the subregion
                ROI(c).Nvox = tempSubProps.Volume(volSortInd(s));
                ROI(c).Nsub = 1; % vestigial
                ROI(c).box_x = [tempSubProps.BoundingBox(volSortInd(s),1), tempSubProps.BoundingBox(volSortInd(s),1)+tempSubProps.BoundingBox(volSortInd(s),4)];%tempSubRegion.BoundingBox(volSortInd(s),:)
                ROI(c).box_y = [tempSubProps.BoundingBox(volSortInd(s),2), tempSubProps.BoundingBox(volSortInd(s),2)+tempSubProps.BoundingBox(volSortInd(s),5)];
                ROI(c).box_z = [tempSubProps.BoundingBox(volSortInd(s),3), tempSubProps.BoundingBox(volSortInd(s),3)+tempSubProps.BoundingBox(volSortInd(s),6)];
                ROI(c).cent = tempSubProps.Centroid(volSortInd(s),:);
                ROI(c).corrMax = tempSubProps.MaxIntensity(volSortInd(s));
                ROI(c).corrMean = tempSubProps.MeanIntensity(volSortInd(s));
                ROI(c).PrinAxLength = tempSubProps.PrincipalAxisLength(volSortInd(s),:);
                ROI(c).footprintInd = unique( sub2ind([expt.Nx, expt.Ny], tempSubProps.VoxelList{volSortInd(s)}(:,2), tempSubProps.VoxelList{volSortInd(s)}(:,1)) );
                ROI(c).footprint = numel(ROI(c).footprintInd); 
                [tempRow, tempCol] = ind2sub( [expt.Nx,expt.Ny], ROI(c).footprintInd );
                ROI(c).footprintXY = [tempRow, tempCol];
                tempSubBinary = false( expt.Nx, expt.Ny ); tempSubBinary(ROI(c).footprintInd) = true; 
                tempSubBinary = imdilate(tempSubBinary, strel('rectangle',[3,3]) ); % imshow( tempSubBinary )
                tempSubEdge = edge(tempSubBinary);
                [ROI(c).edge(:,1), ROI(c).edge(:,2)] = ind2sub( [expt.Nx, expt.Ny], find( tempSubEdge ) );
                ROI(c).crop = [floor(ROI(c).box_x(1)), expt.Ny-floor(ROI(c).box_x(2)), floor(ROI(c).box_y(1)), expt.Nx-floor(ROI(c).box_y(2))];
            
                % How many voxels in each plane?
                ROI(c).zProfile = zeros(1, expt.Nplane);
                for z = 1:expt.Nplane
                     tempPlane = ROI(c).mask_3d(:,:,z);
                     ROI(c).zProfile(z) = sum( tempPlane(:) ); 
                end
                ROI(c).zUse = find(ROI(c).zProfile); %find( squeeze(sum(sum(ROI(c).mask_3d, 1), 2)) ); % which planes were used for this ROI?

                ROI(c).maxProj = max(projStack(:,:,ROI(c).zUse),[],3); 
                %close all;
                %imshow( imadjust(ROI(c).maxProj), [] ); hold on; % 
                %plot( ROI(c).edge(:,2), ROI(c).edge(:,1), '.b' ); hold on;

                % Get cropped summary projection images
                cropProj = projStack(ROI(c).crop(3)+1:end-ROI(c).crop(4), ROI(c).crop(1)+1:end-ROI(c).crop(2), ROI(c).zUse);
                ROI(c).max = max( max( cropProj, [], 3), [], 4 ); % imshow( ROI(c).max, [] )
                ROI(c).mean = mean( mean( cropProj, 3), 4 ); % imshow( ROI(c).mean, [] )
            end
            %ecdf( [ROI(c).sub.footprint] )
        end
    end
    % Remove ROI with minimal footprint
    rFoot = find([ROI.footprint]<minFoot);
    ROI([ROI.footprint]<minFoot) = [];
    fprintf('Erased %i ROI (footprint)', numel(rFoot));
    
    % Combine data across masks
    expt.Nroi = numel( ROI );
    allMask = any( cat(4, ROI(:).mask_3d), 4 ); % Combined mask from all surviving ROI
    allInd = find( allMask(:) );
    zExclude = find(squeeze(~any(allMask,[1,2])));  % find(squeeze(~any(mask_3d_all, [1,2,3])))'; % determine which planes were excluded from analysis
    
    % Construct neuropil regions, excluding any ROI's voxels
    w = waitbar(0,'Generating neuropil ROI...');
    for r = 1:expt.Nroi
        tempSubNP = imdilate(ROI(r).mask_3d,strel('cuboid',[15,15,2])) & ~imdilate(ROI(r).mask_3d,strel('cuboid',[8,8,1])); % imshow( max(tempSubNP,[],3), [] );
        tempSubNP(:,:,zExclude) = false; % suppress planes that were excluded from segmentation due to poor registration
        tempNPind = find(tempSubNP);
        [~,tempNPbad, ~] = intersect(tempNPind, allInd );
        tempNPind(tempNPbad) = []; 
        ROI(r).np = tempNPind;
        % Construct labeled projection image
        labelMat = zeros(expt.Nx, expt.Ny, expt.Nplane, 'uint8');
        labelMat(ROI(r).np) = 1; labelMat(ROI(r).ind) = 2;
        ROI(r).label = max(labelMat, [] , 3);   
        ROI(r).labelCrop = ROI(r).label( ROI(r).crop(3)+1:end-ROI(r).crop(4), ROI(r).crop(1)+1:end-ROI(r).crop(2) );
        waitbar(r/expt.Nroi);
    end
    delete(w);
    
    fprintf('\nExtracting traces');
    Froi = nan(expt.totScan, expt.Nroi); Fnp = nan(expt.totScan, expt.Nroi); Fsub = cell(1, expt.Nroi); FsubNP = cell(1, expt.Nroi);
    for r = 1:expt.Nroi, Fsub{r} = nan( expt.totScan, ROI(r).Nsub ); FsubNP{r} = nan( expt.totScan, ROI(r).Nsub ); end
    w = waitbar(0,'Extracting 3D traces of ROI and subROI...'); % w = 
    for s = 1:expt.totScan-1
        tempVol = double( readSbx( expt.sbx, expt.Nplane*(s-1)+1, expt.Nplane, 1, [] ) ); %double(pipe.io.sbxRead(expt.sbx, expt.Nplane*(t-1)+1, expt.Nplane, 1,[]));
        for r = 1:expt.Nroi
            Froi(s,r) = mean( tempVol(ROI(r).ind) );
            Fnp(s,r) = mean( tempVol(ROI(r).np) );
        end
        waitbar(s/(expt.totScan-1));
    end
    delete(w);
    toc
    F = Froi - Fnp + mean(Fnp,1,'omitnan'); % plot( F(:,1) );

    % Censor scans where each ROI's planes had a bad deformation value
    transAPcat = cat(1, deform.transAP);
    for r = 1:expt.Nroi
        nan_scan = find( isnan(sum(transAPcat(:,ROI(r).zUse),2) ) )';
        F(nan_scan,r) = NaN;
        fprintf('\nROI %i: censored %i scans with bad deformation values', r, numel(nan_scan) );
    end
    % Use the censored, neuropil-subtracted signal to normalize
    tic
    Fo = MovingPercentile(F, basePrct, windowSize, 'pre'); % movprctile(F, basePrct, windowSize, 1); %plot([F(:,1), Fbase(:,1)] );
    dFF = (F-Fo)./Fo;
    % Deconvolve dF/F
    activity = nan(size(dFF,1),size(dFF,2));
    %{
    fprintf('\nDeconvolving censored dF/F');
    rBad = [];
    for r = 1:expt.Nroi
        try
            goodFrames = find(~isnan(dFF(:,r)))';
            dFFcens = dFF(goodFrames,r);
            [~, actCens, ~] = deconvolveCa( dFFcens, 'exp2', 'foopsi', 'optimize_pars', 'optimize_smin', 'maxIter',20, 'window',300 );  %ar1
            activity(goodFrames,r) = actCens;
        catch
            fprintf('\nROI %i:  Deconvolution failed', r);
            rBad = [rBad, r];
        end
    end
    toc
    % Final quality control - eliminate ROIs with failed deconvolution
    for r = 1:expt.Nroi
        if any(abs(dFF(:,r))> 20) 
            fprintf('\nr = %i: bad dF/F',r);
            rBad = [rBad, r];
        elseif median( F(:,r), 'omitnan' ) < 0
            rBad = [rBad, r];
            fprintf('\nr = %i: negative F',r);
        elseif sum(activity(:,r), 'omitnan') == 0%all(activity(:,r) == 0 ) 
            fprintf('\nr = %i: zero activity',r);
            rBad = [rBad, r];
        end
    end
    if ~isempty(rBad)
        fprintf('\nEliminating %i bad ROI', numel(rBad));
        ROI(rBad) = [];
        Froi(:,rBad) = []; Fnp(:,rBad) = []; 
        F(:,rBad) = []; Fo(:,rBad) = []; dFF(:,rBad) = []; activity(:,rBad) = []; %Zdff(:,rBad) = [];
        Fsub(rBad) = []; FsubNP(rBad) = [];
        expt.Nroi = expt.Nroi - numel(rBad);
    end
    %}

    % Find ROIs' brightest scans and use them to make a projection image
    %figure('WindowState','maximized');
    % {
    Npeaks = 30;
    for r = 1:expt.Nroi
        [~,brightScans] = sort( dFF(:,r), 'descend', 'MissingPlacement','last'); brightScans = brightScans(1:Npeaks);
        for p = flip(1:Npeaks)
            peakVol(:,:,:,p) = pipe.io.sbxRead(expt.sbx, expt.Nplane*(brightScans(p)-1)+1, expt.Nplane, 1, []);
        end
        meanPeakVol = mean(peakVol, 4);
        ROI(r).maxProj = max( meanPeakVol(:,:,[ceil(ROI(r).box_z(1)), floor( ROI(r).box_z(2) )]), [], 3); 
        %cropProj = max(meanPeakVol(ROI(r).box_y(1):ROI(r).box_y(2), ROI(r).box_x(1):ROI(r).box_x(2), [ceil(ROI(r).sub(1).box_z(1)), floor( ROI(r).sub(1).box_z(2) )]), [], 3); 
        %{
        cropPrctile = prctile( cropProj(:), [10, 99.5] );
        subplot(2,1,1); cla;
        plot( dFF(:,r) ); hold on;
        plot( brightScans, dFF(brightScans,r), 'ko' ); hold off;
        axis tight;
        subplot(2,1,2); cla;
        imshow( cropProj, cropPrctile ); hold on; % 
        for s = 1:ROI(r).Nsub
            plot( ROI(r).edge(:,2)-ROI(r).box_x(1)+1, ROI(r).edge(:,1)-ROI(r).box_y(1)+1, '.', 'MarkerSize',2 ); %, 'Color',tempColors(s,:)
        end
        impixelinfo;
        pause;
        %}
    end
    %}

    % Break fluor signals up by run
    fprintf('\nBreaking fluor signals up by run...  ')
    Nscan = cellfun(@numel, T);
    scanLims = [0,cumsum(Nscan)];
    Zdff = nanzscore(dFF, 1); %zscore( dFF, 0, 1 );
    for runs = 1:expt.Nruns
        fluor(runs).Froi.ROI = Froi(scanLims(runs)+1:scanLims(runs+1),:); 
        fluor(runs).Fnp.ROI = Fnp(scanLims(runs)+1:scanLims(runs+1),:);       
        fluor(runs).F.ROI = F(scanLims(runs)+1:scanLims(runs+1),:); 
        fluor(runs).Fo.ROI = Fo(scanLims(runs)+1:scanLims(runs+1),:);
        fluor(runs).dFF.ROI = dFF(scanLims(runs)+1:scanLims(runs+1),:);
        fluor(runs).z.ROI = Zdff(scanLims(runs)+1:scanLims(runs+1),:);
        fluor(runs).act.ROI = activity(scanLims(runs)+1:scanLims(runs+1),:);
        % Get subROI raw data
        for roi = 1:expt.Nroi
            fluor(runs).Froi.sub{roi} = Fsub{roi}(scanLims(runs)+1:scanLims(runs+1),:);
            fluor(runs).Fnp.sub{roi} = FsubNP{roi}(scanLims(runs)+1:scanLims(runs+1),:);
        end
        % Correlation between ROIs for this run
        fluor(runs).Froi.corr = corr( fluor(runs).Froi.ROI, 'Rows','complete' );
        fluor(runs).Fnp.corr = corr( fluor(runs).Fnp.ROI, 'Rows','complete' );
        fluor(runs).F.corr = corr( fluor(runs).F.ROI, 'Rows','complete' );
        fluor(runs).Fo.corr = corr( fluor(runs).Fo.ROI, 'Rows','complete' );
        fluor(runs).dFF.corr = corr( fluor(runs).dFF.ROI, 'Rows','complete' );
        fluor(runs).act.corr = corr( fluor(runs).act.ROI, 'Rows','complete' );
        %{
        for r = 1:expt.Nroi
            subplot(3,2,[1,3,5]);
            imshow( ROI(r).mask_2d, [] ); hold on;
            title( sprintf('ROI %i',r) );
            impixelinfo
            sp(1) = subplot(3,2,2); cla;
            plot( fluor(runs).F.ROI(:,r) ); hold on;
            plot( fluor(runs).Fo.ROI(:,r) ); ylabel('F');
            title(sprintf('[run, ROI] = [%i, %i]', runs, r));
            sp(2) = subplot(3,2,4); cla;
            plot( fluor(runs).dFF.ROI(:,r) ); ylabel('dF/F');
            sp(3) = subplot(3,2,6); cla;
            plot( fluor(runs).act.ROI(:,r) ); ylabel('Activity');
            linkaxes(sp,'x'); xlim([-Inf,Inf]);
            pause;
        end
        %}
    end
    toc
    %ViewROI3D(expt, ROI, fluor, 'save',false)
    % Save the results
    fprintf('\nSaving %s', savePath);
    save(savePath, '-v7.3'); % ,'ROIlabel' , 'fluor','ROI','expt','actCorr','Froi','Fnp','dFF','F','Fo','activity','Zdff','IP'
else
    fprintf('\nLoading %s...  ', savePath);
    load(savePath);  %, 'fluor', 'ROI', 'expt', 'actCorr'
end
toc
