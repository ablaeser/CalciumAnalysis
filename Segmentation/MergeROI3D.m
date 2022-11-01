function [axon, expt, stillCosSim, stillCorr, rSilent] = MergeROI3D(expt, T, loco, ROI, fluor, stillEventRaster, varargin) % sbxPath,  , axonLabelFoot, axonLabelVol , stillEpoch
% Use the hierarchical clustering procedure outlined in Liang,... Andermann 2018 to find ROI sets that likely belong to the same axon
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'T', @iscell )
addRequired( IP, 'loco', @isstruct )
addRequired( IP, 'ROI', @isstruct )
addRequired( IP, 'fluor', @isstruct )
addRequired( IP, 'stillEventRaster', @islogical )
%addRequired( IP, 'stillEpoch', @isstruct )
addParameter( IP, 'method', '', @ischar ) % 'cluster'
addParameter( IP, 'zThresh', 2.5, @isnumeric)
addParameter( IP, 'mergeThresh', 0.15, @isnumeric)
addParameter( IP, 'show', false, @islogical)
%addParameter( IP, 'write', true, @islogical)
parse( IP, expt, T, loco, ROI, fluor, stillEventRaster, varargin{:} ); % , stillEpochstillEvent, 
mergeMethod = IP.Results.method;
show = IP.Results.show;
%writeTif = IP.Results.write;
zThresh = IP.Results.zThresh;
mergeThresh = IP.Results.mergeThresh;

%{
% Concatenate all periods of stillness (pre-CSD)
if ~isnan(expt.csd)
    tempFluor = [fluor(1:expt.csd-1).z];
    tempLoco = loco(1:expt.csd-1);
else
    tempFluor = [fluor.z];
    tempLoco = loco;
end
tempFluorCat = vertcat( tempFluor.ROI ); 
tempStateCat = vertcat( tempLoco.stateDown  );
fluorStillPre = tempFluorCat(tempStateCat==1,:);
zStillPre = normalize(fluorStillPre, 1); %zscore(fluorStillPre, 1);

% Detect 2.5 sigma events for each ROI, and calculate correlation with all other ROI during those events
threshStillPre = zStillPre > zThresh;
stillEventRaster = imdilate(threshStillPre, strel('rectangle',[3,1]) );
%}

% Concatenate fluor data
tempFluor = [fluor.z];
fluorCat = vertcat( tempFluor.ROI ); 
Tcat = vertcat(T{:});

% Calculate asymmetric, thresholded correlations during ROI-specific still epoch events
%Nstill = numel(stillEpoch);

%{
epochRasters = cell(1,Nstill);
for s = 1:Nstill
    epochRaster{s} = false(stillEpoch(s).Nscan, expt.Nroi);
    for roi = 1:expt.Nroi %find(stillSumm.fluor.Nevent(s,:) > 0)
        for e = 1:stillEvent(s,roi).Nevent
            epochRaster{s}(stillEvent(s,roi).scan{e},roi) = true;
        end
    end
end
catRaster = cat(1, epochRaster{:});
rSilent = find(sum(catRaster,1) == 0);
%}

stillCorr = zeros(expt.Nroi); rSilent = [];
for roi = 1:expt.Nroi
    tempEventScans = find(stillEventRaster(:,roi)); %find(catRaster) %find( stillEventRaster(:,roi) );
    if ~isempty(tempEventScans)
        stillCorr(:,roi) = corr( fluorCat(tempEventScans,roi), fluorCat(tempEventScans,:), 'rows','complete' );
        stillCorr(roi,roi) = NaN;
        tempCorrThresh = min(0.7, mean(stillCorr(:,roi),'omitnan')+zThresh*std(stillCorr(:,roi),'omitnan'));
        stillCorr(stillCorr(:,roi) < tempCorrThresh,roi) = 0;
    else
        fprintf('\nNo still event scans for ROI %i', roi );
        rSilent = [rSilent, roi];
    end
end
stillCorr(isnan(stillCorr)) = 1;
stillCosDist = squareform( pdist( stillCorr, 'cosine' ) );
stillCosDist(isnan(stillCosDist)) = 1; % otherwise, NaN values from eventless ROIs mess things up
stillLinkage = linkage(stillCosDist, 'weighted'); % , 'correlation'
stillCosSim = 1 - stillCosDist;

% Use chosen method to merge ROIs
if strcmpi(mergeMethod, 'cluster')
    fprintf('\nClustering ROIs into putative afferents using cosine dissimilarity')
    axonClust = cluster(stillLinkage, 'cutoff',mergeThresh, 'criterion','distance');
    % {
    testThresh = 0.1:0.1:3.5; % 0.01,
    Ntest = numel(testThresh); testClust = cell(1,Ntest); 
    for t = 1:Ntest
        testClust{t} = cluster(stillLinkage, 'cutoff',testThresh(t), 'criterion','distance');
    end
    NtestAxon = cellfun(@numel, cellfun(@unique, testClust, 'UniformOutput',false) );
    %plot(testThresh, NtestAxon);
    %xlabel('Threshold'); ylabel('# of putative axons');
    %}
    expt.Naxon = max(axonClust);
    fprintf('\n%i ROI -> %i axons\n', expt.Nroi, expt.Naxon) 
    for a = 1:expt.Naxon
        axon(a).ROI = find(axonClust == a)';
        axon(a).Nroi = numel(axon(a).ROI);
    end
else
    % Use PCA/ICA component from which ROI was originally derived
    fprintf('\nMerging ROIs using PCA/ICA component from which ROI was originally derived')
    preROI = unique( [ROI.preROI] );
    expt.Naxon = numel(preROI);
    fprintf('\n%i ROI -> %i axons\n', expt.Nroi, expt.Naxon) 
    for a = 1:expt.Naxon
        axon(a).ROI = find([ROI.preROI] == preROI(a));
        axon(a).Nroi = numel(axon(a).ROI);
    end
end
%[~,aSort] = sort([axon.Nroi], 'descend' );
%axon = axon(aSort);

% Gather other useful information about the axon
axonLabelFoot = zeros(expt.Nrow, expt.Ncol); % axonLabelVol = zeros(expt.Nrow, expt.Ncol, expt.Nplane);
for a = 1:expt.Naxon
    axon(a).Nvox = sum([ROI(axon(a).ROI).Nvox]);
    axon(a).ind = vertcat( ROI(axon(a).ROI).ind );
    axon(a).neuropil = vertcat( ROI(axon(a).ROI).neuropil );
    if ~isempty( intersect(axon(a).ROI, rSilent) ) 
        axon(a).silent = true;
    else
        axon(a).silent = false;
    end
    tempCrop = vertcat( ROI(axon(a).ROI).crop );
    axon(a).crop = [min(tempCrop(:,1)), max(tempCrop(:,2)), min(tempCrop(:,3)), max(tempCrop(:,4))];
    axon(a).width = axon(a).crop(2) - axon(a).crop(1); 
    axon(a).height = axon(a).crop(4) - axon(a).crop(3);
    axon(a).zUse = unique( vertcat(ROI(axon(a).ROI).zUse) ); % horzcat(ROI(axon(a).ROI).zUse)
    axon(a).similarity = stillCosSim(axon(a).ROI, axon(a).ROI); %nan(axon(a).Nroi); %
    axon(a).correlation = stillCorr(axon(a).ROI, axon(a).ROI); %nan(axon(a).Nroi); %
    tempCorr = axon(a).correlation(:); tempCorr(tempCorr == 1) = [];
    axon(a).medCorr = median(tempCorr);
    axon(a).labelFoot = zeros(expt.Nrow, expt.Ncol);
    axon(a).labelVol = zeros(expt.Nrow, expt.Ncol, expt.Nplane);
    for roi = axon(a).ROI
        axonLabelFoot(ROI(roi).footprintInd) = a;
        %axonLabelVol(ROI(roi).ind) = a;
        axon(a).labelFoot(ROI(roi).footprintInd) = roi;
        axon(a).labelVol(ROI(roi).ind) = roi;
    end
    if axon(a).Nroi > 1
        tempHull = bwconvhull(axon(a).labelFoot );
        %subplot(1,2,1); imshow(axon(a).labelFoot); title('Axon');
        %subplot(1,2,2); imshow(tempHull ); title('Convex Hull'); pause;
        tempHullProps = regionprops(tempHull, 'Orientation', 'Eccentricity', 'MajoraxisLength', 'MinorAxisLength'); % 
        axon(a).orientation = tempHullProps.Orientation;
        axon(a).eccentricity = tempHullProps.Eccentricity;
        axon(a).MajorAxisLength = tempHullProps.MajorAxisLength;
        axon(a).MinorAxisLength = tempHullProps.MinorAxisLength;
        %axon(a).PrinAxLength = [tempHullProps.MajorAxisLength, tempHullProps.MinorAxisLength];
    else
        axon(a).orientation = ROI(axon(a).ROI).orientation; %median( [] );
        axon(a).eccentricity = ROI(axon(a).ROI).eccentricity;
        axon(a).MajorAxisLength = ROI(axon(a).ROI).PrinAxLength(1);
        axon(a).MinorAxisLength = ROI(axon(a).ROI).PrinAxLength(2);
        %axon(a).PrinAxLength = ROI(axon(a).ROI).PrinAxLength;%(1:2)
    end
end
[~,aSort] = sort([axon.medCorr], 'descend', 'MissingPlacement','last');
axon = axon(aSort);

% Figure illustrating the process
if show
    % Big picture
    close all;
    opt = {[0.04,0.03], [0.07,0.04], [0.05,0.03]}; % {[vert, horz], [bottom, top], [left, right] }
    figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
    sp(1) = subtightplot(2,3,1,opt{:});
    imagesc(fluorCat' ); colorbar;
    caxis([-1,3]);
    set(gca,'Xtick',[]);
    ylabel('ROI'); title('Fluor (z score)'); 
    sp(2) = subtightplot(2,3,4,opt{:});
    imagesc(stillEventRaster'); colorbar;
    xlabel('Scan'); ylabel('ROI'); title('Events');  
    linkaxes(sp,'xy');

    subtightplot(2,3, 2,opt{:})
    imagesc(stillCorr); axis square; 
    caxis([-1,1]);
    colormap(bluewhitered);
    colorbar;
    set(gca,'Xtick',[]);
    
    ylabel('ROI'); title('Thresholded, spontaneous correlation');

    subtightplot(2,3, 5,opt{:})
    imagesc(stillCosSim); 
    axis square; 
    caxis([-1,1]);
    colormap(bluewhitered);
    colorbar; 
    xlabel('ROI'); ylabel('ROI'); title('Cosine similarity');
    impixelinfo;

    subtightplot(2,3,3,opt{:});
    dendrogram( stillLinkage, expt.Nroi );
    if strcmpi(mergeMethod, 'cluster')
        for t = 1:Ntest
            text(expt.Nroi, testThresh(t), sprintf('%1.2f->%i',testThresh(t),NtestAxon(t)), 'fontSize',8, 'VerticalAlignment','middle','HorizontalAlignment','left');% pause; sprintf('n=%i',NtestAxon(t))
        end
        line([0,expt.Nroi], mergeThresh*[1,1], 'color','r','lineStyle','--'); 
    end
    set(gca,'FontSize',8); xtickangle(30);
    xlabel('ROI', 'FontSize',12); ylabel('Linkage');
    title(sprintf('%i ROI -> %i axons', expt.Nroi, expt.Naxon ) );

    subtightplot(2,3,6,opt{:});
    imshow(label2rgb(axonLabelFoot)); hold on;
    for a = 1:expt.Naxon
        for roi = axon(a).ROI
            text( ROI(roi).cent(1), ROI(roi).cent(2), sprintf('%i-%i',roi,a), 'FontSize',6, 'HorizontalAlignment','center' );
        end 
    end
    impixelinfo;
    pause;
    % Individual axons
    clf;
    for a = 1:expt.Naxon
        subtightplot(1,3,1,opt{:}); cla;
        tempFull = fluorCat(:,axon(a).ROI);
        tempFullSpread = CalSpread( {1:size(tempFull,1)'}, {tempFull}, 'show',false, 'sepPct',90);
        plot( tempFullSpread{1} ); 
        axis tight;
        set(gca,'Ytick',[]);
        legend( sprintfc('%i', axon(a).ROI) )
        title('Fluorescence');

        subtightplot(1,3,2,opt{:}); cla;
        imagesc(axon(a).correlation); colorbar; %caxis;
        caxis([-1,1]);
        colormap(bluewhitered);
        axis square;
        set(gca, 'Xtick',1:axon(a).Nroi, 'Ytick',1:axon(a).Nroi, 'XtickLabel',sprintfc('%i', axon(a).ROI), 'YtickLabel',sprintfc('%i', axon(a).ROI) )
        title( sprintf('Still Event Correlation: median = %2.2f', axon(a).medCorr));
        
        subtightplot(1,3,3,opt{:}); cla;
        imshow( label2rgb(axon(a).labelFoot ), [] );
        for roi = axon(a).ROI
            text( ROI(roi).cent(1), ROI(roi).cent(2), sprintf('%i',roi), 'FontSize',8, 'HorizontalAlignment','center' );
        end 
        title( sprintf('Axon %i / %i', a, expt.Naxon) );
        impixelinfo;
        pause%(1);
    end
end
%{
if writeTif
    %fprintf('\nWriting summary tifs');
    VisualizeSegmentation( expt, ROI, 'overwrite',true); % axon, 
    
    RGBOpt = struct('overwrite',true, 'message',false, 'append',false, 'big',true, 'color',true );
    
    axonFootPath = strcat(expt.dir, expt.name,'_axonFoot.tif');
    axonCmap = distinguishable_colors(expt.Naxon); 
    if expt.Nplane > 1
        axonVolPath = strcat(expt.dir, expt.name,'_axonVol.tif');
        axonRGBvol = label2rgb3d(axonLabelVol,axonCmap,'w');
        axonRGBvol = permute(axonRGBvol, [1,2,4,3]);
        fprintf('\nWriting %s\n', axonVolPath);
        saveastiff( uint16(axonRGBvol), axonVolPath, RGBOpt );
        
        axonRGBfoot = cat(3, axon.labelFoot);%max( axonRGBvol,[], 4 );
        axonRGBfoot = label2rgb(max( axonRGBfoot, [], 3)); % , 'w' , axonCmap
    else
        axonRGBfoot = label2rgb(axonLabelFoot,axonCmap,'w');
    end

    fprintf('\nWriting %s\n', axonFootPath);
    saveastiff( uint16(axonRGBfoot), axonFootPath, RGBOpt );

end
%}
end