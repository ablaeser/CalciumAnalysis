%% check optotune limits
c = 0;
for x = xPresent
    c = c + 1;
    line( c*[1,1], runInfo{x}(1).otparam(1:2), 'color', 'k' ); hold on;
    %plot( c, runInfo{x}(1).otparam(1:2) ); 
    plot( c, runInfo{x}(1).otparam(1:2), '.' );
end
xlim([0,c+1]);
set(gca, 'Xtick',1:c, 'XtickLabel',{expt(xPresent).name}, 'TickLabelInterpreter','none', 'TickDir','out'); xtickangle(30);

%% z position of ROI
close all;
figure;
zROI = cell(1,Nexpt);
for x = xPresent
    zROI{x} = vertcat(ROI{x}.cent); 
    zROI{x} = zROI{x}(:,3);
    zDensity{x} = histcounts( zROI{x}, 1:expt(x).Nplane, 'Normalization','probability' );
    barh((1:expt(x).Nplane-1)+0.5, zDensity{x});
    ylabel('< Pial    Optotune Plane     Dural >');
    xlabel('Fraction of ROI');
    title( expt(x).name, 'Interpreter','none' );
    %pause;
    
end

%% Heatmaps/rasters sorted by layer: dural (high z) to pial (low z)
for x = xPresent
    %fluorZ = fluor{x}(expt(x).preRuns).z;
    %fluorZcat = vertcat( fluorZ.ROI );
    rasterCat = cat(1, eventRaster{x}{expt(x).preRuns});
    
    tempCent = vertcat( ROI{x}.cent );
    [~,rSortZ] = sort( tempCent(:,3), 'descend' );
    
    imagesc( rasterCat(:,rSortZ)' )
    %imagesc( fluorZcat(:,rSortZ)' );
    title( sprintf('x = %i', x) );
    pause;
end



%% threshold crossings as a function of depth - obsolete
LW = 0.5; LS = '-'; FS = 10; FW = 'bold';
crossing = repmat(struct('thresh',3, 'N',[], 'dur',[], 'totDur',[], 'meanDur',[]), 1, Nexpt);
close all; figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
sp(2) = subplot(2,1,2); sp(1) = subplot(2,1,1); 
linkaxes(sp,'x');
for x = xPresent
    dT = (1/expt(x).scanRate);
    if ismember(x, xCSD)
        preCSDruns = 1:expt(x).csd-1;
    else
        preCSDruns = expt(x).runs;
    end
    tempFluor = [fluor{x}(preCSDruns).z];
    catZroi = vertcat( tempFluor.ROI );
    for r = flip(1:expt(x).Nroi)
        tempThreshCross = catZroi(:,r) > crossing(x).thresh;
        crossingEvents = regionprops( tempThreshCross, 'Area', 'PixelIdxList' );
        crossing(x).N(r) = numel(crossingEvents);
        crossing(x).dur{r} = dT*[crossingEvents.Area];
        crossing(x).totDur(r) = sum(crossing(x).dur{r});
        crossing(x).meanDur(r) = crossing(x).totDur(r)/crossing(x).N(r);
        %{
        plot(catZroi(:,r)); hold on;
        line([1,size(catZroi,1)], crossing(x).thresh*[1,1], 'color','r', 'lineStyle','--');
        title(sprintf('r = %i:  %i crossings, mean duration = %2.1f s', r, crossing(x).N(r), crossing(x).meanDur(r)));
        pause(1);
        cla;
        %}
    end
    
    roiCentNorm = (vertcat(ROI{x}.cent)-1)/expt(x).Nplane;
    % Does activity depend on depth?
    NcrossMdl = fitlm( roiCentNorm(:,3), crossing(x).N ); 
    NcrossMdlCoeff = NcrossMdl.Coefficients; 
    NcrossSlope(x) = NcrossMdlCoeff.Estimate(2);
    crossDurMdl = fitlm( roiCentNorm(:,3), crossing(x).meanDur ); 
    crossDurMdlCoeff = crossDurMdl.Coefficients;
    crossDurSlope(x) = crossDurMdlCoeff.Estimate(2);
    
    subplot(sp(1)); cla;
    plot(roiCentNorm(:,3), crossing(x).N, '.' ); hold on;
    plot(roiCentNorm(:,3), predict(NcrossMdl, roiCentNorm(:,3)), 'k', 'LineWidth',LW, 'lineStyle',LS);
    xlim([1,expt(x).Nplane]);
    text(0.2, 0.95, sprintf('Slope = %2.1f (p = %2.3f)', NcrossSlope(x), NcrossMdlCoeff.pValue(2)), 'Units','normalized', 'FontSize',FS, 'HorizontalAlignment','center', 'FontWeight',FW );
    ylabel( sprintf('# of threshold crossings (z > %2.1f)',  crossing(x).thresh) );
    title( expt(x).name, 'Interpreter','none' );
    box off;
    
    subplot(sp(2)); cla;
    plot(roiCentNorm(:,3), crossing(x).meanDur, '.' ); hold on;
    plot(roiCentNorm(:,3), predict(crossDurMdl, roiCentNorm(:,3)), 'k', 'LineWidth',LW, 'lineStyle',LS);
    text(0.2, 0.95, sprintf('Slope = %2.1f (p = %2.3f)', crossDurSlope(x), crossDurMdlCoeff.pValue(2)), 'Units','normalized', 'FontSize',FS, 'HorizontalAlignment','center', 'FontWeight',FW );
    xlim([0,1]); %xlim([1,expt(x).Nplane]);
    ylabel('Mean duration of crossing'); xlabel('< Pial -  Normalized Depth  - Dural  > ');
    box off;
    pause;
end
%close all;

%%  - obsolete
plot( NcrossSlope, crossDurSlope, '.' );
set(gca, 'XaxisLocation','origin', 'YaxisLocation','origin' );
xlabel('Frequency Slope'); ylabel('Duration Slope');

%% Get cross-correlation of measures of deformation across planes
for x = 
% concatenate pre-CSD deformation


