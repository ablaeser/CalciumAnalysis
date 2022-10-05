% Estimate CSD speed through neuropil
csdWave = cell(1,Nexpt); 
for x = intersect(xFiber, xSlope) %xSlope %intersect(xFiber, xCSD)
    %[csdBout{x}(expt(x).csd), csdStat{x}, csdParam(x)] = PeriCSD(runInfo{x}(expt(x).csd), Tscan{x}{expt(x).csd}, loco{x}(expt(x).csd), deform{x}(expt(x).csd), fluor{x}(expt(x).csd), defVars, 'show',false);
    try
        % Define neuropil volume
        %{
        npVol = true( expt(x).Nrow, expt(x).Ncol, expt(x).Nplane );
        npVol(:,1:segParams{x}.edges(1),:) = false;
        npVol(:,end-segParams{x}.edges(2)+1:end,:) = false;
        npVol(1:segParams{x}.edges(3),:,:) = false;
        npVol(end-segParams{x}.edges(4)+1:end,:,:) = false;
        npVol(:,:,setdiff(1:expt(x).Nplane, segParams{x}.zProj)) = false;
        for r = 1:expt(x).Nroi,  npVol( imdilate(ROI{x}(r).mask_3d,strel('cuboid',[8,8,1])) ) = false;  end
        npVolInd(x) = {find(npVol)};
        %}
        %Tshift = 0; %if x == 37, Tshift = 12.5888; else, Tshift = 0; end
        csdWave{x} = GetNeuropilWaveSpeed(expt(x), catInfo(x), csdBout{x}(expt(x).csd),  {[]}, 'window',[-6,14], 'show',true, 'save','CSD_np_redo', 'overwrite',false); %  npVolInd(x) , 'bin',20, 'Tshift',Tshift
        pause;
    catch
        fprintf('\nx = %i: CSD wave speed failed', x);
    end
end
xWave = find(~cellfun(@isempty, csdWave ));
xWave = setdiff(xWave, [8,11,14]); % these experiments have problems with CSD registration
csdWave([8,11,14]) = cell(1,3);


%% Plot CSD neuropil velocity 
close all; clearvars h;
figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
k = 0;
for x = xWave
    k = k+1;
    h(k) = polarplot( csdWave{x}.angle*[1,1], abs(csdWave{x}.speed)*[0,1], 'LineWidth',2 ); hold on;
    waveName{k} = expt(x).name;
    pause;
end
%thetaticks(0:45:315)
set(gca, 'ThetaTick',[-90,0,90,180], 'ThetaTickLabel',{'Med','Ant','Lat','Post'})
legend(h, waveName, 'Interpreter','none', 'FontSize',14 );
title('CSD Wave Velocity (um/s)', 'FontSize',14 );

%% Get fiber ROI onset of CSD
xSlope = [5,9,30,31,37,41,28];
slopeTablePath = 'D:\MATLAB\Dura\3D\LocoVsCSDslopeTable.mat';
if exist(slopeTablePath, 'file')
    fprintf('Loading %s', slopeTablePath); 
    load(slopeTablePath, 'slopeTable');
else
    for x = xSlope % , 36 30 % [5,9,30,31,37,41, 28, 36] %xWave
        fprintf('\nx = %i, %s', x, expt(x).name);
        fiber{x} = GetFiberWaveSpeed(expt(x), ROI{x}, fiber{x}, csdBout{x}(expt(x).csd), csdWave(x), 'show',false, 'minR2',0.7); % , 'CSD'
        fiber{x} = GetFiberWaveSpeed(expt(x), ROI{x}, fiber{x}, periBout{x}(expt(x).preRuns), {[]}, 'show',true, 'minEffect',0, 'minR2',0.7); %false
    end
    
    
    slopeTable = []; % [x, f, spont, csd, orientation]
    k = 0;
    for x = xSlope %intersect(xWave, xPresent) , 36]
        for f = 1:Nfiber(x)
            k = k+1;
            slopeTable(k,5) = fiber{x}(f).orientation;
            slopeTable(k,1:4) = NaN;
            slopeTable(k,1) = x;
            slopeTable(k,2) = f;
            if isfield(fiber{x}(f), 'spont') && ~isempty(fiber{x}(f).spont) && ~isnan(fiber{x}(f).spont.sepDelaySlope)
                fprintf('\n%s fiber %i: spontaneous', expt(x).name, f);
                slopeTable(k,3) = 100*fiber{x}(f).spont.sepDelaySlope;
            end
            if isfield(fiber{x}(f), 'CSD') && ~isempty(fiber{x}(f).CSD) && ~isnan(fiber{x}(f).CSD.sepDelaySlope)
                fprintf('\n%s fiber %i: CSD', expt(x).name, f);
                slopeTable(k,4) = 100*fiber{x}(f).CSD.sepDelaySlope;
            end
        end
    end
    slopeTable(:,6) = rad2deg(slopeTable(:,5));
    save(slopeTablePath, 'slopeTable');
end

SpontVsCSDspeed = figure;
JitterPlot( abs(slopeTable(:,3:4)), 0.45);
[~,pSlope] = ttest2( abs(slopeTable(:,3)), abs(slopeTable(:,4)));
title(sprintf('p = %2.5f', pSlope)); box off;
set(gca,'Xtick',[1,2], 'XtickLabel',{sprintf('Loco (n=%i)', sum(~isnan(slopeTable(:,3)))),sprintf('CSD (n=%i)', sum(~isnan(slopeTable(:,4))))});
ylabel('Delay per Separation (s/100 um)');
axis square;

figDir = 'C:\Users\ablaeser\Documents\Afferent Paper\';
figPath = strcat(figDir, 'SpontVsCSDspeed.tif'); 
%print(SpontVsCSDspeed, figPath, '-dtiff'); fprintf('\nSaved %s', figPath); 

%% Summarize each fiber's loco vs CSD results

for row = 1:size(slopeTable,1)
    if slopeTable(row,4) == 0 % slopeTable(row,3) ~= 0 && ~isnan(slopeTable(row,3)) %  ~isnan(slopeTable(row,3)) && ~isnan(slopeTable(row,4))
        X = slopeTable(row,1);
        F = slopeTable(row,2);
        % Show the fiber ROIs
        subplot(2,2,[1,3]); cla;
        imshow(fiber{X}(F).labelFoot); title(sprintf('x = %i, fiber %i. orientation = %2.2f  (%2.2f)', X, F, fiber{X}(F).orientation, rad2deg(fiber{X}(F).orientation))); 
        % Show the CSD results
        subplot(2,2,2); cla;
        plot( fiber{X}(F).CSD.sepDelay(:,1), fiber{X}(F).CSD.sepDelay(:,2), '.' ); hold on;
        plot( fiber{X}(F).CSD.sepDelay(:,1), fiber{X}(F).CSD.sepDelayPred );
        title( sprintf('CSD: %2.1f ', 100*fiber{X}(F).CSD.sepDelaySlope )); 
        % Show the loco results
        subplot(2,2,4); cla;
        plot( fiber{X}(F).spont.sepDelay(:,1), fiber{X}(F).spont.sepDelay(:,2), '.' ); hold on;
        plot( fiber{X}(F).spont.sepDelay(:,1), fiber{X}(F).spont.sepDelayPred );
        title( sprintf('Loco: %2.1f ', 100*fiber{X}(F).spont.sepDelaySlope )); 
        
        pause;
    end
end

%% Plot FiberWaveSpeed results - CSD wave onset vs ROI onset
close all; clearvars h; k = 0;
CSDonset_NPvsROI = figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
line([-1,10], [-1,10],'color','r', 'lineStyle','--'); hold on;
axis square;
set(gca,'XaxisLocation','origin','YaxisLocation','origin')
xlabel('Neuropil Onset (s)'); ylabel('ROI Onset (s)');
fiberName = cell(1,0);
fiberColor = distinguishable_colors(30);
FS = 16;
for x = intersect(xFiber, xCSD)
    try
        for f = find(~cellfun(@isempty, {fiber{x}.CSD}))
            k = k+1;
            fiberName{k} = sprintf('%s_fiber%i', expt(x).name, f);
            h(k) = plot( fiber{x}(f).CSD.binOnset, fiber{x}(f).CSD.Tonset, '.', 'MarkerSize',10, 'color',fiberColor(k,:) ); hold all;
            pause;
        end
    catch
        fprintf('%s failed', expt(x).name);
    end
end
set(gca, 'FontSize',FS)
legend(h, fiberName, 'Location','EastOutside', 'Interpreter','none', 'FontSize',FS);
figPath = 'C:\Users\ablaeser\Documents\Afferent Paper\CSDonset_NPvsROI.tif';
print(CSDonset_NPvsROI, figPath, '-dtiff'); fprintf('\nSaved %s', figPath); 

%% Plot FiberWaveSpeed results
close all;
opt = {[0.08,0.08], [0.08,0.06], [0.08,0.02]};  % {[vert, horz], [bottom, top], [left, right] }
figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
for x = 37 
    for f = 1
        subtightplot(2,3,[1,4],opt{:});
        imshow( label2rgb(fiber{x}(f).labelFoot ), [] ); %hold on;
        
        subtightplot(2,3,2,opt{:});
        imagesc( periBout{x}(1).fluor{1}(:,fiber{x}(f).ROI)' );
        title('Locomotion-Associated');
        caxis([-1,3]);
        colorbar;
        
        subtightplot(2,3,5,opt{:});
        imagesc( csdBout{x}(3).fluor{1}(csdWave{x}.scans,fiber{x}(f).ROI)' ); 
        title('Spreading Depression');
        caxis([-1,3]);
        colorbar;
        
        if isfield(fiber{x}(f), 'spont')
            subtightplot(2,3,3,opt{:});
            plot( fiber{x}(f).spont.sepDelay(:,1), fiber{x}(f).spont.sepDelay(:,2), '.' ); hold on;
            plot( fiber{x}(f).spont.sepDelay(:,1), fiber{x}(f).spont.sepDelayPred, 'r--' );
            axis square;
            xlim([0,500]); ylim([-3,3]);
            xlabel('Pairwise Separation (um)'); ylabel('Pairwise Lag (s)'); 
        end
        
        if isfield(fiber{x}(f), 'CSD')
            subtightplot(2,3,6,opt{:});
            plot( fiber{x}(f).CSD.sepDelay(:,1), fiber{x}(f).CSD.sepDelay(:,2), '.' ); hold on;
            plot( fiber{x}(f).CSD.sepDelay(:,1), fiber{x}(f).CSD.sepDelayPred, 'r--' );
            axis square;
            xlim([0,500]); ylim([-3,3]);
            xlabel('Pairwise Separation (um)'); ylabel('Pairwise Lag (s)'); 
        end
    end
end

%% Plot FiberWaveSpeed results - spont only
minLength = 200;
close all;
opt = {[0.08,0.08], [0.08,0.06], [0.08,0.02]};  % {[vert, horz], [bottom, top], [left, right] }
figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
k = 0;
exptColor = distinguishable_colors(numel(xFiber));
for x = intersect(xFiber, xPresent)
    k = k+1;
    [fiberLengthSort, fSort] = sort([fiber{x}.length], 'descend');
    fiberColor = distinguishable_colors(Nfiber(x));
    % Show that there is negligible delay between distant ROIs from the same fiber for spontaneous activity
    for f = fSort(fiberLengthSort > minLength) %1:Nfiber(x)
        roiColor = distinguishable_colors(fiber{x}(f).Nroi);
        % Show fiber ROIs
        subplot(2,3,1:3);
        roiLabel = zeros(size(fiber{x}(f).labelFoot));
        for r = 1:fiber{x}(f).Nroi
            roiLabel(ROI{x}(fiber{x}(f).ROI(r)).footprintInd) = r; %roiLabel(roiLabel==fiber{x}(f).ROI(r)) = r; 
        end
        imshow(label2rgb(roiLabel, roiColor));
        title(sprintf('%s fiber %i', expt(x).name, f), 'Interpreter','none' ); 

        if isfield(fiber{x}(f), 'spont')
            subplot(2,3,4:5); cla;
            for b = find(sum(fiber{x}(f).spont.goodFit,2)>1)' %1:numel(fiber{x}(f).spont.Twindow)
                for r = 1:fiber{x}(f).Nroi %
                    if ~isempty(fiber{x}(f).spont.trace{b,r})
                        plot( fiber{x}(f).spont.Twindow{b}, fiber{x}(f).spont.trace{b,r}, 'color',roiColor(r,:) ); hold on;
                        if fiber{x}(f).spont.goodFit(b,r)
                            plot( fiber{x}(f).spont.Twindow{b}, predict(fiber{x}(f).spont.traceFit{b,r}, fiber{x}(f).spont.Twindow{b}), '--', 'color',roiColor(r,:));
                        end
                    end
                end
                xlim([-Inf,Inf]); ylim([-Inf,Inf]);
                xlabel('Peri-Run Time'); ylabel('Normalized Fluor');
            end

            subplot(2,3,6);
            if size(fiber{x}(f).spont.sepDelay,1) > 1
                plot(fiber{x}(f).spont.sepDelay(:,1), fiber{x}(f).spont.sepDelay(:,2), '.', 'color',fiberColor(f,:)); hold on;
                if ~isempty(fiber{x}(f).spont.sepDelayFit)
                    plot(fiber{x}(f).spont.sepDelay(:,1), polyval([fiber{x}(f).spont.sepDelaySlope, fiber{x}(f).spont.sepDelayConst], fiber{x}(f).spont.sepDelay(:,1) ), ...
                        '--', 'color',fiberColor(f,:) )
                    %plot(fiber{x}(f).spont.sepDelay(:,1), predict( fiber{x}(f).spont.sepDelayFit, fiber{x}(f).spont.sepDelay(:,1)), '--', 'color',fiberColor(f,:) );
                end
            end
            hold off;
            xlabel('Pairwise Separation (um)'); ylabel('Relative Delay (s)'); title('All well-fit bouts');
            axis square; 
            pause%(0.1);
        end
    end
end

%% Plot FiberWaveSpeed results - CSD only
minLength = 200;
close all;
opt = {[0.08,0.08], [0.08,0.06], [0.08,0.02]};  % {[vert, horz], [bottom, top], [left, right] }
figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
k = 0;
exptColor = distinguishable_colors(numel(xFiber));
for x = intersect(xCSD, intersect(xFiber, xPresent))
    k = k+1;
    [fiberLengthSort, fSort] = sort([fiber{x}.length], 'descend');
    fiberColor = distinguishable_colors(Nfiber(x));
    % Show that there is negligible delay between distant ROIs from the same fiber for spontaneous activity
    %tempFiber = fiber{x}; tempFiber = [tempFiber.CSD]
    for f = fSort(fiberLengthSort > minLength) %1:Nfiber(x) find(~cellfun(@isempty, {fiber{x}.CSD}))
        roiColor = distinguishable_colors(fiber{x}(f).Nroi);
        % Show fiber ROIs
        subplot(2,3,1:3);
        roiLabel = zeros(size(fiber{x}(f).labelFoot));
        for r = 1:fiber{x}(f).Nroi
            roiLabel(ROI{x}(fiber{x}(f).ROI(r)).footprintInd) = r; %roiLabel(roiLabel==fiber{x}(f).ROI(r)) = r; 
        end
        imshow(label2rgb(roiLabel, roiColor));
        title(sprintf('%s fiber %i', expt(x).name, f), 'Interpreter','none' ); 
        
        subplot(2,3,4:5); cla;
        for r = 1:fiber{x}(f).Nroi %
            if ~isempty(fiber{x}(f).CSD.trace)
                plot( fiber{x}(f).CSD.Twindow, fiber{x}(f).CSD.trace(:,r), 'color',roiColor(r,:) ); hold on;
                if fiber{x}(f).CSD.goodFit(1,r)
                    plot( fiber{x}(f).CSD.Twindow, predict(fiber{x}(f).CSD.traceFit{r}, fiber{x}(f).CSD.Twindow), '--', 'color',roiColor(r,:));
                end
            end
        end
        xlim([-Inf,Inf]); ylim([-Inf,Inf]);
        xlabel('Peri-Run Time'); ylabel('Normalized Fluor');
        
        subplot(2,3,6);  cla;
        if size(fiber{x}(f).CSD.sepDelay,1) > 1
            plot(fiber{x}(f).CSD.sepDelay(:,1), fiber{x}(f).CSD.sepDelay(:,2), '.', 'color',fiberColor(f,:)); hold on;
            if ~isempty(fiber{x}(f).CSD.sepDelayFit)
                plot(fiber{x}(f).CSD.sepDelay(:,1), polyval([fiber{x}(f).CSD.sepDelaySlope, fiber{x}(f).CSD.sepDelayConst],fiber{x}(f).CSD.sepDelay(:,1) ), '--', 'color',fiberColor(f,:) )
                %plot(fiber{x}(f).CSD.sepDelay(:,1), predict( fiber{x}(f).CSD.sepDelayFit, fiber{x}(f).CSD.sepDelay(:,1)), '--', 'color',fiberColor(f,:) );
            end
        end
        xlabel('Pairwise Separation (um)'); ylabel('Relative Delay (s)'); title('All well-fit bouts');
        axis square; hold off;
        pause%(0.1);
    end
end

%% CSD speed vs fiber orientation
close all;
figure;
for t = find(~isnan(slopeTable(:,4)))'
    x = slopeTable(t,1);
    f = slopeTable(t,2);
    tempAngSep = AngularSeparation(fiber{x}(f).orientation, csdWave{x}.angle, pi); % , false
    polarplot( tempAngSep*[1,1], (100/abs(slopeTable(t,4)))*[0,1], 'LineWidth',2); hold on;
end

%% Calculate CSD cross-correlation within each fiber
opt = {[0.07,0.07], [0.08,0.06], [0.05,0.05]};  % {[vert, horz], [bottom, top], [left, right] }
maxDelay = 10;
minXC = 0.5;
minSep = 150; % (2/expt(x).scanRate)*csdWave{x}.speed;
onsetWindow = [-10, 60]; 
FS = 10;
close all; clearvars sp;  figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
for x = xWave %37 %intersect(xPresent, find(Nfiber>0))   % 37 %
    gaussFilt = MakeGaussFilt( 4, 0, 2/expt(x).scanRate, expt(x).scanRate, false ); 
    dT = 1/expt(x).scanRate;
    tempPreFluorZ = fluor{x}(1:expt(x).csd-1).z;
    tempPreFluorZ = vertcat( tempPreFluorZ.ROI );
    for f = 1:Nfiber(x) 
        % PRE-CSD spontaneous activity
        % Get pre-CSD fiber fluor signals
        fiberPreFluor = tempPreFluorZ(:,fiber{x}(f).ROI); %fiberFluor = fluor{x}(expt(x).csd).z.ROI(onsetWindowScans, fiber{x}(f).ROI); % csdBout{x}(expt(x).csd).scan{1}
        fiberPreFluor(find(isnan(sum(fiberPreFluor,2))),:) = []; % remove missing data scans
        fiberPreFluor = filtfilt(gaussFilt, 1, fiberPreFluor);
        % Calculate pairwise cross-corr and peak delays
        [tempXcorr, tempLagScans] = xcorr( fiberPreFluor, round(maxDelay*expt(x).scanRate), 'normalized' );
        tempLag = dT*((1:numel(tempLagScans))  - median(1:numel(tempLagScans)));
        XCpeak_pre = nan(1, fiber{x}(f).Nroi^2); XCpeakScan_pre = nan(1, fiber{x}(f).Nroi^2);
        for c = flip(1:size(tempXcorr,2)),   [XCpeak_pre(c), XCpeakScan_pre(c)] = max( tempXcorr(:,c) );   end
        XCpeak_pre = reshape(XCpeak_pre, fiber{x}(f).Nroi, fiber{x}(f).Nroi );
        XClagPeak_pre = reshape(tempLag(XCpeakScan_pre), fiber{x}(f).Nroi, fiber{x}(f).Nroi );
        % Estimate spontaneous conduction velocity
        tempPairInd = find( triu(fiber{x}(f).footSep,1) > minSep & XCpeak_pre > minXC ); %use only well-separated, highly correlated pairs of ROI
        if ~isempty(tempPairInd)
            clearvars tempPair;
            [tempPair(:,1), tempPair(:,2)] = ind2sub( [fiber{x}(f).Nroi, fiber{x}(f).Nroi], tempPairInd ); 
            tempPairDisp = fiber{x}(f).cent( tempPair(:,2),1:2) - fiber{x}(f).cent( tempPair(:,1),1:2); % displacement vectors between pairs
            %tempPairDisp(abs(tempPairDisp) < minSep) = NaN; % suppress short dimensions
            fiber{x}(f).spont.pairLag = XClagPeak_pre(tempPairInd);
            fiber{x}(f).spont.pairVelocity = tempPairDisp./repmat(fiber{x}(f).spont.pairLag,1,2);
            fiber{x}(f).spont.velocity = median(fiber{x}(f).spont.pairVelocity, 1, 'omitnan');
            [fiber{x}(f).spont.angle, fiber{x}(f).spont.speed] = cart2pol(fiber{x}(f).spont.velocity(1), -fiber{x}(f).spont.velocity(2));
        else
            fiber{x}(f).spont.pairLag = NaN;
            fiber{x}(f).spont.pairVelocity = nan(1,2);
            fiber{x}(f).spont.velocity = nan(1,2);
            fiber{x}(f).spont.speed = NaN;
            fiber{x}(f).spont.angle = NaN;
            fiber{x}(f).spont.speed = NaN;
        end

        % PLOT SPONTANEOUS PROCESS
        %{
        subtightplot(2,4,[1,2], opt{:}); 
        imshow( label2rgb(fiber{x}(f).labelFoot)); hold on;
        for r = fiber{x}(f).ROI
            text( ROI{x}(r).cent(1), ROI{x}(r).cent(2), sprintf('%i', r), 'HorizontalAlignment','center', 'FontSize',8 );
        end
        title( sprintf('[x,f] = [%i, %i]  Orientation = %2.1f,  Length = %2.1f',x, f, fiber{x}(f).orientation, max(fiber{x}(f).footSep(:))) ); % fiber{x}(f).MajorAxisLength

        subtightplot(2,4,[3,4], opt{:});  
        %fiberFluorSpread = CalSpread( {1:size(fiberFluor,1)}, {fiberFluor}, 'sepPct',50, 'show',false ); %plot( fiberFluorSpread{1}  ); %
        imagesc( fiberPreFluor'); cb = colorbar; cb.Label.String = 'z-score';
        xlabel('Pre-CSD Time (s)');  ylabel('ROI');
        set(gca, 'Ytick',1:fiber{x}(f).Nroi, 'YtickLabel',fiber{x}(f).ROI ); % 'Xtick',XtickOnset, 'XtickLabel',sprintfc('%2.1f', Tonset(XtickOnset)), 
        xlim([-Inf,Inf]);
        title( sprintf('Estimated conduction velocity = [%2.1f,  %2.1f] um/s  (AP, ML)',fiber{x}(f).spont.velocity ) ); % ,  %2.1f

        subtightplot(2,4,5, opt{:});  
        imagesc( fiber{x}(f).footSep ); axis square;
        colorbar;
        title( sprintf('Pairwise Distance (min = %2.2f um)', minSep )); %title('Separation'); 
        set(gca,'Xtick',1:fiber{x}(f).Nroi, 'XtickLabel',fiber{x}(f).ROI, 'Ytick',1:fiber{x}(f).Nroi, 'YtickLabel',fiber{x}(f).ROI); 
        xtickangle(30);

        subtightplot(2,4,6, opt{:});  
        imagesc( XCpeak_pre ); axis square;
        colorbar; caxis([0.5,1]);
        set(gca,'Xtick',1:fiber{x}(f).Nroi, 'XtickLabel',fiber{x}(f).ROI, 'Ytick',1:fiber{x}(f).Nroi, 'YtickLabel',fiber{x}(f).ROI); 
        xtickangle(30);
        title( sprintf('Peak Cross-Correlation (min = %2.2f)', minXC ));  %title('Peak Cross-Correlation'); 
        
        subtightplot(2,4,7, opt{:})
        imagesc( XClagPeak_pre ); axis square; hold on;
        if ~isempty(tempPairInd), plot( tempPair(:,2), tempPair(:,1) ,'k.' ); end
        
        
        CB = colorbar; CB.Label.String = 'Delay (s)';
        title(sprintf('Cross-Correlation Peak Latency (max = %2.1f)', maxDelay)); %title('Cross-Correlation Peak Latency');
        set(gca,'Xtick',1:fiber{x}(f).Nroi, 'XtickLabel',fiber{x}(f).ROI, 'Ytick',1:fiber{x}(f).Nroi, 'YtickLabel',fiber{x}(f).ROI); 
        xtickangle(30);
        impixelinfo;
        
        subtightplot(2,4,8, opt{:});
        tempPairVelocity = fiber{x}(f).spont.pairVelocity; tempPairVelocity(isnan(tempPairVelocity))= 0;
        [pairTheta, pairMag] = cart2pol( tempPairVelocity(:,1), tempPairVelocity(:,2) );
        for p = 1:length(pairTheta)
             polarplot( pairTheta(p)*[1,1], pairMag(p)*[0,1], 'LineWidth',1 ); hold on;
        end
        tempVelocity = fiber{x}(f).spont.velocity; tempVelocity(isnan(tempVelocity))= 0;
        [medTheta, medMag] = cart2pol( tempVelocity(:,1), tempVelocity(:,2) );
        polarplot( medTheta*[1,1], medMag*[0,1], 'LineWidth',1.5, 'Color','k' ); hold off;
        pause;
        %}
        
        % CSD
        % Get peri-CSD fiber fluor signals
        Tperi = Tscan{x}{expt(x).csd} - csdBout{x}(expt(x).csd).Tstart; % Tcat
        if isempty( onsetWindow ), onsetWindow = [0, 0.67*csdBout{x}(expt(x).csd).dur(1)]; end
        onsetWindowScans = find(Tperi >= onsetWindow(1) & Tperi <= onsetWindow(2))';
        Tonset = Tscan{x}{expt(x).csd}(onsetWindowScans) - csdBout{x}(expt(x).csd).Tstart;
        XtickZero = find(Tonset == 0);
        XtickOnset = 1:XtickZero-1:numel(Tonset);  
        fiberCSDfluor = fluor{x}(expt(x).csd).z.ROI(onsetWindowScans, fiber{x}(f).ROI); % csdBout{x}(expt(x).csd).scan{1}
        fiberCSDfluor(find(isnan(sum(fiberCSDfluor,2))),:) = []; % remove missing data scans
        fiberCSDfluor = filtfilt(gaussFilt, 1, fiberCSDfluor);
        % Calculate pairwise cross-corr and peak delays
        [tempXcorr, tempLagScans] = xcorr( fiberCSDfluor, round(maxDelay*expt(x).scanRate), 'normalized' );
        tempLag = dT*((1:numel(tempLagScans))  - median(1:numel(tempLagScans)));
        XCpeak_CSD = nan(1, fiber{x}(f).Nroi^2); XCpeakScan = nan(1, fiber{x}(f).Nroi^2);
        for c = flip(1:size(tempXcorr,2))
            [XCpeak_CSD(c), XCpeakScan(c)] = max( tempXcorr(:,c) );
        end
        XCpeak_CSD = reshape(XCpeak_CSD, fiber{x}(f).Nroi, fiber{x}(f).Nroi );
        XCpeakLag_CSD = reshape(tempLag(XCpeakScan), fiber{x}(f).Nroi, fiber{x}(f).Nroi );
        
        % Estimate CSD conduction velocity
        tempPairInd = find( triu(fiber{x}(f).footSep,1) > minSep & XCpeak_CSD > minXC ); %use only well-separated, highly correlated pairs of ROI
        if ~isempty(tempPairInd)
            clearvars tempPair;
            [tempPair(:,1), tempPair(:,2)] = ind2sub( [fiber{x}(f).Nroi, fiber{x}(f).Nroi], tempPairInd ); 
            tempPairDisp = fiber{x}(f).cent(tempPair(:,2),1:2) - fiber{x}(f).cent(tempPair(:,1),1:2); % displacement vectors between pairs
            %tempPairDisp(abs(tempPairDisp) < minSep) = NaN; % suppress short dimensions
            fiber{x}(f).CSD.pairLag = XCpeakLag_CSD(tempPairInd);
            fiber{x}(f).CSD.pairVelocity = tempPairDisp./repmat(fiber{x}(f).CSD.pairLag,1,2);
            fiber{x}(f).CSD.velocity = median(fiber{x}(f).CSD.pairVelocity, 1, 'omitnan');
            [fiber{x}(f).CSD.angle, fiber{x}(f).CSD.speed] = cart2pol(fiber{x}(f).CSD.velocity(1), -fiber{x}(f).CSD.velocity(2));
        else
            fiber{x}(f).CSD.pairLag = NaN;
            fiber{x}(f).CSD.pairVelocity = nan(1,2);
            fiber{x}(f).CSD.velocity = nan(1,2);
            fiber{x}(f).CSD.speed = NaN;
            fiber{x}(f).CSD.angle = NaN;
            fiber{x}(f).CSD.speed = NaN;
        end
        
        %{ 
        close all;
        figure;
        subtightplot(1,2,1, opt{:}); 
        imshow( label2rgb(fiber{x}(f).labelFoot)); hold on;
        for r = fiber{x}(f).ROI
            text( ROI{x}(r).cent(1), ROI{x}(r).cent(2), sprintf('%i', r), 'HorizontalAlignment','center', 'FontSize',8 );
        end
        DrawArrow(expt(x).Ncol/2, expt(x).Nrow/2, csdWave{x}.velocity(1)/expt(x).scanRate, csdWave{x}.velocity(2)/expt(x).scanRate); % Draw an arrow representing the CSD velocity      
        %title( sprintf('[x,f] = [%i, %i]  Orientation = %2.1f,  Length = %2.1f',x, f, fiber{x}(f).orientation, max(fiber{x}(f).footSep(:))) ); 
        clearvars tempPairInd tempPair tempPairDisp tempPairLag tempPairSepCSD;
        tempPairInd = find( triu(fiber{x}(f).footSep,1) > minSep & XCpeak_CSD > minXC ); %use only well-separated, highly correlated pairs of ROI
        [tempPair(:,1), tempPair(:,2)] = ind2sub( [fiber{x}(f).Nroi, fiber{x}(f).Nroi], tempPairInd ); 
        for p = flip(1:size(tempPair,1))
            tempPairDisp{p} = fiber{x}(f).cent(tempPair(p,2),1:2) - fiber{x}(f).cent(tempPair(p,1),1:2);
            tempPairLag(p) = XCpeakLag_CSD(tempPairInd(p));
            tempPairSepCSD(p) = dot( tempPairDisp{p}, csdWave{x}.velocity/norm(csdWave{x}.velocity) );
            
            DrawArrow(fiber{x}(f).cent(tempPair(p,1),1), fiber{x}(f).cent(tempPair(p,1),2), tempPairDisp{p}(1), tempPairDisp{p}(2));
            %pause;
        end
        
        tempLM = fitlm( tempPairSepCSD, tempPairLag, 'RobustOpts',false );
        tempSlope = tempLM.Coefficients.Estimate(2); % s/um
        subtightplot(1,2,2, opt{:});
        plot(tempPairSepCSD, tempPairLag, '.' ); hold on;
        plot( tempPairSepCSD, predict(tempLM, tempPairSepCSD'), 'lineStyle','-');
        axis square;
        tempXrange = get(gca,'Xlim');
        plot( tempXrange, (1/csdWave{x}.speed)*tempXrange, 'lineStyle','--' )

        %}
        
        

        % PLOT CSD PROCESS
        %{
        subtightplot(2,4,[1,2], opt{:}); 
        imshow( label2rgb(fiber{x}(f).labelFoot)); hold on;
        for r = fiber{x}(f).ROI
            text( ROI{x}(r).cent(1), ROI{x}(r).cent(2), sprintf('%i', r), 'HorizontalAlignment','center', 'FontSize',8 );
        end
        title( sprintf('[x,f] = [%i, %i]  Orientation = %2.1f,  Length = %2.1f',x, f, fiber{x}(f).orientation, max(fiber{x}(f).footSep(:))) ); % fiber{x}(f).MajorAxisLength

        subtightplot(2,4,[3,4], opt{:});  
        %fiberFluorSpread = CalSpread( {1:size(fiberFluor,1)}, {fiberFluor}, 'sepPct',50, 'show',false ); %plot( fiberFluorSpread{1}  ); %
        imagesc( fiberCSDfluor'); cb = colorbar; cb.Label.String = 'z-score';
        xlabel('Peri-CSD Time (s)');  ylabel('ROI');
        set(gca, 'Xtick',XtickOnset, 'XtickLabel',sprintfc('%2.1f', Tonset(XtickOnset)), 'Ytick',1:fiber{x}(f).Nroi, 'YtickLabel',fiber{x}(f).ROI );
        xlim([-Inf,Inf]);
        title( sprintf('Estimated conduction velocity = [%2.1f,  %2.1f] um/s  (AP, ML)',fiber{x}(f).CSD.velocity ) ); % ,  %2.1f

        subtightplot(2,4,5, opt{:});  
        imagesc( fiber{x}(f).footSep ); axis square;
        colorbar;
        title( sprintf('Pairwise Distance (min = %2.2f um)', minSep )); %title('Separation'); 
        set(gca,'Xtick',1:fiber{x}(f).Nroi, 'XtickLabel',fiber{x}(f).ROI, 'Ytick',1:fiber{x}(f).Nroi, 'YtickLabel',fiber{x}(f).ROI, 'FontSize',FS); 
        xtickangle(30);

        subtightplot(2,4,6, opt{:});  
        imagesc( XCpeak_CSD ); axis square;
        colorbar; caxis([0.5,1]);
        set(gca,'Xtick',1:fiber{x}(f).Nroi, 'XtickLabel',fiber{x}(f).ROI, 'Ytick',1:fiber{x}(f).Nroi, 'YtickLabel',fiber{x}(f).ROI, 'FontSize',FS); 
        xtickangle(30);
        title( sprintf('Peak Cross-Correlation (min = %2.2f)', minXC ));  %title('Peak Cross-Correlation'); 
        
        subtightplot(2,4,7, opt{:})
        imagesc( XCpeakLag_CSD ); hold on;
        if ~isempty(tempPairInd), plot( tempPair(:,2), tempPair(:,1) ,'k.' ); end
        axis square;
        CB = colorbar; CB.Label.String = 'Delay (s)';
        title(sprintf('Cross-Correlation Peak Latency (max = %2.1f)', maxDelay)); %title('Cross-Correlation Peak Latency');
        set(gca,'Xtick',1:fiber{x}(f).Nroi, 'XtickLabel',fiber{x}(f).ROI, 'Ytick',1:fiber{x}(f).Nroi, 'YtickLabel',fiber{x}(f).ROI, 'FontSize',FS); 
        xtickangle(30);
        impixelinfo;
        
        subtightplot(2,4,8, opt{:});
        tempPairVelocity = fiber{x}(f).CSD.pairVelocity; tempPairVelocity(isnan(tempPairVelocity))= 0;
        [pairTheta, pairMag] = cart2pol( tempPairVelocity(:,1), tempPairVelocity(:,2) );
        for p = 1:length(pairTheta)
             polarplot( pairTheta(p)*[1,1], pairMag(p)*[0,1], 'LineWidth',1 ); hold on;
        end
        tempVelocity = fiber{x}(f).CSD.velocity; tempVelocity(isnan(tempVelocity))= 0;
        [medTheta, medMag] = cart2pol( tempVelocity(:,1), tempVelocity(:,2) );
        polarplot( medTheta*[1,1], medMag*[0,1], 'LineWidth',1.5, 'Color','k' ); hold off;
        pause; 
        %}
    end
end


%% Summarize INVERSE speed
%opt = {[0.1,0.07], [0.1,0.04], [0.05,0.05]};  % {[vert, horz], [bottom, top], [left, right] }
close all; clearvars h sp; 
ConductionSpeedFig = figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
MS = 10; BW = 0.2;
colorPresent = distinguishable_colors(numel(xPresent));
k = 0;
for x = xPresent
    k = k+1;
    tempAPfiber = abs([fiber{x}.orientation]) < 30;
    tempCSD = [fiber{x}(tempAPfiber).CSD];
    tempCSDspd = abs( 1./vertcat( tempCSD.velocity ) );
    tempCSDvelMean = mean(tempCSDspd, 1, 'omitnan');
    tempCSDvelSEM = std(tempCSDspd, 0, 1, 'omitnan')./sqrt(sum(~isnan(tempCSDspd) ) );
    tempSpont = [fiber{x}.spont]; % (tempAPfiber)
    tempSpontSpd = abs( 1./vertcat( tempSpont.velocity ) );
    tempSpontSpdMean = mean(tempSpontSpd, 1, 'omitnan');
    tempSpontSpdSEM = std(tempSpontSpd, 0, 1, 'omitnan')./sqrt(sum(~isnan(tempSpontSpd) ) );
    
    h(1) = bar( k, 1/abs(csdWave{x}(1,1,1)), BW, 'FaceColor','k' ); hold on;
    h(2) = bar( k+BW, tempCSDvelMean(1), BW, 'FaceColor','c');
    errorbar( k+BW, tempCSDvelMean(1), tempCSDvelSEM(1), 'Color','k' );
    h(3) = bar( k+2*BW, tempSpontSpdMean(1), BW, 'FaceColor','m');
    errorbar( k+2*BW, tempSpontSpdMean(1), tempSpontSpdSEM(1), 'Color','k' );
    
    %pause;
end
ylabel('Inverse A-P Conduction Speed (s/um)');
legend(h, {'Neuropil', 'Axon-CSD', 'Axon-Spontaneous'});
set(gca, 'Xtick', (1:k)+BW, 'XtickLabel',{expt(xPresent).name}, 'TickLabelInterpreter','none' )

%figDir = 'D:\MATLAB\LevyLab\Figures\BoutMovies\CSD\';
%figPath = strcat(figDir, 'InverseConductionSpeed'); % sprintf('InverseConductionSpeed', figDir);
%print(ConductionSpeedFig, figPath, '-dtiff'); fprintf('\nSaved %s', figPath); 

%% Among qualified fibers, how many finite conduction speeds were found?

for x = xPresent
    tempSpont = [fiber{x}.spont];
    vertcat(tempSpont.velocity)
    %Nqual(x) = 
    %Ninf(x) = 
end



%% Summarize conduction speed
%opt = {[0.1,0.07], [0.1,0.04], [0.05,0.05]};  % {[vert, horz], [bottom, top], [left, right] }
close all; clearvars h sp; 
ConductionSpeedFig = figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
MS = 10; BW = 0.2;
colorPresent = distinguishable_colors(numel(xPresent));
k = 0;
for x = xPresent
    tempAPfiber = abs([fiber{x}.orientation]) < 30;
    tempCSD = [fiber{x}(tempAPfiber).CSD];
    tempCSDvel = vertcat( tempCSD.velocity );
    tempCSDvelMean = mean(tempCSDvel, 1, 'omitnan');
    tempCSDvelSEM = std(tempCSDvel, 0, 1, 'omitnan')./sqrt(sum(~isnan(tempCSDvel) ) );
    tempCSDspd = vertcat( tempCSD.speed ); %abs( vertcat( tempCSD.velocity ) );
    tempCSDspdMean = mean(tempCSDspd, 1, 'omitnan');
    tempCSDspdSEM = std(tempCSDspd, 0, 1, 'omitnan')./sqrt(sum(~isnan(tempCSDspd) ) );   
    tempSpont = [fiber{x}.spont]; % (tempAPfiber)
    tempSpontSpd = tempSpont.speed; %abs( vertcat( tempSpont.velocity ) );
    tempSpontSpdMean = mean(tempSpontSpd, 1, 'omitnan');
    tempSpontSpdSEM = std(tempSpontSpd, 0, 1, 'omitnan')./sqrt(sum(~isnan(tempSpontSpd) ) );
    
    k = k+1;
    h(1) = bar( k,csdWave{x}.speed, BW, 'FaceColor','k' ); hold on;
    h(2) = bar( k+BW, tempCSDspdMean(1), BW, 'FaceColor','c');
    errorbar( k+BW, tempCSDspdMean(1), tempCSDvelSEM(1), 'Color','k' );
    h(3) = bar( k+2*BW, tempSpontSpdMean(1), BW, 'FaceColor','m');
    errorbar( k+2*BW, tempSpontSpdMean(1), tempSpontSpdSEM(1), 'Color','k' );
    
    %pause;
end
ylabel('A-P Conduction Speed (s/um)');
legend(h, {'Neuropil', 'Axon-CSD', 'Axon-Spontaneous'});
set(gca, 'Xtick', (1:k)+BW, 'XtickLabel',{expt(xPresent).name}, 'TickLabelInterpreter','none' )

%figDir = 'D:\MATLAB\LevyLab\Figures\BoutMovies\CSD\';
%figPath = strcat(figDir, 'InverseConductionSpeed'); % sprintf('InverseConductionSpeed', figDir);
%print(ConductionSpeedFig, figPath, '-dtiff'); fprintf('\nSaved %s', figPath); 

%% Show fiber orientation vs CSD wave
clearvars sp
figure; 
for x = xPresent
    %sp(1) = subtightplot(1,2,1, opt{:}); 
    polarplot( csdWave{x}.angle*[1,1], csdWave{x}.speed*[0,1], 'k', 'LineWidth',2 ); hold on;

    for f = 1:Nfiber(x)

        polarplot( fiber{x}(f).orientation*[1,1], fiber{x}(f).CSD.speed*[1,1]*[0,1] ); %hold off;
        %title( sprintf('x = %i:  CSD angle = %2.1f',x, csdWave{x}.angle ) );
%{
        subtightplot(1,2,2, opt{:}); 
        imshow( label2rgb(fiber{x}(f).labelFoot)); hold on;
        for r = fiber{x}(f).ROI
            text( ROI{x}(r).cent(1), ROI{x}(r).cent(2), sprintf('%i', r), 'HorizontalAlignment','center', 'FontSize',8 );
        end
        title( sprintf('[x,f] = [%i, %i]  Orientation = %2.1f,  Length = %2.1f',x, f, fiber{x}(f).orientation, max(fiber{x}(f).footSep(:))) );
        %}
       % pause;
        %cla;
    end
    pause;
    cla;
end

    thetaRotData = thetaDirData - thetaRef; % rotate thetaIn so that thetaRef = 0
    thetaRotData(thetaRotData < -pi) = thetaRotData(thetaRotData < -pi) + 2*pi; 
    thetaRotData(thetaRotData >= pi) = thetaRotData(thetaRotData >= pi) - 2*pi; % direction: [-pi, pi)

    subplot(1,2,1);
    polarplot( thetaRef, 1, 'k.', 'MarkerSize',10 ); hold on;
    polarplot( thetaDirData, ones(numel(thetaDirData),1), '.' )
    subplot(1,2,2);
    polarplot( thetaRotData, ones(numel(thetaRotData),1) ,'.' ) 
 %%   
close all;
figure;
for x = xPresent
    for f = 1:Nfiber(x)  
        tempAngSep = AngularSeparation(fiber{x}(f).orientation, csdWave{x}.angle, pi); % , false
        tempSpeedDiff = abs(csdWave{x}.speed - fiber{x}(f).CSD.speed);
        if isinf(tempSpeedDiff), tempSpeedDiff = NaN; end
        polarplot( tempAngSep*[1,1], tempSpeedDiff*[0,1], 'LineWidth',2); hold on;
    end
    pause;
end

%%
allWaves = [csdWave{xPresent}];
allWaveSpeeds = [allWaves.speed];
allAlignedSpeeds = []; allOrthogSpeeds = []; allSpontSpeeds = []; allFiberCSDspeeds = []; allSpontSpeeds = [];
for x = xPresent
    for f = 1:Nfiber(x)
        if AngularSeparation(fiber{x}(f).orientation, csdWave{x}.angle, pi, false) < pi/6
            allAlignedSpeeds = [allAlignedSpeeds, fiber{x}(f).CSD.speed];
        elseif AngularSeparation(fiber{x}(f).orientation, csdWave{x}.angle, pi, false) > 2*pi/6
            allOrthogSpeeds = [allOrthogSpeeds, fiber{x}(f).CSD.speed];
        end      
        allSpontSpeeds = [allSpontSpeeds, fiber{x}(f).spont.speed];
        allFiberCSDspeeds = [allFiberCSDspeeds, fiber{x}(f).CSD.speed];
    end
end
allWaveSpeeds = allWaveSpeeds(~isnan(allWaveSpeeds)); % isfinite(allWaveSpeeds) & 
allAlignedSpeeds = allAlignedSpeeds(~isnan(allAlignedSpeeds)); % isfinite(allAlignedSpeeds) & 
allOrthogSpeeds = allOrthogSpeeds(~isnan(allOrthogSpeeds)); % isfinite(allOrthogSpeeds) & 
allSpontSpeeds = allSpontSpeeds(~isnan(allSpontSpeeds));
allFiberCSDspeeds = allFiberCSDspeeds(~isnan(allFiberCSDspeeds));

infFrac = nan(1,4);
infFrac(1) = sum(~isfinite(allWaveSpeeds))/numel(allWaveSpeeds);
infFrac(2) = sum(~isfinite(allAlignedSpeeds))/numel(allAlignedSpeeds);
infFrac(3) = sum(~isfinite(allOrthogSpeeds))/numel(allOrthogSpeeds);
infFrac(4) = sum(~isfinite(allSpontSpeeds))/numel(allSpontSpeeds);
%infFrac(5) = sum(~isfinite(allFiberCSDspeeds))/numel(allFiberCSDspeeds);

close all; clearvars h sp; 
ConductionSpeedFig = figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'Color','w');
MS = 10; BW = 0.2;
subplot(1,2,1);
JitterPlot( {allWaveSpeeds(isfinite(allWaveSpeeds)), allAlignedSpeeds(isfinite(allAlignedSpeeds)), allOrthogSpeeds(isfinite(allOrthogSpeeds)), allSpontSpeeds(isfinite(allSpontSpeeds))}, 0.2 )
ylabel('Conduction Speed (um/s)'); %ylabel('Inverse Conduction Speed (s/um)'); %
set(gca, 'Xtick', 1:4, 'XtickLabel',{'CSD Wave','Fiber (CSD, Aligned)','Fiber (CSD, Orthogonal)', 'Fiber, Spontaneous'}, 'TickLabelInterpreter','none', 'TickDir','out','box','off'); xtickangle(30);
axis square;

subplot(1,2,2);
bar( infFrac ); 
ylabel('Fraction of Fibers with Infinite Conduction Speed');
set(gca, 'Xtick', 1:5, 'XtickLabel',{'CSD Wave','Fiber (CSD, Aligned)','Fiber (CSD, Orthogonal)', 'Fiber (Spontaneous)'}, 'TickLabelInterpreter','none', 'TickDir','out','box','off' );
xtickangle(30);
axis square;

figDir = 'D:\MATLAB\LevyLab\Figures\BoutMovies\CSD\';
figPath = strcat(figDir, 'ConductionSpeed'); % sprintf('InverseConductionSpeed', figDir);
print(ConductionSpeedFig, figPath, '-dtiff'); fprintf('\nSaved %s', figPath); 


%% Revised figure

% A - still frames of wave progress, with time shown, some fiber ROIs shown, and scalebar
X = 37; % DL75
% Generate movies of fibers during registered CSD bouts
movParam.dir = 'D:\MATLAB\LevyLab\Figures\BoutMovies\CSD\Fibers\'; % 'D:\MATLAB\LevyLab\Figures\CSD\Movies\';
mkdir(movParam.dir)
movParam.fmtSpec = '%2.1f';
movParam.binT = 1;
movParam.displayPct = [2,99.7];
movParam.sbx = [];
movParam.edges = []; % [80,80,20,20]; %[80,90,60,110]; %
movParam.aviRate = 10; % frames per second
movParam.level = 'ind';
movParam.edges = segParams{X}.edges;
if expt(X).Nplane > 1
    movParam.sourceSbx = 'sbx_interp';  
else
    movParam.sourceSbx = 'sbx_affine'; %'sbxz'; %
end
movParam.scalebar = MakeScaleBar( round(expt(X).umPerPixel*[50,50]), {[0,expt(X).Ncol-segParams{X}.edges(1)-segParams{X}.edges(2)]+0.5, [0,expt(X).Nrow-segParams{X}.edges(3)-segParams{X}.edges(4)]+0.5},...
    [0.1,0.95], [0,0], 'label',false, 'color','w', 'show',false );
movParam.zProj = unique(round(fiber{X}(1).cent(:,3)))';
% Write CSD movie
movParam.regType = sprintf('_aff_ExampleFibers_CSD'); %'
movParam.Toffset = 14.5;
movParam.Tperi = [0,0];
movParam.boutType = 'csd';
[~,TcsdMov,~,csdFrames] = WriteBoutMovies(expt(X), catInfo(X), Tscan{X}, loco{X}, csdBout{X}, movParam, ROI{X}(fiber{X}(1).ROI) );
%[~,TcsdMov,~,csdFrames] = WriteBoutMovies(expt(X), catInfo(X), Tscan{X}, loco{X}, csdBout{X}, movParam, ROI{X}, fiber{X}(1)); % , 'overwrite',true (1)  , ROI{X}  , 'scalebar',
% Write example loco bout movie
movParam.boutType = 'loco';
movParam.Tperi = [4,0];
movParam.regType = sprintf('_aff_ExampleFibers_spont'); %'
movParam.Toffset = 0;
[~,TspontMov,~,spontFrames] = WriteBoutMovies(expt(X), catInfo(X), Tscan{X}, loco{X}, periBout{X}, movParam, ROI{X}(fiber{X}(1).ROI), 'run',1, 'bout',1); % , fiber{X}(1), 'overwrite',true (1)  , ROI{X}  , 'scalebar',

%%
close all; clearvars h sp; 
plotOpt = {[0.08,0.02], [0.07, 0.001], [0.04, 0.06]};  % {[vert, horz], [bottom, top], [left, right] }
frameOpt = {[0.01,0.01], [0.07, 0.001], [0.04, 0.06]};  % {[vert, horz], [bottom, top], [left, right] }
CSDspreadFig = figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'Color','w');
subtightplot(3,4,1, frameOpt{:});  
imshow( spontFrames(3).cdata, [] ); axis image;
ylabel('Spontaneous');

subtightplot(3,4,2, frameOpt{:});  
imshow( spontFrames(7).cdata, [] ); axis image;

subtightplot(3,4,3, frameOpt{:});  
imshow( spontFrames(11).cdata, [] ); axis image;

subtightplot(3,4,4, plotOpt{:}); 
imagesc([fiber{X}(1).spont.trace{1,:}]' );
CB = colorbar;
CB.Label.String = 'Normalized Fluorescence';
XtickZero = find(fiber{X}(1).spont.Twindow{1} == 0);
XtickOnset = 1:XtickZero-1:numel(fiber{X}(1).spont.Twindow{1});
set(gca,'Xtick',XtickOnset, 'XtickLabel',sprintfc('%2.1f', fiber{X}(1).spont.Twindow{1}(XtickOnset)), 'Ytick',[]);
xlabel('Peri-Locomotion Time (s)'); ylabel('ROI');

subtightplot(3,4,5, frameOpt{:});  
imshow( csdFrames(19).cdata, [] ); axis image;
ylabel('Spreading Depression');

subtightplot(3,4,6, frameOpt{:});  
imshow( csdFrames(22).cdata, [] ); axis image;

subtightplot(3,4,7, frameOpt{:});  
imshow( csdFrames(25).cdata, [] ); axis image;

subtightplot(3,4,8, plotOpt{:}); 
imagesc([fiber{X}(1).CSD.trace]' );
CB = colorbar;
CB.Label.String = 'Normalized Fluorescence';
XtickZero = 7; %find(fiber{X}(1).CSD.Twindow == 0);
XtickOnset = 1:XtickZero-1:numel(fiber{X}(1).CSD.Twindow);
set(gca,'Xtick',XtickOnset, 'XtickLabel',sprintfc('%2.1f', fiber{X}(1).CSD.Twindow(XtickOnset)), 'Ytick',[]);
xlabel('Peri-CSD Time (s)'); ylabel('ROI');

subtightplot(3,4,9, plotOpt{:});
plot(fiber{X}(1).spont.sepDelay(:,1), fiber{X}(1).spont.sepDelay(:,2), '.', 'color','k'); hold on;
plot(fiber{X}(1).spont.sepDelay(:,1), polyval([fiber{X}(1).spont.sepDelaySlope, fiber{X}(1).spont.sepDelayConst],fiber{X}(1).spont.sepDelay(:,1) ), '--', 'color','k' )
hold off;
xlabel('Pairwise Separation (\mum)'); ylabel('Relative Delay (s)'); title('Spontaneous');
xlim([0,450]); ylim([-5,2]); 
axis square; 

subtightplot(3,4,10,plotOpt{:});
plot(fiber{X}(1).CSD.sepDelay(:,1), fiber{X}(1).CSD.sepDelay(:,2), '.', 'color','k'); hold on;
if ~isempty(fiber{X}(1).CSD.sepDelayFit)
    plot(fiber{X}(1).CSD.sepDelay(:,1), polyval([fiber{X}(1).CSD.sepDelaySlope, fiber{X}(1).CSD.sepDelayConst],fiber{X}(1).CSD.sepDelay(:,1) ), '--', 'color','k' )
    %plot(fiber{X}(1).CSD.sepDelay(:,1), predict( fiber{X}(1).CSD.sepDelayFit, fiber{X}(1).CSD.sepDelay(:,1)), '--', 'color',fiberColor(1,:) );
    %pause;
end
xlim([0,450]); ylim([-5,2]); 
hold off;
xlabel('Pairwise Separation (\mum)'); ylabel('Relative Delay (s)'); title('CSD');
axis square; 

subtightplot(3,4,11,plotOpt{:})
JitterPlot( abs(slopeTable(:,3:4)), 0.45);
[~,pSlope] = ttest2( abs(slopeTable(:,3)), abs(slopeTable(:,4)));
title(sprintf('p = %2.5f', pSlope)); box off;
set(gca,'Xtick',[1,2], 'XtickLabel',{sprintf('Loco (n=%i)', sum(~isnan(slopeTable(:,3)))),sprintf('CSD (n=%i)', sum(~isnan(slopeTable(:,4))))});
ylabel('Delay per Separation (s/100 \mum)');
xtickangle(25);
axis square;

subtightplot(3,4,12,plotOpt{:});
clearvars h; k = 0;
line([-1,10], [-1,10],'color','r', 'lineStyle','--'); hold on;
axis square;
%set(gca,'XaxisLocation','origin','YaxisLocation','origin')
xlabel('Neuropil Onset (s)'); ylabel('ROI Onset (s)');
fiberName = cell(1,0);
fiberColor = distinguishable_colors(30);
%FS = 16;
for x = intersect(xFiber, xCSD)
    try
        for f = find(~cellfun(@isempty, {fiber{x}.CSD}))
            k = k+1;
            fiberName{k} = sprintf('%s_fiber%i', expt(x).name, f);
            h(k) = plot( fiber{x}(f).CSD.binOnset, fiber{x}(f).CSD.Tonset, '.', 'MarkerSize',10, 'color',fiberColor(k,:) ); hold all;
            %pause;
        end
    catch
        fprintf('%s failed', expt(x).name);
    end
end
%set(gca, 'FontSize',FS)
%legend(h, fiberName, 'Location','EastOutside', 'Interpreter','none', 'FontSize',FS);
figDir = 'C:\Users\ablaeser\Documents\Afferent Paper\';
figPath = strcat(figDir, 'csdSpread_Draft'); % sprintf('InverseConductionSpeed', figDir);
%print(CSDspreadFig, figPath, '-dtiff'); fprintf('\nSaved %s', figPath); 

%% Across focal planes
close all;
figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
for x = xPresent
    for run = expt(x).csd
        TperiCSD = csdBout{x}(run).T{1} - csdBout{x}(run).Tstart;
        imagesc( fluor{x}(expt(x).csd).dFF.plane(csdBout{x}(run).scan{1},:)' );
        CB = colorbar('EastOutside'); CB.Label.String = 'dF/F'; CB.Label.FontWeight = 'bold';
        
        XtickZero = find(TperiCSD == 0);
        XtickCSD = 1:XtickZero-1:numel(TperiCSD);      
        set(gca, 'Xtick',XtickCSD, 'XtickLabel',sprintfc('%2.1f', TperiCSD(XtickCSD)) );       
        title( sprintf('%s', expt(x).name ), 'Interpreter','none');
    end
    pause;
end

%% Plot peri-CSD deformation
close all;
for x = xWave
    Tcsd = Tscan{x}{expt(x).csd} - csdBout{x}(expt(x).csd).Tstart(1);
    if expt(x).Nplane > 1
        plot( Tcsd, mean(deform{x}(expt(x).csd).scaleMag(:,segParams{x}.zProj), 2, 'omitnan') ); hold on;
    else
        plot( Tcsd, deform{x}(expt(x).csd).scaleMag ); %hold on;
    end
    xlim([-30,180]);
    pause;
    %csdBout{x}(expt(x).csd)
end


%% Compare deformation during locomotion vs CSD
scaleLim = [0,8];
close all;
figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
for x = xWave
    % Show the most deformational bout
    bMax = []; maxDur = [];
    for run = 1:expt(x).csd-1 %b = 
        [maxDur(run), bMax(run)] = max( periBout{x}(run).dur ); % find the most deformational bout
        %[maxDur(run), bMax(run)] = max( mean( periBout{x}(run).stat.scaleMag.effect(:,segParams{x}.zProj,1), 2, 'omitnan' ) ); % find the most deformational bout
    end
    [~,maxRun] = max(maxDur);
    
    subplot(1,2,1);
    imagesc( periBout{x}(maxRun).scaleMag{bMax(maxRun)}' );
    caxis(scaleLim)
    Tloco = periBout{x}(maxRun).T{bMax(maxRun)} - periBout{x}(maxRun).Tstart(bMax(maxRun));
    find(Tloco==0)
    set(gca,'Xtick', 1:find(Tloco==0)-1:numel(Tloco), 'XtickLabel', sprintfc('%2.0f', Tloco(1:find(Tloco==0)-1:numel(Tloco))) );
    
    subplot(1,2,2);
    imagesc( csdBout{x}(expt(x).csd).scaleMag{1}' ); 
    caxis(scaleLim)
    Tcsd = csdBout{x}(expt(x).csd).T{1} - csdBout{x}(expt(x).csd).Tstart(1);
    set(gca,'Xtick', 1:8:numel(Tcsd), 'XtickLabel', sprintfc('%2.0f', Tcsd(1:8:numel(Tcsd))) );
    
    impixelinfo
    pause;
end

%%
close all;
binRate = 15.49/30;
scaleEvent = cell(1,Nexpt); NscaleEvent = cell(1,Nexpt); 
for x = xAcute
    scaleEvent{x} = cell(1, expt(x).Nruns);
    for run = 1:expt(x).Nruns 
        % Calculate scaling magnitude, and normalize using still periods
        tempScaleMag = mean(deform{x}(run).scaleMag(:,segParams{x}.zProj), 2, 'omitnan');
        tempStillScan = [stillEpoch{x}([stillEpoch{x}.run] == run).scan];
        tempMean = mean(tempScaleMag(tempStillScan), 'omitnan');
        tempStd = std(tempScaleMag(tempStillScan), 'omitnan');
        tempScaleZ = (tempScaleMag-tempMean)/tempStd;
        % Downsample to equalize scanrates
        tempScaleT = BinDownMean( Tscan{x}{run}, expt(x).scanRate/binRate );
        tempScaleMag = BinDownMean( tempScaleMag, expt(x).scanRate/binRate );
        tempScaleZ = BinDownMean( tempScaleZ, expt(x).scanRate/binRate );
        % Detect events
        [scaleEvent{x}{run}, NscaleEvent{x}(run), ~, ~] = DetectEvents(tempScaleT, tempScaleZ, tempScaleMag, 'minMag',0, 'minDur',2, 'show',false); % NscaleEvent{x}(run,:)  expt(x), 
    end
end

% Find scaleEvent associated with stillness and locomotion
windowPad = [-2,2]; padDur = diff(windowPad);
csdPad = [-2,180];
SASE = cell(1,Nexpt); Nsase = cell(1,Nexpt); LASE = cell(1,Nexpt); Nlase = cell(1,Nexpt); CASE = cell(1,Nexpt); Ncase = cell(1,Nexpt);
tic;
for x = xAcute
    SASE{x} = cell(1,size(NscaleEvent{x},2));  LASE{x} = cell(1,size(NscaleEvent{x},2)); CASE{x} = cell(1,size(NscaleEvent{x},2));
    for run = 1:expt(x).Nruns           
        for roi = 1 %find(NscaleEvent{x}(run,:) > 0) %find(cellfun(@numel, scaleEvent{x}{run} ) > 0) %1:expt(x).Nroi % 1:expt(x).Nruns
            tempTstart = [scaleEvent{x}{run}{roi}.Tstart];
            % Find scale events within epochs of stillness
            for s = find([stillEpoch{x}.run] == run) %stillSumm{x}.sPre %1:stillSumm{x}.Nepoch
                eStill = find(tempTstart >= stillEpoch{x}(s).Tstart & tempTstart <= stillEpoch{x}(s).Tstop);
                if ~isempty(eStill)
                    SASE{x}{roi} = [SASE{x}{roi}, scaleEvent{x}{run}{roi}(eStill)];
                end
            end
            % Find scale events around locomotive bouts
            if periBout{x}(run).Nbout > 0
                periBout{x}(run).bIso = find(periBout{x}(run).iso(:,1) > minIso)';
                for b = periBout{x}(run).bIso %1:periBout{x}(run).Nbout
                    periBoutWindow = [periBout{x}(run).T{b}(periBout{x}(run).boutScan{b}(1)),  periBout{x}(run).T{b}(periBout{x}(run).boutScan{b}(end))] + windowPad;
                    eLoco = find(tempTstart >= periBoutWindow(1) & tempTstart <= periBoutWindow(2));
                    if ~isempty(eLoco),  LASE{x}{roi} = [LASE{x}{roi}, scaleEvent{x}{run}{roi}(eLoco)];  end
                end
            end
            % Find scale events around CSD bout
            if run == expt(x).csd
                csdWindow = [csdBout{x}(run).Tstart, csdBout{x}(run).Tstop] + csdPad;
                eCSD = find(tempTstart >= csdWindow(1) & tempTstart <= csdWindow(2));
                CASE{x}{roi} = [CASE{x}{roi}, scaleEvent{x}{run}{roi}(eCSD)];
            end
        end
    end
    Nsase{x} = cellfun(@numel, SASE{x});  Nlase{x} = cellfun(@numel, LASE{x}); Ncase{x} = cellfun(@numel, CASE{x});
    toc
end

%% Compare locomotion-associated scaling events pre/acute/post-CSD
xAcute = [5, 9, 28, 31, 37, 41]; % experiments where CSD wave was well-registered
SASEpre = cell(1,Nexpt); SASEpost = cell(1,Nexpt); LASEpre = cell(1,Nexpt); LASEpost = cell(1,Nexpt);
SASEpreAUC = cell(1,Nexpt); SASEpostAUC = cell(1,Nexpt); LASEpreAUC = cell(1,Nexpt); LASEpostAUC = cell(1,Nexpt); CASEauc = cell(1,Nexpt);
for x = xAcute
    % Define time intervals pre/acute/post-CSD
    preLims = [0, Tscan{x}{expt(x).csd}(1)];
    postLims = [csdBout{x}(expt(x).csd).Tstop+csdPad(2), csdBout{x}(expt(x).csd).Tstop+csdPad(2)+3600];    
    for r = 1:numel(SASE{x})
        if ~isempty(LASE{x}{r})
            tempTlase = [LASE{x}{r}.Tstart];
            LASEpre{x} = LASE{x}{r}(tempTlase < preLims(2));
            LASEpost{x} = LASE{x}{r}(tempTlase > postLims(1));
        end
        if ~isempty(SASE{x}{r})
            tempTsase = [SASE{x}{r}.Tstart];
            SASEpre{x} = SASE{x}{r}(tempTsase < preLims(2));
            SASEpost{x} = SASE{x}{r}(tempTsase > postLims(1));
        end


        if ~isempty(SASEpre{x}), SASEpreAUC{x} = [SASEpre{x}.magAUC]';  else,  SASEpreAUC{x} = NaN;  end
        if ~isempty(SASEpost{x}), SASEpostAUC{x} = [SASEpost{x}.magAUC]'; else, SASEpostAUC{x} = NaN; end
        if ~isempty(LASEpre{x}), LASEpreAUC{x} = [LASEpre{x}.magAUC]'; else, LASEpreAUC{x} = NaN; end
        if ~isempty(LASEpost{x}), LASEpostAUC{x} = [LASEpost{x}.magAUC]'; else, LASEpostAUC{x} = NaN; end
        if ~isempty(CASE{x}), CASEauc{x} = [CASE{x}{1}.magAUC]'; else, CASEauc{x} = NaN; end
        %sumMat = cell2padmat({SASEpreAUC{x}, LASEpreAUC{x}, CASEsum{x}, SASEpostAUC{x}, LASEpostAUC{x}});
        %{
        JitterPlot(sumMat, 0.5);
        set(gca,'Xtick',1:5, 'XtickLabel',{'Still Pre-CSD','Loco Pre-CSD','CSD','Still Post-CSD','Loco Post-CSD'});
        xtickangle(30);
        axis square;
        title(sprintf('x = %i',x))
        pause; 
        cla;
        %}
        %{
        [Xscale, Fscale] = ecdfCol( cell2padmat({SASEpreAUC{x}, SASEpostAUC{x}, LASEpreAUC{x}, LASEpostAUC{x}, CASEauc{x}}) );  
        plot( Xscale, Fscale ); colormap(distinguishable_colors(size(Xscale,2) ) )
        legend({'Still-Pre','Still-Post','Loco-Pre','Loco-Post','CSD'}, 'Location','SouthEast')
        title(sprintf('x = %i',x))
        pause;
        %}
    end
end

%normMat = [cellfun(@max, SASEpreAUC(xAcute))', cellfun(@max, LASEpreAUC(xAcute))', cellfun(@max, CASEauc(xAcute))', cellfun(@max, SASEpostAUC(xAcute))'];
%normMat = normMat./normMat(:,1);
%plot(1:4,  normMat)


[Xscale, Fscale] = ecdfCol( cell2padmat({vertcat(SASEpreAUC{xAcute}), vertcat(SASEpostAUC{xAcute}), vertcat(LASEpreAUC{xAcute}), vertcat(LASEpostAUC{xAcute}), vertcat(CASEauc{xAcute})}) );  
plot( Xscale, Fscale );
colormap(distinguishable_colors(size(Xscale,2) ) )
legend({'Still-Pre','Still-Post','Loco-Pre','Loco-Post','CSD'}, 'Location','SouthEast');
set(gca, 'Xscale','log', 'TickDir','out')
ylabel('Cumulative Distribution'); xlabel('Scale Event AUC (um*s)');

%%
close all;
figure;
for x = xAcute
    Tsase = ([SASE{x}{1}.Tstart]-csdBout{x}(expt(x).csd).Tstart)/60;
    Tlase = ([LASE{x}{1}.Tstart]-csdBout{x}(expt(x).csd).Tstart)/60;
    Tcase = ([CASE{x}{1}.Tstart]-csdBout{x}(expt(x).csd).Tstart)/60;
    cla;
    plot(Tcase, [CASE{x}{1}.magAUC], 'b.' ); hold on; 
    plot(Tcase, [CASE{x}{1}.magAUC], 'b' )
    plot(Tlase, [LASE{x}{1}.magAUC], 'go' ); %hold on; 
    plot(Tlase, [LASE{x}{1}.magAUC], 'g' )
    plot(Tsase, [SASE{x}{1}.magAUC], 'rx' ); % hold on; 
    plot(Tsase, [SASE{x}{1}.magAUC], 'r' )
    xlim([-70,70]);
    set(gca,'XTick',-60:20:60, 'box','off');
    xlabel('Peri-CSD time (min)'); 
    ylabel('Scale Event AUC (um*s)');
    pause;
end
