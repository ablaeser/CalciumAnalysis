%%
tuningDataPath = 'D:\MATLAB\Dura\tuningAnalysis_corrected.mat'; % _220404
if exist(tuningDataPath, 'file')
    fprintf('\nLoading %s', tuningDataPath); 
    load(tuningDataPath);
else
    dTheta = pi/6;  dPct = 5; testPct = [60,95]; Nboot = 500;  minNorm = 0.5; maxWidth = 60;
    scalePC = cell(1,Nexpt); PCangle = nan(1,Nexpt); scaleAngStats = cell(1,Nexpt); %testPctInd = nan(Nexpt,2);
    testStruct = repmat( struct('pct',testPct, 'pctInd',[], 'range',[], 'scan',[], 'dir',[], 'or',[], 'fluor',[], 'scale',[]), 1, Nexpt );
    rMech = cell(1,Nexpt); Nmech = nan(1,Nexpt);
    dirResponse = cell(1,Nexpt);  tuneResult = cell(1,Nexpt); Ntuned = nan(1,Nexpt); tunedFrac = nan(1,Nexpt); %Nor = nan(1,Nexpt); orFrac = nan(1,Nexpt);
    for x = x2D
        fprintf('\nx = %i: %s', x, expt(x).name)
        preScaleAP = vertcat(deform{x}(expt(x).preRuns).scaleAP);
        preScaleML = vertcat(deform{x}(expt(x).preRuns).scaleML);
        preScaleMag = vertcat(deform{x}(expt(x).preRuns).scaleMag);
        preScaleAngle = vertcat(deform{x}(expt(x).preRuns).scaleAngle);
        preFluor = [fluor{x}(expt(x).preRuns).z];
        preFluorZ = vertcat( preFluor.ROI );
        preState = vertcat( loco{x}(expt(x).preRuns).stateDown );
        [ ~, preScaleOr, ~, thetaBin, ~, ~, scaleAngStats{x} ] = AngularDist(preScaleAP, preScaleML, 'dPct', dPct, 'dThetaDir', dTheta, 'shift',dTheta/2, 'show',false ); % get scale mag percentiles within 30 deg angular bins
        % Find the principal component of scale and its angle (by convention, AP is considered the X-axis)
        tempPC = pca( [preScaleAP, preScaleML] );
        scalePC{x} = [tempPC(1,1), tempPC(2,1)];  % scalePC{x} = [1,1];
        [PCangle(x), ~] = cart2pol( scalePC{x}(1), scalePC{x}(2) );

        % Define test range: what magnitude of scale is consistently observed across all angles?
        [~,testStruct(x).pctInd(2)] = min( abs(scaleAngStats{x}.dir.pctRange - testPct(2)) );
        [~,testStruct(x).pctInd(1)] = min( abs(scaleAngStats{x}.dir.pctRange - testPct(1)) );
        testStruct(x).range = [min( scaleAngStats{x}.dir.pct(:,testStruct(x).pctInd(1)) ), min( scaleAngStats{x}.dir.pct(:,testStruct(x).pctInd(2)) )];  % testRange{x}
        testStruct(x).scan = find( preScaleMag >= testStruct(x).range(1) & preScaleMag <= testStruct(x).range(2) );
        testStruct(x).dir = preScaleAngle(testStruct(x).scan);
        testStruct(x).or = preScaleOr(testStruct(x).scan);
        testStruct(x).scale = preScaleMag(testStruct(x).scan);
        
        % DEFINE AND SUBTRACT BASELINE (STILL) FLUOR...
        baseInd = find( preScaleMag < testStruct(x).range(1) );
        testStruct(x).fluorBase = mean(preFluorZ(baseInd,:), 'omitnan');   %mean(preFluorZ(preState == 1,:), 'omitnan');
        preFluorZsub = preFluorZ - testStruct(x).fluorBase;
        %preFluorZsub(preFluorZsub < 0) = 0; % suppress negative values
        %tic
        
        % ...AND LOOK FOR DIRECTIONAL TUNING OF TEST RANGE DATA
        rMech{x} = sort([preCSD_2D_summary{x}.rDeform, preCSD_2D_summary{x}.rMixed], 'descend');
        Nmech(x) = numel(rMech{x});
        testStruct(x).fluor = nan(numel(testStruct(x).scan), Nmech(x)); % expt(x).Nroi
        k = 0;
        for r = rMech{x} % flip(1:expt(x).Nroi) %find(tempAngle > 0 & tempAngle < pi/2) % flip(1:expt(x).Nroi) %[preCSD_2D_summary{x}.rDeform, preCSD_2D_summary{x}.rMixed] %
            k = k+1;
            testStruct(x).fluor(:,k) = preFluorZsub(testStruct(x).scan, r);
            [ dirResponse{x}(k), ~, tuneResult{x}(k) ] = AngularTuning( thetaBin, testStruct(x).dir, testStruct(x).or, testStruct(x).fluor(:,k), Nboot, false ); % , testCI{x}(k)  orResponse{x}(k)
            tuneResult{x}(k).dir.ROI = r;
            tuneResult{x}(k).dir.tuned = tuneResult{x}(k).dir.angleCIwidth <= maxWidth && tuneResult{x}(k).dir.norm >= minNorm;
            tuneResult{x}(k).dir.tunePCdiff = AngularSeparation( tuneResult{x}(k).dir.angle, PCangle(x) );
            tuneResult{x}(k).dir.tuneElongDiff = AngularSeparation( tuneResult{x}(k).dir.angle, (pi/180)*ROI{x}(r).orientation, pi );
        end
        tempDir = [tuneResult{x}.dir];  
        Ntuned(x) = sum([tempDir.tuned]);
        tunedFrac(x) = Ntuned(x)/Nmech(x); %expt(x).Nroi;
    end
    fprintf('\nSaving %s', tuningDataPath)
    save(tuningDataPath, 'thetaBin','dTheta','dPct','testPct','Nboot','minNorm','maxWidth','scaleAngStats','scalePC','PCangle','rMech','Nmech','dirResponse','testStruct','tuneResult','Ntuned','tunedFrac'); % 'orResponse'
end


%% 
scaleBinLims = -3:0.1:3;
scaleBinTicks = linspace(1, numel(scaleBinLims), 7); % scaleBinLims(scaleBinTicks)
zeroTick = find(scaleBinLims ==0);
Tacute = 5;
close all;
figure('WindowState','maximized');
for x = x2Dcsd %xPresent
    % Post scale
    postScaleAP = vertcat(deform{x}(expt(x).preRuns).scaleAP);
    preScaleML = vertcat(deform{x}(expt(x).preRuns).scaleML);

    preScaleDensity = histcounts2(preScaleAP, preScaleML, scaleBinLims, scaleBinLims, 'normalization','probability'); %,'XBinLimits',scaleBinLims, 'YBinLimits',scaleBinLims,
    subplot(1,3,1); % figure%
    imagesc( flip(preScaleDensity, 1)); hold on;
    colorbar;
    line(scaleBinTicks([1,end]), zeroTick*[1,1],'color','k');
    line(zeroTick*[1,1], scaleBinTicks([1,end]), 'color','k');
    axis square;
    set(gca, 'Xtick',scaleBinTicks, 'XtickLabel',scaleBinLims(scaleBinTicks), 'Ytick',scaleBinTicks, 'YtickLabel',flip(scaleBinLims(scaleBinTicks)));
    xlabel('A-P Scaling (um)'); ylabel('M-L Scaling (um)');
    colormap(bluewhitered);
    title( sprintf('x = %i: Pre-CSD', x) );
    
    if ~isnan(expt(x).csd)
        postScaleAP = vertcat(deform{x}(expt(x).csd:expt(x).csd+1).scaleAP);
        postScaleML = vertcat(deform{x}(expt(x).csd:expt(x).csd+1).scaleML);
        % Break post-CSD data into acute and post-CSD periods
        acuteInd = find(Tscan{x}{expt(x).csd} - Tscan{x}{expt(x).csd}(1) <= Tacute*60);
        %postInd = find(Tscan{x}{expt(x).csd} - Tscan{x}{expt(x).csd}(1) > Tacute*60);
        acuteScaleAP = postScaleAP(acuteInd);
        acuteScaleML = postScaleML(acuteInd);
        postScaleAP(acuteInd) = [];
        postScaleML(acuteInd) = [];
        acuteScaleDensity = histcounts2(acuteScaleAP, acuteScaleML, scaleBinLims, scaleBinLims, 'normalization','probability'); %,'XBinLimits',scaleBinLims, 'YBinLimits',scaleBinLims,
        postScaleDensity = histcounts2(postScaleAP, postScaleML, scaleBinLims, scaleBinLims, 'normalization','probability'); %,'XBinLimits',scaleBinLims, 'YBinLimits',scaleBinLims,
        
        subplot(1,3,2);
        imagesc( flip(acuteScaleDensity, 1)); hold on;
        colorbar;
        line(scaleBinTicks([1,end]), zeroTick*[1,1],'color','k');
        line(zeroTick*[1,1], scaleBinTicks([1,end]), 'color','k');
        axis square;
        set(gca, 'Xtick',scaleBinTicks, 'XtickLabel',scaleBinLims(scaleBinTicks), 'Ytick',scaleBinTicks, 'YtickLabel',flip(scaleBinLims(scaleBinTicks)));
        xlabel('A-P Scaling (um)'); ylabel('M-L Scaling (um)');
        colormap(bluewhitered);
        title( sprintf('Acute CSD (T < %2.1f min)', Tacute) );
        
        subplot(1,3,3);
        imagesc( flip(postScaleDensity, 1)); hold on;
        colorbar;
        line(scaleBinTicks([1,end]), zeroTick*[1,1],'color','k');
        line(zeroTick*[1,1], scaleBinTicks([1,end]), 'color','k');
        axis square;
        set(gca, 'Xtick',scaleBinTicks, 'XtickLabel',scaleBinLims(scaleBinTicks), 'Ytick',scaleBinTicks, 'YtickLabel',flip(scaleBinLims(scaleBinTicks)));
        xlabel('A-P Scaling (um)'); ylabel('M-L Scaling (um)');
        colormap(bluewhitered);
        title( 'Post-CSD' );
    end
    pause;
end


%% Plot angular response to test range data for all ROI
Ndisp = 12;
Ncol = round(Ndisp/3); 
Nrow = round(Ndisp/4); 
opt = {[0.05,0.05], [0.06,0.06], [0.07, 0.04]};  % {[vert, horz], [bottom, top], [left, right] }
FS = 14;
close all; 
DirTuning_all = figure('WindowState','maximized', 'color','w');
for x = 43 % xPresent
    %figPath = sprintf('%sDirTuning_all_%s_%s_%i.pdf', figDir, metadata{x,1}.mouse, metadata{x,1}.date, metadata{x,1}.run); 
    tempResult = [tuneResult{x}.dir];
    [~,rSort] = sort( [tempResult.norm], 'descend', 'MissingPlacement','last' );
    for p = 1:ceil(expt(x).Nroi/Ndisp)
        rTemp = (1:Ndisp) + (p-1)*Ndisp; rTemp = rTemp( rTemp <= expt(x).Nroi );
        for c = 1:min(Ndisp, numel(rTemp))
            maxResp = max(dirResponse{x}(rSort(rTemp(c))).mean);
            subtightplot(Nrow, Ncol, c, opt{:} )
            polarplot( thetaBin.dir.cent, [dirResponse{x}(rSort(rTemp(c))).mean, dirResponse{x}(rSort(rTemp(c))).mean(1)], 'k', 'LineWidth',2 ); hold on;
            %polarplot( thetaBin.dir.cent, dirResponse{x}(rSort(rTemp(c))).mean-dirResponse{x}(rSort(rTemp(c))).sem, 'k--', 'LineWidth',1 ); 
            %polarplot( thetaBin.dir.cent, dirResponse{x}(rSort(rTemp(c))).mean+dirResponse{x}(rSort(rTemp(c))).sem, 'k:', 'LineWidth',1 ); 
            polarplot( tuneResult{x}(rSort(rTemp(c))).dir.angle*[1,1], maxResp*[0,1], 'r', 'LineWidth',1.5 );
            polarplot( tuneResult{x}(rSort(rTemp(c))).dir.angleCI(1)*[1,1], maxResp*[0,1], 'r--', 'LineWidth',1 );
            polarplot( tuneResult{x}(rSort(rTemp(c))).dir.angleCI(2)*[1,1], maxResp*[0,1], 'r--', 'LineWidth',1 );
            titleStr = sprintf('r = %i: DSI = %2.2f, width = %2.2f', rSort(rTemp(c)), tuneResult{x}(rSort(rTemp(c))).dir.norm, tuneResult{x}(rSort(rTemp(c))).dir.width ); % sprintf('r = %i: Direction sensitivity index = %1.3f. CI = [%1.3f, %1.3f]', rSort(rTemp(c)), tuneResult{x}(rSort(rTemp(c))).dir.norm, testCI{x}(rSort(rTemp(c))).dir.norm(1), testCI{x}(rSort(rTemp(c))).dir.norm(2) );
            title( titleStr );
            set(gca,'FontSize',FS, 'ThetaTick',[], 'Rtick',[]) % , 'ThetaTickLabels',{'', 'AP-Expansion','', 'AP-Compression'} , 'ThetaTick',[0, 90, 180, 270]
        end
        pause%(1); %pause; 
        %export_fig( figPath, '-pdf', '-painters','-q101', '-append', DirTuning_all );
        clf;
    end
end



%% Plot 

for x = 10
    for r = 1:expt(x).Nroi
        polarplot( thetaBin.dir.cent, [dirResponse{x}(r).mean, dirResponse{x}(r).mean(1)] ); hold on;
        polarplot( thetaBin.dir.cent(1:thetaBin.dir.Nbin), dirResponse{x}(r).mean, '.' );  %hold on;
        polarplot( tuneResult{x}(r).dir.angle*[1,1], tuneResult{x}(r).dir.mag*[0,1], 'color','k' );
        pause;
        cla;
    end
    cla;
end



