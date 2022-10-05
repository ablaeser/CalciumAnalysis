%% Calculate CSD cross-correlation within each fiber
for x = intersect(xPresent, find(Nfiber>0))   % 37 %
    if x == 37 
        csdOnsetWindow = [13, 27]; 
    else
        csdOnsetWindow = [-2, 12]; 
    end

    gaussFilt = MakeGaussFilt( 4, 0, 1/expt(x).scanRate, expt(x).scanRate, false );  % 8*(1/expt(x).scanRate)
    % Get pre-CSD fiber fluor signals
    tempPreFluorZ = fluor{x}(1:expt(x).csd-1).z;
    tempPreFluorZ = vertcat( tempPreFluorZ.ROI );
    % Get peri-CSD fiber fluor signals
    csdFluor =  fluor{x}(expt(x).csd).z.ROI; % figure; imagesc(csdFluor')
    for r = find(sum( isnan(csdFluor), 1)) % 1:expt(x).Nroi
        firstGoodScan = find(~isnan(csdFluor(:,r)), 1, 'first');
        nanScan = find(isnan(csdFluor(:,r)) );
        nanScan(nanScan < firstGoodScan) = [];
        goodScan = 1:size(csdFluor,1); goodScan(nanScan) = [];
        preRepairFluor = csdFluor(:,r);
        repairFluor = preRepairFluor;
        repairFluor(nanScan) = interp1( goodScan, preRepairFluor(goodScan), nanScan, 'spline' );
        csdFluor(:,r) = repairFluor;
        %plot( 1:size(csdFluor,1), [repairFluor, preRepairFluor ]);
        %pause;
        %cla;
    end   
    csdFluor(find(~isnan(sum(csdFluor,2))),:) = filtfilt(gaussFilt, 1, csdFluor(find(~isnan(sum(csdFluor,2))),:)); %[];
    %csdFluor(find(isnan(sum(csdFluor,2))),:) = []; % remove missing data scans
    %csdFluor = filtfilt(gaussFilt, 1, csdFluor);
    Tperi = Tscan{x}{expt(x).csd} - csdBout{x}(expt(x).csd).Tstart; % Tcat
    if isempty( csdOnsetWindow ), csdOnsetWindow = [0, 0.67*csdBout{x}(expt(x).csd).dur(1)]; end
    csdWaveScans = find(Tperi >= csdOnsetWindow(1) & Tperi <= csdOnsetWindow(2))';   
    
    csdFluor = csdFluor(csdWaveScans,:); % fluor{x}(expt(x).csd).dFF.ROI(csdScans, :);
    %{
    tempCent = vertcat(ROI{x}.cent);
    [~,rSortCent] = sort(tempCent(:,1), 'ascend');
    imagesc( csdFluor(:,rSortCent)' );
    title( sprintf('%s', expt(x).name), 'Interpreter','none');
    pause;
    cla;
    %}
    % {
    for f = 1:Nfiber(x) 
        % Get pre-CSD fiber fluor signals
        fiber{x}(f).spont.fluor = tempPreFluorZ(:,fiber{x}(f).ROI); %fiberFluor = fluor{x}(expt(x).csd).z.ROI(onsetWindowScans, fiber{x}(f).ROI); % csdBout{x}(expt(x).csd).scan{1}
        %fiber{x}(f).spont.fluor(find(isnan(sum(fiber{x}(f).spont.fluor,2))),:) = []; % remove missing data scans
        %fiber{x}(f).spont.fluor = filtfilt(gaussFilt, 1, fiber{x}(f).spont.fluor);        
        % peri-CSD
        fiber{x}(f).CSD.fluor = csdFluor(:, fiber{x}(f).ROI); % csdBout(expt.csd).scan{1}
        %fiber{x}(f).CSD.fluor(find(isnan(sum(fiber{x}(f).CSD.fluor,2))),:) = []; % remove missing data scans
        %fiber{x}(f).CSD.fluor = filtfilt(gaussFilt, 1, fiber{x}(f).CSD.fluor);
    end
    fiber{x} = GetFiberWaveSpeed(expt(x), ROI{x}, fiber{x}, 'CSD', 'minSep',max([20, csdWave{x}.speed/expt(x).scanRate]), 'show',true);
    fiber{x} = GetFiberWaveSpeed(expt(x), ROI{x}, fiber{x}, 'spont', 'minSep',max([20, csdWave{x}.speed/expt(x).scanRate]), 'minXC',0.5 , 'show',true); % true   
    %}
end


%%

angVsLPS = [];  k = 0;
%close all; figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
for x = intersect(xPresent, find(Nfiber>0)) %28 
    csdLagPerSep = 1000*(1/csdWave{x}.speed);
    %{
    subplot(3,2,2);
    polarplot( csdWave{x}.angle*[1,1], csdWave{x}.speed*[0,1], 'LineWidth',2 ); hold on;
    title('CSD Wave');
    %}
    for f = 1:Nfiber(x)
        k = k+1;
        tempAngSep = AngularSeparation(fiber{x}(f).orientation+pi, csdWave{x}.angle+pi, pi, false);
        angVsLPS(k,1) = tempAngSep;
        angVsLPS(k,2) = fiber{x}(f).CSD.medLagPerSep/csdLagPerSep;
        angVsLPS(k,3) = fiber{x}(f).spont.medLagPerSep/csdLagPerSep;
        %{
        subplot(3,2,1); cla;
        imshow( label2rgb(fiber{x}(f).labelFoot)); hold on;
        for r = fiber{x}(f).ROI
            text( ROI{x}(r).cent(1), ROI{x}(r).cent(2), sprintf('%i', r), 'HorizontalAlignment','center', 'FontSize',8 );
        end
        %DrawArrow(expt.Ncol/2, expt.Nrow/2, csdWave.velocity(1)/expt.scanRate, csdWave.velocity(2)/expt.scanRate); % Draw an arrow representing the CSD velocity      
        title( sprintf('%s, fiber %i', expt(x).name, f), 'Interpreter','none' ); 
        
        subplot(3,2,3); cla;
        for p = 1:numel(fiber{x}(f).CSD.pairAngle)
            polarplot( fiber{x}(f).CSD.pairAngle(p)*[1,1], fiber{x}(f).CSD.pairLagPerSep(p)/csdLagPerSep*[0,1], 'LineWidth',1 ); hold on;
        end
        title('Paired CSD Delays (normalized to CSD wave)');
        rlim([0,1.2]);

        subplot(3,2,4); cla;
        polarplot( tempAngSep*[1,1], fiber{x}(f).CSD.medLagPerSep/csdLagPerSep*[0,1], 'LineWidth',1 ); %fiber{x}(f).orientation
        title('Median CSD Delay per Separation (normalized, relative to CSD wave)');
        rlim([0,1.2]);
        
        subplot(3,2,5); cla;
        polarplot(0,NaN); hold on;
        for p = 1:numel(fiber{x}(f).spont.pairAngle)
            polarplot( fiber{x}(f).spont.pairAngle(p)*[1,1], fiber{x}(f).spont.pairLagPerSep(p)/csdLagPerSep*[0,1], 'LineWidth',1 ); hold on;
        end
        title('Paired Spontaneous Delays (normalized to CSD wave)');
        rlim([0,1.2]);
            
        subplot(3,2,6); cla;
        polarplot(0,NaN); hold on;
        polarplot( tempAngSep*[1,1], fiber{x}(f).spont.medLagPerSep/csdLagPerSep*[0,1], 'LineWidth',1 ); %fiber{x}(f).orientation
        title('Median Spontaneous Delay per Separation (normalized, relative to CSD wave)');
        rlim([0,1.2]);

        pause;
        %}
    end
end

%% Bin by angle and calculate mean lag per separation by fiber orientation, relative to CSD wave speed and direction
dTheta = pi/12;
thetaLim = 0:dTheta:pi/2;
NthetaBin = numel(thetaLim)-1;
LPS_CSD = cell(1,NthetaBin); LPS_spont = cell(1,NthetaBin); thetaBinCent = zeros(1,NthetaBin);
for t = 1:NthetaBin
    binInd = find( angVsLPS(:,1) >= thetaLim(t) & angVsLPS(:,1) < thetaLim(t+1) );
    thetaBinCent(t) = (thetaLim(t)+thetaLim(t+1))/2;
    LPS_CSD{t} = angVsLPS(binInd,2);
    LPS_spont{t} = angVsLPS(binInd,3);
end
cellfun(@nanmean, LPS_CSD)
% {
close all; figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
%subplot(1,2,1);
h(1) = polarplot(thetaBinCent, cellfun(@nanmean, LPS_CSD) ); hold on;
h(2) = polarplot(thetaBinCent, cellfun(@nanmean, LPS_spont) );
thetalim([0,90]);
rlim([0,1.1]);
set(gca,'ThetaTick', rad2deg(thetaLim), 'Rtick',[0,0.5,1]);
legend(h, {'CSD','Spontaneous'});
tempAx = get(gca);
tempAx.RAxis.Label.String = 'Mean Normalized Delay (Fraction of CSD Wave)';
tempAx.ThetaAxis.Label.String = 'Alignment with CSD Wave';
pause(1);
tempAx.ThetaAxis.Label.Rotation = -45;
tempAx.ThetaAxis.Label.Position(2) = 1.15;
tempAx.ThetaAxis.Label.Position(1) = 45;

