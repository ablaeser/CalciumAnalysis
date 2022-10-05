figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'color','w');
for x = x3Dcsd %xPresent
    preScans = 1:expt(x).scanLims(expt(x).csd); % Pre-CSD scans
    csdScans = (expt(x).scanLims(expt(x).csd) + csdBout{x}(expt(x).csd).scan{1}(1)):(expt(x).scanLims(expt(x).csd) + csdBout{x}(expt(x).csd).scan{1}(end));
    shortScans = (expt(x).scanLims(expt(x).csd) + csdBout{x}(expt(x).csd).scan{1}(end)+100):expt(x).scanLims(expt(x).csd+1); % Scans just after CSD
    longScans = expt(x).scanLims(expt(x).csd+1)+1:expt(x).scanLims(end); % Scans long after CSD
    epochColor = distinguishable_colors(4);
    for v = 5%1:NdefVars
        if isempty(strfind(defVars{v}, 'shift'))
            tempPlanes = 1:expt(x).Nplane;   
        else
            tempPlanes = 4:expt(x).Nplane-3;  
        end
        deepPlanes = intersect( 1:floor(expt(x).Nplane/2), tempPlanes );
        shallowPlanes = intersect( ceil(expt(x).Nplane/2):expt(x).Nplane, tempPlanes );

        tempDefCat = vertcat( deform{x}.(defVars{v}) ); %deformCat{x}.(defVars{v}); %
        tempDefStd = std(tempDefCat(:,tempPlanes), 0, 2, 'omitnan');
        tempDeepStd = std(tempDefCat(:,deepPlanes), 0, 2, 'omitnan');
        tempShallowStd = std(tempDefCat(:,shallowPlanes), 0, 2, 'omitnan');

        tempMeanPreDef = mean(tempDefCat(preScans,tempPlanes), 1, 'omitnan');
        tempMeanCsdDef = mean(tempDefCat(csdScans,tempPlanes), 1, 'omitnan');
        tempMeanShortDef = mean(tempDefCat(shortScans,tempPlanes), 1, 'omitnan');
        tempMeanLongDef = mean(tempDefCat(longScans,tempPlanes), 1, 'omitnan');        
        [~, tempMaxPlane] = max( tempDefCat(:,tempPlanes), [], 2 );
        tempMaxPlane = tempMaxPlane + tempPlanes(1) - 1;
        tempMaxInd = sub2ind( size(tempDefCat), (1:size(tempDefCat,1))', tempMaxPlane);
        tempMaxDef = tempDefCat(tempMaxInd);
        [~, tempMinPlane] = min( tempDefCat(:,tempPlanes), [], 2 );
        tempMinPlane = tempMinPlane + tempPlanes(1) - 1;
        tempMinInd = sub2ind( size(tempDefCat), (1:size(tempDefCat,1))', tempMinPlane);
        tempMinDef = tempDefCat(tempMinInd);

        sp(1) = subplot(4,2,1);
        imagesc( tempDefCat' );
        set(gca,'Xtick', expt(x).scanLims, 'XtickLabel',[], 'box','off', 'TickDir','out', 'Ytick',tempPlanes([1,end]) ); % 
        title( sprintf('%s: %s', expt(x).name, defVars{v} ), 'Interpreter','none' );
        ylabel('Plane');
        
        sp(2) = subplot(4,2,3);
        plot( [tempDefStd, tempShallowStd, tempDeepStd]  ); hold on;
        Ytemp = get(gca,'Ylim');  Ytemp = Ytemp(2); % - 0.05;
        line( preScans([1,end]), Ytemp*[1,1], 'color',epochColor(1,:), 'LineWidth',2 )
        line( csdScans([1,end]), Ytemp*[1,1], 'color',epochColor(2,:), 'LineWidth',2 )
        line( shortScans([1,end]), Ytemp*[1,1], 'color',epochColor(3,:), 'LineWidth',2 )
        line( longScans([1,end]), Ytemp*[1,1], 'color',epochColor(4,:), 'LineWidth',2 )        
        set(gca,'Xtick', expt(x).scanLims, 'XtickLabel',[] , 'box','off', 'TickDir','out' );
        ylabel('Std Dev'); 
        
        sp(3) = subplot(4,2,5);
        plot( [tempMaxDef, tempMinDef] ); hold on;
        ylabel('Max and Min'); 
        
        sp(4) = subplot(4,2,7);
        plot( [tempMaxPlane, tempMinPlane] ); hold on;
        ylabel('Plane'); 
        
        %sp(3) = subplot(3,2,5);
        %plot( vertcat(loco{x}.speedDown) );
        %set(gca,'Xtick', expt(x).scanLims, 'box','off', 'TickDir','out' );
        %ylabel('Locomotion State'); xlabel('Volume Scan');

        linkaxes(sp, 'x');
        xlim([-Inf,Inf]);
        impixelinfo;
        
        subplot(4,2,[2,4]);
        plot(tempPlanes, [tempMeanPreDef', tempMeanCsdDef', tempMeanShortDef', tempMeanLongDef'] );
        axis square;
        colororder( gca, epochColor );
        legend('Pre-CSD', 'CSD', 'Just after', 'Long after' )
        xlabel('Plane'); ylabel( defVars{v} );
        
        subplot(4,2,[6,8]);
        plot( tempShallowStd, tempDeepStd, '.' ); hold on;
        Ytemp = get(gca, 'Ylim');
        line([0,Ytemp], [0,Ytemp], 'color','k', 'lineStyle','--'); 
        xlabel('Shallow Planes Std Dev'); ylabel('Deep Planes Std Dev');
        axis square
        
        pause;
        clf;
    end
end

%%
close all;
figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'color','w');
tiledlayout(NdefVars, 1, 'TileSpacing','compact')
for x = 5 %x3Dcsd %xPresent
    preScans = 1:expt(x).scanLims(expt(x).csd); % Pre-CSD scans
    csdScans = (expt(x).scanLims(expt(x).csd) + csdBout{x}(expt(x).csd).scan{1}(1)):(expt(x).scanLims(expt(x).csd) + csdBout{x}(expt(x).csd).scan{1}(end));
    shortScans = (expt(x).scanLims(expt(x).csd) + csdBout{x}(expt(x).csd).scan{1}(end)+100):expt(x).scanLims(expt(x).csd+1); % Scans just after CSD
    longScans = expt(x).scanLims(expt(x).csd+1)+1:expt(x).scanLims(end); % Scans long after CSD
    epochColor = distinguishable_colors(4);
    for v = 1:NdefVars
        tempPlanes = 4:expt(x).Nplane-3;  
        tempDefCat = deformCat{x}.(defVars{v});
        
        nexttile();
        imagesc( tempDefCat' ); hold on;
        set(gca,'Xtick', expt(x).scanLims, 'XtickLabel',[], 'box','off', 'TickDir','out', 'Ytick',tempPlanes([1,end]) ); hold on;
        if v == 1
            Yline = -0.5*[1,1];
            line( preScans([1,end]), Yline, 'color',epochColor(1,:), 'LineWidth',2 )
            line( csdScans([1,end]), Yline, 'color',epochColor(2,:), 'LineWidth',2 )
            line( shortScans([1,end]), Yline, 'color',epochColor(3,:), 'LineWidth',2 )
            line( longScans([1,end]), Yline, 'color',epochColor(4,:), 'LineWidth',2 )  
            title( sprintf('%s: %s', expt(x).name, defVars{v} ), 'Interpreter','none' );
            ylim([Yline(1), Inf]);
            %break
        else
            title( sprintf('%s', defVars{v} ), 'Interpreter','none' );
        end
        ylabel('Plane');
        hold off;
        %pause; 
    end
end
%%
defVars = {'transAP', 'transML', 'scaleAP', 'scaleML', 'stretchAP', 'stretchML', 'shearAP', 'shearML', 'shiftZ'}; %, 'DshiftZ'
NdefVars = numel( defVars ); %, 'dShiftZ'
close all;
figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'color','w');
N3Dcsd= numel(x3Dcsd);
tiledlayout(NdefVars, N3Dcsd, 'TileSpacing','compact', 'Padding','compact')
epochColor = distinguishable_colors(4);
k = 0;
for x = x3Dcsd %xPresent
    preScans = 1:expt(x).scanLims(expt(x).csd); % Pre-CSD scans
    csdScans = (expt(x).scanLims(expt(x).csd) + csdBout{x}(expt(x).csd).scan{1}(1)):(expt(x).scanLims(expt(x).csd) + csdBout{x}(expt(x).csd).scan{1}(end));
    shortScans = (expt(x).scanLims(expt(x).csd) + csdBout{x}(expt(x).csd).scan{1}(end)+100):expt(x).scanLims(expt(x).csd+1); % Scans just after CSD
    longScans = expt(x).scanLims(expt(x).csd+1)+1:expt(x).scanLims(end); % Scans long after CSD
    k = k+1;
    for v = 1:NdefVars
        tempPlanes = 4:expt(x).Nplane-3;  
        deepPlanes = intersect( 1:floor(expt(x).Nplane/2), tempPlanes );
        shallowPlanes = intersect( ceil(expt(x).Nplane/2):expt(x).Nplane, tempPlanes );

        tempDefCat = deformCat{x}.(defVars{v}); %vertcat( deform{x}.(defVars{v}) );
        tempDefStd = std(tempDefCat(:,tempPlanes), 0, 2, 'omitnan');
        tempDeepStd = std(tempDefCat(:,deepPlanes), 0, 2, 'omitnan');
        tempShallowStd = std(tempDefCat(:,shallowPlanes), 0, 2, 'omitnan');

        tempMedPreDef = median(tempDefCat(preScans,tempPlanes), 1, 'omitnan');
        tempMedCsdDef = median(tempDefCat(csdScans,tempPlanes), 1, 'omitnan');
        tempMedShortDef = median(tempDefCat(shortScans,tempPlanes), 1, 'omitnan');
        tempMedLongDef = median(tempDefCat(longScans,tempPlanes), 1, 'omitnan');        
 
         
        nexttile((v-1)*N3Dcsd+k); %subtightplot(Nrow, NdefVars, ) % v+(k-1)*N3Dcsd
        plot(tempPlanes, [tempMedPreDef', tempMedCsdDef', tempMedShortDef', tempMedLongDef'] );
        xlim([1,expt(x).Nplane])
        axis square;
        colororder( gca, epochColor );
        if v == NdefVars
            xlabel('Plane'); 
            set(gca, 'Xtick',[1,expt(x).Nplane]);
        else
            set(gca, 'Xtick',[]);
        end
        if k == 1, ylabel( defVars{v} ); end
        if v == 1, title( expt(x).name, 'Interpreter','none', 'FontSize',8 ); end
        if v == 1 && k == N3Dcsd, legend('Pre-CSD', 'CSD', 'Just after', 'Long after', 'Location','EastOutside' ); end
        %pause
    end
end


%% Plot standard deviance of translation, scale, shear as a function of 
k = 0;
defConst = zeros(N3Dcsd, NdefVars); defSlope = zeros(N3Dcsd, NdefVars);
figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'color','w');
tiledlayout(NdefVars, N3Dcsd, 'TileSpacing','compact', 'Padding','compact')
for x = x3Dcsd
    k = k+1;
    if ~isnan(expt(x).csd),  preCSDruns = 1:expt(x).csd-1;  else,  preCSDruns = 1:expt(x).Nruns;  end
    stillScanCat = find( vertcat( loco{x}(preCSDruns).stateDown ) == 1 );
    tempPlanes = 4:expt(x).Nplane-3;  
    for v = 1:NdefVars
        tempDefCat = deformCat{x}.(defVars{v});
        tempMedPreDef = median(tempDefCat(stillScanCat,tempPlanes), 1, 'omitnan');
        tempFit = fitlm( tempPlanes', tempMedPreDef', 'RobustOpts',true );
        if tempFit.Coefficients.pValue(1) < 0.05
            defConst(k,v) = tempFit.Coefficients.Estimate(1);
        end
        if tempFit.Coefficients.pValue(2) < 0.05
            defSlope(k,v) = tempFit.Coefficients.Estimate(2);
        end
        nexttile((v-1)*N3Dcsd+k);
        plot( tempPlanes, tempMedPreDef, 'b'); hold on;
        plot( tempPlanes, polyval([defSlope(k,v), defConst(k,v)], tempPlanes'), 'r--');
        if v == NdefVars
            xlabel('Plane'); 
            set(gca, 'Xtick',[1,expt(x).Nplane]);
        else
            set(gca, 'Xtick',[]);
        end
        if k == 1, ylabel( defVars{v} ); end
        if v == 1, title( expt(x).name, 'Interpreter','none', 'FontSize',8 ); end
        xlim([1,expt(x).Nplane]);
        %if v == 1 && k == N3Dcsd, legend('Pre-CSD', 'CSD', 'Just after', 'Long after', 'Location','EastOutside' ); end
    end
    %{
    transCat = sqrt(deformCat{x}.transAP.^2 + deformCat{x}.transML.^2); % vertcat( deform{x}(preCSDruns).transMag );
    scaleCat = sqrt(deformCat{x}.scaleAP.^2 + deformCat{x}.scaleML.^2);  %vertcat( deform{x}(preCSDruns).scaleMag );
    shearCat = sqrt(deformCat{x}.shearAP.^2 + deformCat{x}.shearML.^2);  %vertcat( deform{x}(preCSDruns).shearMag );
    relDepth = (1:expt(x).Nplane)./expt(x).Nplane;

    transMean = mean(transCat(stillScanCat,:), 1, 'omitnan');
    scaleMean = mean(scaleCat(stillScanCat,:), 1, 'omitnan');
    shearMean = mean(shearCat(stillScanCat,:), 1, 'omitnan');
    transStd = std( transCat(stillScanCat,:), 0, 1, 'omitnan'); % mad(transCat(stillScanCat,:), 1, 1); %
    scaleStd = std( scaleCat(stillScanCat,:), 0, 1, 'omitnan'); % mad(scaleCat(stillScanCat,:), 1, 1); %
    shearStd = std( shearCat(stillScanCat,:), 0, 1, 'omitnan'); % mad(shearCat(stillScanCat,:), 1, 1); %
    badPlane = unique([find(relDepth<=2/30), find( sum(isnan(transCat), 1)/size(transCat,1) > 0.1), find(relDepth>=29/30)]);
    goodPlane = 4:expt(x).Nplane-3;
 
    transDepthMdl = fitlm( relDepth(goodPlane), transMean(goodPlane), 'RobustOpts',true ); 
    transDepthMdlCoeff = transDepthMdl.Coefficients;
    transDepthSlope(x) = transDepthMdlCoeff.Estimate(2);
    scaleDepthMdl = fitlm( relDepth(goodPlane), scaleMean(goodPlane), 'RobustOpts',true ); 
    scaleDepthMdlCoeff = scaleDepthMdl.Coefficients;
    scaleDepthSlope(x) = scaleDepthMdlCoeff.Estimate(2);
    shearDepthMdl = fitlm( relDepth(goodPlane), shearMean(goodPlane), 'RobustOpts',true ); 
    shearDepthMdlCoeff = shearDepthMdl.Coefficients;
    shearDepthSlope(x) = shearDepthMdlCoeff.Estimate(2);
    
    
    subplot(1,3,1);
    errorbar( relDepth, transMean, transStd, 'Marker','.'); hold on;
    plot( relDepth(goodPlane)', predict(transDepthMdl,relDepth(goodPlane)') ) 
    axis square;
    set(gca,'Xtick', [0,1], 'XtickLabel',{'Deep','Superficial'});
    xlim([-0.05, 1.05]); ylim([-Inf,Inf]);
    xlabel('Depth'); ylabel('Translation'); title(sprintf('x = %i: %s',x, expt(x).name),'Interpreter','none'); % %ylabel('MAD'); 
    
    subplot(1,3,2);
    errorbar( relDepth, scaleMean, scaleStd, 'Marker','.'); hold on;
    plot( relDepth(goodPlane)', predict(scaleDepthMdl,relDepth(goodPlane)') ) 
    axis square;
    set(gca,'Xtick', [0,1], 'XtickLabel',{'Deep','Superficial'});    
    xlim([0, 1.05]); ylim([-Inf,Inf]);
    xlabel('Depth'); ylabel('Scale'); %ylabel('MAD'); title(sprintf('Scale'));
    
    subplot(1,3,3);
    errorbar( relDepth, shearMean, shearStd, 'Marker','.'); hold on;
    plot( relDepth(goodPlane)', predict(shearDepthMdl,relDepth(goodPlane)') ) 
    axis square;
    set(gca,'Xtick', [0,1], 'XtickLabel',{'Deep','Superficial'});
    xlim([0, 1.05]); ylim([-Inf,Inf]);
    xlabel('Depth'); ylabel('Shear'); % ylabel('MAD'); title(sprintf('Shear'));
    %}
    %pause;
    %clf;
end
%%
figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'color','w');
tiledlayout(1,9, 'TileSpacing','compact')
for v = 1:NdefVars
    nexttile()
    bar(defSlope(:,v));
    title( defVars{v} );
    ylabel('Slope (\Deltadeformation per plane)'); 
    set(gca,'Xtick',1:N3Dcsd, 'XtickLabel',{expt(x3Dcsd).name}, 'TickLabelInterpreter','none');
    axis square;
    xtickangle(60);
    xlim([0,N3Dcsd+1]);
end

%%
figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'color','w');
subplot(1,3,1);
bar( transDepthSlope(x3Dcsd) ) 
axis square;
ylabel('Slope (mean |deformation|/depth)');
title('Translation');
set(gca,'Xtick',1:numel(x3Dcsd), 'XtickLabel',{expt(x3Dcsd).name}, 'TickLabelInterpreter','none') %
xtickangle(45);
xlim([0,numel(x3Dcsd)+1]);

subplot(1,3,2);
bar( scaleDepthSlope(x3Dcsd) )
axis square;
%ylabel('Slope (mean deform magnitude/plane)');
title('Scale');
set(gca,'Xtick',1:numel(x3Dcsd), 'XtickLabel',{expt(x3Dcsd).name}, 'TickLabelInterpreter','none') %
xtickangle(45);
xlim([0,numel(x3Dcsd)+1]);

subplot(1,3,3);
bar( shearDepthSlope(x3Dcsd) );
axis square;
title('Shear');
set(gca,'Xtick',1:numel(x3Dcsd), 'XtickLabel',{expt(x3Dcsd).name}, 'TickLabelInterpreter','none') %
xtickangle(45);
xlim([0,numel(x3Dcsd)+1]);
%JitterPlot( [transDepthSlope', scaleDepthSlope', shearDepthSlope'], 0.45)

%% Get deformation by 5 minute increments
deformInt = 5; % minutes

Tint = cell(1,Nexpt); scaleAPint = cell(1,Nexpt); scaleMLint = cell(1,Nexpt); Tcsd = nan(1,Nexpt);
for x = x3Dcsd
    k = 0;
    APpix = expt(x).Ncol-segParams{x}.edges(1)-segParams{x}.edges(2);
    MLpix = expt(x).Nrow-segParams{x}.edges(3)-segParams{x}.edges(4);
    NintScan = round(60*deformInt*expt(x).scanRate);
    Tcat = vertcat( Tscan{x}{:} );
    stateCat = vertcat( loco{x}.stateDown );
    for run = 1:expt(x).Nruns
        firstScan = 1;
        for i = 1:ceil(expt(x).Nscan(run)/NintScan)
            k = k+1;
            currRunScans = firstScan:firstScan+NintScan-1;
            currRunScans(currRunScans > expt(x).Nscan(run)) = [];
            currRunScans(loco{x}(run).stateDown(currRunScans) == 2) = [];
            Tcurr = Tscan{x}{run}(currRunScans);
            Tint{x}(k,1) = mean(Tcurr)/60;
            [~,catScans] = ismember(Tcurr, Tcat);
            scaleAPcurr = deformCat{x}.scaleAP(catScans,:);
            scaleAPcurr(scaleAPcurr < deformLim.scale(1) | scaleAPcurr > deformLim.scale(2)) = NaN;
            scaleAPcurr = expt(x).umPerPixel*APpix./scaleAPcurr - expt(x).umPerPixel*APpix;
            scaleAPint{x}(k,:) = median(scaleAPcurr, 1, 'omitnan');

            scaleMLcurr = deformCat{x}.scaleML(catScans,:);
            scaleMLcurr(scaleMLcurr < deformLim.scale(1) | scaleMLcurr > deformLim.scale(2)) = NaN;
            scaleMLcurr = expt(x).umPerPixel*MLpix./scaleMLcurr - expt(x).umPerPixel*MLpix; %deformCat{x}.scaleML(catScans,:);
            scaleMLint{x}(k,:) = median(scaleMLcurr, 1, 'omitnan');
            firstScan = firstScan + NintScan;
        end
    end
end
%%
figure;
for x = x3Dcsd
    cla
    Tcsd = median(csdBout{x}(expt(x).csd).T{1}/60);
    tempPlanes = 4:expt(x).Nplane-3; 
    plot(Tint{x}-Tcsd,  mean(scaleAPint{x}(:,tempPlanes), 2,'omitnan'), 'b'); hold on;
    h(1) = plot(Tint{x}-Tcsd,  mean(scaleAPint{x}(:,tempPlanes), 2,'omitnan'), 'b.'); hold on;
    plot(Tint{x}-Tcsd,  mean(scaleMLint{x}(:,tempPlanes), 2,'omitnan'), 'r');
    h(2) = plot(Tint{x}-Tcsd,  mean(scaleMLint{x}(:,tempPlanes), 2,'omitnan'), 'r.');
    Ylims = get(gca, 'Ylim');
    h(3) = line(0*[1,1], Ylims, 'color','k', 'lineStyle','--');
    
    legend(h, {'AP','ML', 'CSD'}, 'Location','northeastoutside'); 
    xlabel('Peri-CSD Time (min)'); ylabel('Scale (um)'); title(sprintf('%s', expt(x).name), 'interpreter','none');
    
    pause;
end

