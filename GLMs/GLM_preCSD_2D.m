%% Use GLM to assess contribution of different variables
preCSD_2D_pred = cell(1,Nexpt); preCSD_2D_resp = cell(1,Nexpt); preCSD_2D_opts = cell(1,Nexpt); preCSD_2D_result = cell(1,Nexpt); preCSD_2D_summary = cell(1,Nexpt);
%figDir = 'D:\MATLAB\LevyLab\Figures\3D\GLM\'; mkdir( figDir )
exptColor = distinguishable_colors(Npresent);
GLMname = 'preCSD_2D_fluor_TScShDtStDshZVAL'; % XC
GLMrate = 15.49/3;
for x = x2D %xPresent %intersect( find(cellfun(@isempty, preCSD_2D_result)), xPresent ) %x3Dcsd  (xPresent>21)
    % GLMparallel options
    preCSD_2D_opts{x}.name = sprintf('%s_%s', expt{x}.name, GLMname); %strcat(expt{x}.name, , '_preCSDglm');
    %preCSD_2D_opts{x}.show = true;
    preCSD_2D_opts{x}.rShow = NaN;
    preCSD_2D_opts{x}.figDir = ''; % figDir;
    preCSD_2D_opts{x}.alpha = 0.01;  % The regularization parameter, default is 0.01
    preCSD_2D_opts{x}.standardize = true; 
    preCSD_2D_opts{x}.trainFrac = 0.75; % 1; %
    preCSD_2D_opts{x}.Ncycle = 20;
    preCSD_2D_opts{x}.distribution = 'gaussian'; % 'poisson'; %  
    preCSD_2D_opts{x}.CVfold = 10;
    preCSD.opts(x).nlamda = 1000;
    preCSD.opts(x).maxit = 5*10^5;
    preCSD_2D_opts{x}.minDev = 0.05; 
    preCSD_2D_opts{x}.minDevFrac = 0.1;
    preCSD_2D_opts{x}.maxP = 0.05;
    preCSD_2D_opts{x}.Nshuff = 0; % Nshuff = preCSD_2D_opts{x}.Nshuff; 
    preCSD_2D_opts{x}.minShuff = 15; 
    preCSD_2D_opts{x}.window = [-4,4]; % [0,0]; % [-0.5, 0.5]; % 
    preCSD_2D_opts{x}.lopo = true; %false; %
    preCSD_2D_opts{x}.frameRate = GLMrate; 
    preCSD_2D_opts{x}.binSize = expt{x}.scanRate/GLMrate;
    preCSD_2D_opts{x}.minShuffFrame = round( preCSD_2D_opts{x}.frameRate*preCSD_2D_opts{x}.minShuff );
    windowFrame = round(preCSD_2D_opts{x}.window*preCSD_2D_opts{x}.frameRate);
    preCSD_2D_opts{x}.shiftFrame = windowFrame(1):windowFrame(2);
    preCSD_2D_opts{x}.maxShift = max( abs(windowFrame) );
    preCSD_2D_opts{x}.Nshift = numel( preCSD_2D_opts{x}.shiftFrame );  %Nshift = preCSDOpts(x).Nshift;
    preCSD_2D_opts{x}.lags = preCSD_2D_opts{x}.shiftFrame/preCSD_2D_opts{x}.frameRate;
    preCSD_2D_opts{x}.xVar = 'Time';
    % Define predictors (Concatenate pre-CSD data)
    %{
    tempTrans = BinDownMean( mean( vertcat(deform{x}(expt{x}.preRuns).transMag), 2, 'omitnan'), preCSD_2D_opts{x}.binSize );
    tempTransSpd = BinDownMean( mean( vertcat(deform{x}(expt{x}.preRuns).DtransMag), 2, 'omitnan'), preCSD_2D_opts{x}.binSize );
    tempScaleMag = BinDownMean( mean( vertcat(deform{x}(expt{x}.preRuns).scaleMag), 2, 'omitnan'), preCSD_2D_opts{x}.binSize );
    tempStretchMag = BinDownMean( mean( vertcat(deform{x}(expt{x}.preRuns).stretchMag), 2, 'omitnan'), preCSD_2D_opts{x}.binSize );
    tempShearMag = BinDownMean( mean( vertcat(deform{x}(expt{x}.preRuns).shearMag), 2, 'omitnan'), preCSD_2D_opts{x}.binSize );
    tempDshearMag = BinDownMean( mean( vertcat(deform{x}(expt{x}.preRuns).DshearMag), 2, 'omitnan'), preCSD_2D_opts{x}.binSize ); %circshift(tempAccelCat, 5); %   
    tempVelocityCat = BinDownMean( vertcat(loco{x}(expt{x}.preRuns).Vdown), preCSD_2D_opts{x}.binSize );
    tempAccelCat = abs(BinDownMean( vertcat(loco{x}(expt{x}.preRuns).Adown), preCSD_2D_opts{x}.binSize )); 
    tempStateCat = BinDownMean( vertcat(loco{x}(expt{x}.preRuns).stateDown), preCSD_2D_opts{x}.binSize );
    preCSD_2D_pred{x} = struct('data',[], 'name',[], 'N',NaN, 'TB',[], 'lopo',[], 'fam',[]); 
    preCSD_2D_pred{x}.data = [tempTrans, tempScaleMag, tempShearMag, tempTransSpd, tempStretchMag, tempDshearMag, tempVelocityCat, tempAccelCat, tempStateCat]; % 
    preCSD_2D_pred{x}.name = {'|Translation|', '|Scale|', '|Shear|','TransSpeed', '|Stretch|', '|Shear| Rate', 'Velocity', '|Accel|', 'Locomotive State'}; % 
    preCSD_2D_pred{x}.N = size(preCSD_2D_pred{x}.data,2);
    for p = flip(1:preCSD_2D_pred{x}.N), preCSD_2D_pred{x}.lopo.name{p} = ['No ',preCSD_2D_pred{x}.name{p}]; end
    % Set up leave-one-family-out
    preCSD_2D_pred{x}.fam.col = {1:3, 4:6, 1:6, 7:9}; 
    preCSD_2D_pred{x}.fam.N = numel(preCSD_2D_pred{x}.fam.col); 
    preCSD_2D_pred{x}.fam.name = {'Reg', 'Deriv', 'Deformation', 'Kinematics'};  % 'Loco. State'

    % Define responses
    tempFluor = [fluor{x}(expt{x}.preRuns).z]; % [fluor{x}(expt{x}.preRuns).act]; % 
    tempFluorCat = BinDownMean( vertcat(tempFluor.ROI), preCSD_2D_opts{x}.binSize ); 
    preCSD_2D_resp{x}.data = tempFluorCat; % tempScaleMag; %
    preCSD_2D_resp{x}.N = size(preCSD_2D_resp{x}.data, 2); 
    preCSD_2D_resp{x}.name = sprintfc('Fluor %i', 1:preCSD_2D_resp{x}.N);

    % Remove scans with missing data 
    nanFrame = find(any(isnan([preCSD_2D_pred{x}.data, preCSD_2D_resp{x}.data]),2));
    fprintf('\nRemoving %i NaN-containing frames', numel(nanFrame));
    preCSD_2D_pred{x}.data(nanFrame,:) = [];  preCSD_2D_resp{x}.data(nanFrame,:) = [];
    %}
    % Run the GLM
    preCSD_2D_opts{x}.load = true; % false; %     
    preCSD_2D_opts{x}.saveRoot = expt{x}.dir; %sprintf('%s%s', expt{x}.dir  ); %''; % , preCSD_2D_opts{x}.name
    [ preCSD_2D_result{x}, preCSD_2D_summary{x}, preCSD_2D_opts{x}, preCSD_2D_pred{x}, preCSD_2D_resp{x} ] = GLMparallel( preCSD_2D_pred{x}, preCSD_2D_resp{x}, preCSD_2D_opts{x} ); 
    preCSD_2D_summary{x} = SummarizeGLM( preCSD_2D_result{x}, preCSD_2D_pred{x}, preCSD_2D_resp{x}, preCSD_2D_opts{x} );

    %preCSD_2D_opts{x}.rShow = 1:preCSD_2D_resp{x}.N;
    %ViewGLM(preCSD_2D_pred{x}, preCSD_2D_resp{x}, preCSD_2D_opts{x}, preCSD_2D_result{x}, preCSD_2D_summary{x});
    %pause;
end

%% View GLM results
for x = find(~cellfun(@isempty, preCSD_2D_result)) %  xPresent
    preCSD_2D_opts{x}.rShow = preCSD_2D_summary{x}.rGood; % 2; %1:7; %1:LocoDeform_resp{x}.N; %NaN; % 1:LocoDeform_resp{x}.N; %NaN;
    preCSD_2D_opts{x}.xVar = 'Time';
    ViewGLM(preCSD_2D_pred{x}, preCSD_2D_resp{x}, preCSD_2D_opts{x}, preCSD_2D_result{x}, preCSD_2D_summary{x}); %GLMresultFig = 
end

%% Compare GLM to data for each experiment
close all; clearvars sp SP;
PreGLMresults = figure('WindowState','maximized', 'color','w');
opt = {[0.02,0.07], [0.08,0.03], [0.04,0.02]};  % {[vert, horz], [bottom, top], [left, right] }\
rightOpt = {[0.1,0.07], [0.1,0.03], [0.04,0.02]};  % {[vert, horz], [bottom, top], [left, right] }\
jitterWidth = 0.45;
xAngle = 30;
Nrow = preCSD_2D_pred{xPresent(1)}.N+1;   Ncol = 3;
spGrid = reshape( 1:Nrow*Ncol, Ncol, Nrow )';
for x = find(~cellfun(@isempty, preCSD_2D_result)) % xPresent xPresent
    sp(Nrow) = subtightplot(Nrow, 3, 1:2, opt{:});
    [~,sortInd] = sort( preCSD_2D_summary{x}.dev, 'descend');
    imagesc( preCSD_2D_resp{x}.data(:,sortInd)' ); hold on;
    caxis([-5,5]);
    colormap(bluewhitered);
    ylabel('Fluor', 'Interpreter','none');
    title( sprintf('%s', expt{x}.name), 'Interpreter','none');
    set(gca,'TickDir','out', 'TickLength',[0.003,0], 'box','off', 'Xtick',[]); % , 'Ytick',onStruct(x).fluor.responder
    line( size(preCSD_2D_resp{x}.data,1)*[0,1], (preCSD_2D_summary{x}.Ngood+0.5)*[1,1], 'color','k');
    %text( repmat(size(preCSD_2D_resp{x}(r).data,1)+1, preCSD_2D_summary{x}.Ngood, 1), preCSD_2D_summary{x}.rGood+0.5, '*', 'VerticalAlignment','middle', 'FontSize',8);
    impixelinfo;
    
    for v = 1:preCSD_2D_pred{x}(r).N
        sp(v) = subtightplot(Nrow, Ncol, spGrid(v+1,1:2), opt{:});
        plot( preCSD_2D_pred{x}(r).data(:,v) ); hold on;
        ylabel(preCSD_2D_pred{x}(r).name{v}, 'Interpreter','none');
        xlim([-Inf,Inf]);
        if v < preCSD_2D_pred{x}(r).N
            set(gca,'TickDir','out', 'TickLength',[0.003,0], 'box','off', 'XtickLabel',[]);
        else
            set(gca,'TickDir','out', 'TickLength',[0.003,0], 'box','off');
        end
    end
    xlabel('Scan');
    linkaxes(sp,'x');
    xlim([-Inf,Inf]);
    
    subtightplot(2,3,3, rightOpt{:});
    bar([preCSD_2D_summary{x}.Ngood]/expt{x}.Nroi ); % numel(onStruct(x).fluor.responder),   , numel(rLocoPreFit{x})
    set(gca,'Xtick',1, 'XtickLabel',{'Fit'}, 'box','off'); % 'Loco','Fit','Both'  :3
    ylabel('Fraction of ROI');
    ylim([0,1]);
    
    subtightplot(2,3,6, rightOpt{:});
    JitterPlot( vertcat(preCSD_2D_summary{x}.lopo.devFrac(:,preCSD_2D_summary{x}.rGood), preCSD_2D_summary{x}.lofo.devFrac(:,preCSD_2D_summary{x}.rGood))', jitterWidth ); hold on;
    set(gca, 'Xtick',1:preCSD_2D_pred{x}(r).N+preCSD_2D_pred{x}(r).fam.N,  'XtickLabel', [preCSD_2D_summary{x}.lopo.name, preCSD_2D_summary{x}.lofo.name], 'TickDir','out', 'TickLength',[0.003,0], 'TickLabelInterpreter','none', 'box','off' ); 
    line([0,preCSD_2D_pred{x}(r).N+preCSD_2D_pred{x}(r).fam.N+1], [1,1], 'color','k', 'lineStyle','--');
    xlim([0,preCSD_2D_pred{x}(r).N+preCSD_2D_pred{x}(r).fam.N+1]); 
    ylim([0,Inf]); 
    ylabel('Fraction of total deviance'); title('Leave One Predictor Out (well-fit units only)');
    xtickangle(xAngle);
    
    %{
    figPath = sprintf('%s%s_Deviance.tif', figDir, GLMname);
    if exist(figPath,'file'), delete(figPath); end
    fprintf('\nSaving %s', figPath);
    export_fig( figPath, '-pdf', '-painters','-q101', '-append', LocoSensitivePrePost ); pause(1);
    print(PreGLMresults, figPath, '-dtiff' ); 
    pause%(1);   
    clf;
    %}
    pause; clf;
end

%% Show single examples of each subtype
for x = xPresent
    preCSD_2D_opts{x}.xVar = 'Time';
    %{
    if preCSD_2D_summary{x}.nMixed > 0
        preCSD_2D_opts{x}.rShow = preCSD_2D_summary{x}.rMixed;
        ViewGLM(preCSD_2D_pred{x}, preCSD_2D_resp{x}, preCSD_2D_opts{x}, preCSD_2D_result{x}, preCSD_2D_summary{x});
    end
    %}
    if preCSD_2D_summary{x}.nDeform > 0
        preCSD_2D_opts{x}.rShow = preCSD_2D_summary{x}.rDeform;
        ViewGLM(preCSD_2D_pred{x}, preCSD_2D_resp{x}, preCSD_2D_opts{x}, preCSD_2D_result{x}, preCSD_2D_summary{x});
    end
    %{
    if preCSD_2D_summary{x}.nLoco > 0
        preCSD_2D_opts{x}.rShow = preCSD_2D_summary{x}.rLoco;
        ViewGLM(preCSD_2D_pred{x}, preCSD_2D_resp{x}, preCSD_2D_opts{x}, preCSD_2D_result{x}, preCSD_2D_summary{x});
    end
    %}
end


%% Pool results across experiments
preCSDdevPool = []; goodDevPool = [];
preCSDdevFracPool = []; %lofoDevFracPool = [];
preCSD_2D_Nsubtype = []; k = 0;
for x = find(~cellfun(@isempty, preCSD_2D_result)) % xPresentxPresent%find(~cellfun(@isempty, rLocoPreFit))
    k = k+1;
    preCSD_2D_Nsubtype(k,:) = [preCSD_2D_summary{x}.nIns, preCSD_2D_summary{x}.nMixed, preCSD_2D_summary{x}.nDeform, preCSD_2D_summary{x}.nLoco]; % /expt{x}.Nroi
    if ~isempty( preCSD_2D_summary{x}.rGood )
        preCSDdevPool = [preCSDdevPool, preCSD_2D_summary{x}.dev]; 
        goodDevPool = [goodDevPool, preCSD_2D_summary{x}.dev( preCSD_2D_summary{x}.rGood )];
        preCSDdevFracPool = [preCSDdevFracPool, vertcat(preCSD_2D_summary{x}.lopo.devFrac(:,preCSD_2D_summary{x}.rGood), preCSD_2D_summary{x}.lofo.devFrac(:,preCSD_2D_summary{x}.rGood) )]; % rLocoPreFit{x}, :
    end
end
preCSD_2D_Nsubtype(k+1,:) = sum(preCSD_2D_Nsubtype, 1);
preCSD_2D_subtypeFrac = preCSD_2D_Nsubtype./repmat( sum(preCSD_2D_Nsubtype,2), 1, 4);

figure; bar(preCSD_2D_subtypeFrac,'stacked')


%% Summarize deviance explained
opt = {[0.02,0.07], [0.1,0.07], [0.09,0.09]};  % {[vert, horz], [bottom, top], [left, right] }
DevianceFig = figure('WindowState','maximized', 'color','w');
k = 1; clearvars h;
subtightplot(1,3,1,opt{:});
for x = find(~cellfun(@isempty, preCSD_2D_result)) % xPresent 
    [Ftemp, Xtemp] = ecdf( preCSD_2D_summary{x}.dev ); hold on;
    h(k) = plot(Xtemp, Ftemp, 'color',exptColor(k,:) );
    k = k + 1;
end
[Fdev, Xdev] = ecdf( preCSDdevPool );
h(k) = plot( Xdev, Fdev, 'color','k', 'LineWidth',2 ); 
axis square;
legend(h, {expt(xPresent).name, 'Pooled', 'Threshold'}, 'Location','southEast', 'Interpreter','none', 'AutoUpdate',false );
xlim([0, 0.4]);
line(preCSD_2D_opts{xPresent(1)}.minDev*[1,1], [0,1], 'Color','r', 'LineStyle','--'); % h(k+1) = 
xlabel('Deviance Explained'); ylabel('Fraction of Units');
title( sprintf('%s Fit Results', GLMname), 'Interpreter','none' );

subtightplot(1,3,2,opt{:});
JitterPlot( 1 - preCSDdevFracPool', 0.45, 'ErrorCap',10, 'monochrome',0.6); hold on;
line([0,size(preCSDdevFracPool,1)+1], [0,0], 'color','k');
axis square;
ylabel('Relative Explanatory Value'); %ylabel('Cumulative Distribution');
ylim([-1,1]);
set(gca,'Xtick', 1:size(preCSDdevFracPool,1), 'XtickLabel',[preCSD_2D_pred{xPresent(1)}.name, preCSD_2D_pred{xPresent(1)}.fam.name]);
xtickangle(30);
title('Well-Fit Units');

subtightplot(1,3,3,opt{:});
bar(preCSD_2D_subtypeFrac,'stacked');
ylabel('Fraction of Units');
barPos = get(gca, 'Position');
xlim([0, size(preCSD_2D_subtypeFrac,1)+1]);
ylim([0,1]);
axis square;
title( sprintf('Subtype Breakdown (min deviance explained = %2.2f)', preCSD_2D_opts{x}.minDev ));
legend('Insensitive','Mixed','Deform-only','Loco-only', 'Location','EastOutside');
set(gca, 'Xtick',1:size(preCSD_2D_subtypeFrac,1), 'XtickLabel',{expt(xPresent).name, 'Pooled'}, 'TickLabelInterpreter','none', 'FontSize',10, 'Position',barPos );
xtickangle(30);

% Save the figure
figPath = sprintf('%s%sDevianceResults.tif', figDir, GLMname);
if exist(figPath,'file'), delete(figPath); end
fprintf('\nSaving %s', figPath);
print(DevianceFig, figPath, '-dtiff', '-r300'); %pause(1); clf;

%% Show mean coeff by type

typeColor = getColorSet(preCSD_2D_pred{xPresent(1)}.N,10,true);  %distinguishable_colors(preCSD_2D_pred{xPresent(1)}.N);
close all;
SubtypeCoeffFig = figure('WindowState','maximized', 'color','w');
opt = {[0.02,0.07], [0.1,0.1], [0.1,0.06]};  % {[vert, horz], [bottom, top], [left, right] }
LW = 1.5;
colororder( typeColor ) 
for x = find(~cellfun(@isempty, preCSD_2D_result)) % xPresent
    tempInsCoeff = cat(3, preCSD_2D_result{x}( preCSD_2D_summary{x}.rIns ).coeff );
    meanInsCoeff = mean(tempInsCoeff, 3, 'omitnan' );
    
    tempMixedCoeff = cat(3, preCSD_2D_result{x}( preCSD_2D_summary{x}.rMixed ).coeff );
    meanMixedCoeff = mean(tempMixedCoeff, 3, 'omitnan' );
    
    tempDeformCoeff = cat(3, preCSD_2D_result{x}( preCSD_2D_summary{x}.rDeform ).coeff );
    meanDeformCoeff = mean(tempDeformCoeff, 3, 'omitnan' );
    
    tempLocoCoeff = cat(3, preCSD_2D_result{x}( preCSD_2D_summary{x}.rLoco ).coeff );
    meanLocoCoeff = mean(tempLocoCoeff, 3, 'omitnan' );
    meanLocoCoeff = meanLocoCoeff./max(meanLocoCoeff);
    meanLocoCoeff(isinf(meanLocoCoeff) | isnan(meanLocoCoeff)) = 0;
    
    sp(1) = subtightplot(1,4,1,opt{:});
    if preCSD_2D_summary{x}.nIns > 0
        plot(preCSD_2D_opts{x}.lags,  meanInsCoeff, 'LineWidth',LW );
    end
    axis square;
    tempPos = get(gca,'Position');
    xlabel('Lag (s)'); ylabel('Coefficient'); 
    title( sprintf('Insensitive (n = %i)', preCSD_2D_summary{x}.nIns) );
    legend(preCSD_2D_pred{x}.name, 'Location','NorthWest', 'AutoUpdate',false)
    set(gca,'Position',tempPos);
    
    sp(2) = subtightplot(1,4,2,opt{:});
    if preCSD_2D_summary{x}.nMixed > 0
        plot(preCSD_2D_opts{x}.lags,  meanMixedCoeff, 'LineWidth',LW ); 
    end
    axis square;
    title( sprintf('Mixed (n = %i)', preCSD_2D_summary{x}.nMixed) );
    xlabel('Lag (s)'); 
    
    sp(3) = subtightplot(1,4,3,opt{:});
    if preCSD_2D_summary{x}.nDeform > 0
        plot(preCSD_2D_opts{x}.lags,  meanDeformCoeff, 'LineWidth',LW ); 
    end
    axis square;
    title( sprintf('Deformation-dependent (n = %i)', preCSD_2D_summary{x}.nDeform) );
    xlabel('Lag (s)'); %title('Deformation-dependent'); % ylabel('Coefficient');
    
    sp(4) = subtightplot(1,4,4,opt{:});
    if preCSD_2D_summary{x}.nLoco > 0
        plot(preCSD_2D_opts{x}.lags,  meanLocoCoeff, 'LineWidth',LW );
    end
    axis square;
    xlabel('Lag (s)'); 
    title(sprintf('Locomotion-dependent (n = %i)', preCSD_2D_summary{x}.nLoco)); % ylabel('Coefficient'); 
    linkaxes(sp,'xy');
    
    pause;
    
    % Save the figure
    %{
    figPath = sprintf('%s%s_%s_SubtypeCoeff.tif', figDir, GLMname, expt{x}.name);
    if exist(figPath,'file'), delete(figPath); end
    print(SubtypeCoeffFig, figPath, '-dtiff', '-r300' );  fprintf('\nSaved %s\n', figPath);
    %}
    clf;
end

%%
devSummMat = cell(1,Nexpt); devSummPool = [];
opt = {[0.02,0.05], [0.02,0.02], [0.06,0.03]};  % {[vert, horz], [bottom, top], [left, right] }
GLMresultFig = figure('WindowState','maximized', 'color','w');

for x = xPresent
    cla;
    devSummMat{x} = [preCSD_2D_summary{x}.dev; preCSD_2D_summary{x}.lopo.dev; preCSD_2D_summary{x}.lofo.dev];
    %devSummPool = [devSummPool, devSummMat{x}];
    % {
    subtightplot(1,1,1,opt{:});
    imagesc( devSummMat{x} ); %imagesc( devSummPool ); %
    axis image;
    set(gca, 'Ytick', 1:size(devSummMat{x}, 1), 'YtickLabel', [{'All'}; preCSD_2D_summary{x}.lopo.name(:); preCSD_2D_summary{x}.lofo.name(:)], ...
        'Xtick',1:10:preCSD_2D_resp{x}.N, 'TickDir','out', 'TickLength',[0.003,0], 'FontSize',8); % 'XtickLabel',preCSD_2D_resp{x}.name
    title(sprintf('%s: GLM Fit Summary', preCSD_2D_opts{x}.name), 'Interpreter','none');
    %xtickangle(30);
    
    CB = colorbar; CB.Label.String = 'Deviance Explained';
    impixelinfo;
    pause;
    %}
end

devSummCat = cat(2, devSummMat{:});
devSummMed = median(devSummCat, 2, 'omitnan' );

GLMdevFig = figure('WindowState','maximized', 'color','w');
imagesc( devSummMed' );
set(gca,'XtickLabel',  [{'All'}; preCSD_2D_summary{x}.lopo.name(:); preCSD_2D_summary{x}.lofo.name(:)], ...
    'YtickLabel',preCSD_2D_resp{x}.name, 'TickDir','out', 'TickLength',[0.003,0]);
title('GLM Fit Summary (All 3D Data)', 'Interpreter','none');
xtickangle(30);
axis image;
CB = colorbar; CB.Label.String = 'Median Deviance Explained';
%figPath = sprintf('%s%s_DevSummary.tif', figDir, GLMname );
%fprintf('\nSaving %s\n', figPath);
%print( GLMdevFig, figPath, '-dtiff');
impixelinfo;

%% Plot coefficient values for each predictor, unit and experiment
close all; clearvars sp SP;
PreGLMcoeff = figure('WindowState','maximized', 'color','w');
opt = {[0.03,0.04], [0.09,0.05], [0.04,0.02]};  % {[vert, horz], [bottom, top], [left, right] }
for v = 1:preCSD_2D_pred{x}(r).N
    c = 0;
    for x = xPresent %find(~cellfun(@isempty, rLocoPreFit)) %
        c = c+1;
        if preCSD_2D_summary{x}.Ngood > 0
            subtightplot(1, preCSD_2D_pred{x}(r).N, v, opt{:});
            plot(c, preCSD_2D_summary{x}.peakCoeff(preCSD_2D_summary{x}.rGood,v), 'k.' ); hold on; %rLocoPreFit{x}
        end
    end
    line([0,c+1], [0,0], 'color','k');
    title( preCSD_2D_pred{x}(r).name{v}, 'Interpreter','none' );
    if v == 1, ylabel('Peak Coefficient'); end
    xlim([0,c+1]);
    set(gca, 'Xtick', 1:c, 'XtickLabel', {expt(xPresent).name}, 'TickLabelInterpreter','none' );
    xtickangle(45);
    axis square;
end

%% Plot peak coefficient vs latency values for each predictor, good unit, and experiment
close all; clearvars sp SP;
PreGLMcoeff = figure('WindowState','maximized', 'color','w');
opt = {[0.03,0.04], [0.09,0.05], [0.04,0.02]};  % {[vert, horz], [bottom, top], [left, right] }

for v = 1:preCSD_2D_pred{x}(r).N
    c = 0;
    for x = xPresent %find(~cellfun(@isempty, rLocoPreFit)) %
        c = c+1;
        if preCSD_2D_summary{x}.Ngood > 0
            subtightplot(1, preCSD_2D_pred{x}(r).N, v, opt{:});
            plot(preCSD_2D_summary{x}.peakLag(preCSD_2D_summary{x}.rGood,v), preCSD_2D_summary{x}.peakCoeff(preCSD_2D_summary{x}.rGood,v), '.', 'color',exptColor(c,:) ); hold on; %rLocoPreFit{x}
        end
    end
    title( preCSD_2D_pred{xPresent(1)}.name{v}, 'Interpreter','none' );
    if v == 1, ylabel('Coefficient'); end
    xlabel('Lag (s)');
    xlim([-6,6]);
    axis square;
end

%% Check mechanosensitive units for sigmoidal stimulus response curve to various forms of deformation PRE-CSD

savePath = strcat(figDir,'sigmoidsAnalysis.mat');
if exist(savePath, 'file')
    load(savePath); fprintf('\nLoading %s', savePath); 
else
    responderMat = cell(1,Nexpt);
    sigmResp_pre = repmat( struct('speed',[], 'Nspeed',NaN, 'speedFrac',NaN, 'trans',[], 'Ntrans',NaN, 'transFrac',NaN, 'scale',[], 'Nscale',NaN, 'scaleFrac',NaN, 'stretch',[], 'Nstretch',NaN, 'stretchFrac',NaN, ...
        'shear',[], 'Nshear',NaN, 'shearFrac',NaN, 'shearRate',[], 'NshearRate',NaN, 'shearRateFrac',NaN, 'poly',[], 'Npoly',NaN, 'polyFrac',NaN), 1, Nexpt);
    tic
    for x = x2D
        fluorPre = [fluor{x}(expt{x}.preRuns).dFF];
        fluorPre = vertcat(fluorPre.ROI);
        speedPre = loco{x}(expt{x}.preRuns);
        speedPre = vertcat(speedPre.speedDown);
        deformPre = deform{x}(expt{x}.preRuns);
        transMagPre = vertcat(deformPre.transMag);
        scaleMagPre = vertcat(deformPre.scaleMag);
        stretchMagPre = vertcat(deformPre.stretchMag);
        shearMagPre = vertcat(deformPre.shearMag);
        shearRateMagPre = vertcat(deformPre.DshearMag);

        rCheck = [preCSD_2D_summary{x}.rDeform, preCSD_2D_summary{x}.rMixed];
        Ncheck = numel(rCheck);
        responderMat{x} = nan(5, expt{x}.Nroi);
        responderMat{x}(:,rCheck) = 0;

        % Scaling
        [scaleStim{x}, scaleResp{x}, scaleResult{x}, scaleSumm{x}, scaleFit{x}] = StimResponse(scaleMagPre, fluorPre, 0, 10, 0, 'fit',true, 'show',false); % , 'show',false  dffResp{x}
        sigmResp_pre(x).scale = intersect(rCheck, scaleSumm{x}.rGood);
        sigmResp_pre(x).Nscale = numel(sigmResp_pre(x).scale);
        sigmResp_pre(x).scaleFrac = sigmResp_pre(x).Nscale/Ncheck;
        responderMat{x}(2,sigmResp_pre(x).scale) = 1;
        sigmResp_pre(x).scaleThresh = [scaleResult{x}(sigmResp_pre(x).scale).thresh];
        sigmResp_pre(x).scaleSlope = [scaleResult{x}(sigmResp_pre(x).scale).rate];
        %{
        % Speed
        [speedStim{x}, speedResp{x}, speedResult{x}, speedSumm{x}, speedFit{x}] = StimResponse(speedPre, fluorPre, 0, 10, 0, 'fit',true); %   dffResp{x} , 'show',true
        sigmResp_pre(x).speed = intersect([preCSD_2D_summary{x}.rLoco, preCSD_2D_summary{x}.rMixed], speedSumm{x}.rGood);
        sigmResp_pre(x).Nspeed = numel(sigmResp_pre(x).speed);
        sigmResp_pre(x).speedFrac = sigmResp_pre(x).Nspeed/numel([preCSD_2D_summary{x}.rLoco, preCSD_2D_summary{x}.rMixed]);
        sigmResp_pre(x).speedThresh = [speedResult{x}(sigmResp_pre(x).speed).thresh];
        sigmResp_pre(x).speedSlope = [speedResult{x}(sigmResp_pre(x).speed).rate];
        % Translation
        [transStim{x}, transResp{x}, transResult{x}, transSumm{x}, transFit{x}] = StimResponse(transMagPre, fluorPre, 0, 10, 0, 'fit',true); % , 'show',false  dffResp{x}
        sigmResp_pre(x).trans = intersect(rCheck, transSumm{x}.rGood);
        sigmResp_pre(x).Ntrans = numel(sigmResp_pre(x).trans);
        sigmResp_pre(x).transFrac = sigmResp_pre(x).Ntrans/Ncheck;
        responderMat{x}(1,sigmResp_pre(x).trans) = 1;
        sigmResp_pre(x).transThresh = [transResult{x}(sigmResp_pre(x).trans).thresh];
        sigmResp_pre(x).transSlope = [transResult{x}(sigmResp_pre(x).trans).rate];

        % Stretch
        [stretchStim{x}, stretchResp{x}, stretchResult{x}, stretchSumm{x}, stretchFit{x}] = StimResponse(stretchMagPre, fluorPre, 0, 10, 0, 'fit',true); % , 'show',false
        sigmResp_pre(x).stretch = intersect(rCheck, stretchSumm{x}.rGood);
        sigmResp_pre(x).Nstretch = numel(sigmResp_pre(x).stretch);
        sigmResp_pre(x).stretchFrac = sigmResp_pre(x).Nstretch/Ncheck;
        responderMat{x}(3,sigmResp_pre(x).stretch) = 1;
        sigmResp_pre(x).stretchThresh = [stretchResult{x}(sigmResp_pre(x).stretch).thresh];
        sigmResp_pre(x).stretchSlope = [stretchResult{x}(sigmResp_pre(x).stretch).rate];
        % Shear
        [shearStim{x}, shearResp{x}, shearResult{x}, shearSumm{x}, shearFit{x}] = StimResponse(shearMagPre, fluorPre, 0, 10, 0, 'fit',true); % , 'show',false
        sigmResp_pre(x).shear = intersect(rCheck, shearSumm{x}.rGood);
        sigmResp_pre(x).Nshear = numel(sigmResp_pre(x).shear);
        sigmResp_pre(x).shearFrac = sigmResp_pre(x).Nshear/Ncheck;
        responderMat{x}(4,sigmResp_pre(x).shear) = 1;
        sigmResp_pre(x).shearThresh = [shearResult{x}(sigmResp_pre(x).shear).thresh];
        sigmResp_pre(x).shearSlope = [shearResult{x}(sigmResp_pre(x).shear).rate];
        % Shear rate
        [shearRateStim{x}, shearRateResp{x}, shearRateResult{x}, shearRateSumm{x}, shearRateFit{x}] = StimResponse(shearRateMagPre, fluorPre, 0, 10, 0, 'fit',true); % , 'show',false
        sigmResp_pre(x).shearRate = intersect(rCheck, shearRateSumm{x}.rGood);
        sigmResp_pre(x).NshearRate = numel(sigmResp_pre(x).shearRate);
        sigmResp_pre(x).shearRateFrac = sigmResp_pre(x).NshearRate/Ncheck;
        responderMat{x}(5,sigmResp_pre(x).shearRate) = 1;
        sigmResp_pre(x).shearRateThresh = [shearRateResult{x}(sigmResp_pre(x).shearRate).thresh];
        sigmResp_pre(x).shearRateSlope = [shearRateResult{x}(sigmResp_pre(x).shearRate).rate];  

        % Poly-responders
        sigmResp_pre(x).poly = find( sum(responderMat{x}, 1) > 1);
        sigmResp_pre(x).Npoly = numel(sigmResp_pre(x).poly);
        sigmResp_pre(x).polyFrac = sigmResp_pre(x).Npoly/Ncheck;
        %}
        toc
        %bar( [deformResponders(x).transFrac, deformResponders(x).scaleFrac, deformResponders(x).stretchFrac, deformResponders(x).shearFrac] )
    end
    
    save(savePath, 'sigmResp_pre','responderMat','scaleStim','scaleResp','scaleResult','scaleSumm','scaleFit', '-v7.3'); 
    fprintf('\nSaving %s', savePath);
end

%% Check mechanosensitive units for sigmoidal stimulus response curve to various forms of deformation POST-CSD
postCSD_cutoff = 15; % how many minutes after CSD to start grabbing data?
responderMat_post = cell(1,Nexpt);
sigmResp_post = repmat( struct('speed',[], 'Nspeed',NaN, 'speedFrac',NaN, 'trans',[], 'Ntrans',NaN, 'transFrac',NaN, 'scale',[], 'Nscale',NaN, 'scaleFrac',NaN, 'stretch',[], 'Nstretch',NaN, 'stretchFrac',NaN, ...
    'shear',[], 'Nshear',NaN, 'shearFrac',NaN, 'shearRate',[], 'NshearRate',NaN, 'shearRateFrac',NaN, 'poly',[], 'Npoly',NaN, 'polyFrac',NaN), 1, Nexpt);
tic
for x = x2D
    postRuns = 3:4;
    Tpost = vertcat(Tscan{x}{postRuns});
    Tpost = Tpost - Tpost(1);
    postScans = find( Tpost > 60*postCSD_cutoff );
    fluorPost = [fluor{x}(postRuns).dFF];
    fluorPost = vertcat(fluorPost.ROI);
    fluorPost = fluorPost(postScans,:);
    speedPost = loco{x}(postRuns);
    speedPost = vertcat(speedPost.speedDown);
    speedPost = speedPost(postScans,:);
    deformPost = deform{x}(postRuns);
    transMagPost = vertcat(deformPost.transMag);
    transMagPost = transMagPost(postScans,:);
    scaleMagPost = vertcat(deformPost.scaleMag);
    scaleMagPost = scaleMagPost(postScans,:);
    stretchMagPost = vertcat(deformPost.stretchMag);
    stretchMagPost = stretchMagPost(postScans,:);
    shearMagPost = vertcat(deformPost.shearMag);
    shearMagPost = shearMagPost(postScans,:);
    shearRateMagPost = vertcat(deformPost.DshearMag);
    shearRateMagPost = shearRateMagPost(postScans,:);

    rCheck = [preCSD_2D_summary{x}.rDeform, preCSD_2D_summary{x}.rMixed];
    Ncheck = numel(rCheck);
    responderMat_post{x} = nan(5, expt{x}.Nroi);
    responderMat_post{x}(:,rCheck) = 0;
    
    % Speed
    [speedStim{x}, speedResp{x}, speedResult{x}, speedSumm{x}, speedFit{x}] = StimResponse(speedPost, fluorPost, 0, 10, 0, 'fit',true); %   dffResp{x} , 'show',true
    sigmResp_post(x).speed = intersect([preCSD_2D_summary{x}.rLoco, preCSD_2D_summary{x}.rMixed], speedSumm{x}.rGood);
    sigmResp_post(x).Nspeed = numel(sigmResp_post(x).speed);
    sigmResp_post(x).speedFrac = sigmResp_post(x).Nspeed/numel([preCSD_2D_summary{x}.rLoco, preCSD_2D_summary{x}.rMixed]);
    sigmResp_post(x).speedThresh = [speedResult{x}(sigmResp_post(x).speed).thresh];
    sigmResp_post(x).speedSlope = [speedResult{x}(sigmResp_post(x).speed).rate];
    % Translation
    [transStim{x}, transResp{x}, transResult{x}, transSumm{x}, transFit{x}] = StimResponse(transMagPost, fluorPost, 0, 10, 0, 'fit',true); % , 'show',false  dffResp{x}
    sigmResp_post(x).trans = intersect(rCheck, transSumm{x}.rGood);
    sigmResp_post(x).Ntrans = numel(sigmResp_post(x).trans);
    sigmResp_post(x).transFrac = sigmResp_post(x).Ntrans/Ncheck;
    responderMat_post{x}(1,sigmResp_post(x).trans) = 1;
    sigmResp_post(x).transThresh = [transResult{x}(sigmResp_post(x).trans).thresh];
    sigmResp_post(x).transSlope = [transResult{x}(sigmResp_post(x).trans).rate];
    % Scaling
    [scaleStim{x}, scaleResp{x}, scaleResult{x}, scaleSumm{x}, scaleFit{x}] = StimResponse(scaleMagPost, fluorPost, 0, 10, 0, 'fit',true, 'show',false); % , 'show',false  dffResp{x}
    sigmResp_post(x).scale = intersect(rCheck, scaleSumm{x}.rGood);
    sigmResp_post(x).Nscale = numel(sigmResp_post(x).scale);
    sigmResp_post(x).scaleFrac = sigmResp_post(x).Nscale/Ncheck;
    responderMat_post{x}(2,sigmResp_post(x).scale) = 1;
    sigmResp_post(x).scaleThresh = [scaleResult{x}(sigmResp_post(x).scale).thresh];
    sigmResp_post(x).scaleSlope = [scaleResult{x}(sigmResp_post(x).scale).rate];
    % Stretch
    [stretchStim{x}, stretchResp{x}, stretchResult{x}, stretchSumm{x}, stretchFit{x}] = StimResponse(stretchMagPost, fluorPost, 0, 10, 0, 'fit',true); % , 'show',false
    sigmResp_post(x).stretch = intersect(rCheck, stretchSumm{x}.rGood);
    sigmResp_post(x).Nstretch = numel(sigmResp_post(x).stretch);
    sigmResp_post(x).stretchFrac = sigmResp_post(x).Nstretch/Ncheck;
    responderMat_post{x}(3,sigmResp_post(x).stretch) = 1;
    sigmResp_post(x).stretchThresh = [stretchResult{x}(sigmResp_post(x).stretch).thresh];
    sigmResp_post(x).stretchSlope = [stretchResult{x}(sigmResp_post(x).stretch).rate];
    % Shear
    [shearStim{x}, shearResp{x}, shearResult{x}, shearSumm{x}, shearFit{x}] = StimResponse(shearMagPost, fluorPost, 0, 10, 0, 'fit',true); % , 'show',false
    sigmResp_post(x).shear = intersect(rCheck, shearSumm{x}.rGood);
    sigmResp_post(x).Nshear = numel(sigmResp_post(x).shear);
    sigmResp_post(x).shearFrac = sigmResp_post(x).Nshear/Ncheck;
    responderMat_post{x}(4,sigmResp_post(x).shear) = 1;
    sigmResp_post(x).shearThresh = [shearResult{x}(sigmResp_post(x).shear).thresh];
    sigmResp_post(x).shearSlope = [shearResult{x}(sigmResp_post(x).shear).rate];
    % Shear rate
    [shearRateStim{x}, shearRateResp{x}, shearRateResult{x}, shearRateSumm{x}, shearRateFit{x}] = StimResponse(shearRateMagPost, fluorPost, 0, 10, 0, 'fit',true); % , 'show',false
    sigmResp_post(x).shearRate = intersect(rCheck, shearRateSumm{x}.rGood);
    sigmResp_post(x).NshearRate = numel(sigmResp_post(x).shearRate);
    sigmResp_post(x).shearRateFrac = sigmResp_post(x).NshearRate/Ncheck;
    responderMat_post{x}(5,sigmResp_post(x).shearRate) = 1;
    sigmResp_post(x).shearRateThresh = [shearRateResult{x}(sigmResp_post(x).shearRate).thresh];
    sigmResp_post(x).shearRateSlope = [shearRateResult{x}(sigmResp_post(x).shearRate).rate];  
    
    % Poly-responders
    sigmResp_post(x).poly = find( sum(responderMat_post{x}, 1) > 1);
    sigmResp_post(x).Npoly = numel(sigmResp_post(x).poly);
    sigmResp_post(x).polyFrac = sigmResp_post(x).Npoly/Ncheck;
    
    toc
    %bar( [deformResponders(x).transFrac, deformResponders(x).scaleFrac, deformResponders(x).stretchFrac, deformResponders(x).shearFrac] )
end

%% Resummarize the data in terms of simple mechanosensivity, and deviance explained that is attributed to of deformation
close all;
figure('WindowState','maximized');  %('Units','inches', 'Position', [6, 5, 3.5, 2.5], 'Color','w');   %
for x = xPresent
    lofoDeformCol = find(strcmpi(preCSD_2D_pred{x}.fam.name, 'Deformation'));
    lofoOtherCol = setdiff(1:preCSD_2D_pred{x}.fam.N, lofoDeformCol);
    tempResult = [preCSD_2D_result{x}];
    tempDev = [tempResult.dev]';
    tempDev(tempDev < preCSD_2D_opts{x}.minDev) = NaN; % suppress units that were poorly fit overall
    tempLofo = [tempResult.lofo];
    tempLofoDev = vertcat(tempLofo.dev);
    tempDev( tempLofoDev(:,lofoDeformCol) < preCSD_2D_opts{x}.minDev ) = NaN; % suppress units where dropping deformation ruins the fit
    
    famDevMargin = repmat(tempDev,1,preCSD_2D_pred{x}.fam.N) - tempLofoDev; % how much of the total deviance explained is lost after dropping each family? (larger number = that family is more important)
    famDevMargin( famDevMargin(:,lofoDeformCol) )

    [sortDev, sortDevInd] = sort(tempDev, 'descend', 'MissingPlacement','last');
    sortDevInd(find(isnan(sortDev), 1, 'first'):end) = [];
    tempMargin = famDevMargin(sortDevInd,lofoOtherCol);
    
    mechanoMat = preCSD_2D_resp{X}.data(:,sortDevInd); 
    subplot(1,2,1);
    imagesc(mechanoMat');
    xlim(scanLims);
    caxis([-3,3]);
    colormap(bluewhitered)
    
    subplot(1,2,2);
    bar(tempMargin, 'stacked'); hold on;
    legend(preCSD_2D_pred{x}.fam.name(lofoOtherCol))
    plot( tempDev(sortDevInd), 'x' )

    pause;
end





