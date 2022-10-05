%% Use GLM to assess contribution of different variables
postCSD_pred = cell(1,Nexpt); postCSD_resp = cell(1,Nexpt); postCSD_opts = cell(1,Nexpt); postCSD_result = cell(1,Nexpt); postCSD_summary = cell(1,Nexpt);
figDir = 'D:\MATLAB\LevyLab\Figures\3D\GLM\'; mkdir( figDir )
GLMname = 'postCSD_3D_GLM_TScShDtStDshZL';
GLMrate = 15.49/30;
for x = x3Dcsd %xPresent
    % GLMparallel options
    postCSD_opts{x}.name = sprintf('%s_%s', expt(x).name, GLMname); %strcat(expt(x).name, , '_postCSDglm');
    postCSD_opts{x}.rShow = NaN;
    postCSD_opts{x}.figDir = ''; % figDir;
    postCSD_opts{x}.alpha = 0.01;  % The regularization parameter, default is 0.01
    postCSD_opts{x}.standardize = true; 
    postCSD_opts{x}.trainFrac = 0.75; % 1; %
    postCSD_opts{x}.Ncycle = 20;
    postCSD_opts{x}.distribution = 'gaussian'; % 'poisson'; %  
    postCSD_opts{x}.CVfold = 10;
    postCSD_opts{x}.nlamda = 1000;
    postCSD_opts{x}.maxit = 5*10^5;
    postCSD_opts{x}.minDev = 0.05; 
    postCSD_opts{x}.minDevFrac = 0.1;
    postCSD_opts{x}.maxP = 0.05;
    postCSD_opts{x}.Nshuff = 0;  
    postCSD_opts{x}.minShuff = 15; 
    postCSD_opts{x}.window = [-4,4]; % [0,0]; % [-0.5, 0.5]; % 
    postCSD_opts{x}.lopo = true; %false; %
    postCSD_opts{x}.frameRate = GLMrate; %expt(x).scanRate; 
    postCSD_opts{x}.binSize = expt(x).scanRate/GLMrate;
    postCSD_opts{x}.minShuffFrame = round( postCSD_opts{x}.frameRate*postCSD_opts{x}.minShuff );
    windowFrame = round(postCSD_opts{x}.window*postCSD_opts{x}.frameRate);
    postCSD_opts{x}.shiftFrame = windowFrame(1):windowFrame(2);
    postCSD_opts{x}.maxShift = max( abs(windowFrame) );
    postCSD_opts{x}.Nshift = numel( postCSD_opts{x}.shiftFrame );  %Nshift = postCSDOpts(x).Nshift;
    postCSD_opts{x}.lags = postCSD_opts{x}.shiftFrame/postCSD_opts{x}.frameRate;
    postCSD_opts{x}.xVar = 'Time';

    % Concatenate input variables post-CSD, excluding acute period
    % Define predictors
    postRuns = expt(x).csd:expt(x).Nruns;
    tempT = BinDownMean( vertcat(Tscan{x}{postRuns}), postCSD_opts{x}.binSize );
    postScans = find(tempT > csdBout{x}(expt(x).csd).Tstop + 180 & tempT < csdBout{x}(expt(x).csd).Tstop + 3600 ); % exclude the acute CSD wave from this analysis
    tempTrans = BinDownMean( mean( vertcat(deform{x}(postRuns).transMag), 2, 'omitnan'), postCSD_opts{x}.binSize );
    tempTransSpd = BinDownMean( mean( vertcat(deform{x}(postRuns).DtransMag), 2, 'omitnan'), postCSD_opts{x}.binSize );
    tempScaleMag = BinDownMean( mean( vertcat(deform{x}(postRuns).scaleMag), 2, 'omitnan'), postCSD_opts{x}.binSize );
    tempStretchMag = BinDownMean( mean( vertcat(deform{x}(postRuns).stretchMag), 2, 'omitnan'), postCSD_opts{x}.binSize );
    tempShearMag = BinDownMean( mean( vertcat(deform{x}(postRuns).shearMag), 2, 'omitnan'), postCSD_opts{x}.binSize );
    tempDshearMag = BinDownMean( mean( vertcat(deform{x}(postRuns).DshearMag), 2, 'omitnan'), postCSD_opts{x}.binSize ); 
    tempShift = vertcat(deform{x}(postRuns).shiftZ);
    tempShiftMean =  BinDownMean( mean( tempShift(:,4:end-3), 2, 'omitnan' ), postCSD_opts{x}.binSize );
    tempStateCat = BinDownMean( vertcat(loco{x}(postRuns).stateDown), postCSD_opts{x}.binSize );

    postCSD_pred{x} = struct('data',[], 'name',[], 'N',NaN, 'TB',[], 'lopo',[], 'fam',[]); 
    postCSD_pred{x}.data = [tempTrans, tempScaleMag, tempShearMag, tempTransSpd, tempStretchMag, tempDshearMag, tempShiftMean, tempStateCat]; % , tempShiftMean tempStretchMag, tempExp, tempComp, 
    postCSD_pred{x}.data = postCSD_pred{x}.data(postScans,:);
    postCSD_pred{x}.name = {'|Translation|', '|Scale|', '|Shear|','TransSpeed', '|Stretch|' ,'|Shear| Rate', 'Z Shift', 'Locomotive State'}; % ,'Speed',  'Str-Exp', 'Str-Comp',
    postCSD_pred{x}.N = size(postCSD_pred{x}.data,2);
    for p = flip(1:postCSD_pred{x}.N), postCSD_pred{x}.lopo.name{p} = ['No ',postCSD_pred{x}.name{p}]; end
    % Set up leave-one-family-out
    postCSD_pred{x}.fam.col = {1:3, 4:6, 1:7, 8}; %{1:2, 3:4, 5:6, 7:8, 9:10, 11:12};  % {1:12};%{1, 2:3, 4:5, 6:7, 8, 9}; 
    postCSD_pred{x}.fam.N = numel(postCSD_pred{x}.fam.col); 
    postCSD_pred{x}.fam.name = {'Reg', 'Deriv', 'Deformation', 'Loco. State'}; %{'All'};%  'Onset Time',

    % Define response
    tempFluor = [fluor{x}(postRuns).z]; %[fluor{x}.z]; %[fluor{x}.act]; % 
    postCSD_resp{x}.data = BinDownMean( vertcat(tempFluor.ROI), postCSD_opts{x}.binSize ); 
    postCSD_resp{x}.data = postCSD_resp{x}.data(postScans,:);
    postCSD_resp{x}.N = size(postCSD_resp{x}.data,2); 
    postCSD_resp{x}.name = sprintfc('Fluor %i', 1:postCSD_resp{x}.N);

    % Remove scans with missing data 
    nanFrame = find(any(isnan([postCSD_pred{x}.data, postCSD_resp{x}.data]),2)); % find( isnan(sum(pred(x).data,2)) ); 
    fprintf('\nRemoving %i NaN-containing frames', numel(nanFrame));
    postCSD_pred{x}.data(nanFrame,:) = []; postCSD_resp{x}.data(nanFrame,:) = [];

    % Run the GLM
    postCSD_opts{x}.load = true; % false; %  
    postCSD_opts{x}.saveRoot = expt(x).dir; %
    [postCSD_result{x}, postCSD_summary{x}, ~, postCSD_pred{x}, postCSD_resp{x}] = GLMparallel(postCSD_pred{x}, postCSD_resp{x}, postCSD_opts{x}); % postCSD_result{x}, postCSD_summary(x), ~, postCSDPred(x), postCSDResp(x)
    postCSD_summary{x} = SummarizeGLM(postCSD_result{x}, postCSD_pred{x}, postCSD_resp{x}, postCSD_opts{x});
end

%% Compare GLM to data for each experiment
close all; clearvars sp SP;
PreGLMresults = figure('WindowState','maximized', 'color','w');
opt = {[0.02,0.07], [0.06,0.03], [0.04,0.02]};  % {[vert, horz], [bottom, top], [left, right] }\
rightOpt = {[0.1,0.07], [0.1,0.03], [0.04,0.02]};  % {[vert, horz], [bottom, top], [left, right] }\
jitterWidth = 0.45;
xAngle = 30;
Nrow = postCSD_pred{x3Dcsd(1)}(1).N+1; Ncol = 3;
spGrid = reshape( 1:Nrow*Ncol, Ncol, Nrow )';
for x = 37 %x3Dcsd
    sp(postCSD_pred{x}.N+1) = subtightplot(postCSD_pred{x}.N+1, 3, 1:2, opt{:});
    imagesc( postCSD_resp{x}.data' );
    ylabel('Fluor', 'Interpreter','none');
    title( sprintf('%s', expt(x).name), 'Interpreter','none');
    set(gca,'TickDir','out', 'TickLength',[0.003,0], 'box','off', 'Xtick',[]); % , 'Ytick',onStruct(x).fluor.responder
    text( repmat(size(postCSD_resp{x}.data,1)+1, postCSD_summary{x}.Ngood, 1), postCSD_summary{x}.rGood+0.5, '*', 'VerticalAlignment','middle', 'FontSize',8);
    impixelinfo;
    
    for v = 1:postCSD_pred{x}.N
        sp(v) = subtightplot(Nrow, Ncol, spGrid(v+1,1:2), opt{:});
        plot( postCSD_pred{x}.data(:,v) ); hold on;
        ylabel(postCSD_pred{x}.name{v}, 'Interpreter','none');
        xlim([-Inf,Inf]);
        if v < postCSD_pred{x}.N
            set(gca,'TickDir','out', 'TickLength',[0.003,0], 'box','off', 'XtickLabel',[]);
        else
            set(gca,'TickDir','out', 'TickLength',[0.003,0], 'box','off');
        end
    end
    xlabel('Scan');
    
    subtightplot(3,3,3, rightOpt{:});
    bar([postCSD_summary{x}.Ngood]/expt(x).Nroi ); % numel(onStruct(x).fluor.responder),   , numel(rLocoPreFit{x})
    set(gca,'Xtick',1, 'XtickLabel',{'Fit'}, 'box','off'); % 'Loco','Fit','Both'  :3
    ylabel('Fraction of ROI');
    ylim([0,1]);
    
    subtightplot(3,3,6, rightOpt{:});
    JitterPlot( postCSD_summary{x}.lopo.devFrac(:,postCSD_summary{x}.rGood)', jitterWidth ); hold on;
    line([0,postCSD_pred{x}.N+1], [1,1], 'color','k', 'lineStyle','--');
    xlim([0,postCSD_pred{x}.N+1]); ylim([0,Inf]); 
    ylabel('Fraction of total deviance'); title('Leave One Predictor Out (well-fit units only)');
    set(gca, 'Xtick',1:postCSD_pred{x}.N,  'XtickLabel', postCSD_summary{x}.lopo.name, 'TickDir','out', 'TickLength',[0.003,0], 'TickLabelInterpreter','none', 'box','off' ); 
    xtickangle(xAngle);
    
    
    subtightplot(3,3,9, rightOpt{:});
    JitterPlot( postCSD_summary{x}.lofo.devFrac(:,postCSD_summary{x}.rGood)', jitterWidth ); hold on;
    line([0,postCSD_pred{x}.fam.N]+0.5, [1,1], 'color','k', 'lineStyle','--');
    xlim([0,postCSD_pred{x}.fam.N]+0.5);
    ylabel('Fraction of total deviance'); title('Leave One Family Out (well-fit units only)');
    set(gca, 'Xtick',1:postCSD_pred{x}.fam.N,  'XtickLabel', postCSD_summary{x}.lofo.name, 'TickDir','out', 'TickLength',[0.003,0], 'TickLabelInterpreter','none', 'box','off' ); 
    xtickangle(xAngle);
    ylim([0,Inf]);
    
    linkaxes(sp,'x');
    % {
    figPath = sprintf('%s%s_Deviance.tif', figDir, GLMname);
    if exist(figPath,'file'), delete(figPath); end
    fprintf('\nSaving %s', figPath);
    %export_fig( figPath, '-pdf', '-painters','-q101', '-append', LocoSensitivePrePost ); pause(1);
    %print(PreGLMresults, figPath, '-dtiff' ); 
    pause%(1);   
    clf;
    %}
    %pause; clf;
end

%% Divide units into Insensitiveensitive, Mixed, loco-only or deformation-only units
postCSD_Nsubtype = []; k = 0;
for x = intersect( find(~cellfun(@isempty, postCSD_result)), x3Dcsd ) %x3Dcsd
   
    k = k+1;
    postCSD_Nsubtype(k,:) = [postCSD_summary{x}.nInsensitive, postCSD_summary{x}.nMixed, postCSD_summary{x}.nDeform, postCSD_summary{x}.nLoco]; % /expt(x).Nroi
end
postCSD_Nsubtype(k+1,:) = sum(postCSD_Nsubtype, 1);
postCSD_subtypeFrac = postCSD_Nsubtype./repmat( sum(postCSD_Nsubtype,2), 1, 4);

bar(postCSD_subtypeFrac,'stacked')

%% Show single examples of each subtype
for x = x3Dcsd
    postCSD_opts{x}.xVar = 'Time';
    %{
    if postCSD_summary{x}.nMixed > 0
        postCSD_opts{x}.rShow = postCSD_summary{x}.rMixed;
        ViewGLM(postCSD_pred{x}, postCSD_resp{x}, postCSD_opts{x}, postCSD_result{x}, postCSD_summary{x});
    end
    %}
    if postCSD_summary{x}.nDeform > 0
        postCSD_opts{x}.rShow = postCSD_summary{x}.rDeform;
        ViewGLM(postCSD_pred{x}, postCSD_resp{x}, postCSD_opts{x}, postCSD_result{x}, postCSD_summary{x});
    end
    %{
    if postCSD_summary{x}.nLoco > 0
        postCSD_opts{x}.rShow = postCSD_summary{x}.rLoco;
        ViewGLM(postCSD_pred{x}, postCSD_resp{x}, postCSD_opts{x}, postCSD_result{x}, postCSD_summary{x});
    end
    %}
end

%% Pool results across experiments
postCSDdevPool = []; goodDevPool = [];
postCSDdevFracPool = []; %lofoDevFracPool = [];
for x = x3Dcsd%find(~cellfun(@isempty, rLocoPreFit))
    if ~isempty( postCSD_summary{x}.rGood )
        postCSDdevPool = [postCSDdevPool, postCSD_summary{x}.dev]; 
        goodDevPool = [goodDevPool, postCSD_summary{x}.dev( postCSD_summary{x}.rGood )];
        postCSDdevFracPool = [postCSDdevFracPool, vertcat(postCSD_summary{x}.lopo.devFrac(:,postCSD_summary{x}.rGood), postCSD_summary{x}.lofo.devFrac(:,postCSD_summary{x}.rGood) )]; % rLocoPreFit{x}, :
    end
end
%% Summarize deviance explained
opt = {[0.02,0.07], [0.1,0.07], [0.09,0.09]};  % {[vert, horz], [bottom, top], [left, right] }
DevianceFig = figure('WindowState','maximized', 'color','w');
k = 1; clearvars h;
subtightplot(1,3,1,opt{:});
for x = x3Dcsd
    [Ftemp, Xtemp] = ecdf( postCSD_summary{x}.dev ); hold on;
    h(k) = plot(Xtemp, Ftemp, 'color',0.7*[1,1,1] );
    k = k + 1;
end
[Fdev, Xdev] = ecdf( postCSDdevPool );
h(k) = plot( Xdev, Fdev, 'color','k', 'LineWidth',2 ); 
axis square;
legend(h, {expt(x3Dcsd).name, 'Pooled', 'Threshold'}, 'Location','southEast', 'Interpreter','none', 'AutoUpdate',false );
xlim([0, 0.6]);
line(postCSD_opts{x3Dcsd(1)}.minDev*[1,1], [0,1], 'Color','r', 'LineStyle','--'); % h(k+1) = 
xlabel('Deviance Explained'); ylabel('Fraction of Units');
title( sprintf('%s Fit Results', GLMname), 'Interpreter','none' );

subtightplot(1,3,2,opt{:});
JitterPlot( 1 - postCSDdevFracPool', 0.5, 'ErrorCap',10, 'monochrome',0.6); hold on;
line([0,size(postCSDdevFracPool,1)+1], [0,0], 'color','k');
axis square;
ylabel('Relative Explanatory Value'); %ylabel('Cumulative Distribution');
ylim([-1,1]);
set(gca,'Xtick', 1:size(postCSDdevFracPool,1), 'XtickLabel',[postCSD_pred{x3Dcsd(1)}.name, postCSD_pred{x3Dcsd(1)}.fam.name]);
xtickangle(30);
title('Well-Fit Units');

subtightplot(1,3,3,opt{:});
bar(postCSD_subtypeFrac,'stacked');
ylabel('Fraction of Units');
barPos = get(gca, 'Position');
xlim([0, size(postCSD_subtypeFrac,1)+1]);
axis square;
title('Subtype Breakdown');
legend('Insensitive','Mixed','Deform-only','Loco-only', 'Location','EastOutside');
set(gca, 'Xtick',1:size(postCSD_subtypeFrac,1), 'XtickLabel',{expt(x3Dcsd).name, 'Pooled'}, 'TickLabelInterpreter','none', 'FontSize',10, 'Position',barPos );
xtickangle(30);

% Save the figure
figPath = sprintf('%s%s_DevianceResults.tif', figDir, GLMname);
if exist(figPath,'file'), delete(figPath); end
fprintf('\nSaving %s', figPath);
print(DevianceFig, figPath, '-dtiff', '-r300'); %pause(1); clf;

%% Show mean coeff by type
Ntype = 4;
subtypeColor = distinguishable_colors(Ntype);
close all;
SubtypeCoeffFig = figure('WindowState','maximized', 'color','w');
opt = {[0.02,0.07], [0.1,0.1], [0.1,0.06]};  % {[vert, horz], [bottom, top], [left, right] }
LW = 1.5;
colororder( subtypeColor ) 
for x = x3Dcsd
    tempInsensitiveCoeff = cat(3, postCSD_result{x}( postCSD_summary{x}.rIns ).coeff );
    meanInsensitiveCoeff = mean(tempInsensitiveCoeff, 3, 'omitnan' );
    
    tempMixedCoeff = cat(3, postCSD_result{x}( postCSD_summary{x}.rMixed ).coeff );
    meanMixedCoeff = mean(tempMixedCoeff, 3, 'omitnan' );
    
    tempDeformCoeff = cat(3, postCSD_result{x}( postCSD_summary{x}.rDeform ).coeff );
    meanDeformCoeff = mean(tempDeformCoeff, 3, 'omitnan' );
    
    tempLocoCoeff = cat(3, postCSD_result{x}( postCSD_summary{x}.rLoco ).coeff );
    meanLocoCoeff = mean(tempLocoCoeff, 3, 'omitnan' );
    
    sp(1) = subtightplot(1,4,1,opt{:});
    if postCSD_summary{x}.nIns > 0
        plot(postCSD_opts{x}.lags,  meanInsCoeff, 'LineWidth',LW );
    end
    axis square;
    tempPos = get(gca,'Position');
    xlabel('Lag (s)'); ylabel('Coefficient'); 
    title( sprintf('Insensitive (n = %i)', postCSD_summary{x}.nIns) );
    legend(postCSD_pred{x}.name, 'Location','NorthWest', 'AutoUpdate',false)
    set(gca,'Position',tempPos);
    
    sp(2) = subtightplot(1,4,2,opt{:});
    if postCSD_summary{x}.nMixed > 0
        plot(postCSD_opts{x}.lags,  meanMixedCoeff, 'LineWidth',LW ); 
    end
    axis square;
    title( sprintf('Mixed (n = %i)', postCSD_summary{x}.nMixed) );
    xlabel('Lag (s)'); 
    
    sp(3) = subtightplot(1,4,3,opt{:});
    if postCSD_summary{x}.nDeform > 0
        plot(postCSD_opts{x}.lags,  meanDeformCoeff, 'LineWidth',LW ); 
    end
    axis square;
    title( sprintf('Deformation-dependent (n = %i)', postCSD_summary{x}.nDeform) );
    xlabel('Lag (s)'); %title('Deformation-dependent'); % ylabel('Coefficient');
    
    sp(4) = subtightplot(1,4,4,opt{:});
    if postCSD_summary{x}.nLoco > 0
        plot(postCSD_opts{x}.lags,  meanLocoCoeff, 'LineWidth',LW );
    end
    axis square;
    xlabel('Lag (s)'); 
    title(sprintf('Locomotion-dependent (n = %i)', postCSD_summary{x}.nLoco)); % ylabel('Coefficient'); 
    linkaxes(sp,'xy');
    
    pause;
    
    % Save the figure
    figPath = sprintf('%s%s_%s_SubtypeCoeff.tif', figDir, GLMname, expt(x).name);
    if exist(figPath,'file'), delete(figPath); end
    %print(SubtypeCoeffFig, figPath, '-dtiff', '-r300' );  fprintf('\nSaved %s\n', figPath);

    clf;
end

%% Show  coeff by subtype
subtype = {'Insensitive', 'Mixed', 'Deform', 'Loco'};
Nsubtype = 4;
pooledCoeff = struct('Insensitive',[], 'Mixed',[], 'Deform',[], 'Loco',[]);
for x = x3Dcsd
    pooledCoeff.Insensitive = cat(3, pooledCoeff.Insensitive, postCSD_result{x}( postCSD_summary{x}.rIns ).coeff );
    pooledCoeff.Mixed = cat(3, pooledCoeff.Mixed, postCSD_result{x}( postCSD_summary{x}.rMixed).coeff );
    pooledCoeff.Deform = cat(3, pooledCoeff.Deform, postCSD_result{x}( postCSD_summary{x}.rDeform ).coeff );
    pooledCoeff.Loco = cat(3, pooledCoeff.Loco, postCSD_result{x}( postCSD_summary{x}.rLoco ).coeff );
    
end
pooledLags = postCSD_opts{x3Dcsd(1)}.lags;

Ncol = postCSD_pred{x3Dcsd(1)}.N;
k = 0;
close all;
opt = {[0.03,0.03], [0.07,0.05], [0.05,0.03]};  % {[vert, horz], [bottom, top], [left, right] }
SubtypeCoeffFig = figure('WindowState','maximized', 'color','w');
for row = 1:Ntype
    for col = 1:postCSD_pred{x3Dcsd(1)}.N
        k = k + 1;
        subtightplot(Ntype, Ncol, k, opt{:});
        plot( pooledLags, squeeze(pooledCoeff.(subtype{row})(:,col,:) ), 'color',[0,0,0,0.1] );
        axis square;
        if row == 1, title( postCSD_pred{x3Dcsd(1)}.name{col} ); end
        if col == 1, ylabel( sprintf('%s coeff', subtype{row} )); end
        if row == Nsubtype, xlabel('Lag (s)'); end
    end
    %pause;
end
% Save the figure
figPath = sprintf('%s%s_SubtypeCoeff.tif', figDir, GLMname);
if exist(figPath,'file'), delete(figPath); end
print(SubtypeCoeffFig, figPath, '-dtiff', '-r300' );  fprintf('\nSaved %s\n', figPath);

%%
devSummMat = cell(1,Nexpt); devSummPool = [];
opt = {[0.02,0.05], [0.02,0.02], [0.06,0.03]};  % {[vert, horz], [bottom, top], [left, right] }
GLMresultFig = figure('WindowState','maximized', 'color','w');

for x = x3Dcsd
    cla;
    devSummMat{x} = [postCSD_summary{x}.dev; postCSD_summary{x}.lopo.dev; postCSD_summary{x}.lofo.dev];
    %devSummPool = [devSummPool, devSummMat{x}];
    % {
    subtightplot(1,1,1,opt{:});
    imagesc( devSummMat{x} ); %imagesc( devSummPool ); %
    axis image;
    set(gca, 'Ytick', 1:size(devSummMat{x}, 1), 'YtickLabel', [{'All'}; postCSD_summary{x}.lopo.name(:); postCSD_summary{x}.lofo.name(:)], ...
        'Xtick',1:10:postCSD_resp{x}.N, 'TickDir','out', 'TickLength',[0.003,0], 'FontSize',8); % 'XtickLabel',postCSD_resp{x}.name
    title(sprintf('%s: GLM Fit Summary', postCSD_opts{x}.name), 'Interpreter','none');
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
set(gca,'XtickLabel',  [{'All'}; postCSD_summary{x}.lopo.name(:); postCSD_summary{x}.lofo.name(:)], ...
    'YtickLabel',postCSD_resp{x}.name, 'TickDir','out', 'TickLength',[0.003,0]);
title('GLM Fit Summary (All 3D Data)', 'Interpreter','none');
xtickangle(30);
axis image;
CB = colorbar; CB.Label.String = 'Median Deviance Explained';
figPath = sprintf('%s%s_DevSummary.tif', figDir, GLMname );
fprintf('\nSaving %s\n', figPath);
print( GLMdevFig, figPath, '-dtiff');
impixelinfo;

%% Plot coefficient values for each predictor, unit and experiment
close all; clearvars sp SP;
PreGLMcoeff = figure('WindowState','maximized', 'color','w');
opt = {[0.03,0.04], [0.09,0.05], [0.04,0.02]};  % {[vert, horz], [bottom, top], [left, right] }
for v = 1:postCSD_pred{x}.N
    c = 0;
    for x = x3Dcsd %find(~cellfun(@isempty, rLocoPreFit)) %
        c = c+1;
        if postCSD_summary{x}.Ngood > 0
            subtightplot(1, postCSD_pred{x}.N, v, opt{:});
            plot(c, postCSD_summary{x}.peakCoeff(postCSD_summary{x}.rGood,v), 'k.' ); hold on; %rLocoPreFit{x}
        end
    end
    line([0,c+1], [0,0], 'color','k');
    title( postCSD_pred{x}.name{v}, 'Interpreter','none' );
    if v == 1, ylabel('Peak Coefficient'); end
    xlim([0,c+1]);
    set(gca, 'Xtick', 1:c, 'XtickLabel', {expt(x3Dcsd).name}, 'TickLabelInterpreter','none' );
    xtickangle(45);
    axis square;
end

%% Plot peak coefficient vs latency values for each predictor, good unit, and experiment
close all; clearvars sp SP;
PreGLMcoeff = figure('WindowState','maximized', 'color','w');
opt = {[0.03,0.04], [0.09,0.05], [0.04,0.02]};  % {[vert, horz], [bottom, top], [left, right] }
exptColor = distinguishable_colors(numel(x3Dcsd));
for v = 1:postCSD_pred{x}.N
    c = 0;
    for x = x3Dcsd %find(~cellfun(@isempty, rLocoPreFit)) %
        c = c+1;
        if postCSD_summary{x}.Ngood > 0
            subtightplot(1, postCSD_pred{x}.N, v, opt{:});
            plot(postCSD_summary{x}.peakLag(postCSD_summary{x}.rGood,v), postCSD_summary{x}.peakCoeff(postCSD_summary{x}.rGood,v), '.', 'color',exptColor(c,:) ); hold on; %rLocoPreFit{x}
        end
    end
    title( postCSD_pred{x3Dcsd(1)}.name{v}, 'Interpreter','none' );
    if v == 1, ylabel('Coefficient'); end
    xlabel('Lag (s)');
    xlim([-6,6]);
    axis square;
end

%% Identify the best-fit units from each experiment
for x = x3Dcsd
    [devSort, rDevSort] = sort( [postCSD_result{x}.dev], 'descend' );
    %[maxDevExp, rMaxDevExp] = max( [postCSD_result{x}.dev] );
    WriteROIproj(expt(x), ROI{x}, 'edges',segParams{x}.edges, 'overwrite',true, 'rSet',rDevSort(1:10), 'buffer',20*[1,1,1,1]); % ROI{x} =   
    %ViewResults3D( expt(x), Tscan{x}, deform{x}, loco{x}, fluor{x}, allVars, ROI{x}, 'cat',true, 'limits', viewLims ); 
end


