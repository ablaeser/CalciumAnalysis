%% Use GLM to assess contribution of different variables
onPreResult = cell(1,Nexpt);
figDir = 'D:\MATLAB\LevyLab\Figures\3D\Loco\';
clearvars onPrePred onPreResp onPreResult onPreSumm
rLocoPreFit = cell(1,Nexpt);
for x = x3Dcsd
    % GLMparallel options
    onPreOpts(x).show = true;
    onPreOpts(x).rShow = [];
    onPreOpts(x).figDir = figDir;
    onPreOpts(x).alpha = 0.01;  % The regularization parameter, default is 0.01
    onPreOpts(x).standardize = true; 
    onPreOpts(x).trainFrac = 0.75; % 1; %
    onPreOpts(x).Ncycle = 20;
    onPreOpts(x).distribution = 'gaussian'; % 'poisson';
    onPreOpts(x).CVfold = 10;
    onPreOpts(x).minDev = 0.1;  minDev = onPreOpts(x).minDev;
    onPreOpts(x).minDevFrac = 0.1;
    onPreOpts(x).maxP = 0.05;
    onPreOpts(x).Nshuff = 0;  Nshuff = onPreOpts(x).Nshuff;
    onPreOpts(x).minShuff = 15; 
    onPreOpts(x).window = [0, 0]; % [0,0]; % [-0.5, 0.5]; % 
    onPreOpts(x).lopo = true; %false; %
    onPreOpts(x).frameRate = expt(x).scanRate; 
    onPreOpts(x).minShuffFrame = round( onPreOpts(x).frameRate*onPreOpts(x).minShuff );
    windowFrame = round(onPreOpts(x).window*onPreOpts(x).frameRate);
    onPreOpts(x).shiftFrame = windowFrame(1):windowFrame(2);
    onPreOpts(x).maxShift = max( abs(windowFrame) );
    onPreOpts(x).Nshift = numel( onPreOpts(x).shiftFrame );  %Nshift = onPreOpts(x).Nshift;
    onPreOpts(x).lags = onPreOpts(x).shiftFrame/onPreOpts(x).frameRate;
    onPreOpts(x).name = strcat(expt(x).name, '_onPreGLM');

    % Concatenate pre-CSD bouts
    onSubCat = [];
    for v = 1:NallVars
        subData = onStruct(x).(allVars{v}).data(:,:,onStruct(x).csd.pre.bout) - onStruct(x).(allVars{v}).pre(:,:,onStruct(x).csd.pre.bout); % subtract each bout's pre-onset mean
        onSubCat.(allVars{v}) = reshape( permute( subData, [1,3,2]), ...
            size(onStruct(x).(allVars{v}).data,1)*onStruct(x).csd.pre.Nbout, size(onStruct(x).(allVars{v}).data,2) );       
    end
    onSubCat.T = repmat( onStruct(x).T, onStruct(x).csd.pre.Nbout, 1);
    
    % Define predictors
    onPrePred(x) = struct('data',[], 'name',[], 'N',NaN, 'TB',[], 'lopo',[], 'fam',[]); % onPrePred(x).data = []; %onPrePred(x).data = onSubCat.T; onPrePred(x).name = {'Time'};
    for v = 2:NallVars
        onPrePred(x).data = [onPrePred(x).data, mean(onSubCat.(allVars{v}), 2, 'omitnan')]; %rand(1000,3); %
        onPrePred(x).name = [onPrePred(x).name, allVars(v)];
    end
    onPrePred(x).N = size(onPrePred(x).data,2);
    for p = flip(1:onPrePred(x).N), onPrePred(x).lopo.name{p} = ['No ',onPrePred(x).name{p}]; end
    onPrePred(x).fam.col = {1:2, 3:4, 5:6, 7:8, 9, 10};  %{1, 2:3, 4:5, 6:7, 8, 9}; 
    onPrePred(x).fam.N = numel(onPrePred(x).fam.col); 
    onPrePred(x).fam.name = {'Translation', 'Scale', 'Stretch', 'Shear', 'Z Shift', 'Speed'}; % 'Onset Time',
    
    % Define responses
    onPreResp(x).data = onSubCat.fluor; %[merged(x).scaleap-medAP, merged(x).scaleml-medML]; % , merged(x).dScAP, merged(x).dScML
    onPreResp(x).N = size(onPreResp(x).data,2); 
    onPreResp(x).name = repmat({'Fluor'}, 1, onPreResp(x).N ); %{'Scale-AP','Scale-ML'}; % , 'dScAP', 'dScML
    nanFrame = find(any(isnan([onPrePred(x).data, onPreResp(x).data]),2)); % find( isnan(sum(pred(x).data,2)) ); 
    fprintf('\nRemoved %i NaN-containing frames', numel(nanFrame));
    onPrePred(x).data(nanFrame,:) = []; onPreResp(x).data(nanFrame,:) = [];
    
    % Run the GLM
    onPreOpts(x).load = false; %  true; %
    onPreOpts(x).saveRoot = sprintf('%s%s', expt(x).dir, onPreOpts(x).name  ); %''; %
    [onPreResult{x}, onPreSumm(x), ~, onPrePred(x), onPreResp(x)] = GLMparallel(onPrePred(x), onPreResp(x), onPreOpts(x)); % onPreResult{x}, onPreSumm(x), ~, onPrePred(x), onPreResp(x)
    rLocoPreFit{x} = intersect(onStruct(x).fluor.responder, onPreSumm(x).rFit);
    nLocoPreFit(x) = numel( rLocoPreFit{x} );
end

%% Summarize GLM results
close all; clearvars sp SP;
OnPreGLMresults = figure('WindowState','maximized', 'color','w');
leftOpt = {[0.02,0.07], [0.06,0.03], [0.04,0.02]};  % {[vert, horz], [bottom, top], [left, right] }\
rightOpt = {[0.1,0.07], [0.1,0.03], [0.04,0.02]};  % {[vert, horz], [bottom, top], [left, right] }\
jitterWidth = 0.45;
xAngle = 30;
Nrow = onPrePred(x).N+1; Ncol = 3;
spGrid = reshape( 1:Nrow*Ncol, Ncol, Nrow )';
for x = xPresent
    subtightplot(onPrePred(x).N+1, 3, 1:2, leftOpt{:});
    imagesc( onPreResp(x).data' );
    ylabel('Fluor', 'Interpreter','none');
    set(gca,'TickDir','out', 'TickLength',[0.003,0], 'box','off', 'Ytick',onStruct(x).fluor.responder, 'Xtick',[]);
    text( repmat(size(onPreResp(x).data,1)+1, onPreSumm(x).Nfit, 1), onPreSumm(x).rFit+0.5, '*', 'VerticalAlignment','middle', 'FontSize',8);
    impixelinfo;
    
    for v = 1:onPrePred(x).N
        subtightplot(Nrow, Ncol, spGrid(v+1,1:2), leftOpt{:});
        plot( onPrePred(x).data(:,v) ); hold on;
        ylabel(onPrePred(x).name{v}, 'Interpreter','none');
        xlim([-Inf,Inf]);
        if v < onPrePred(x).N
            set(gca,'TickDir','out', 'TickLength',[0.003,0], 'box','off', 'XtickLabel',[]);
        else
            set(gca,'TickDir','out', 'TickLength',[0.003,0], 'box','off');
        end
    end
    xlabel('Scan');
    
    subtightplot(3,3,3, rightOpt{:});
    bar([numel(onStruct(x).fluor.responder), onPreSumm(x).Nfit, numel(rLocoPreFit{x})]/expt(x).Nroi );
    set(gca,'Xtick',1:3, 'XtickLabel',{'Loco','Fit','Both'}, 'box','off');
    ylabel('Fraction of ROI');
    ylim([0,1]);
    
    subtightplot(3,3,6, rightOpt{:});
    JitterPlot( onPreSumm(x).lopo.devFrac(onPreSumm(x).rFit,:), jitterWidth ); hold on;
    line([0,onPrePred(x).N], [1,1], 'color','k', 'lineStyle','--');
    xlim([0,onPrePred(x).N]);
    ylabel('Fraction of total deviance'); title('Leave One Predictor Out (well-fit units only)');
    set(gca, 'Xtick',1:onPrePred(x).N,  'XtickLabel', onPreSumm(x).lopo.name, 'TickDir','out', 'TickLength',[0.003,0], 'TickLabelInterpreter','none', 'box','off' ); 
    xtickangle(xAngle);
    ylim([0,Inf]);
    
    subtightplot(3,3,9, rightOpt{:});
    JitterPlot( onPreSumm(x).lofo.devFrac(onPreSumm(x).rFit,:), jitterWidth ); hold on;
    line([0,onPrePred(x).fam.N]+0.5, [1,1], 'color','k', 'lineStyle','--');
    xlim([0,onPrePred(x).fam.N]+0.5);
    ylabel('Fraction of total deviance'); title('Leave One Family Out (well-fit units only)');
    set(gca, 'Xtick',1:onPrePred(x).fam.N,  'XtickLabel', onPreSumm(x).lofo.name, 'TickDir','out', 'TickLength',[0.003,0], 'TickLabelInterpreter','none', 'box','off' ); 
    xtickangle(xAngle);
    ylim([0,Inf]);
    
    % {
    figPath = sprintf('%s%s.tif', figDir, onPreOpts(x).name);
    if exist(figPath,'file'), delete(figPath); end
    fprintf('\nSaving %s', figPath);
    %export_fig( figPath, '-pdf', '-painters','-q101', '-append', LocoSensitivePrePost ); pause(1);
    print(OnPreGLMresults, figPath, '-dtiff' ); pause(1);   
    clf;
    %}
    %pause; clf;
end

%% Across experiments, summarize fraction of units locomotion sensitive, well-fit and both

preSensFitFrac = nan(0,3);
c = 0;
for x = xPresent
    c = c+1;
    preSensFitFrac(c,:) = [numel(rLocoPre{x}), onPreSumm(x).Nfit, numel(intersect( rLocoPre{x}, onPreSumm(x).rFit ))]/expt(x).Nroi;
end

GLMlocoFitExpt = figure('WindowState','maximized', 'color','w');
bar(100*preSensFitFrac);
ylabel('% of ROI');
legend('Locomotion-driven', 'Well-fit', 'Both', 'Location','NorthWest');
set(gca,'Xtick', 1:numel(xPresent), 'XtickLabel', {expt(xPresent).name}, 'TickLabelInterpreter','none', 'box','off', 'TickDir','out' ); % ,
xtickangle(30);
figPath = sprintf('%sGLMlocoFitExpt.tif', figDir);
if exist(figPath,'file'), delete(figPath); end
fprintf('\nSaving %s', figPath);
print(GLMlocoFitExpt, figPath, '-dtiff' );

%% Plot coefficient values for each predictor, unit and experiment
close all; clearvars sp SP;
OnPreGLMcoeff = figure('WindowState','maximized', 'color','w');
opt = {[0.03,0.04], [0.09,0.05], [0.04,0.02]};  % {[vert, horz], [bottom, top], [left, right] }
for v = 1:onPrePred(x).N
    c = 0;
    for x = find(~cellfun(@isempty, rLocoPreFit)) %xPresent
        c = c+1;
        subtightplot(1, onPrePred(x).N, v, opt{:});
        plot(c, onPreSumm(x).peakCoeff(rLocoPreFit{x},v), 'k.' ); hold on;
    end
    line([0,c+1], [0,0], 'color','k');
    title( onPrePred(x).name{v}, 'Interpreter','none' );
    if v == 1, ylabel('Coefficient'); end
    xlim([0,c+1]);
    set(gca, 'Xtick', 1:c, 'XtickLabel', {expt(xPresent).name}, 'TickLabelInterpreter','none' );
    xtickangle(45);
    axis square;
end

%% LOFO/LOPO results, pooling well-fit, locomotion-sensitive ROI across experiments
lopoPool = []; lofoPool = [];
for x = find(~cellfun(@isempty, rLocoPreFit))
    lopoPool = vertcat( lopoPool, onPreSumm(x).lopo.devFrac(rLocoPreFit{x}, :) );
    lofoPool = vertcat( lofoPool, onPreSumm(x).lofo.devFrac(rLocoPreFit{x}, :) );
end
close all; clearvars sp SP;
OnPreLeaveOneOutResults = figure('WindowState','maximized', 'color','w');
opt = {[0.03,0.07], [0.09,0.05], [0.08,0.07]};  % {[vert, horz], [bottom, top], [left, right] }

subtightplot(1,2,1,opt{:});
JitterPlot( lopoPool, 0.45 );
line([0,onPrePred(1).N]+0.5, [1,1], 'color','k', 'lineStyle','--');
xlim([0,onPrePred(1).N]+0.5);
xtickangle(xAngle);
set(gca,'Xtick',1:onPrePred(1).N, 'XtickLabel',onPrePred(1).name, 'TickLabelInterpreter','none' );
title('Leave One Predictor Out (locomotion-sensitive, well-fit units only)'); ylabel('Fraction of total deviance');
axis square;

subtightplot(1,2,2,opt{:});
JitterPlot( lofoPool, 0.45 );
line([0,onPrePred(1).fam.N]+0.5, [1,1], 'color','k', 'lineStyle','--');
xlim([0,onPrePred(1).fam.N]+0.5);
set(gca,'Xtick',1:onPrePred(1).fam.N, 'XtickLabel',onPrePred(1).fam.name, 'TickLabelInterpreter','none' );
title('Leave One Family Out'); ylabel('Fraction of total deviance');
axis square;

figPath = sprintf('%sLeaveOneOutResults_onsetPreCSD.tif', figDir);
if exist(figPath,'file'), delete(figPath); end
fprintf('\nSaving %s', figPath);
print(OnPreLeaveOneOutResults, figPath, '-dtiff' ); %pause(1); clf;