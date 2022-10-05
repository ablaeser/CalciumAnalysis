%% Use GLM to assess contribution of different variables
locoDeform_pre_pred = cell(1,Nexpt); locoDeform_pre_resp = cell(1,Nexpt); locoDeform_pre_opts = cell(1,Nexpt); locoDeform_pre_result = cell(1,Nexpt); locoDeform_pre_summary = cell(1,Nexpt);
figDir = 'D:\MATLAB\LevyLab\Figures\3D\GLM\'; mkdir( figDir )
%rLocoPreFit = cell(1,Nexpt);
GLMname = 'LocoDeform_3D_preCSD';
GLMrate = 15.49/30;
for x = xPresent % x3Dcsd % 
    % GLMparallel options
    locoDeform_pre_opts{x}.name = sprintf('%s_%s', expt(x).name, GLMname); %strcat(expt(x).name, , '_locoDeform_preglm');
    %locoDeform_pre_opts{x}.show = true;
    locoDeform_pre_opts{x}.rShow = NaN;
    locoDeform_pre_opts{x}.figDir = ''; % figDir;
    locoDeform_pre_opts{x}.alpha = 0.01;  % The regularization parameter, default is 0.01
    locoDeform_pre_opts{x}.standardize = true; 
    locoDeform_pre_opts{x}.trainFrac = 0.75; % 1; %
    locoDeform_pre_opts{x}.Ncycle = 20;
    locoDeform_pre_opts{x}.distribution = 'gaussian'; % 'poisson'; %  
    locoDeform_pre_opts{x}.CVfold = 10;
    locoDeform_pre_opts{x}.nlamda = 1000;
    locoDeform_pre_opts{x}.maxit = 5*10^5;
    locoDeform_pre_opts{x}.minDev = 0.05;  
    locoDeform_pre_opts{x}.minDevFrac = 0.1;
    locoDeform_pre_opts{x}.maxP = 0.05;
    locoDeform_pre_opts{x}.Nshuff = 0;  Nshuff = locoDeform_pre_opts{x}.Nshuff;
    locoDeform_pre_opts{x}.minShuff = 15; 
    locoDeform_pre_opts{x}.window = [-60,60]; % [0,0]; % [-0.5, 0.5]; % 
    locoDeform_pre_opts{x}.lopo = true; %false; %
    locoDeform_pre_opts{x}.frameRate = GLMrate; %expt(x).scanRate; 
    locoDeform_pre_opts{x}.binSize = expt(x).scanRate/GLMrate;
    locoDeform_pre_opts{x}.minShuffFrame = round( locoDeform_pre_opts{x}.frameRate*locoDeform_pre_opts{x}.minShuff );
    windowFrame = [ceil(locoDeform_pre_opts{x}.window(1)*locoDeform_pre_opts{x}.frameRate), floor(locoDeform_pre_opts{x}.window(2)*locoDeform_pre_opts{x}.frameRate)]; % floor(locoDeform_pre_opts{x}.window*locoDeform_pre_opts{x}.frameRate);
    locoDeform_pre_opts{x}.shiftFrame = windowFrame(1):windowFrame(2);
    locoDeform_pre_opts{x}.maxShift = max( abs(windowFrame) );
    locoDeform_pre_opts{x}.Nshift = numel( locoDeform_pre_opts{x}.shiftFrame );  %Nshift = locoDeform_preOpts(x).Nshift;
    locoDeform_pre_opts{x}.lags = locoDeform_pre_opts{x}.shiftFrame/locoDeform_pre_opts{x}.frameRate;

    % Define predictors
    tempVelocityCat = BinDownMean( vertcat(loco{x}(expt(x).preRuns).Vdown), locoDeform_pre_opts{x}.binSize );
    tempAccelCat = BinDownMean( vertcat(loco{x}(expt(x).preRuns).Adown), locoDeform_pre_opts{x}.binSize ); %[NaN; diff(tempVelocityCat)];
    tempSpeedCat = BinDownMean( vertcat(loco{x}(expt(x).preRuns).speedDown), locoDeform_pre_opts{x}.binSize );
    tempStateCat = BinDownMean( vertcat(loco{x}(expt(x).preRuns).stateDown), locoDeform_pre_opts{x}.binSize );
    locoDeform_pre_pred{x} = struct('data',[], 'name',[], 'N',NaN, 'TB',[], 'lopo',[], 'fam',[]); % LocoDeformPred(x).data = []; %LocoDeformPred(x).data = LocoDeformCat.T; LocoDeformPred(x).name = {'Time'};
    locoDeform_pre_pred{x}.data = [tempVelocityCat, tempAccelCat, tempStateCat]; % tempSpeedCat,, 
    locoDeform_pre_pred{x}.name = {'Velocity','Acceleration','State'}; % ,'Speed',
    locoDeform_pre_pred{x}.N = size(locoDeform_pre_pred{x}.data,2);
    for p = flip(1:locoDeform_pre_pred{x}.N), locoDeform_pre_pred{x}.lopo.name{p} = ['No ',locoDeform_pre_pred{x}.name{p}]; end
    % Set up leave-one-family-out
    locoDeform_pre_pred{x}.fam.col = {1:2, 3}; %{1:2, 3:4, 5:6, 7:8, 9:10, 11:12};  % {1:12};%{1, 2:3, 4:5, 6:7, 8, 9}; 
    locoDeform_pre_pred{x}.fam.N = numel(locoDeform_pre_pred{x}.fam.col); 
    locoDeform_pre_pred{x}.fam.name = {'Kinematics', 'State'}; %{'All'};%  'Onset Time',
    
    % Define responses
    tempTrans = vertcat(deform{x}(expt(x).preRuns).transMag);
    tempTrans = BinDownMean( mean( tempTrans(:,segParams{x}.zProj), 2, 'omitnan'), locoDeform_pre_opts{x}.binSize ); % (:,segParams{x}.zProj)
    tempTransSpd = vertcat(deform{x}(expt(x).preRuns).DtransMag);
    tempTransSpd = BinDownMean( mean( tempTransSpd(:,segParams{x}.zProj), 2, 'omitnan'), locoDeform_pre_opts{x}.binSize ); % vertcat(deform{x}(expt(x).preRuns).DtransMag)
    tempScaleMag = vertcat(deform{x}(expt(x).preRuns).scaleMag);
    tempScaleMag = BinDownMean( mean( tempScaleMag(:,segParams{x}.zProj), 2, 'omitnan'), locoDeform_pre_opts{x}.binSize );
    tempStretchMag = vertcat(deform{x}(expt(x).preRuns).stretchMag);
    tempStretchMag = BinDownMean( mean( tempStretchMag(:,segParams{x}.zProj), 2, 'omitnan'), locoDeform_pre_opts{x}.binSize );
    tempShearMag = vertcat(deform{x}(expt(x).preRuns).shearMag);
    tempShearMag = BinDownMean( mean( tempShearMag(:,segParams{x}.zProj), 2, 'omitnan'), locoDeform_pre_opts{x}.binSize );
    tempShearRate = vertcat(deform{x}(expt(x).preRuns).DshearMag);
    tempShearRate = BinDownMean( mean( tempShearRate(:,segParams{x}.zProj), 2, 'omitnan'), locoDeform_pre_opts{x}.binSize );
    tempShift = vertcat(deform{x}(expt(x).preRuns).shiftZ);
    tempShiftMean =  BinDownMean( mean( tempShift(:,4:end-3), 2, 'omitnan' ), locoDeform_pre_opts{x}.binSize );

    locoDeform_pre_resp{x}.data = [tempTrans, tempScaleMag, tempShearMag, tempTransSpd, tempStretchMag, tempShearRate, tempShiftMean]; % tempScaleMag; %  , tempDshearMag tempExp, tempComp,  , tempCompressMean
    tempNormData = normalize(locoDeform_pre_resp{x}.data, 1); % locoDeform_post_resp{x}.data
    locoDeform_pre_resp{x}.data(abs(tempNormData) > 5) = NaN; % suppress extreme outliers
    locoDeform_pre_resp{x}.name = {'|Translation|', '|Scale|', '|Shear|', 'Trans. Speed', '|Stretch|', 'Shear Rate', 'Z Shift'}; % , 'Z Shift', 'Compression' , 'Z Compression' 
    locoDeform_pre_resp{x}.N = size(locoDeform_pre_resp{x}.data,2); 
    
    % Remove scans with missing data 
    nanFrame = find(any(isnan([locoDeform_pre_pred{x}.data, locoDeform_pre_resp{x}.data]),2)); % find( isnan(sum(pred(x).data,2)) ); 
    fprintf('\nRemoving %i NaN-containing frames', numel(nanFrame));
    locoDeform_pre_pred{x}.data(nanFrame,:) = []; locoDeform_pre_resp{x}.data(nanFrame,:) = [];

    % Run the GLM
    locoDeform_pre_opts{x}.load = true; %  false; % 
    locoDeform_pre_opts{x}.saveRoot = expt(x).dir; % sprintf('%s%s', expt(x).dir, locoDeform_pre_opts{x}.name  ); % ''; %  
    [locoDeform_pre_result{x}, locoDeform_pre_summary{x}, locoDeform_pre_opts{x}, locoDeform_pre_pred{x}, locoDeform_pre_resp{x}] = GLMparallel(locoDeform_pre_pred{x}, locoDeform_pre_resp{x}, locoDeform_pre_opts{x}); % LocoDeform_result{x}, LocoDeform_summary(x), ~, LocoDeformPred(x), LocoDeformResp(x)
    %ViewGLM(LocoDeform_pred{x}, LocoDeform_resp{x}, locoDeform_pre_opts{x}, LocoDeform_result{x}, LocoDeform_summary{x}); %GLMresultFig = 
    %pause;
end

%% View GLM results
for x = 30 %xPresent
    locoDeform_pre_opts{x}.rShow = 1:locoDeform_pre_resp{x}.N; % 3; % 1:LocoDeform_resp{x}.N; %NaN;
    locoDeform_pre_opts{x}.xVar = 'Time';
    ViewGLM(locoDeform_pre_pred{x}, locoDeform_pre_resp{x}, locoDeform_pre_opts{x}, locoDeform_pre_result{x}, locoDeform_pre_summary{x}); %GLMresultFig = 
end

%%
devSummMat = cell(1,Nexpt);
GLMresultFig = figure('WindowState','maximized', 'color','w');
for x = xPresent
    cla;
    devSummMat{x} = [locoDeform_pre_summary{x}.dev; locoDeform_pre_summary{x}.lopo.dev; locoDeform_pre_summary{x}.lofo.dev];
    % {
    imagesc( devSummMat{x}' );
    set(gca,'XtickLabel',  [{'All'}; locoDeform_pre_summary{xPresent(1)}.lopo.name(:); locoDeform_pre_summary{xPresent(1)}.lofo.name(:)], ...
        'Ytick',1:locoDeform_pre_resp{xPresent(1)}.N , 'YtickLabel',locoDeform_pre_resp{xPresent(1)}.name, 'TickDir','out', 'TickLength',[0.003,0]);
    title(sprintf('%s: GLM Fit Summary', locoDeform_pre_opts{x}.name), 'Interpreter','none');
    xtickangle(30);
    axis image;
    CB = colorbar; CB.Label.String = 'Deviance Explained';
    caxis([0,0.2]);
    impixelinfo;
    pause;
    %}
end

devSummCat = cat(3, devSummMat{:});
devSummMed = median(devSummCat, 3, 'omitnan' );

locoDeform_pre_GLMdevFig = figure('WindowState','maximized', 'color','w');
imagesc( devSummMed' );
set(gca,'XtickLabel',  [{'All'}; locoDeform_pre_summary{x}.lopo.name(:); locoDeform_pre_summary{x}.lofo.name(:)], ...
    'YtickLabel',locoDeform_pre_resp{x}.name, 'TickDir','out', 'TickLength',[0.003,0]);
title('GLM Fit Summary (All 3D Data)', 'Interpreter','none');
xtickangle(30);
axis image;
CB = colorbar; CB.Label.String = 'Median Deviance Explained';
figPath = sprintf('%s%s_DevSummary.tif', figDir, GLMname );
fprintf('\nSaving %s\n', figPath);
%print( locoDeform_pre_GLMdevFig, figPath, '-dtiff');
impixelinfo;

%% Show relative explanatory value from submodels
exptColor = distinguishable_colors(Npresent);
locoDeform_pre_relExpFig = figure('WindowState','maximized', 'color','w');
Ncol = locoDeform_pre_resp{xPresent(1)}.N;  
%tiledlayout(1, Ncol);
opt = {[0.05,0.04], [0.06,0.02], [0.04,0.02]};
Ndrops = locoDeform_pre_pred{xPresent(1)}.N+locoDeform_pre_pred{xPresent(1)}.fam.N;
for col = 1:Ncol
    tempDevFracMat = nan( Npresent, Ndrops );
    k = 0;
    for x = xPresent
        k = k + 1;
        if locoDeform_pre_result{x}(col).dev > locoDeform_pre_opts{x}.minDev
            tempDevFracMat(k,:) = [locoDeform_pre_result{x}(col).lopo.devFrac, locoDeform_pre_result{x}(col).lofo.devFrac];
        end
    end
    subtightplot(1,Ncol,col,opt{:}); %nexttile
    JitterPlot( 1-tempDevFracMat, 0.45 ); hold on;
    line(0.5+[0,Ndrops], [0,0], 'color','k'); 
    set(gca, 'Xtick',1:Ndrops); %, 'Xticklabel'
    ylim([-0.4, 1.1]); xlim(0.5+[0,Ndrops])
    %colororder(exptColor)
    axis square;
    if col == 1, ylabel('Relative Explanatory Value'); end
    title( locoDeform_pre_resp{xPresent(1)}.name{col} );
    set(gca, 'XtickLabel', [locoDeform_pre_pred{xPresent(1)}.name, locoDeform_pre_pred{xPresent(1)}.fam.name])  % [locoDeform_pre_summary{x}.lopo.name(:); locoDeform_pre_summary{x}.lofo.name(:)]
    xtickangle(45)
end
figPath = sprintf('%s%s_ExpValSummary.tif', figDir, GLMname );
fprintf('\nSaving %s\n', figPath);
print( locoDeform_pre_relExpFig, figPath, '-dtiff');

%% Show time filters per response
exptColor = distinguishable_colors(Npresent);
locoDeform_pre_GLMcoeffFig = figure('WindowState','maximized', 'color','w');
Ncol = locoDeform_pre_resp{xPresent(1)}.N;   Nrow = locoDeform_pre_pred{xPresent(1)}.N;
tiledlayout(Nrow, Ncol);
k = 0;
for row = 1:Nrow
    for col = 1:Ncol
        nexttile
        tempCoeffMat = zeros( numel(locoDeform_pre_opts{xPresent(1)}.lags), Npresent );
        k = 0;
        for x = xPresent
            k = k + 1;
            if locoDeform_pre_result{x}(col).dev > locoDeform_pre_opts{x}.minDev
                tempCoeffMat(:,k) = locoDeform_pre_result{x}(col).coeff(:,row);
            end
        end
        plot( locoDeform_pre_opts{x}.lags, tempCoeffMat ); hold on;
        colororder(exptColor)
        axis square;
        if col == 1, ylabel(sprintf('%s Coefficient', locoDeform_pre_pred{xPresent(1)}.name{row} )); end
        if row == 1, title( locoDeform_pre_resp{xPresent(1)}.name{col} ); end
        if row == Nrow, xlabel('Response Lag (s)'); end
        %pause;
    end
end

figPath = sprintf('%s%s_CoeffSummary.tif', figDir, GLMname );
fprintf('\nSaving %s\n', figPath);
print( locoDeform_pre_GLMcoeffFig, figPath, '-dtiff');
