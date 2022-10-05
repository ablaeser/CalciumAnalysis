%% Use GLM to assess contribution of different variables
locoDeform_csd_pred = cell(1,Nexpt); locoDeform_csd_resp = cell(1,Nexpt); locoDeform_csd_opts = cell(1,Nexpt); locoDeform_csd_result = cell(1,Nexpt); locoDeform_csd_summary = cell(1,Nexpt);
figDir = 'D:\MATLAB\LevyLab\Figures\3D\GLM\'; mkdir( figDir )
GLMname = 'LocoDeform_CSD';
GLMrate = 15.49/30;
for x = [5, 9, 28, 31, 37, 41] %xWave%x3Dcsd %xPresent %
    % GLMparallel options
    locoDeform_csd_opts{x}.name = sprintf('%s_%s', expt(x).name, GLMname); %strcat(expt(x).name, , '_locoDeform_csdglm');
    %locoDeform_csd_opts{x}.show = true;
    locoDeform_csd_opts{x}.rShow = NaN;
    locoDeform_csd_opts{x}.figDir = ''; % figDir;
    locoDeform_csd_opts{x}.alpha = 0.01;  % The regularization parameter, default is 0.01
    locoDeform_csd_opts{x}.standardize = true; 
    locoDeform_csd_opts{x}.trainFrac = 0.75; % 1; %
    locoDeform_csd_opts{x}.Ncycle = 20;
    locoDeform_csd_opts{x}.distribution = 'gaussian'; % 'poisson'; %  
    locoDeform_csd_opts{x}.CVfold = 10;
    locoDeform_csd_opts{x}.nlamda = 1000;
    locoDeform_csd_opts{x}.maxit = 5*10^5;
    locoDeform_csd_opts{x}.minDev = 0.05;  
    locoDeform_csd_opts{x}.minDevFrac = 0.1;
    locoDeform_csd_opts{x}.maxP = 0.05;
    locoDeform_csd_opts{x}.Nshuff = 0;  Nshuff = locoDeform_csd_opts{x}.Nshuff;
    locoDeform_csd_opts{x}.minShuff = 15; 
    locoDeform_csd_opts{x}.window = [-60,60]; % [0,0]; % [-0.5, 0.5]; % 
    locoDeform_csd_opts{x}.lopo = true; %false; %
    locoDeform_csd_opts{x}.frameRate = GLMrate; %expt(x).scanRate; 
    locoDeform_csd_opts{x}.binSize = expt(x).scanRate/GLMrate;
    locoDeform_csd_opts{x}.minShuffFrame = round( locoDeform_csd_opts{x}.frameRate*locoDeform_csd_opts{x}.minShuff );
    windowFrame = [ceil(locoDeform_csd_opts{x}.window(1)*locoDeform_csd_opts{x}.frameRate), floor(locoDeform_csd_opts{x}.window(2)*locoDeform_csd_opts{x}.frameRate)]; % floor(locoDeform_csd_opts{x}.window*locoDeform_csd_opts{x}.frameRate);
    locoDeform_csd_opts{x}.shiftFrame = windowFrame(1):windowFrame(2);
    locoDeform_csd_opts{x}.maxShift = max( abs(windowFrame) );
    locoDeform_csd_opts{x}.Nshift = numel( locoDeform_csd_opts{x}.shiftFrame );  %Nshift = locoDeform_csdOpts(x).Nshift;
    locoDeform_csd_opts{x}.lags = locoDeform_csd_opts{x}.shiftFrame/locoDeform_csd_opts{x}.frameRate;

    % Concatenate and bin input variables 
    tempVelocityCat = BinDownMean( vertcat(loco{x}(expt(x).csd).Vdown), locoDeform_csd_opts{x}.binSize );
    tempAccelCat = BinDownMean( vertcat(loco{x}(expt(x).csd).Adown), locoDeform_csd_opts{x}.binSize ); %[NaN; diff(tempVelocityCat)];
    tempSpeedCat = BinDownMean( vertcat(loco{x}(expt(x).csd).speedDown), locoDeform_csd_opts{x}.binSize );
    tempStateCat = BinDownMean( vertcat(loco{x}(expt(x).csd).stateDown), locoDeform_csd_opts{x}.binSize );
    tempCSDcat = ones(expt(x).Nscan(expt(x).csd),1);
    tempCSDcat( csdBout{x}(expt(x).csd).boutScan{1} ) = 2; 
    tempCSDcat = BinDownMean( tempCSDcat, locoDeform_csd_opts{x}.binSize );
    
    % Define predictors
    locoDeform_csd_pred{x} = struct('data',[], 'name',[], 'N',NaN, 'TB',[], 'lopo',[], 'fam',[]); % LocoDeformPred(x).data = []; %LocoDeformPred(x).data = LocoDeformCat.T; LocoDeformPred(x).name = {'Time'};
    locoDeform_csd_pred{x}.data = [tempVelocityCat, tempAccelCat, tempStateCat, tempCSDcat]; % tempSpeedCat,, 
    locoDeform_csd_pred{x}.name = {'Velocity','Acceleration','Loco','CSD'}; % ,'Speed',
    locoDeform_csd_pred{x}.N = size(locoDeform_csd_pred{x}.data,2);
    for p = flip(1:locoDeform_csd_pred{x}.N), locoDeform_csd_pred{x}.lopo.name{p} = ['No ',locoDeform_csd_pred{x}.name{p}]; end
    % Set up leave-one-family-out
    locoDeform_csd_pred{x}.fam.col = {1:2, 3, 4}; %{1:2, 3:4, 5:6, 7:8, 9:10, 11:12};  % {1:12};%{1, 2:3, 4:5, 6:7, 8, 9}; 
    locoDeform_csd_pred{x}.fam.N = numel(locoDeform_csd_pred{x}.fam.col); 
    locoDeform_csd_pred{x}.fam.name = {'Kinematics','Loco','CSD'}; %{'All'};%  'Onset Time',

    % Define responses
    tempTrans = vertcat(deform{x}(expt(x).csd).transMag);
    tempTrans = BinDownMean( mean( tempTrans(:,segParams{x}.zProj), 2, 'omitnan'), locoDeform_csd_opts{x}.binSize ); % (:,segParams{x}.zProj)
    tempTransSpd = vertcat(deform{x}(expt(x).csd).DtransMag);
    tempTransSpd = BinDownMean( mean( tempTransSpd(:,segParams{x}.zProj), 2, 'omitnan'), locoDeform_csd_opts{x}.binSize ); % vertcat(deform{x}(expt(x).csd).DtransMag)
    tempScaleMag = vertcat(deform{x}(expt(x).csd).scaleMag);
    tempScaleMag = BinDownMean( mean( tempScaleMag(:,segParams{x}.zProj), 2, 'omitnan'), locoDeform_csd_opts{x}.binSize );
    tempStretchMag = vertcat(deform{x}(expt(x).csd).stretchMag);
    tempStretchMag = BinDownMean( mean( tempStretchMag(:,segParams{x}.zProj), 2, 'omitnan'), locoDeform_csd_opts{x}.binSize );
    tempShearMag = vertcat(deform{x}(expt(x).csd).shearMag);
    tempShearMag = BinDownMean( mean( tempShearMag(:,segParams{x}.zProj), 2, 'omitnan'), locoDeform_csd_opts{x}.binSize );
    tempShearRate = vertcat(deform{x}(expt(x).csd).DshearMag);
    tempShearRate = BinDownMean( mean( tempShearRate(:,segParams{x}.zProj), 2, 'omitnan'), locoDeform_csd_opts{x}.binSize );
    if expt(x).Nplane > 1
        tempShift = vertcat(deform{x}(expt(x).csd).shiftZ);
        tempShiftMean =  BinDownMean( mean( tempShift(:,4:end-3), 2, 'omitnan' ), locoDeform_csd_opts{x}.binSize );
        locoDeform_csd_resp{x}.data = [tempTrans, tempTransSpd, tempScaleMag, tempStretchMag, tempShearMag, tempShearRate, tempShiftMean]; 
        locoDeform_csd_resp{x}.name = {'|Translation|','Trans. Speed','|Scale|','|Stretch|','|Shear|', 'Shear Rate','Z Shift'};
    else
        %tempShiftMean = zeros(size(tempScaleMag,1), 1);
        locoDeform_csd_resp{x}.data = [tempTrans, tempTransSpd, tempScaleMag, tempStretchMag, tempShearMag, tempShearRate]; 
        locoDeform_csd_resp{x}.name = {'|Translation|','Trans. Speed','|Scale|','|Stretch|','|Shear|', 'Shear Rate'};
    end
    locoDeform_csd_resp{x}.data = normalize(locoDeform_csd_resp{x}.data, 1);
    locoDeform_csd_resp{x}.data(abs(locoDeform_csd_resp{x}.data) > 5) = NaN; % suppress extreme outliers
    locoDeform_csd_resp{x}.N = size(locoDeform_csd_resp{x}.data,2); 
    % Remove scans with missing data 
    nanFrame = find(any(isnan([locoDeform_csd_pred{x}.data, locoDeform_csd_resp{x}.data]),2)); % find( isnan(sum(pred(x).data,2)) ); 
    fprintf('\nRemoving %i NaN-containing frames', numel(nanFrame));
    locoDeform_csd_pred{x}.data(nanFrame,:) = []; locoDeform_csd_resp{x}.data(nanFrame,:) = [];

    % Run the GLM
    locoDeform_csd_opts{x}.load = true; %false; %  
    locoDeform_csd_opts{x}.saveRoot = sprintf('%s%s', expt(x).dir, locoDeform_csd_opts{x}.name  ); % ''; %  
    [locoDeform_csd_result{x}, locoDeform_csd_summary{x}, locoDeform_csd_opts{x}, locoDeform_csd_pred{x}, locoDeform_csd_resp{x}] = ...
        GLMparallel(locoDeform_csd_pred{x}, locoDeform_csd_resp{x}, locoDeform_csd_opts{x}); 
end

%% View GLM results
for x = find(~cellfun(@isempty, locoDeform_csd_summary))
    %locoDeform_csd_opts{x}.rShow = 2; %1:locoDeform_csd_resp{x}.N; % 1:LocoDeform_resp{x}.N; %NaN;
    locoDeform_csd_opts{x}.xVar = 'Time';
    ViewGLM(locoDeform_csd_pred{x}, locoDeform_csd_resp{x}, locoDeform_csd_opts{x}, locoDeform_csd_result{x}, locoDeform_csd_summary{x}); %GLMresultFig = 
end

%%
devSummMat = cell(1,Nexpt);
GLMresultFig = figure('WindowState','maximized', 'color','w');
for x = x3Dcsd
    cla;
    devSummMat{x} = [locoDeform_csd_summary{x}.dev; locoDeform_csd_summary{x}.lopo.dev; locoDeform_csd_summary{x}.lofo.dev];
    % {
    imagesc( devSummMat{x}' );
    set(gca,'XtickLabel',  [{'All'}; locoDeform_csd_summary{x3Dcsd(1)}.lopo.name(:); locoDeform_csd_summary{x3Dcsd(1)}.lofo.name(:)], ...
        'Ytick',1:locoDeform_csd_resp{x3Dcsd(1)}.N , 'YtickLabel',locoDeform_csd_resp{x3Dcsd(1)}.name, 'TickDir','out', 'TickLength',[0.003,0]);
    title(sprintf('%s: GLM Fit Summary', locoDeform_csd_opts{x}.name), 'Interpreter','none');
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

locoDeform_csd_GLMdevFig = figure('WindowState','maximized', 'color','w');
imagesc( devSummMed' );
set(gca,'XtickLabel',  [{'All'}; locoDeform_csd_summary{x}.lopo.name(:); locoDeform_csd_summary{x}.lofo.name(:)], ...
    'YtickLabel',locoDeform_csd_resp{x}.name, 'TickDir','out', 'TickLength',[0.003,0]);
title('GLM Fit Summary (All 3D Data)', 'Interpreter','none');
xtickangle(30);
axis image;
CB = colorbar; CB.Label.String = 'Median Deviance Explained';
figPath = sprintf('%s%s_DevSummary.tif', figDir, GLMname );
fprintf('\nSaving %s\n', figPath);
%print( locoDeform_csd_GLMdevFig, figPath, '-dtiff');
impixelinfo;

%% Show relative explanatory value from submodels
exptColor = distinguishable_colors(Npresent);
locoDeform_csd_relExpFig = figure('WindowState','maximized', 'color','w');
Ncol = locoDeform_csd_resp{x3Dcsd(1)}.N;  
%tiledlayout(1, Ncol);
opt = {[0.05,0.04], [0.06,0.02], [0.04,0.02]};
Ndrops = locoDeform_csd_pred{x3Dcsd(1)}.N+locoDeform_csd_pred{x3Dcsd(1)}.fam.N;
for col = 1:Ncol
    tempDevFracMat = nan( Npresent, Ndrops );
    k = 0;
    for x = x3Dcsd
        k = k + 1;
        if locoDeform_csd_result{x}(col).dev > locoDeform_csd_opts{x}.minDev
            tempDevFracMat(k,:) = [locoDeform_csd_result{x}(col).lopo.devFrac, locoDeform_csd_result{x}(col).lofo.devFrac];
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
    title( locoDeform_csd_resp{x3Dcsd(1)}.name{col} );
    set(gca, 'XtickLabel', [locoDeform_csd_pred{x3Dcsd(1)}.name, locoDeform_csd_pred{x3Dcsd(1)}.fam.name])  % [locoDeform_csd_summary{x}.lopo.name(:); locoDeform_csd_summary{x}.lofo.name(:)]
    xtickangle(45)
end
figPath = sprintf('%s%s_ExpValSummary.tif', figDir, GLMname );
fprintf('\nSaving %s\n', figPath);
print( locoDeform_csd_relExpFig, figPath, '-dtiff');

%% Show time filters per response
exptColor = distinguishable_colors(Npresent);
locoDeform_csd_GLMcoeffFig = figure('WindowState','maximized', 'color','w');
Ncol = locoDeform_csd_resp{x3Dcsd(1)}.N;   Nrow = locoDeform_csd_pred{x3Dcsd(1)}.N;
tiledlayout(Nrow, Ncol);
k = 0;
for row = 1:Nrow
    for col = 1:Ncol
        nexttile
        tempCoeffMat = zeros( numel(locoDeform_csd_opts{x3Dcsd(1)}.lags), Npresent );
        k = 0;
        for x = x3Dcsd
            k = k + 1;
            if locoDeform_csd_result{x}(col).dev > locoDeform_csd_opts{x}.minDev
                tempCoeffMat(:,k) = locoDeform_csd_result{x}(col).coeff(:,row);
            end
        end
        plot( locoDeform_csd_opts{x}.lags, tempCoeffMat ); hold on;
        colororder(exptColor)
        axis square;
        if col == 1, ylabel(sprintf('%s Coefficient', locoDeform_csd_pred{x3Dcsd(1)}.name{row} )); end
        if row == 1, title( locoDeform_csd_resp{x3Dcsd(1)}.name{col} ); end
        if row == Nrow, xlabel('Response Lag (s)'); end
        %pause;
    end
end

figPath = sprintf('%s%s_CoeffSummary.tif', figDir, GLMname );
fprintf('\nSaving %s\n', figPath);
print( locoDeform_csd_GLMcoeffFig, figPath, '-dtiff');
