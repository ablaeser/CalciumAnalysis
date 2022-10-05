%% Use GLM to assess contribution of different variables
LKF_pred = cell(1,Nexpt); LKF_resp = cell(1,Nexpt); LKF_opts = cell(1,Nexpt); LKF_result = cell(1,Nexpt); LKF_summary = cell(1,Nexpt);
%figDir = 'D:\MATLAB\LevyLab\Figures\3D\GLM\'; mkdir( figDir )
GLMname = 'LocoKinFluor_ctrl_accelMag_60s'; % 'LocoKinFluor_2D_VelAcc'; %  'LocoKinFluor_2D' 
GLMrate = 15.49/30;
for x = xCtrl %xPresent %x3Dcsd
    % GLMparallel options
    LKF_opts{x}.name = sprintf('%s_%s', expt{x}.name, GLMname); %strcat(expt{x}.name, , '_LKFglm');
    LKF_opts{x}.rShow = NaN;
    LKF_opts{x}.figDir = ''; % figDir;
    LKF_opts{x}.alpha = 0.01;  % The regularization parameter, default is 0.01
    LKF_opts{x}.standardize = true; 
    LKF_opts{x}.trainFrac = 0.75; % 1; %
    LKF_opts{x}.Ncycle = 20;
    LKF_opts{x}.distribution = 'gaussian'; % 'poisson'; %  
    LKF_opts{x}.CVfold = 10;
    LocoDeform.opts(x).nlamda = 1000;
    LocoDeform.opts(x).maxit = 5*10^5;
    LKF_opts{x}.minDev = 0.05;  minDev = LKF_opts{x}.minDev;
    LKF_opts{x}.minDevFrac = 0.1;
    LKF_opts{x}.maxP = 0.05;
    LKF_opts{x}.Nshuff = 0;  Nshuff = LKF_opts{x}.Nshuff;
    LKF_opts{x}.minShuff = 15; 
    LKF_opts{x}.window = [-60,60]; % [0,0]; % [-0.5, 0.5]; % 
    LKF_opts{x}.lopo = true; %false; %
    LKF_opts{x}.frameRate = GLMrate; %expt{x}.scanRate; 
    LKF_opts{x}.binSize = round(expt{x}.scanRate/GLMrate);
    LKF_opts{x}.minShuffFrame = round( LKF_opts{x}.frameRate*LKF_opts{x}.minShuff );
    windowFrame = [ceil(LKF_opts{x}.window(1)*LKF_opts{x}.frameRate), floor(LKF_opts{x}.window(2)*LKF_opts{x}.frameRate)]; % floor(LKF_opts{x}.window*LKF_opts{x}.frameRate);
    LKF_opts{x}.shiftFrame = windowFrame(1):windowFrame(2);
    LKF_opts{x}.maxShift = max( abs(windowFrame) );
    LKF_opts{x}.Nshift = numel( LKF_opts{x}.shiftFrame );  %Nshift = LKFOpts(x).Nshift;
    LKF_opts{x}.lags = LKF_opts{x}.shiftFrame/LKF_opts{x}.frameRate;

    % Concatenate and bin input variables 
    if ~isnan(expt{x}.csd),  preCSDruns = 1:expt{x}.csd-1;  else,  preCSDruns = 1:expt{x}.Nruns;  end

    % Define predictors
    tempVelocityCat = BinDownMean( vertcat(loco{x}(preCSDruns).Vdown), LKF_opts{x}.binSize );
    tempAccelCat = BinDownMean( abs(vertcat(loco{x}(preCSDruns).Adown)), LKF_opts{x}.binSize ); %[NaN; diff(tempVelocityCat)];
    %tempSpeedCat = BinDownMean( vertcat(loco{x}(preCSDruns).speedDown), LKF_opts{x}.binSize );
    tempStateCat = BinDownMean( vertcat(loco{x}(preCSDruns).stateDown), LKF_opts{x}.binSize );
    LKF_pred{x} = struct('data',[], 'name',[], 'N',NaN, 'TB',[], 'lopo',[], 'fam',[]); 
    LKF_pred{x}.data = [tempVelocityCat, tempAccelCat, tempStateCat]; % ; % tempSpeedCat,, , tempStateCat 
    LKF_pred{x}.name = {'Velocity','|Accel|','State'}; % ,'Speed'  {}; %,  ,'State'
    LKF_pred{x}.N = size(LKF_pred{x}.data,2);
    for p = flip(1:LKF_pred{x}.N), LKF_pred{x}.lopo.name{p} = ['No ',LKF_pred{x}.name{p}]; end
    % Set up leave-one-family-out
    LKF_pred{x}.fam.col = {}; % {1:2, 3}; % 1:2, 3 %{1:2, 3:4, 5:6, 7:8, 9:10, 11:12};  % {1:12};%{1, 2:3, 4:5, 6:7, 8, 9}; 
    LKF_pred{x}.fam.N = numel(LKF_pred{x}.fam.col); 
    LKF_pred{x}.fam.name = {}; %{'Kinematics', 'State'}; % {}; %{'All'};%  'Onset Time',

    % Define responses
    tempFluor = [fluor{x}(preCSDruns).z]; % [fluor{x}(preCSDruns).act]; % 
    tempFluorCat = BinDownMean( vertcat(tempFluor.ROI), LKF_opts{x}.binSize ); 
    %tempFluorCat(tempFluorCat < 0) = 0;
    LKF_resp{x}.data = tempFluorCat; % tempScaleMag; %
    LKF_resp{x}.N = size(LKF_resp{x}.data, 2); 
    LKF_resp{x}.name = sprintfc('Fluor %i', 1:LKF_resp{x}.N);
    
    % Remove scans with missing data 
    nanFrame = find(any(isnan([LKF_pred{x}.data, LKF_resp{x}.data]),2)); % find( isnan(sum(pred(x).data,2)) ); 
    fprintf('\nRemoving %i NaN-containing frames', numel(nanFrame));
    LKF_pred{x}.data(nanFrame,:) = []; LKF_resp{x}.data(nanFrame,:) = [];

    % Run the GLM
    LKF_opts{x}.load = true; %  false; %    
    LKF_opts{x}.saveRoot = expt{x}.dir; %; sprintf('%s%s', expt{x}.dir, LKF_opts{x}.name  );  
    [LKF_result{x}, LKF_summary{x}, LKF_opts{x}, LKF_pred{x}, LKF_resp{x}] = GLMparallel(LKF_pred{x}, LKF_resp{x}, LKF_opts{x}); 

end

%% View GLM results
for x = xCtrl %xPresent
    LKF_opts{x}.rShow = boutResponse(x).pre.inh; % 1:expt{x}.Nroi; %LKF_summary{x}.rGood; %1:LocoDeform_resp{x}.N; %NaN; % 1:LocoDeform_resp{x}.N; %NaN;
    LKF_opts{x}.xVar = 'Time';
    ViewGLM(LKF_pred{x}, LKF_resp{x}, LKF_opts{x}, LKF_result{x}, LKF_summary{x}); %GLMresultFig = 
end

%% Pool results across experiments
LKFdevPool = []; goodDevPool = [];
LKFdevFracPool = []; %lofoDevFracPool = [];
LKF_Nsubtype = []; k = 0;
for x = xPresent
    k = k+1;
    if ~isempty( LKF_summary{x}.rGood )
        LKFdevPool = [LKFdevPool, LKF_summary{x}.dev]; 
        goodDevPool = [goodDevPool, LKF_summary{x}.dev( LKF_summary{x}.rGood )];
        LKFdevFracPool = [LKFdevFracPool, vertcat(LKF_summary{x}.lopo.devFrac(:,LKF_summary{x}.rGood) , LKF_summary{x}.lofo.devFrac(:,LKF_summary{x}.rGood) )]; %
    end
end

% Summarize deviance explained
opt = {[0.02,0.07], [0.1,0.07], [0.09,0.09]};  % {[vert, horz], [bottom, top], [left, right] }
DevianceFig = figure('WindowState','maximized', 'color','w');
k = 1; clearvars h;
subtightplot(1,2,1,opt{:});
for x = xPresent
    [Ftemp, Xtemp] = ecdf( LKF_summary{x}.dev ); hold on;
    h(k) = plot(Xtemp, Ftemp, 'color',0.7*[1,1,1] ); % exptColor(k,:)
    k = k + 1;
end
[Fdev, Xdev] = ecdf( LKFdevPool );
h(k) = plot( Xdev, Fdev, 'color','k', 'LineWidth',2 ); 
axis square;
legend(h, {expt(xPresent).name, 'Pooled', 'Threshold'}, 'Location','southEast', 'Interpreter','none', 'AutoUpdate',false );
xlim([0, 0.2]);
line(LKF_opts{xPresent(1)}.minDev*[1,1], [0,1], 'Color','r', 'LineStyle','--'); % h(k+1) = 
xlabel('Deviance Explained'); ylabel('Fraction of Units');
title( sprintf('%s Fit Results', GLMname), 'Interpreter','none' );

subtightplot(1,2,2,opt{:});
JitterPlot( 1 - LKFdevFracPool', 0.45, 'ErrorCap',10, 'monochrome',0.6); hold on;
line([0,size(LKFdevFracPool,1)+1], [0,0], 'color','k');
axis square;
ylabel('Relative Explanatory Value'); %ylabel('Cumulative Distribution');
ylim([-1,1]);
set(gca,'Xtick', 1:size(LKFdevFracPool,1), 'XtickLabel',[LKF_pred{xPresent(1)}.name, LKF_pred{xPresent(1)}.fam.name]);
xtickangle(30);
title('Well-Fit Units');




%% Show time filters per response
exptColor = distinguishable_colors(Npresent);
LKF_GLMcoeffFig = figure('WindowState','maximized', 'color','w');
Ncol = LKF_resp{xPresent(1)}.N;   Nrow = LKF_pred{xPresent(1)}.N;
tiledlayout(Nrow, Ncol);
k = 0;
for row = 1:Nrow
    for col = 1:Ncol
        nexttile
        try
            tempCoeffMat = zeros( numel(LKF_opts{xPresent(1)}.lags), Npresent );
            k = 0;
            for x = xPresent
                k = k + 1;
                if LKF_result{x}(col).dev > LKF_opts{x}.minDev
                    tempCoeffMat(:,k) = LKF_result{x}(col).coeff(:,row);
                end
            end
            plot( LKF_opts{x}.lags, tempCoeffMat ); hold on;
            colororder(exptColor)
            axis square;
            if col == 1, ylabel(sprintf('%s Coefficient', LKF_pred{xPresent(1)}.name{row} )); end
            if row == 1, title( LKF_resp{xPresent(1)}.name{col} ); end
            if row == Nrow, xlabel('Response Lag (s)'); end
        end
        %pause;
    end
end

figPath = sprintf('%s%s_CoeffSummary.tif', figDir, GLMname );
fprintf('\nSaving %s\n', figPath);
%print( LKF_GLMcoeffFig, figPath, '-dtiff');

%% Break out coefficients by bout response subtype

for x = xPresent %x2D
    excCoeffCat = cat(3, LKF_result{x}( boutResponse(x).pre.exc ).coeff);
    neutCoeffCat = cat(3, LKF_result{x}( boutResponse(x).pre.neut ).coeff);
    inhCoeffCat = cat(3, LKF_result{x}( boutResponse(x).pre.inh ).coeff);
    
    LKF_lags = linspace(LKF_opts{x}.window(1), LKF_opts{x}.window(2), LKF_opts{x}.Nshift);
    
    subplot(1,3,1);
    plot( mean(excCoeffCat(:,1,:), 3), 'b' ); hold on;
    plot( mean(neutCoeffCat(:,1,:), 3), 'k' ); hold on;
    plot( mean(inhCoeffCat(:,1,:), 3), 'r' ); hold on;
    axis square;
    title('dF/F per velocity')
    
    subplot(1,3,2);
    plot( mean(excCoeffCat(:,2,:), 3), 'b' ); hold on;
    plot( mean(neutCoeffCat(:,2,:), 3), 'k' ); hold on;
    plot( mean(inhCoeffCat(:,2,:), 3), 'r' ); hold on;
    axis square;
    title('dF/F per accel')
    
    subplot(1,3,3);
    plot( mean(excCoeffCat(:,3,:), 3), 'b' ); hold on;
    plot( mean(neutCoeffCat(:,3,:), 3), 'k' ); hold on;
    plot( mean(inhCoeffCat(:,3,:), 3), 'r' ); hold on;
    axis square;
    title('dF/F per loco state')
    pause;
end

%% Breakout coefficients by bout response subtype - heatmaps
Nshift = LKF_opts{x2D(1)}.Nshift;
typeCoeff = {nan(Nshift, 3, 0), nan(Nshift, 3, 0), nan(Nshift, 3, 0)};
typeVal = {nan(0,3), nan(0,3), nan(0,3)};
typeDev = cell(1,3);
for x = find(~cellfun(@isempty, LKF_summary))
    goodExc = intersect(LKF_summary{x}.rGood, boutResponse(x).pre.exc);
    goodNeut = intersect(LKF_summary{x}.rGood, boutResponse(x).pre.neut);
    goodInh = intersect(LKF_summary{x}.rGood, boutResponse(x).pre.inh);

    typeDev{1} = cat(2, typeDev{1}, [LKF_result{x}(goodExc).dev]);
    typeDev{2} = cat(2, typeDev{2}, [LKF_result{x}(goodNeut).dev]);
    typeDev{3} = cat(2, typeDev{3}, [LKF_result{x}(goodInh).dev]);
    
    typeCoeff{1} = cat(3, typeCoeff{1}, LKF_result{x}(goodExc).coeff);
    typeCoeff{2} = cat(3, typeCoeff{2}, LKF_result{x}(goodNeut).coeff);
    typeCoeff{3} = cat(3, typeCoeff{3}, LKF_result{x}(goodInh).coeff);
    
    typeVal{1} = vertcat(typeVal{1}, 1-LKF_summary{x}.lopo.devFrac(:,goodExc)');
    typeVal{2} = vertcat(typeVal{2}, 1-LKF_summary{x}.lopo.devFrac(:,goodNeut)');
    typeVal{3} = vertcat(typeVal{3}, 1-LKF_summary{x}.lopo.devFrac(:,goodInh)');
end
typeName = ["Bout-Activated","Bout-Neutral","Bout-Suppressed"];


Tshift = linspace(-60,60, LKF_opts{x}.Nshift);
Ttick = linspace(-60,60, 5);
coeffSubtypeFig_shaded = figure('WindowState','maximized', 'color','w');
opt = {[0.05,0.05], [0.1,0.07], [0.09,0.09]};  % {[vert, horz], [bottom, top], [left, right] }
for p = 1:3 % bout response type
    for t = 1:3 %  predictor
        subtightplot(3,3, sub2ind([3,3], t,p), opt{:});
        tempCoeff = squeeze(typeCoeff{t}(:,p,:));
        meanCoeff = mean(tempCoeff, 2);
        semCoeff = std(tempCoeff,0,2); % SEM(tempCoeff, 2);
        coeffShadeMat = [meanCoeff-semCoeff, meanCoeff, meanCoeff+semCoeff]';
        line([-60, 60], [0,0], 'color','k');  hold on;
        plotshaded(Tshift, coeffShadeMat, '-');
        axis tight; axis square; 
        set(gca,'Xtick',Ttick);
        
        if t == 1, ylabel([LKF_pred{x}.name{p},' coeff']); end
        if p == 3; xlabel('Delay (s)'); end
        if p == 1, title(typeName(t)); end
        
        %text(0.9, 0.9, sprintf('Val = %2.2f +/ %2.2f', mean(typeVal{t}(:,p)), SEM(typeVal{t}(:,p))  ), 'Units','normalized' )
        %pause;
    end
end
figPath = sprintf('%scoeffSubtypeFig_shaded.pdf', figDir);
exportgraphics(coeffSubtypeFig_shaded, figPath, 'Resolution',300); fprintf('\nSaved %s', figPath);


devSubtypeFig = figure('WindowState','maximized', 'color','w');
subplot(1,4,1);
JitterPlot(cell2padmat(typeDev), 'ErrorBar',3, 'ErrorCap',20, 'new',false);
axis square;
ylabel('Deviance explained');
set(gca,'Xtick',1:3, 'XtickLabel',typeName);
xtickangle(30);

for t = 1:3
    subplot(1,4,t+1);
    bar(mean(typeVal{t})); hold on;
    errorbar(1:3, mean(typeVal{t}), SEM(typeVal{t}), 'linestyle','none', 'color','k');
    title(typeName{t})
    if t == 1, ylabel('Rel. exp. value'); end
    %JitterPlot(typeVal{t}, 0.2, 'paired',true, 'new',false);
    %axis tight; 
    xlim([0.5,3.5]);
    axis square;
    set(gca,'Xtick',[1,2,3], 'XtickLabel',{'Velocity','|Accel|','State'});
end
figPath = sprintf('%sdevSubtypeFig.pdf', figDir);
exportgraphics(devSubtypeFig, figPath, 'Resolution',300); fprintf('\nSaved %s', figPath);

valMean = cellfun(@mean, typeVal, 'UniformOutput',false);
valMean = vertcat(valMean{[1,3]}); % exclude bout-neutral ROI
valSEM = cellfun(@SEM, typeVal, 'UniformOutput',false);
valSEM = vertcat(valSEM{[1,3]}); % % exclude bout-neutral ROI
GroupedBarError(valMean', valSEM', 'ErrorBar',1.5, 'ErrorCap',3)


%{
coeffSubtypeFig_heatmaps = figure('WindowState','maximized', 'color','w');
for t = 1:3 % bout response type
    for p = 1:3 %  predictor
        subtightplot(3,3, sub2ind([3,3], t,p));
        tempCoeff = squeeze(typeCoeff{t}(:,p,:))';
        imagesc(tempCoeff); 
        colormap(bluewhitered);
        caxis(max(tempCoeff(:))*[-1,1]);
        colormap(bluewhitered)
        colorbar();
        pause;
    end
end
%}

