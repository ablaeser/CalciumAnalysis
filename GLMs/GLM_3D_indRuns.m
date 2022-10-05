%% Use GLM to assess contribution of different variables
runGLM_pred = cell(1,Nexpt); runGLM_resp = cell(1,Nexpt); runGLM_opts = cell(1,Nexpt); runGLM_result = cell(1,Nexpt); runGLM_summary = cell(1,Nexpt);
figDir = 'D:\MATLAB\LevyLab\Figures\3D\GLM\'; mkdir( figDir )
GLMname = 'runGLM_TScShDtStDshZVAL_80pct';
GLMrate = 15.49/30;
for x = x3D %xPresent %
    for run = flip(1:expt(x).Nruns)
        % GLMparallel options
        runGLM_opts{x}(run).name = sprintf('%s_%s_run%i', expt(x).name, GLMname, run); %strcat(expt(x).name, , '_preCSDglm');
        runGLM_opts{x}(run).rShow = NaN;
        runGLM_opts{x}(run).figDir = ''; % figDir;
        runGLM_opts{x}(run).alpha = 0.01;  % The regularization parameter, default is 0.01
        runGLM_opts{x}(run).standardize = true; 
        runGLM_opts{x}(run).trainFrac = 0.75; % 1; %
        runGLM_opts{x}(run).Ncycle = 20;
        runGLM_opts{x}(run).distribution = 'gaussian'; % 'poisson'; %  
        runGLM_opts{x}(run).CVfold = 10;
        runGLM_opts{x}(run).nlamda = 1000;
        runGLM_opts{x}(run).maxit = 5*10^5;
        runGLM_opts{x}(run).minDev = 0.05; 
        runGLM_opts{x}(run).minDevFrac = 0.1;
        runGLM_opts{x}(run).maxP = 0.05;
        runGLM_opts{x}(run).Nshuff = 0;  
        runGLM_opts{x}(run).minShuff = 15; 
        runGLM_opts{x}(run).window = [-4,4]; % [0,0]; % [-0.5, 0.5]; % 
        runGLM_opts{x}(run).lopo = true; %false; %
        runGLM_opts{x}(run).frameRate = GLMrate; %expt(x).scanRate; 
        runGLM_opts{x}(run).binSize = expt(x).scanRate/GLMrate;
        runGLM_opts{x}(run).minShuffFrame = round( runGLM_opts{x}(run).frameRate*runGLM_opts{x}(run).minShuff );
        windowFrame = round(runGLM_opts{x}(run).window*runGLM_opts{x}(run).frameRate);
        runGLM_opts{x}(run).shiftFrame = windowFrame(1):windowFrame(2);
        runGLM_opts{x}(run).maxShift = max( abs(windowFrame) );
        runGLM_opts{x}(run).Nshift = numel( runGLM_opts{x}(run).shiftFrame );  %Nshift = preCSDOpts(x).Nshift;
        runGLM_opts{x}(run).lags = runGLM_opts{x}(run).shiftFrame/runGLM_opts{x}(run).frameRate;
        runGLM_opts{x}(run).xVar = 'Time';

        % Concatenate input variables pre-CSD
        % Define predictors
        transMag = deform{x}(run).transMag;
        transPred = prctile(transMag(:,segParams{x}.zProj),80,2);
        transPred(normalize(transPred) > 5) = NaN;

        scaleMag = deform{x}(run).scaleMag;
        scalePred = prctile(scaleMag(:,segParams{x}.zProj),80,2);
        scalePred(normalize(scalePred) > 5) = NaN;
        %{
        sp(1) = subplot(2,1,1); imagesc( catScaleMag' ); set(gca,'Ytick',segParams{x}(runs).zProj, 'TickDir','out');
        title(sprintf('%s', expt(x).name ), 'Interpreter','none');
        impixelinfo;
        sp(2) = subplot(2,1,2);
        plot( scalePred ); hold on; axis tight;
        %plot( max(catScaleMag(:,segParams{x}(runs).zProj),[],2) );  
        %plot( prctile(catScaleMag(:,segParams{x}(runs).zProj),80,2) ); hold on; axis tight;
        %plot( prctile(catScaleMag(:,segParams{x}(runs).zProj),90,2) ); hold on; axis tight;
        %plot( mean(catScaleMag(:,segParams{x}(runs).zProj),2, 'omitnan') );
        %plot( median(catScaleMag(:,segParams{x}(runs).zProj),2, 'omitnan'))
        linkaxes(sp,'x');
        pause; clf;
        %}
        shearMag = deform{x}(run).shearMag;
        shearPred = prctile(shearMag(:,segParams{x}.zProj),80,2);
        shearPred(normalize(shearPred) > 5) = NaN;

        transSpdMag = deform{x}(run).DtransMag;
        transSpdPred = prctile(transSpdMag(:,segParams{x}.zProj),80,2);
        transSpdPred(normalize(transSpdPred) > 5) = NaN;

        stretchMag = deform{x}(run).stretchMag;
        stretchPred = prctile(stretchMag(:,segParams{x}.zProj),80,2);
        stretchPred(normalize(stretchPred) > 5) = NaN;

        shearRateMag = deform{x}(run).DshearMag;
        shearRatePred = prctile(shearRateMag(:,segParams{x}.zProj),80,2);
        shearRatePred(normalize(shearRatePred) > 5) = NaN;

        tempShift = deform{x}(run).shiftZ;
        shiftPred = prctile(tempShift(:,segParams{x}.zProj),80,2);
        shiftPred(normalize(shiftPred) > 5) = NaN;

        %{
        sp(1) = subplot(2,1,1); imagesc( tempShift' ); set(gca,'Ytick',intersect(segParams{x}(runs).zProj, 4:expt(x).Nplane-3), 'TickDir','out');
        title(sprintf('%s', expt(x).name ), 'Interpreter','none');
        impixelinfo;
        sp(2) = subplot(2,1,2);
        plot( shiftPred ); hold on; axis tight;

        linkaxes(sp,'x');
        pause; clf;
        %}

        runGLM_pred{x}(run) = struct('data',[], 'name',[], 'N',NaN, 'TB',[], 'lopo',[], 'fam',[]); 
        runGLM_pred{x}(run).data = [transPred, scalePred, shearPred, transSpdPred, stretchPred, shearRatePred, shiftPred, loco{x}(run).Vdown, abs(loco{x}(run).Adown), loco{x}(run).stateDown]; 
        runGLM_pred{x}(run).data = BinDownMean( runGLM_pred{x}(run).data, runGLM_opts{x}(run).binSize );
        runGLM_pred{x}(run).name = {'|Translation|', '|Scale|', '|Shear|', 'TransSpd', 'Stretch', 'Shear Rate','Z Shift', 'Velocity', '|Accel|', 'Locomotive State'}; 
        runGLM_pred{x}(run).N = size(runGLM_pred{x}(run).data,2);
        for p = flip(1:runGLM_pred{x}(run).N), runGLM_pred{x}(run).lopo.name{p} = ['No ',runGLM_pred{x}(run).name{p}]; end
        % Set up leave-one-family-out
        runGLM_pred{x}(run).fam.col = {1:7,8:10}; %{1:4, 5:7}; %{1:2, 3:4, 5:6, 7:8, 9:10, 11:12};  % {1:12};%{1, 2:3, 4:5, 6:7, 8, 9}; 
        runGLM_pred{x}(run).fam.N = numel(runGLM_pred{x}(run).fam.col); 
        runGLM_pred{x}(run).fam.name = {'Deformation', 'Locomotion'}; %{'All'};%  'Onset Time',

        % Define response
        runGLM_resp{x}(run).data = BinDownMean( fluor{x}(run).z.ROI, runGLM_opts{x}(run).binSize ); 
        runGLM_resp{x}(run).N = size(runGLM_resp{x}(run).data, 2); 
        runGLM_resp{x}(run).name = sprintfc('Fluor %i', 1:runGLM_resp{x}(run).N);

        % Remove scans with any missing data 
        nanFrame = find(any(isnan([runGLM_pred{x}(run).data, runGLM_resp{x}(run).data]),2)); % find( isnan(sum(pred(x).data,2)) ); 
        fprintf('\nRemoving %i NaN-containing frames', numel(nanFrame));
        runGLM_pred{x}(run).data(nanFrame,:) = []; runGLM_resp{x}(run).data(nanFrame,:) = [];

        % Run the GLM
        runGLM_opts{x}(run).load = true; % false; % 
        runGLM_opts{x}(run).saveRoot = expt(x).dir; %''; %
        try
            [runGLM_result{x}{run}, runGLM_summary{x}{run}, ~, runGLM_pred{x}(run), runGLM_resp{x}(run)] = GLMparallel(runGLM_pred{x}(run), runGLM_resp{x}(run), runGLM_opts{x}(run));
            runGLM_summary{x}{run} = SummarizeGLM(runGLM_result{x}{run}, runGLM_pred{x}(run), runGLM_resp{x}(run), runGLM_opts{x}(run));
        catch
            fprintf('\n[x,run] = [%i, %i] GLM failed\n', x, run) 
        end
    end
end

%% Does the fraction of mechanosensitive units increase post-CSD?
k = 0;
mechFrac = nan(Npresent, 6);
for x = xCSD % xPresent
    k = k+1;
    for run = 1:expt(x).Nruns
        try
            mechFrac(k,run) = (runGLM_summary{x}{run}.nMixed + runGLM_summary{x}{run}.nDeform)/expt(x).Nroi;
        catch
            fprintf('\n[x,run] = [%i, %i] GLM failed\n', x, run)
        end
    end
end
mechFrac(:,3) = []; % exclude the CSD run itself
JitterPlot( mechFrac, 'dependent',true );
ylabel('% ROI mechanosensitive');
set(gca,'Xtick',1:5, 'XtickLabel',{'Pre1','Pre2','Post1','Post2','Post3'});

%% Does the AUC of scaling coefficients increase post-CSD?
close all;
scaleAUC = cell(1,Nexpt);
k = 0;
for x = 40 %x3Dcsd % %xPresent %xCSD % 31 %40
    scaleAUC{x} = nan(expt(x).Nroi, expt(x).Nruns);
    for run = 1:expt(x).Nruns
        posLags = find(runGLM_opts{x}(run).lags >= 0);
        scaleCol = find(contains(runGLM_pred{x}(run).name, 'Scale'));
        try
            for roi = [runGLM_summary{x}{run}.rMixed, runGLM_summary{x}{run}.rDeform]
                scaleAUC{x}(roi,run) = sum(runGLM_result{x}{run}(roi).coeff(posLags,scaleCol), 1);
            end
        end
    end
    % Plot the data
    % {
    JitterPlot( scaleAUC{x}, 'dependent',true );
    title( sprintf('x = %i: %s', x, expt(x).name), 'Interpreter','none' ); 
    set(gca, 'Xtick',1:expt(x).Nruns); xlabel('Run'); ylabel('Scaling coefficient AUC');
    %}
    % repeated measures 2-way ANOVA
    % {
    clearvars dataTable WithinStructure rm ranovaTable runComp csdComp
    dataTable = array2table(scaleAUC{x}, 'VariableNames',sprintfc('Run%i',1:expt(x).Nruns) );
    %WithinStructure = table((1:expt(x).Nruns)', [ones(1, expt(x).csd-1), 2*ones(1, 1+expt(x).Nruns-expt(x).csd)]','VariableNames',{'Run','CSD'});
    WithinStructure = table((1:expt(x).Nruns)', [1,1,2,2]','VariableNames',{'Run','CSD'});
    WithinStructure.CSD = categorical(WithinStructure.CSD);
    wilxNote = sprintf('Run1-Run%i~1', expt(x).Nruns);
    try
        rm = fitrm(dataTable, wilxNote, 'WithinDesign',WithinStructure ); % 'Run1-Run6~1'
        runComp = multcompare( rm, 'Run');
        csdComp = multcompare( rm, 'CSD');
        ranovaTable = ranova(rm,'WithinModel','Run*CSD'); % 'Run+CSD+Run*CSD' gives the same result
        title( sprintf('x = %i: %s. pCSD = %2.4f', x, expt(x).name, csdComp.pValue(1)), 'Interpreter','none' ); 
    catch
        fprintf('\nx = %i repeated measures model failed', x);
    end
    pause; clf;
    %}
end

%% Pool AUC across experiments and run 2-way anova
scaleAUCpool = [];
for x = x3Dcsd
    % pad array to standard format: 2 pre-CSD columns, 4 post CSD
    tempAUC = [nan(size(scaleAUC{x},1),3-expt(x).csd), scaleAUC{x}]; 
    tempAUC = [tempAUC, nan( size(tempAUC,1), 6-size(tempAUC,2) )];
    scaleAUCpool = [scaleAUCpool; tempAUC];
end
scaleAUCpool(sum(isnan(scaleAUCpool), 2)>=5,:) = [];
JitterPlot( scaleAUCpool, 'dependent',true );
title( sprintf('x = %i: %s', x, expt(x).name), 'Interpreter','none' );
set(gca, 'Xtick',1:expt(x).Nruns); xlabel('Run'); ylabel('Scaling coefficient AUC');
% Run repeated measures 2-way ANOVA on the pooled data
clearvars dataTable WithinStructure rm ranovaTable runComp csdComp
dataTable = array2table(scaleAUCpool, 'VariableNames',sprintfc('Run%i',1:expt(x).Nruns) );
WithinStructure = table((1:expt(x).Nruns)', [1,1,2,2,2,2]','VariableNames',{'Run','CSD'}); % [ones(1, expt(x).csd-1), 2*ones(1, 1+expt(x).Nruns-expt(x).csd)]'
WithinStructure.CSD = categorical(WithinStructure.CSD);
rm = fitrm(dataTable, 'Run1-Run6~1', 'WithinDesign',WithinStructure ); %
runComp = multcompare( rm, 'Run');
csdComp = multcompare( rm, 'CSD');
ranovaTable = ranova(rm,'WithinModel','Run*CSD'); % 'Run+CSD+Run*CSD' gives the same result
title( sprintf('x = %i: %s. pCSD = %2.4f', x, expt(x).name, csdComp.pValue(1)), 'Interpreter','none' );

%% use ANOVA to checx for run and pre-post CSD effects
%{
% two-way anova without repeated measures
sampleData = scaleAUC{x}(:);
runGrouping = repmat(1:size(scaleAUC{x},2), size(scaleAUC{x},1), 1 ); runGrouping = runGrouping(:);
preGrouping = repmat(["Pre","Pre","Post","Post","Post","Post"], size(scaleAUC{x},1), 1 ); preGrouping = preGrouping(:);
%anovan(sampleData, {runGrouping, preGrouping}, 'continuous',1, 'model','interaction', 'varnames',["Run","CSD"] )
%}
% repeated measures 2-way ANOVA
dataTable = array2table(scaleAUC{x}, 'VariableNames',sprintfc('Run%i',1:6) );
WithinStructure = table([1:size(scaleAUC{x},2)]', [1,1,2,2,2,2]','VariableNames',{'Run','CSD'});
WithinStructure.CSD = categorical(WithinStructure.CSD);
rm = fitrm(dataTable, 'Run1-Run6~1', 'WithinDesign',WithinStructure );
ranovaTable = ranova(rm,'WithinModel','Run*CSD'); % 'Run+CSD+Run*CSD' gives the same result
multcompare( rm, 'CSD')
%% two-factor repeated measures anova (https://www.mathworxs.com/matlabcentral/answers/467222-ranova-with-two-within-factors?s_tid=srchtitle)
datatable = array2table(rand([20,4])); %cell2table(struct2cell([pre.A, pre.B, post.A, post.B]'));
datatable.Properties.VariableNames = {'pre_A','pre_B','post_A','post_B'};
% When you have more than one repeated-measures factor, you must set up a table
% to indicate the levels on each factor for each of your different variables.
% Here is the command you need for this case:
WithinStructure = table([1 1 2 2]',[1 2 1 2]','VariableNames',{'PrePost','TreatAB'});
% The 4 different rows of the WithinStructure table correspond to the 4 different
% columns, 'pre_A','pre_B','post_A','post_B', respectively, in your data table.
% Each 'pre_A','pre_B','post_A','post_B' column is coded as 1/2 on the PrePost factor
% and as 1/2 on the TreatAB factor.
% Indicate that the 1's and 2's of WithinStructure are category labels
% rather than regression-type numerical covariates:
WithinStructure.PrePost = categorical(WithinStructure.PrePost);
WithinStructure.TreatAB = categorical(WithinStructure.TreatAB);
% Now pass the WithinStructure table to fitrm so that it xnows how the different
% columns correspond to the different levels of the repeated-measures factors.
rm = fitrm(datatable, 'pre_A,pre_B,post_A,post_B~1','WithinDesign',WithinStructure);
% Finally, you need to specify the repeated-measures factors again when you call ranova, lixe this:
ranovatable = ranova(rm,'WithinModel','PrePost*TreatAB');

%% repeated measures anova example (https://www.mathworxs.com/help/stats/repeatedmeasuresmodel.ranova.html)
load('longitudinalData.mat');
Gender = ['F' 'F' 'F' 'F' 'F' 'F' 'F' 'F' 'M' 'M' 'M' 'M' 'M' 'M' 'M' 'M']';
t = table(Gender,Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5), 'VariableNames',{'Gender','t0','t2','t4','t6','t8'});
Time = [0 2 4 6 8]';
rm = fitrm(t,'t0-t8 ~ Gender','WithinDesign',Time);
ranovatbl = ranova(rm)

%% https://www.mathworxs.com/help/stats/repeatedmeasuresmodel.ranova.html#bt42zie-A
load repeatedmeas
rm = fitrm(between,'y1-y8 ~ Group*Gender + Age + IQ','WithinDesign',within);
ranovatbl = ranova(rm); % Perform repeated measures analysis of variance.
[ranovatbl,A,C,D] = ranova(rm,'WithinModel','w1+w2'); %Specify the model for the within-subject factors. Also display the matrices used in the hypothesis test.
T = multcompare(rm,'Group'); %'w1'