%% Revised figure

% A - still frames of wave progress, with time shown, some fiber ROIs shown, and scalebar
X = 37; % DL75
% Generate movies of fibers during registered CSD bouts
movParam.dir = 'D:\MATLAB\LevyLab\Figures\BoutMovies\CSD\Fibers\'; % 'D:\MATLAB\LevyLab\Figures\CSD\Movies\';
mkdir(movParam.dir)
movParam.fmtSpec = '%2.1f';
movParam.Tperi = [0,0]; %[15,10]; % time before/after bout to show
movParam.zProj = [];
movParam.binT = 1;
movParam.displayPct = [1,99.9];
movParam.sbx = [];
movParam.boutType = 'csd';
movParam.edges = []; % [80,80,20,20]; %[80,90,60,110]; %
movParam.aviRate = 10; % frames per second
movParam.level = 'ind';
movParam.edges = segParams{X}.edges;
if expt(X).Nplane > 1
    movParam.sourceSbx = 'sbx_interp';  
else
    movParam.sourceSbx = 'sbx_affine'; %'sbxz'; %
end

movParam.zProj = unique(round(fiber{x}(f).cent(:,3)))';
movParam.regType = sprintf('_aff_ExampleFibers'); %'
movParam.scalebar = MakeScaleBar( round(expt(X).umPerPixel*[50,50]), {[0,expt(X).Ncol-segParams{X}.edges(1)-segParams{X}.edges(2)]+0.5, [0,expt(X).Nrow-segParams{X}.edges(3)-segParams{X}.edges(4)]+0.5},...
    [0.1,0.95], [0,0], 'label',false, 'color','w', 'show',false );
movParam.Toffset = 14.5;
[~,Tmov,~,movFrames] = WriteBoutMovies(expt(X), catInfo(X), Tscan{X}, loco{X}, csdBout{X}, movParam, ROI{X}, fiber{X}([1,3])); % , 'overwrite',true (1)  , ROI{X}  , 'scalebar',
%%
close all; clearvars h sp; 
topOpt = {[0.02,0.02], [0.08, 0.01], [0.06, 0.06]};  % {[vert, horz], [bottom, top], [left, right] }
botOpt = {[0.05,0.06], [0.08, 0.01], [0.06, 0.06]};  % {[vert, horz], [bottom, top], [left, right] }
CSDspreadFig = figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'Color','w');
subtightplot(3,3,1, topOpt{:});  
imshow( movFrames(19).cdata, [] ); axis image;

subtightplot(3,3,2, topOpt{:});  
imshow( movFrames(22).cdata, [] ); axis image;

subtightplot(3,3,3, topOpt{:});  
imshow( movFrames(25).cdata, [] ); axis image;

%%

%% Plot FiberWaveSpeed results
close all;
opt = {[0.08,0.08], [0.08,0.06], [0.08,0.02]};  % {[vert, horz], [bottom, top], [left, right] }
figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
for x = 37 
    for f = 1
        subtightplot(3,2,1:2,opt{:});
        imshow( label2rgb(fiber{x}(f).labelFoot ), [] ); %hold on;
        
        subtightplot(3,2,3,opt{:});
        imagesc( periBout{x}(1).fluor{1}(:,fiber{x}(f).ROI)' );
        title('Locomotion-Associated');
        caxis([-1,3]);
        colorbar;
        
        subtightplot(3,2,5,opt{:});
        imagesc( csdBout{x}(3).fluor{1}(csdWave{x}.scans,fiber{x}(f).ROI)' ); 
        title('Spreading Depression');
        caxis([-1,3]);
        colorbar;
        
        if isfield(fiber{x}(f), 'spont')
            subtightplot(3,2,4,opt{:});
            plot( fiber{x}(f).spont.sepDelay(:,1), fiber{x}(f).spont.sepDelay(:,2), '.' ); hold on;
            plot( fiber{x}(f).spont.sepDelay(:,1), fiber{x}(f).spont.sepDelayPred, 'r--' );
            axis square;
            xlim([0,500]); ylim([-3,3]);
            xlabel('Pairwise Separation (um)'); ylabel('Pairwise Lag (s)'); 
        end
        
        if isfield(fiber{x}(f), 'CSD')
            subtightplot(3,2,6,opt{:});
            plot( fiber{x}(f).CSD.sepDelay(:,1), fiber{x}(f).CSD.sepDelay(:,2), '.' ); hold on;
            plot( fiber{x}(f).CSD.sepDelay(:,1), fiber{x}(f).CSD.sepDelayPred, 'r--' );
            axis square;
            xlim([0,500]); ylim([-3,3]);
            xlabel('Pairwise Separation (um)'); ylabel('Pairwise Lag (s)'); 
        end
    end
end