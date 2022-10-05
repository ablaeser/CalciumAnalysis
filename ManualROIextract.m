function [ ROI, Nroi, F, dFF, C, S, O, Sbout ] = ManualROIextract( metadata, T, varargin )

IP = inputParser;
addRequired( IP, 'metadata', @isstruct )
addRequired( IP, 'T', @isnumeric )
addOptional( IP, 'duraMovie', NaN, @isnumeric )
%addParameter( IP, 'server', [] )
addParameter( IP, 'pad', 3, @isnumeric )
addParameter( IP, 'gaussWidth', 0.25, @isnumeric )
addParameter( IP, 'gaussSigma', 0.1, @isnumeric )
addParameter( IP, 'show', false, @islogical )
addParameter( IP, 'save', true, @islogical )
addParameter( IP, 'overwrite', false, @islogical )
parse( IP, metadata, T, varargin{:} ); 
duraMovie = IP.Results.duraMovie;
padDist = IP.Results.pad;
gaussWidth = IP.Results.gaussWidth;
gaussSigma = IP.Results.gaussSigma;
show = IP.Results.show;
saveToggle = IP.Results.save;
overwrite = IP.Results.overwrite;

fileRoot = sprintf('%s%s_%s_run%i', metadata.folder, metadata.mouse, metadata.date, metadata.run );
matFile = sprintf('%s-ManualFluor.mat', fileRoot ); 

tic
if ~exist(matFile,'file') || overwrite
    % Get the ROIs
    [ ROI, Nroi ] = MakeManualROI( metadata, padDist );
    if show
        figure('Units','normalized','OuterPosition',[0,0,1,1]);
        imshow( metadata.std, [] ); hold all;
        for r = 1:Nroi-1
            plot( ROI(r).edge(:,1), ROI(r).edge(:,2), '--' );
        end
        plot( ROI(Nroi).edge(:,1), ROI(Nroi).edge(:,2), 'w--' );
        title( sprintf('%s - %s, run %i: Std Dev Projection and %i ROIs', metadata.mouse, metadata.date, metadata.run, Nroi) );
        orient( gcf, 'landscape' )
        figPath = [fileRoot, '-ManualROIs.pdf'];
        print( gcf, figPath, '-dpdf', '-fillpage', '-r0' );
    end
    if nargout > 2 && numel(duraMovie) > 1
        % Extract signals 
        F = ExtractSignal( metadata, duraMovie, ROI );
        %if show, CalSpread( {T}, {F}, 'title','Raw Data'  ); end  % , Tshade
        if nargout > 3
            % Filter the traces (Gaussian)
            gaussFilt = MakeGaussFilt( gaussWidth, 0, gaussSigma, metadata.rate, false ); % 
            Flow =  filtfilt( gaussFilt, 1, F ); % double(  )
            fLow = Flow(:,1:Nroi-1) - repmat( Flow(:,Nroi), 1, Nroi-1 ); % Background subtraction
            % convert to dF/Fo
            fprintf('\n\nConverting to dF/F...\n');
            dFF = zeros( metadata.Ngood, Nroi );   brightFrames = cell(1,Nroi); %Fbase = zeros(1,Nroi);
            Npeak = 3; PeakDist = round(5*metadata.rate);
            for r = 1:Nroi-1
                [~, ~, Otemp] = deconvolveCa( fLow(:,r), 'ar1', 'constrained', 'optimize_b', 'optimize_pars' ); % use deconvolveCa to estimate Fo
                ROI(r).base = Otemp.b;  %Fbase(r) = Otemp.b;
                dFF(:,r) = (fLow(:,r) - ROI(r).base)/ROI(r).base;
                [~,brightFrames{r}] = findpeaks( dFF(:,r), 'MinPeakDistance', PeakDist, 'Npeaks',Npeak, 'SortStr','descend' );
            end
            [~, ~, Oback] = deconvolveCa( Flow(:,Nroi), 'ar1', 'constrained', 'optimize_b', 'optimize_pars' ); % use deconvolveCa to estimate Fo
            ROI(Nroi).base = Oback.b; %Fbase(Nroi) = Oback.b;
            dFF(:,Nroi) = (Flow(:,Nroi) - Oback.b)/Oback.b;
            % Create summary images for each ROI based on brightFrames
            [~,brightFrames{Nroi}] = findpeaks( dFF(:,Nroi), 'MinPeakDistance', PeakDist, 'Npeaks',Npeak, 'SortStr','descend' );
            for r = 1:Nroi
                ROI(r).bright = mean( duraMovie(:,:,brightFrames{r}), 3); 
            end
            % Deconvolve the dF/F signals
            if nargout > 4
                tic
                fprintf('\nStarting deconvolution... '); 
                C = zeros(metadata.Ngood, Nroi); S = zeros(metadata.Ngood, Nroi); Sbout = cell(1,Nroi);
                for r = Nroi:-1:1 
                    fprintf('\n%i / %i', r, Nroi );
                    [Ctemp, Stemp, Otemp] = deconvolveCa( dFF(:,r), 'ar1', 'constrained', 'optimize_b', 'optimize_pars' ); % , 'optimize_smin'
                    C(:,r) = Ctemp; S(:,r) = Stemp; O(r) = Otemp;
                    [~, Sbout{r}, ~] = speedBouts(T, S(:,r), 'minFrame', 1, 'minPeak', 0, 'mergeSep', 0, 'finMinFrm', 1, 'show',false, 'pause',false ); 
                end
                %Stot = sum(S,2);
                fprintf('  Done!\n'); toc
                if saveToggle
                    fprintf('\nSaving %s...  ', matFile ); tic;
                    save( matFile, 'metadata', 'ROI', 'Nroi', 'T', 'F', 'fLow', 'dFF', 'C', 'S', 'O', 'Sbout' ); toc
                end
            end
        end
    end
else
    fprintf('\nLoading %s   ', matFile ); 
    load( matFile ); toc
end

if show && nargout > 2
    AllSpread = figure('units','normalized','outerposition',[0 0 1 1],'color','w'); % ,'FontSize',20
    set(AllSpread,'DefaultAxesFontSize',16); % ,'DefaultAxesTickDir','out'
    opt = {[0.05,0.08], [0.09,0.04], [0.06,0.06]};  % {[vert, horz], [bottom, top], [left, right] }
    sp(1) = subtightplot(2,2,1,opt{:});
    Fspread = CalSpread( {T}, {F}, 'show',false );
    plot( T, Fspread{1} ); ylabel('F');  %ylim([-3, Inf]);
    set(gca,'TickDir','out', 'box','off','Xtick',[]); axis tight;
    
    sp(2) = subtightplot(2,2,2,opt{:});
    dFFspread = CalSpread( {T}, {dFF}, 'show',false );
    plot( T, dFFspread{1} ); ylabel('dF/F'); % ylim([-2, Inf]);
    set(gca,'TickDir','out', 'box','off','Xtick',[]); axis tight;
    
    sp(3) = subtightplot(2,2,3,opt{:});
    Cspread = CalSpread( {T}, {C}, 'show',false );
    plot( T, Cspread{1} ); xlabel('Time (s)'); ylabel('C'); % ylim([-2, Inf]);
    set(gca,'TickDir','out', 'box','off'); axis tight;
    
    sp(4) = subtightplot(2,2,4,opt{:});
    Sspread = CalSpread( {T}, {S}, 'show',false );
    plot( T, Sspread{1} ); xlabel('Time (s)'); ylabel('Activity'); %ylim([-2, Inf]);
    set(gca,'TickDir','out', 'box','off'); axis tight;
    linkaxes(sp,'x');
    xlim([T(1), T(end)]);
    
    figPath = sprintf('%s%s_%s_run%i_Fluor.pdf',metadata.folder, metadata.mouse, metadata.date, metadata.run ); % %[metadata.folder, meta 'FluorExtraction.pdf'];
    if saveToggle && (~exist(figPath,'file') || overwrite) 
        orient( AllSpread, 'landscape' )
        fprintf('\nSaving %s\n', figPath );
        print( AllSpread, figPath, '-dpdf', '-fillpage', '-r0' ); 
    end
end
end