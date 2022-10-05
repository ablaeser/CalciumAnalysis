function [mask_out, Nmask_out, overlap] = FilterMasks(mask_in, varargin)
IP = inputParser;
addRequired( IP, 'mask_in', @iscell )
addParameter( IP, 'break', false, @islogical)
addParameter( IP, 'minFoot', 0, @isnumeric )
addParameter( IP, 'edge', [0,0,0,0], @isnumeric )
addParameter( IP, 'show', false, @islogical )
parse( IP, mask_in, varargin{:} );
breakUp = IP.Results.break;
minFoot = IP.Results.minFoot;
edge = IP.Results.edge;
show = IP.Results.show;

% Filter out blobs below minFoot
NmaskIn = numel(mask_in);
maskDim = size(mask_in{1});
mask_out = []; %cell(0); % cell(1,Nchunk); 
if show
    close all;
    figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
end
k = 0;
mask_in_Npix = cellfun(@sum, mask_in, repmat({'all'},1,NmaskIn));
[~,sortInd] = sort( mask_in_Npix, 'descend' );
mask_in = mask_in(sortInd);
for m = 1:NmaskIn % 1:flip(1:size(mask_in,2))
    % Crop, decompose and filter masks into preliminary ROI
    mask_in{m}(:,1:edge(1)) = false;
    mask_in{m}(:,maskDim(2)-edge(2):end) = false;
    mask_in{m}(1:edge(3),:) = false;
    mask_in{m}(maskDim(1)-edge(4):end,:) = false;
    if show
        subplot(1,2,1);
        imshow(mask_in{m}, [])
        title(sprintf('mask %i', m) );
    end
    tempConnComp = bwconncomp( mask_in{m}, 8 );
    tempSubProps = regionprops( tempConnComp, 'Area' ); % , 'Centroid', 'BoundingBox', 'MeanIntensity', 'MaxIntensity', 'VoxelList', 'VoxelIdxList','PrincipalAxisLength'
    goodSub = find([tempSubProps.Area] >= minFoot);
    if breakUp 
        if ~isempty(goodSub)
            for s = goodSub
                k = k+1;
                mask_out(k).mask = false(maskDim); % 
                mask_out(k).mask(tempConnComp.PixelIdxList{s}) = true;
                mask_out(k).ind = find(mask_out(k).mask(:));
                mask_out(k).area = numel(mask_out(k).ind); 
                mask_out(k).pc = m;
                if show
                    subplot(1,2,2);
                    imshow(mask_out(k).mask, []); 
                    title(sprintf('Filtered (minFoot = %i)', minFoot));
                    impixelinfo;
                    pause%(1);
                end
            end
        else
            fprintf('\nmask %i did not survive\n', m)
        end
    else
        
        if ~isempty(goodSub)
            mask_out(m).mask = false(maskDim); %
            for s = goodSub
                mask_out(m).mask(tempConnComp.PixelIdxList{s}) = true;
                mask_out(m).ind = find(mask_out(m).mask(:));
                mask_out(m).area = numel(mask_out(m).ind); 
                mask_out(m).pc = m;
            end
        else
            fprintf('\nmask %i did not survive\n', m)
            mask_out(m).mask = [];
            mask_out(m).ind = [];
            mask_out(m).area = 0; 
            mask_out(m).pc = m;
        end
        %if isempty(mask_out(m).mask), break; end
        if show
            subplot(1,2,2);
            imshow(mask_out(m).mask, []); 
            title(sprintf('Filtered (minFoot = %i)', minFoot));
            impixelinfo;
            pause%(1);
        end
    end
end
mask_out([mask_out.area] < minFoot) = [];
[~,kSort] = sort([mask_out.area], 'descend');
mask_out = mask_out(kSort);
Nmask_out = numel(mask_out);
%show = true;
if show
    allROIim = zeros(maskDim);
    for r = 1:Nmask_out
        allROIim(mask_out(r).ind) = r;
    end
    allROIim = label2rgb(allROIim );
    figure;
    imshow(allROIim, [])
end


%Merge masks that overlap substantially
tic
overlap = nan(Nmask_out);
for M = 1:Nmask_out-1
    for m = M+1:Nmask_out
        commonInd = intersect(mask_out(M).ind, mask_out(m).ind);
        Ncommon = numel(commonInd);
        overlap(M,m) = Ncommon/mask_out(M).area;
        overlap(m,M) = Ncommon/mask_out(m).area; 
    end
end
toc
if show
    %close all; 
    figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
    imagesc( overlap ); 
    ylabel('Mask M'); xlabel('m'); title('ROI Overlap')
    colorbar; caxis([0,1]);
    axis image;
    impixelinfo;
end

minOverlap = 0.6;
if show
    %close all;
    figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
    sp(4) = subplot(2,2,4); sp(3) = subplot(2,2,3);  sp(2) = subplot(2,2,2); sp(1) = subplot(2,2,1); 
    linkaxes(sp,'xy');
end
tic
for M = 1:Nmask_out % flip( )
    mMatch = find( overlap(M,:) > minOverlap ); % unique([find( overlap(M,:) > minOverlap ), find( overlap(M,:) > minOverlap )]);
    if mask_out(M).area > 0 && ~isempty(mMatch) % any(ROI(M).mask,'all')
        % Pre-absorption
        if show
            subplot(sp(1));
            imshow(mask_out(M).mask, [])
            title( sprintf('mask %i', M) );
        end
        for m = mMatch %minOverlap%flip(1:size(mask_final{c},2))
            if show
                subplot(sp(2));
                imshow(mask_out(m).mask, []); 
                title( sprintf('mask %i. Overlap = [%2.1f,  %2.1f] pct', m, 100*overlap(M,m), 100*overlap(m,M) ) );
                %pause;
            end
            % Absorb mask
            fprintf('\nMerging ROI %i and %i', M, m);
            mask_out(M).mask = mask_out(M).mask | mask_out(m).mask;
            mask_out(m).mask = false(maskDim);
            mask_out(m).area = 0;
            
            % Post-absorption
            if show
                subplot(sp(3));
                imshow(mask_out(M).mask, [])
                title( 'Post-absorption' );
                subplot(sp(4));
                imshow(mask_out(m).mask, []); 
                title( 'Post-absorption' );
                pause(1);
            end
        end
    end
    mask_out(M).ind = find(mask_out(M).mask);
    mask_out(M).area = numel(mask_out(M).ind);
    mask_out(M).pc = [mask_out(M).pc, mask_out(m).pc];
    overlap(:,M) = NaN; % prevent this mask from being absorbed later
end
toc

% Purge blanked masks
mask_out([mask_out.area]==0) = [];
Nmask_out = numel(mask_out);
[~,rSort] = sort([mask_out.area], 'descend'); % 
mask_out = mask_out(rSort);

if show
    allROIim = zeros(maskDim);
    for r = 1:Nmask_out
        allROIim(mask_out(r).ind) = r;
    end
    allROIim = label2rgb(allROIim );
    figure;
    imshow(allROIim, [])
end
end

