function [mask, mask_np, varargout] = PCAICA_2D_extract(zproj_interp, varargin) % path_interp, , Nx, Ny, Nt  trace_raw, trace_hp, Nf
p = inputParser;
%addOptional(p,'crop_amt',0.1); %proportion to crop edges
addOptional(p,'axons',false); %axon flag for pipe.extract.pcaica
addOptional(p,'blur',1); %gaussian blur width for zproj trace extraction
addOptional(p,'hp_cutoff',20); %high pass filter cutoff for traces
addOptional(p,'edges',[0,0,0,0]);
addParameter(p,'overwrite',false)
addParameter(p,'save','',@ischar)
parse(p,varargin{:});
p = p.Results;
tic;
% Metadata
[Nrow, Ncol, Nscan] = size( zproj_interp );
%{
[fdir, fname, ~] = fileparts(path_interp);
% PCA/ICA segmentation of interpolated movie
icaPath = strcat(fdir,filesep,fname,'.icadata');
if ~exist(icaPath,'file')
    icaguidata = pipe.extract.pcaica(zproj_interp,'axons',p.axons,'smoothing_width',4,'npcs',2000);
    save(icaPath,'icaguidata');
else
    load(icaPath, '-mat'); % ,'icaguidata'
end
%}
icaguidata = pipe.extract.pcaica(zproj_interp, 'axons',p.axons, 'smoothing_width',4, 'npcs',2000, 'firstpctokeep',2);
toc
%reshape all of the masks
F = [icaguidata.ica.filter];
F = reshape(F,Nrow,Ncol,[]);
F_bw = F>0; %binarize the masks

%set an overlap threshold for the edges
overlap = 100;
%make a mask of the edges of the video

%crop edges of video
Xrange = p.edges(3)+1:Nrow-p.edges(4); %round(p.crop_amt*Nx):round((1-p.crop_amt)*Nx);
Yrange = p.edges(1)+1:Ncol-p.edges(2); % round(p.crop_amt*Ny):round((1-p.crop_amt)*Ny);

edge = zeros(size(F_bw,1),size(F_bw,2));
edge([Xrange(1),Xrange(numel(Xrange))],:) = 1;
edge(:,[Yrange(1),Yrange(numel(Yrange))]) = 1;
%imshow( edge )

good_idx = 1;
for m = 1:size(F_bw,3)
    ROI = F_bw(:,:,m); %get temporary ROI
    %check if that ROI overlaps with the edges too much
    if ROI(:)'*edge(:) > overlap
        continue;
    else
        F_filter(:,:,good_idx) = ROI;
        good_idx = good_idx + 1;
    end
end
Nf = size(F_filter,3); %number of 2D filters
toc

% extract trace estimates from 2d mask and zproj movie
mask = cell(1,Nf); mask_np = cell(1,Nf);  trace_raw = cell(1,Nf); trace_hp = cell(1,Nf);
for m = 1:Nf
    mask{m} = F_filter(:,:,m);
    mask_np{m} = imdilate(mask{m},strel('disk',3)) & ~imdilate(mask{m},strel('disk',1)); %make neuropil mask (2 pixels wide, offset by 2 pixels from mask)
end
if nargout > 2
    fprintf('\nExtracting traces from z-projection... ');
    zproj_interp = imgaussfilt(double(zproj_interp), p.blur);
    for m = 1:Nf % parfor is slower
        %Npx = sum(mask{m}(:)); %number of pixels in the mask
        trace_raw{m} = mask{m}(:)'*reshape(zproj_interp,[],Nscan); %get raw estimate trace from the mean projected movie
        trace_hp{m} = trace_raw{m} - movmedian(trace_raw{m},p.hp_cutoff); %high-pass filter
    end
    varargout{1} = trace_raw;
    varargout{2} = trace_hp;
    varargout{3} = Nf;
end
toc

% Save results
if ~isempty(p.save)
    fprintf('\nSaving %s', p.save);
    save(p.save, 'mask','mask_np','Nf','trace_raw','trace_hp', 'p', '-v7.3');
    toc
end

end