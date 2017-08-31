function dict_video = show_dictionary_video(dictionary, h, w, f)
    % Function to visualize the learnt video dictionary as a gif.
        
    npatches = size(dictionary, 2);
    nimg_row = ceil(sqrt(npatches));
    nimg_col = ceil(npatches/nimg_row);
    dict_video = zeros(nimg_row, nimg_col, f);
    
    atom_idx = 1;
    
    for r = 1:nimg_row
        for c = 1:nimg_col
            vid_patch = reshape(dictionary(:, atom_idx), h, w, f);
            vid_patch = (vid_patch-min(vid_patch(:)))/(max(vid_patch(:))-min(vid_patch(:)));
            dict_video((r-1)*h+1:r*h, (c-1)*w+1:c*w, :) = ...
                vid_patch;
            atom_idx = atom_idx + 1;
			if atom_idx > npatches
				return;
			end
        end
    end
end