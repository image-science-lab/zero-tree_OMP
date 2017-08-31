function atom_im = atom_blow(atom, h_idx, w_idx)
    % Function to blow up a video atom into xy, xt and yt slices
    %
    % Inputs:
    %   atom: 3D video atom
    %   h_idx: Row of the xt-slice
    %   w_idx: Column of the yt-slice
    %
    % Outputs:
    %   atom_im: Image of the atom
    
    [h, w, t] = size(atom);
    atom = mnormalize(atom);
    
    atom_im = ones(round(1.1*h + t), round(1.1*w + t));
    atom_im(1:h, 1:w) = atom(:, :, 1);
    atom_im(h_idx, :) = 1;
    atom_im(:, w_idx) = 1;
    atom_im(1:h, end-t+1:end) = squeeze(atom(:, w_idx, :));
    atom_im(end-t+1:end, 1:w) = squeeze(atom(h_idx, :, :))';
end