function code_block = get_code_block(h, w, t)
    % Function to get a coding block for coding a  video or hyperspectral
    % image, such that there is only one sample across the third axis.
    %
    % Inputs:
    %   h: Height of block
    %   w: Width of block
    %   t: Number of frames.
    %
    % Output:
    %   code_block: hxwxt code block.
    
    code_block = zeros(h, w, t);
    
    % if the number of elements in a frame is less than the number of
    % frames, then we can't code.
    if h*w < t
        fprintf('Number of elements less than number of frames\n');
        code_block = 0;
        return;
    end
    
    % Create random indices.
    rand_idx = randperm(h*w);
    n = floor(h*w/t);
    % First pass, just keep filling in the blocks.
    for idx = 1:t
        tmp = zeros(h, w);
        tmp(rand_idx((idx-1)*n + 1:idx*n)) = 1;
        code_block(:, :, idx) = tmp;
    end
    
    % Second pass, finish the left over.
    n_idx = t*n + 1;
    while n_idx <= h*w
        tmp = zeros(h, w);
        tmp(rand_idx(n_idx)) = 1;
        code_block(:, :, n_idx - t*n) = code_block(:, :, n_idx - t*n) + tmp;
        n_idx = n_idx + 1;
    end
end