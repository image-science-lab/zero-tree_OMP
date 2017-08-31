function val = rsnr(x, xhat)
    % Return reconstruction SNR in dB.
    val = 20*log10(norm(x(:))/norm(x(:) - xhat(:)));
end