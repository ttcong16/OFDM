clearvars; close all;

%% Parameters
n_subcarriers = 128;
n_symbols = 7;
scs = 15e3;
fs = n_subcarriers * scs; %#ok<NASGU>
qam_order = 16;

%% OFDM Grid
qam_symbols = randi(qam_order, n_subcarriers, n_symbols) - 1;

% MATLAB 2015 replacement for:
% ofdm_grid = qammod(qam_symbols, qam_order, 'gray','UnitAveragePower', true);
ofdm_grid = qammod(qam_symbols, qam_order, 0, 'gray');
ref = qammod((0:qam_order-1).', qam_order, 0, 'gray');
scale = sqrt(mean(abs(ref).^2));
ofdm_grid = ofdm_grid ./ scale;

%% Constellation Plot
figure;
plot(real(ofdm_grid), imag(ofdm_grid), '.', ...
    'MarkerSize', 30, 'Color', [1 0.341 0.341]); % '#FF5757' -> RGB
axlim = max(max(abs(ofdm_grid))) + 0.05;
ylim([-axlim, axlim]); xlim([-axlim, axlim]); axis square;
xlabel('In-phase'); ylabel('Quadrature');

%% OFDM Modulation Using Transformation Matrix
k = (0 : n_subcarriers - 1).';
n = (0 : n_subcarriers - 1);
dft_matrix = exp(-2i * pi / n_subcarriers * k * n);
idft_matrix = dft_matrix' / n_subcarriers;

ofdm_grid_shifted = [
    ofdm_grid(n_subcarriers/2 + 1 : end, :); ...
    ofdm_grid(1 : n_subcarriers/2, :)];
ofdm_signal_idft = idft_matrix * ofdm_grid_shifted;
ofdm_signal_idft = reshape(ofdm_signal_idft, [], 1);

%% OFDM Modulation Using IFFT
ofdm_signal_ifft = ifft(ifftshift(ofdm_grid, 1));
ofdm_signal_ifft = reshape(ofdm_signal_ifft, [], 1);

error_idft_ifft = max(abs(ofdm_signal_ifft - ofdm_signal_idft));
fprintf('IDFT vs IFFT error: %g\n', error_idft_ifft);

%% OFDM Modulation Using ofdmmod Function (MATLAB 2015 equivalent)
% Original: ofdmmod(ofdm_grid, n_subcarriers, 0)
% With CP=0, this is exactly IFFT with the same shifting convention used above.
ofdm_signal_ofdmmod = ifft(ifftshift(ofdm_grid, 1));
ofdm_signal_ofdmmod = reshape(ofdm_signal_ofdmmod, [], 1);

error_ifft_ofdmmod = max(abs(ofdm_signal_ofdmmod - ofdm_signal_ifft));
fprintf('IFFT vs OFDMMOD error: %g\n', error_ifft_ofdmmod);

%% OFDM Demodulation Using ofdmdemod Function (MATLAB 2015 equivalent)
% Original: ofdmdemod(ofdm_signal_ofdmmod, n_subcarriers, 0)
% With CP=0: reshape -> FFT -> fftshift (to undo ifftshift)
rx_mat = reshape(ofdm_signal_ofdmmod, n_subcarriers, []);
ofdm_grid_reconstructed = fftshift(fft(rx_mat, n_subcarriers, 1), 1);

error_reconstruction = max(max(abs(ofdm_grid_reconstructed - ofdm_grid)));
fprintf('Reconstruction error: %g\n', error_reconstruction);
