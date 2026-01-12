
%   SISO OFDM + AWGN + Rayleigh (1 tap + 5 tap)
%   Equalizer: ZF + MMSE
%   BER + CONSTELLATION + BER vs SNR

clear; clc; close all;

%% ---THAM S? OFDM ---
N  = 64;          % s? subcarrier
CP = 16;          % cyclic prefix
numSymbols = 100; % s? OFDM symbol

bitsPerSymbol = 2; % QPSK
numBits = N * numSymbols * bitsPerSymbol;

%% --- SINH BIT ---
bits = randi([0 1], numBits, 1);

%% --- QPSK MODULATION (Gray mapping) ---
bitPairs = reshape(bits, 2, []).';
numSyms = length(bitPairs);

symbols = zeros(numSyms,1);
for k = 1:numSyms
    b1 = bitPairs(k,1);
    b2 = bitPairs(k,2);
    if     b1==0 && b2==0, symbols(k) = ( 1 + 1j)/sqrt(2);
    elseif b1==0 && b2==1, symbols(k) = (-1 + 1j)/sqrt(2);
    elseif b1==1 && b2==1, symbols(k) = (-1 - 1j)/sqrt(2);
    else                   symbols(k) = ( 1 - 1j)/sqrt(2);
    end
end

%% --- OFDM IFFT ---
ofdm_frame = reshape(symbols, N, numSymbols);
ifft_out = ifft(ofdm_frame, N);

%% --- ADD CYCLIC PREFIX ---
with_cp = [ifft_out(end-CP+1:end,:); ifft_out];

%% --- GHÉP SERIAL ---
tx_signal = with_cp(:);

disp('OFDM TX hoan tat!');

%   KÊNH AWGN + RAYLEIGH 1 TAP + RAYLEIGH 5 TAP
SNR_dB = 20;

%% --- AWGN ---
rx_awgn = awgn(tx_signal, SNR_dB, 'measured');

%% --- Rayleigh 1 tap ---
h1 = (randn + 1j*randn)/sqrt(2);
rx_flat = h1 * tx_signal;
rx_flat = awgn(rx_flat, SNR_dB, 'measured');

%% --- Rayleigh 5 tap (Multipath) ---
h5 = (randn(1,5)+1j*randn(1,5))/sqrt(2);
rx_5tap = conv(tx_signal, h5);
rx_5tap = rx_5tap(1:length(tx_signal));
rx_5tap = awgn(rx_5tap, SNR_dB, 'measured');


%  RX: FFT + DEMOD (AWGN)
rx_mat = reshape(rx_awgn, N+CP, numSymbols);
rx_noCP = rx_mat(CP+1:end,:);
rx_fft = fft(rx_noCP, N);
rx_syms = rx_fft(:);

dem_awgn = zeros(numBits,1);
for k = 1:numSyms
    re = real(rx_syms(k)); im = imag(rx_syms(k));
    if re>=0 && im>=0, dem_awgn(2*k-1:2*k)=[0;0];
    elseif re<0 && im>=0, dem_awgn(2*k-1:2*k)=[0;1];
    elseif re<0 && im<0, dem_awgn(2*k-1:2*k)=[1;1];
    else, dem_awgn(2*k-1:2*k)=[1;0];
    end
end
BER_awgn = sum(bits~=dem_awgn)/numBits;

%  RX: RAYLEIGH 1 TAP — ZF + MMSE

rx_flat_mat = reshape(rx_flat, N+CP, numSymbols);
rx_flat_noCP = rx_flat_mat(CP+1:end,:);
rx_flat_fft = fft(rx_flat_noCP, N);

H1 = h1 * ones(N, numSymbols);

% ZF Equalizer
rx_flat_ZF = rx_flat_fft ./ H1;

% MMSE Equalizer
sigma2 = 10^(-SNR_dB/10);
rx_flat_MMSE = (conj(H1) ./ (abs(H1).^2 + sigma2)) .* rx_flat_fft;


%% --- DEMOD 1 TAP ZF ---
rx1 = rx_flat_ZF(:);
dem_flat = zeros(numBits,1);
for k = 1:numSyms
    re = real(rx1(k)); im = imag(rx1(k));
    if re>=0 && im>=0, dem_flat(2*k-1:2*k)=[0;0];
    elseif re<0 && im>=0, dem_flat(2*k-1:2*k)=[0;1];
    elseif re<0 && im<0, dem_flat(2*k-1:2*k)=[1;1];
    else, dem_flat(2*k-1:2*k)=[1;0];
    end
end
BER_flat_ZF = sum(bits~=dem_flat)/numBits;

%% --- DEMOD 1 TAP MMSE ---
rx1m = rx_flat_MMSE(:);
dem_flat2 = zeros(numBits,1);
for k = 1:numSyms
    re = real(rx1m(k)); im = imag(rx1m(k));
    if re>=0 && im>=0, dem_flat2(2*k-1:2*k)=[0;0];
    elseif re<0 && im>=0, dem_flat2(2*k-1:2*k)=[0;1];
    elseif re<0 && im<0, dem_flat2(2*k-1:2*k)=[1;1];
    else, dem_flat2(2*k-1:2*k)=[1;0];
    end
end
BER_flat_MMSE = sum(bits~=dem_flat2)/numBits;


%  RX: RAYLEIGH 5 TAP — ZF + MMSE

rx5_mat = reshape(rx_5tap, N+CP, numSymbols);
rx5_noCP = rx5_mat(CP+1:end,:);
rx5_fft = fft(rx5_noCP, N);

H5f = fft(h5, N).';
H5mat = repmat(H5f, 1, numSymbols);

% ZF
rx5_ZF = rx5_fft ./ H5mat;

% MMSE
rx5_MMSE = (conj(H5mat) ./ (abs(H5mat).^2 + sigma2)) .* rx5_fft;


%% --- DEMOD 5 TAP ZF ---
rx5 = rx5_ZF(:);
dem_5tap = zeros(numBits,1);
for k = 1:numSyms
    re = real(rx5(k)); im = imag(rx5(k));
    if re>=0 && im>=0, dem_5tap(2*k-1:2*k)=[0;0];
    elseif re<0 && im>=0, dem_5tap(2*k-1:2*k)=[0;1];
    elseif re<0 && im<0, dem_5tap(2*k-1:2*k)=[1;1];
    else, dem_5tap(2*k-1:2*k)=[1;0];
    end
end
BER_5tap_ZF = sum(bits~=dem_5tap)/numBits;

%% --- DEMOD 5 TAP MMSE ---
rx5m = rx5_MMSE(:);
dem_5tap2 = zeros(numBits,1);
for k = 1:numSyms
    re = real(rx5m(k)); im = imag(rx5m(k));
    if re>=0 && im>=0, dem_5tap2(2*k-1:2*k)=[0;0];
    elseif re<0 && im>=0, dem_5tap2(2*k-1:2*k)=[0;1];
    elseif re<0 && im<0, dem_5tap2(2*k-1:2*k)=[1;1];
    else, dem_5tap2(2*k-1:2*k)=[1;0];
    end
end
BER_5tap_MMSE = sum(bits~=dem_5tap2)/numBits;


%  K?T QU? BER

fprintf('\n===== BER SUMMARY =====\n');
fprintf('AWGN                 : %.6f\n', BER_awgn);
fprintf('Rayleigh 1 tap ZF    : %.6f\n', BER_flat_ZF);
fprintf('Rayleigh 1 tap MMSE  : %.6f\n', BER_flat_MMSE);
fprintf('Rayleigh 5 tap ZF    : %.6f\n', BER_5tap_ZF);
fprintf('Rayleigh 5 tap MMSE  : %.6f\n', BER_5tap_MMSE);


%  CONSTELLATION PLOT

figure;
subplot(1,3,1);
plot(real(symbols), imag(symbols), 'b.');
title('Original QPSK'); axis equal; grid on;

subplot(1,3,2);
plot(real(rx5_fft(:)), imag(rx5_fft(:)), 'r.');
title('5 tap - Before Equalizer'); axis equal; grid on;

subplot(1,3,3);
plot(real(rx5_MMSE(:)), imag(rx5_MMSE(:)), 'g.');
title('5 tap - After MMSE'); axis equal; grid on;


%  BER vs SNR cho AWGN / Rayleigh ZF / Rayleigh MMSE

SNR_range = 0:2:30;
BER_awgn_curve = zeros(size(SNR_range));
BER_r5_ZF = zeros(size(SNR_range));
BER_r5_MMSE = zeros(size(SNR_range));

for idx = 1:length(SNR_range)
    SNR = SNR_range(idx);
    sigma2 = 10^(-SNR/10);

    %% --- AWGN ---
    rx_awgn2 = awgn(tx_signal, SNR, 'measured');
    mat = reshape(rx_awgn2, N+CP, numSymbols);
    sy = fft(mat(CP+1:end,:), N);
    dem = QPSK_demod(sy(:));
    BER_awgn_curve(idx) = sum(bits~=dem)/numBits;

    %% --- 5 tap ---
    h5 = (randn(1,5)+1j*randn(1,5))/sqrt(2);
    r = conv(tx_signal,h5); r = r(1:length(tx_signal));
    r = awgn(r, SNR, 'measured');
    mat = reshape(r, N+CP, numSymbols);
    sy = fft(mat(CP+1:end,:), N);
    H = fft(h5, N).'; H = repmat(H,1,numSymbols);

    % ZF
    eq_ZF = sy./H;
    dem = QPSK_demod(eq_ZF(:));
    BER_r5_ZF(idx) = sum(bits~=dem)/numBits;

    % MMSE
    eq_MMSE = (conj(H)./(abs(H).^2 + sigma2)).*sy;
    dem = QPSK_demod(eq_MMSE(:));
    BER_r5_MMSE(idx) = sum(bits~=dem)/numBits;
end

figure;
semilogy(SNR_range, BER_awgn_curve,'bo-','LineWidth',2); hold on;
semilogy(SNR_range, BER_r5_ZF, 'rs--','LineWidth',2);
semilogy(SNR_range, BER_r5_MMSE,'g^-','LineWidth',2);
grid on; xlabel('SNR (dB)'); ylabel('BER');
legend('AWGN','Rayleigh 5 tap ZF','Rayleigh 5 tap MMSE');
title('BER theo SNR (OFDM SISO)');

