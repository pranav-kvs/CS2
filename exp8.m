clearvars; clc;

nSym = 10^4; % Number of OFDM Symbols to transmit
EbN0dB = 0:2:20; % Bit to noise ratio
MOD_TYPE = 'MQAM'; % Modulation type 'MPSK' or 'MQAM'
M = 64; % Choose modulation order for the chosen MOD_TYPE
%N = 64; % FFT size or total number of subcarriers (used + unused) 64
Ncp = 16; % Number of symbols in the cyclic prefix

k = log2(M); % Number of bits per modulated symbol
EsN0dB = 10*log10(k) + EbN0dB; % Convert to symbol energy to noise ratio
errors = zeros(1, length(EsN0dB)); % To store symbol errors

% Monte Carlo Simulation
for i = 1:length(EsN0dB)
    for j = 1:nSym

        % Transmitter
        d = ceil(M .* rand(1, M)); % Uniform distributed random syms from 1:M
        [X, ref] = modulation_mapper(MOD_TYPE, M, d); 
        x = ifft(X, M); % IDFT
        s = add_cyclic_prefix(x, Ncp); % Add CP

        % Channel
        r = add_awgn_noise(s, EsN0dB(i)); % Add AWGN noise

        % Receiver
        y = remove_cyclic_prefix(r, Ncp, M); % Remove CP
        Y = fft(y, M); % DFT
        [~, dcap] = iqOptDetector(Y, ref); % Demapper using IQ detector

        % Error counter
        d = d(:).';  % Ensure d is a row vector
        dcap = dcap(:).';  % Ensure dcap is a row vector
        numErrors = sum(d ~= dcap); % Count number of symbol errors
        errors(i) = errors(i) + numErrors; % Accumulate symbol errors
    end
end

SER_sim = errors / (nSym * M);

SER_theory = ser_awgn(EbN0dB, MOD_TYPE, M);

% Plot theoretical curves and simulated BER points
plot(EbN0dB, log10(SER_sim), 'ro'); 
hold on;
plot(EbN0dB, log10(SER_theory), 'r'); 
grid on;

title(['Performance of ', num2str(M), ' ', MOD_TYPE, ' OFDM over AWGN channel']);
xlabel('Eb/N0 (dB)');
ylabel('Symbol Error Rate');
%legend('simulated', 'theoretical');                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     