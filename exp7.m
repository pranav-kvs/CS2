%% DSSS

M = 2; % BPSK modulation
numBits = 1000; % Number of bits
chipRate = 10; % Chips per bit
snr = 10; % Signal-to-noise ratio in dB
fs = 1000; % Sampling frequency

% Generate random bits
dataBits = randi([0 1], numBits, 1);
modulatedData = 2*dataBits - 1; % BPSK mapping (0 -> -1, 1 -> 1)

% Generate PN sequence (chips)
pnSequence = randi([0 1], numBits*chipRate, 1);
pnSequence = 2*pnSequence - 1; % Convert to bipolar (-1, 1)

% Spread the signal
spreadSignal = repelem(modulatedData, chipRate) .* pnSequence;

% Transmit over AWGN channel
receivedSignal = awgn(spreadSignal, snr, 'measured');

% Despread the signal
despreadSignal = receivedSignal .* pnSequence;
despreadBits = sum(reshape(despreadSignal, chipRate, numBits), 1)' / chipRate;

% Calculate FFT for plotting spectra
n = length(spreadSignal);
frequencies = (-n/2:n/2-1)*(fs/n);

% Message signal spectrum
messageSpectrum = fftshift(fft(modulatedData, n));
% PN code spectrum
pnSpectrum = fftshift(fft(pnSequence, n));
% Modulated signal spectrum (spread signal)
spreadSignalSpectrum = fftshift(fft(spreadSignal, n));
% Received signal spectrum
receivedSignalSpectrum = fftshift(fft(receivedSignal, n));
% Despread signal spectrum
despreadSignalSpectrum = fftshift(fft(despreadSignal, n));
% Demodulated signal spectrum (before decision device)
demodulatedSignalSpectrum = fftshift(fft(despreadBits, n));

% Plotting
figure;
% Message signal spectrum
subplot(3,2,1);
plot(frequencies, abs(messageSpectrum));
title('Message Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
% PN code spectrum
subplot(3,2,2);
plot(frequencies, abs(pnSpectrum));
title('PN Code Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
% Spread signal spectrum
subplot(3,2,3);
plot(frequencies, abs(spreadSignalSpectrum));
title('Spread Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
% Received signal spectrum
subplot(3,2,4);
plot(frequencies, abs(receivedSignalSpectrum));
title('Received Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
% Despread signal spectrum
subplot(3,2,5);
plot(frequencies, abs(despreadSignalSpectrum));
title('Despread Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

subplot(3,2,6);
plot(frequencies, abs(demodulatedSignalSpectrum));
title('Demodulated Signal Spectrum (Before Decision Device)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
fprintf('DSSS BER: %e\n', BER_DSSS);

%% FHSS

numBits = 20; % Number of bits
numHops = 6; % Number of frequency hops
hopFrequencies = [1, 2, 3, 4, 5, 6] * 1e3; % Frequencies in Hz
hopDuration = 1e-3; % Duration of each hop in seconds
snr = 10; % Signal-to-noise ratio in dB
fs = 20e3; % Sampling frequency in Hz

% Generate random bits
dataBits = randi([0 1], numBits, 1);
modulatedData = 2*dataBits - 1; % BPSK mapping (0 -> -1, 1 -> 1)

bpskSignal = repelem(modulatedData, hopDuration*fs); % Create rectangular wave
t = (0:1/fs:hopDuration-1/fs)';

hopSignal = [];
fhssSignal = [];

for i = 1:numBits
    hopId = mod(i-1, numHops) + 1;
    freq = hopFrequencies(hopId);
    carrier = cos(2*pi*freq*t);
    hopSignal = [hopSignal; carrier];
    fhssSignal = [fhssSignal; bpskSignal((i-1)*length(t)+1:i*length(t)) .* carrier];
end

% Transmit over AWGN channel
receivedSignal = awgn(fhssSignal, snr, 'measured');
receivedBits = zeros(numBits, 1);

% Demodulate FHSS signal
demodulatedSignal = [];

for i = 1:numBits
    hopId = mod(i-1, numHops) + 1;
    freq = hopFrequencies(hopId);
    carrier = cos(2*pi*freq*t);
    segment = receivedSignal((i-1)*length(t)+1:i*length(t));
    demodulated = segment .* carrier;
    demodulatedSignal = [demodulatedSignal; demodulated];
    receivedBits(i) = sum(demodulated) > 0;
end

% Plotting
figure;
% Original 20-bit sequence
subplot(3, 2, 1);
stem(dataBits, 'filled');
title('Original 20-bit Sequence');
xlabel('Bit Index');
ylabel('Bit Value');
grid on;

% BPSK modulated signal (rectangular wave)
subplot(3, 2, 2);
plot(bpskSignal);
title('BPSK Modulated Signal (Rectangular Wave)');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

% Spread signal with 6 frequencies
subplot(3, 2, 3);
plot(hopSignal);
title('Spread Signal with 6 Frequencies');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

% Frequency ho
pped spread spectrum signal at transmitter
subplot(3, 2, 4);
plot(fhssSignal);
title('Frequency Hopped Spread Spectrum Signal');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

% Demodulated BPSK signal from widespread signal
subplot(3, 2, 5);
plot(demodulatedSignal);
title('Demodulated BPSK Signal from Spread Signal');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

% Original transmitted bits sequence at the receiver
subplot(3, 2, 6);
stem(receivedBits, 'filled');
title('Original Transmitted Bits Sequence at Receiver');
xlabel('Bit Index');
ylabel('Bit Value');
grid on;

% Calculate and display BER
%BER_FHSS = sum(dataBits ~= receivedBits) / numBits;
%fprintf('FHSS BER: %e\n', BER_FHSS);
