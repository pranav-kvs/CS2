
N = 10^5; % Number of input symbols
EbN0dB = -4:5:24; % Define EbN0dB range for simulation
M = 64; % M-QAM modulation order
k = log2(M); % Bits per symbol
g = 0.9; phi = 8; dc_i = 1.9; dc_q = 1.7; % Receiver impairments

EsN0dB = EbN0dB + 10*log10(k); % Converting Eb/N0 to Es/N0
SER1 = zeros(length(EsN0dB), 1); % Symbol Error rates (No compensation)
SER2 = SER1; % Symbol Error rates (DC compensation only)
SER3 = SER1; % Symbol Error rates (DC comp & Blind IQ compensation)
SER4 = SER1; % Symbol Error rates (DC comp & Pilot IQ compensation)

d = ceil(M .* rand(1, N)); % Random data symbols drawn from [1,2,..,M]
[s, ref] = mqam_modulator(M, d); % MQAM symbols & reference constellation

for i = 1:length(EsN0dB)
    r = add_awgn_noise(s, EsN0dB(i));
    z = receiver_impairments(r, g, phi, dc_i, dc_q); % Add impairments
    v = dc_compensation(z); % DC compensation
    y3 = blind_iq_compensation(v); % Blind IQ compensation
    [Kest, Pest] = pilot_iq_imb_est(g, phi, dc_i, dc_q); % Pilot based estimation
    y4 = iqImb_compensation(v, Kest, Pest); % IQ comp. using estimated values

    % IQ Detectors
    [tx1, d1] = iqOptDetector(z, ref); % No compensation
    [tx2, d2] = iqOptDetector(v, ref); % DC compensation only
    [tx3, d3] = iqOptDetector(y3, ref); % DC & blind IQ comp.
    [tx4, d4] = iqOptDetector(y4, ref); % DC & pilot IQ comp.

    % Symbol Error Rate Computation
    SER1(i) = sum(d ~= d1) / N; 
    SER2(i) = sum(d ~= d2) / N;
    SER3(i) = sum(d ~= d3) / N; 
    SER4(i) = sum(d ~= d4) / N;
end

theoreticalSER = ser_awgn(EbN0dB, 'MQAM', M); % Theoretical SER
figure; 
semilogy(EbN0dB, SER1, 'r*-'); 
hold on;
semilogy(EbN0dB, SER2, 'bO-'); 
semilogy(EbN0dB, SER3, 'g^-');
semilogy(EbN0dB, SER4, 'm*-'); 
semilogy(EbN0dB, theoreticalSER, 'k');

legend('No compensation', 'DC comp only', 'Sim- DC & blind iq comp', 'Sim- DC & pilot iq comp', 'Theoretical');
xlabel('E_b/N_0 (dB)'); ylabel('Symbol Error Rate (Ps)');
title('Probability of Symbol Error 64-QAM signals');

%%

clearvars; clc;
M = 64; % M-QAM modulation order
N = 1000; % To generate random symbols
d = ceil(M .* rand(N, 1)); % Random data symbol generation
s = mqam_modulator(M, d); % M-QAM modulated symbols (s)

g_1 = 0.8; phi_1 = 0;  dc_i_1 = 0;   dc_q_1 = 0; % Gain mismatch only
g_2 = 1;   phi_2 = 12; dc_i_2 = 0;   dc_q_2 = 0; % Phase mismatch only
g_3 = 1;   phi_3 = 0;  dc_i_3 = 0.5; dc_q_3 = 0.5; % DC offsets only
g_4 = 0.8; phi_4 = 12; dc_i_4 = 0.5; dc_q_4 = 0.5; % All impairments

r1 = receiver_impairments(s, g_1, phi_1, dc_i_1, dc_q_1);
r2 = receiver_impairments(s, g_2, phi_2, dc_i_2, dc_q_2);
r3 = receiver_impairments(s, g_3, phi_3, dc_i_3, dc_q_3);
r4 = receiver_impairments(s, g_4, phi_4, dc_i_4, dc_q_4);

subplot(2,2,1); plot(real(s), imag(s), 'b.'); hold on;
plot(real(r1), imag(r1), 'r.'); title('IQ Gain mismatch only')
subplot(2,2,2); plot(real(s), imag(s), 'b.'); hold on;
plot(real(r2), imag(r2), 'r.'); title('IQ Phase mismatch only')
subplot(2,2,3); plot(real(s), imag(s), 'b.'); hold on;
plot(real(r3), imag(r3), 'r.'); title('DC offsets only')
subplot(2,2,4); plot(real(s), imag(s), 'b.'); hold on;
plot(real(r4), imag(r4), 'r.'); title('IQ impairments & DC offsets');