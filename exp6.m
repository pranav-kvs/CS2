
nSym = 10^5; % Number of symbols to transmit
EbN0dB = -4:2:14; % Eb/N0 range in dB for simulation
MOD_TYPE = 'PSK'; % Set 'PSK' or 'QAM' or 'PAM' or 'FSK'
arrayOfM = [2, 4, 8, 16, 32]; % Array of M values to simulate
%arrayOfM=[4,16,64,256]; % Uncomment this line if MOD_TYPE='QAM'
COHERENCE = 'coherent'; % 'coherent'/'noncoherent'-only for FSK

plotColor = ['b','g','r','c','m','k']; % Plot colors
p = 1; 
%legendString = cell(1, length(arrayOfM) * 2); % For legend entries

for M = arrayOfM
    
    k = log2(M); EsN0dB = EbN0dB + 10*log10(k); % EsN0dB calculation
    SER_sim = zeros(1, length(EbN0dB)); % Simulated Symbol error rates
    d = ceil(M .* rand(1, nSym)); % Uniform random symbols from 1:M
    s = modulate(MOD_TYPE, M, d, COHERENCE); % (Refer Chapter 3)
    
    for i = 1:length(EsN0dB)
        r = add_awgn_noise(s, EsN0dB(i)); % Add AWGN noise
        dCap = demodulate(MOD_TYPE, M, r, COHERENCE); % (Refer Chapter 3)
        SER_sim(i) = sum(d ~= dCap) / nSym; % SER computation
    end
    
    SER_theory = ser_awgn(EbN0dB, MOD_TYPE, M, COHERENCE); % Theory SER
    
    semilogy(EbN0dB, SER_sim, [plotColor(p) '*']); 
    semilogy(EbN0dB, SER_theory, plotColor(p));
    %legendString{2*p-1} = ['Sim ', num2str(M), '-', MOD_TYPE];
    %legendString{2*p} = ['Theory ', num2str(M), '-', MOD_TYPE]; 
    p = p + 1;
    hold on;
end

% Set axis scaling and limits for consistent scaling
%xlim([min(EbN0dB) max(EbN0dB)]);
%ylim([1e-6 1]); % Adjust according to expected SER range

%legend(legendString);
legend('M=2','sim 2','M=4','sim 4','M=8','sim 8','M=16','sim 16','M=32','sim 32');
xlabel('Eb/N0 (dB)'); ylabel('SER (Ps)');
title(['Probability of Symbol Error for M-', MOD_TYPE, ' over AWGN']);
grid on; % Add grid for better visualization
