%% 16-QAM modulation and Output 
M = 16;
x = (0:M-1)';
y = qammod(x,M);
scatterplot(y)


%% 64-QAM modulation using random numbers
M = 64;
x = randi([0 M-1],1000,1);
y = qammod(x,M);
avgPower = mean(abs(y).^2);

%% 64-QAM modulation using bits
M = 64;
k = log2(M);
data = randi([0 1],1000*k,1);
dataInMatrix = reshape(data,length(data)/k,k);
dataSymbolsIn = bi2de(dataInMatrix);
y = qammod(dataSymbolsIn,M);
EbNo = 10;
snr = EbNo + 10*log10(k);
receivedsignal = awgn(y,snr,'measured');
sPlotFig = scatterplot(receivedsignal,1,0,'g.');
hold on
scatterplot(y,1,0,'k*',sPlotFig);





