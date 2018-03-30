clc;
clear all;
%%
N = 10^5;    % No. of symbols for which simulation has to run
Nt = 4;      % No. of transmitting antennas
Nr = 4;      % No. of receiving antennas
M = 4;       % Modulation order for signal constellation
N_ant = log2(Nt); % No. of bits transmiited as antenna index
lut = qammod([0:M-1],M); % Look up table for QAM modulation
Eavg = sum(abs(lut).^2)/M; % Avg. energy of constellation
N_lut = lut/sqrt(Eavg);  % Normalized energy of constellation
EbNo = 0:2:16;  %Range of SNR for simulation

error = zeros(1,length(EbNo));   % Initialization of total error
error_ant = zeros(1,length(EbNo));  % Initialization of antenna error
error_sym = zeros(1,length(EbNo));  % Initialization of symbol error

for itr = 1:1:length(EbNo)  % Iterate for each value of SNR
   for sym = 1:N            % Iterate for each symbol
       H = 1/sqrt(2)*(randn(Nr,Nt)+1i*(randn(Nr,Nt))); % Generate rayleigh fading channel
       n = 1/sqrt(2)*(randn(Nr,1)+1i*(randn(Nr,1)));  % Generate awgn noise
       
       %% Transmission
       input_bits = randi([0,1],log2(Nt)+log2(M),1);  % Generate input bits for SM modulation
       ant_select = input_bits(1:log2(Nt));           % Extract bits for antenna selection
       sym_select = input_bits(log2(Nt)+1:end);       % Extract bits for choosing symbol from constellation
       ant_index = bi2de(ant_select')+1;              % Convert bits to antenna index
       sym_index = bi2de(sym_select');                % Convert bits to symbol index
       x_i = qammod(sym_index,M)/sqrt(Eavg);          % QAM modulation of input bits
       x  = zeros(Nt,1);                              % Initialize all antennas
       x(ant_index) = x_i;                            % Symbol transmiited from respective antenna
       y = H*x;                                        
       r = y+10^(-EbNo(itr)/10)*n;                    % Received signal with awgn noise
       
       
       %% Detection Technique 2
       max = 0.0;
       for i=1:1:Nt                            % Iterate through all antennas
          val = abs(H(:,i)'*r)/norm(H(:,i));   % Compute the value of expression given in slides.
          if val > max                           % finding argument maximum.
              max = val;
              detect_ant_index = i;              % detected antenna index
          end
       end
       
       min = 1000000.0;                         % Initialization of minimum
       for m=1:1:M                              % Iterate through all symbols
           val = norm(H(:,m)*x_i);              % Compute the value of expression given in slides
           val = val - 2*real(H(:,m)'*y*x_i');
           if val < min                         % Computing argument minimum
              min = val;                        
              detected_sym_index = m;           % detected symbol index
           end
       end
       
   end
end


%% Detection Technique 1
       z = zeros(Nr,Nt,M);
       d = zeros(Nr,Nt,M);
       for m=1:M
          z(:,:,m) = H*N_lut(m);
          d(:,:,m) = repmat(r,1,Nt)-z(:,:,m);
          
       end
       
