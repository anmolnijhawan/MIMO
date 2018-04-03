clc;
clear all;
%% Initialization
N = 10^6;    % No. of symbols for which simulation has to run
Nt = 2;      % No. of transmitting antennas
Nr = 4;      % No. of receiving antennas
M = 4;       % Modulation order for signal constellation
lut = qammod([0:M-1],M); % Look up table for QAM modulation
Eavg = sum(abs(lut).^2)/M; % Avg. energy of constellation
N_lut = lut/sqrt(Eavg);  % Normalized energy of constellation
EbNo = 0:2:16;  %Range of SNR for simulation

error = zeros(1,length(EbNo));   % Initialization of total error
for itr = 1:1:length(EbNo)  % Iterate for each value of SNR
   for sym = 1:N            % Iterate for each symbol
       H = 1/sqrt(2)*(randn(Nr,Nt)+1i*(randn(Nr,Nt))); % Generate rayleigh fading channel
       n = 1/sqrt(2)*(randn(Nr,1)+1i*(randn(Nr,1)));  % Generate awgn noise
       
       %% Transmission
       input_bits = randi([0,1],log2(Nt)+log2(M)+log2(Nt),1);  % Generate input bits for SM modulation
       r_ant_select = input_bits(log2(M)+1:log2(Nt)+log2(M));           % Extract bits for antenna selection for real part
       i_ant_select = input_bits(log2(Nt)+log2(M)+1:end);           % Extract bits for antenna selection for imaginary part
       sym_select = input_bits(1:log2(M));       % Extract bits for choosing symbol from constellation
       r_ant_index = bi2de(r_ant_select')+1;              % Convert bits to antenna index
       i_ant_index = bi2de(i_ant_select')+1;              % Convert bits to antenna index
       sym_index = bi2de(sym_select');                % Convert bits to symbol index
       x_i = qammod(sym_index,M)/sqrt(Eavg);          % QAM modulation of input bits
       x  = zeros(Nt,1);                              % Initialize all antennas
       x(r_ant_index) = real(x_i);                            % Symbol transmiited from respective antenna
       x(i_ant_index) = imag(x_i)*i;
       y = H*x;                                        
       r = y+10^(-EbNo(itr)/20)*n;                    % Received signal with awgn noise
       
       
       %% Detection Technique 2
       min = 100000.0;
       for m=0:1:M-1
           for l_r=1:1:Nt
               for l_i=1:1:Nt
                   x_  = qammod(m,M)/sqrt(Eavg);
                   x_r = real(x_);
                   x_i = imag(x_);
                   g = H(:,l_r)*x_r + H(:,l_i)*x_i*1i; 
                   val = norm(g) - 2*real(r'*g);
                   if val<min
                       min = val;
                       detect_ant_index_r = l_r;
                       detect_ant_index_i = l_i;
                       detected_sym_index = m;
                       
                   end
                end
          end
       end
       
       %% Error calculation
       
       error_ant_r = (r_ant_index~=detect_ant_index_r);   % Error in antenna index
       error_ant_i = (i_ant_index~=detect_ant_index_i); 
       error_sym = (detected_sym_index~=sym_index); % Error in symbol index
       error(1,itr) = error(1,itr) + ((error_ant_r+error_ant_i+error_sym)~=0);   % +error_sym
       
   end
end

error_rate = error/N;
semilogy(EbNo,error_rate);
grid on;
hold on;
