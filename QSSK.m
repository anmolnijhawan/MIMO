%% Initialization
N = 10^3;    % No. of symbols for which simulation has to run
Nt = 4;      % No. of transmitting antennas
Nr = 4;      % No. of receiving antennas
N_ant = log2(Nt); % No. of bits transmiited as antenna index

EbNo = 0:2:18;  %Range of SNR for simulation

error = zeros(1,length(EbNo));   % Initialization of total error
error_ant = zeros(1,length(EbNo));  % Initialization of antenna error
for itr = 1:1:length(EbNo)  % Iterate for each value of SNR
   for sym = 1:N            % Iterate for each symbol
       H = 1/sqrt(2)*(randn(Nr,Nt)+1i*(randn(Nr,Nt))); % Generate rayleigh fading channel
       n = 1/sqrt(2)*(randn(Nr,1)+1i*(randn(Nr,1)));  % Generate awgn noise
       
       %% Transmission
       input_bits = randi([0,1],2*log2(Nt),1);          % Generate input bits for SM modulation
       ant_select_r = input_bits(1:log2(Nt));           % Extract bits for antenna selection
       ant_select_i = input_bits(log2(Nt)+1:end);
       ant_index_r = bi2de(ant_select_r')+1;              % Convert bits to antenna index
       ant_index_i = bi2de(ant_select_i')+1;
       x  = zeros(Nt,1);                              % Initialize all antennas
       x(ant_index_r) = 1;                              % Symbol transmiited from respective antenna
       x(ant_index_i) = 1*1i;
       y = H*x;                                        
       r = y+10^(-EbNo(itr)/20)*n;                    % Received signal with awgn noise
       
       min_t = 10000.0;
       for j = 1:1:Nt
           for k = 1:1:Nt
               val = norm(r-H(:,j)-(H(:,k)*1i))^2;
               if val < min_t
                  min_t = val;
                  detect_antenna_r = j;
                  detect_antenna_i = k;
               end
           end
       end
       if( mod(sym,1000)==0)
        ((itr-1)*N+sym)/(N*length(EbNo))*100   
       end
       error_ant_r = (ant_index_r~=detect_antenna_r);   % Error in antenna index
       error_ant_i = (ant_index_i~=detect_antenna_i); 
       error(1,itr) = error(1,itr) + ((error_ant_r+error_ant_i)~=0);   % +error_sym
       
   end
end

error_rate_QSSK = error/N;
fig=semilogy(EbNo,error_rate_QSSK,'y');
grid on;
hold on;
