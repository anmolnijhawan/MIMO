%% Initialization
N = 10^4;    % No. of symbols for which simulation has to run
Nt = 16;      % No. of transmitting antennas
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
       input_bits = randi([0,1],log2(Nt),1);          % Generate input bits for SM modulation
       ant_select = input_bits(1:log2(Nt));           % Extract bits for antenna selection
       ant_index = bi2de(ant_select')+1;              % Convert bits to antenna index
       x  = zeros(Nt,1);                              % Initialize all antennas
       x(ant_index) = 1;                              % Symbol transmiited from respective antenna
       y = H*x;                                        
       r = y+2*10^(-EbNo(itr)/20)*n;                    % Received signal with awgn noise
       
       max_t = 0.0;
       for j=1:1:Nt
          val = real((r-(H(:,j)/2))'*H(:,j));
          if val>max_t
              max_t = val;
              detect_ant_index = j;
          end    
       end
       if( mod(sym,1000)==0)
        ((itr-1)*N+sym)/(N*length(EbNo))*100   
       end
       %% Error calculation
       
       error_ant = (ant_index~=detect_ant_index);   % Error in antenna index
       error(1,itr) = error(1,itr) + error_ant; 
   end       
end

error_rate_SSK = error/N;
fig=semilogy(EbNo,error_rate_SSK,'b');
grid on;
hold on;
