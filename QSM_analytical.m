%clc;
%clear all;
%% Initialization
%N = 10^4;    % No. of symbols for which simulation has to run
Nt = 2;      % No. of transmitting antennas
Nr = 4;      % No. of receiving antennas
M = 4;       % Modulation order for signal constellation
lut = qammod([0:M-1],M); % Look up table for QAM modulation
Eavg = sum(abs(lut).^2)/M; % Avg. energy of constellation
N_lut = lut/sqrt(Eavg);  % Normalized energy of constellation
EbNo = 0:0.5:20;  %Range of SNR for simulation
ABEP=zeros(1,length(EbNo));
asymptotic_BEP=zeros(1,length(EbNo));


       %% Analytical performance
     for itr = 1:1:length(EbNo)
       ABEP_sum=0;
       asymptotic_sum=0;
       for ant_index_r_1=0:1:Nt-1
         for ant_index_i_1=0:1:Nt-1
           for ant_index_r_2=0:1:Nt-1
             for ant_index_i_2=0:1:Nt-1
               for sym_index_1=0:1:M-1
                 x_1  = qammod(sym_index_1,M);%/sqrt(Eavg);
                 x_r_1 = real(x_1);
                 x_i_1 = imag(x_1);
                 for sym_index_2=0:1:M-1
                   x_2  = qammod(sym_index_2,M);%/sqrt(Eavg);
                   x_r_2 = real(x_2);
                   x_i_2 = imag(x_2);
                   if((ant_index_r_1~=ant_index_r_2) && (ant_index_i_1~=ant_index_i_2))
                     exp_mean=10^(EbNo(itr)/10)/2*(abs(x_r_1)^2+abs(x_i_1)^2+abs(x_r_2)^2+abs(x_i_2)^2);   %variance of channel gains =1
                   elseif((ant_index_r_1==ant_index_r_2) && (ant_index_i_1~=ant_index_i_2))
                     exp_mean=10^(EbNo(itr)/10)/2*(abs(x_r_1-x_r_2)^2+abs(x_i_1)^2+abs(x_i_2)^2);
                   elseif((ant_index_r_1~=ant_index_r_2) && (ant_index_i_1==ant_index_i_2))
                     exp_mean=10^(EbNo(itr)/10)/2*(abs(x_r_1)^2+abs(x_i_1-x_i_2)^2+abs(x_r_2)^2);
                   else
                     exp_mean=10^(EbNo(itr)/10)/2*(abs(x_r_1-x_r_2)^2+abs(x_i_1-x_i_2)^2);  
                   end
                   gam=0.5*(1-sqrt(0.5*exp_mean/(1+0.5*exp_mean)));
                   sum_exp=0;
                   for k=0:1:Nr-1
                     sum_exp=sum_exp+nchoosek(Nr-1+k,k)*((1-gam)^k);
                   end  
                   asymptotic_PEP=2^(Nr-1)*gamma(Nr+0.5)*((1/(exp_mean+0.1))^Nr)/sqrt(3.14*factorial(Nr));
                   Avg_PEP=gam^Nr*sum_exp;
                   bit_errors=sum(xor([de2bi(ant_index_r_1,log2(Nt)) de2bi(ant_index_i_1,log2(Nt)) de2bi(sym_index_1,log2(M))],[de2bi(ant_index_r_2,log2(Nt)) de2bi(ant_index_i_2,log2(Nt)) de2bi(sym_index_2,log2(M))])); %hamming distance between the codewords     
                   ABEP_sum=ABEP_sum+Avg_PEP*bit_errors;                
                   asymptotic_sum=asymptotic_sum+asymptotic_PEP*bit_errors;
                 end
               end
             end
           end
         end
       end
       ABEP(itr)=ABEP_sum/(M*Nt*Nt*(log2(M)+2*log2(Nt))); % average bit error probability
       asymptotic_BEP(itr)=asymptotic_sum/(M*Nt*Nt*(log2(M)+2*log2(Nt)));
     end   %% Error calculation 
   %end

semilogy(EbNo,ABEP,'b');

grid on;
hold on;
semilogy(EbNo,asymptotic_BEP,'v');