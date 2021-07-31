close all;clear;clc;
 	
fft_size = 64; 
nDSC = 52; % number of data subcarriers
nbOFDMsymb = 52; %number of bit per OFDM symbols
nSymb = 10000; % number of symbols
 
data=randi([0 1],1,nbOFDMsymb*nSymb);
%modulation
data_mod = 2*data-1;
data_mod = reshape(data_mod,nbOFDMsymb,nSymb)'; 
 
Eb_N0_dB = 0:10;
Es_N0_dB = Eb_N0_dB + 10*log10(nDSC/fft_size) + 10*log10(64/80); 
error=zeros(1,11);
for ii = 1:length(Eb_N0_dB)
 
   tx_F = [zeros(nSymb,6) data_mod(:,(1:nbOFDMsymb/2)) zeros(nSymb,1) data_mod(:,(nbOFDMsymb/2+1:nbOFDMsymb)) zeros(nSymb,5)] ;
   
   tx_ifft = (fft_size/sqrt(nDSC))*ifft(fftshift(tx_F')).';
   tx_cyclic = [tx_ifft(:,(((3/4)*fft_size)+1:fft_size)) tx_ifft];
   tx = reshape(tx_cyclic.',1,nSymb*size(tx_cyclic,2));
   %noise
   noise = 1/sqrt(2)*(randn(1,nSymb*size(tx_cyclic,2))+1i.*randn(1,nSymb*size(tx_cyclic,2)));
 
   % AWGN
   rx = sqrt(80/64)*tx + 10^(-Es_N0_dB(ii)/20)*noise;
 
   % Receiver
   rx = reshape(rx.',size(tx_cyclic,2),nSymb).'; % formatting the received vector into symbols
   rx = rx(:,(((1/4)*fft_size)+1:size(tx_cyclic,2))); % removing cyclic prefix
 
   rx_fft = (sqrt(nDSC)/fft_size)*fftshift(fft(rx.')).'; 
   rx_dem = rx_fft(:,[6+(1:nbOFDMsymb/2) 7+(nbOFDMsymb/2+1:nbOFDMsymb)]); 
 
   % BPSK demodulation
   rx_dem = 2*floor(real(rx_dem/2)) + 1;
   for jj=1:nSymb
       for kk=1:nDSC
           if(rx_dem(jj,kk)>1)
             rx_dem(jj,kk)=1;
           elseif(rx_dem(jj,kk)<-1)
              rx_dem(jj,kk)=-1; 
           end
       end
   end
   
   % converting modulated values into bits
   data_dem = (rx_dem+1)/2;
   data_dem = reshape(data_dem.',nbOFDMsymb*nSymb,1).';
 
   % counting the errors
   error(ii) = size(find(data_dem - data),2);
 
end
 
BER_simulation = error/(nSymb*nbOFDMsymb);
BER_theorical = qfunc(sqrt(2.*(10.^(Eb_N0_dB/10))));
 
figure
semilogy(Eb_N0_dB,BER_theorical,Eb_N0_dB,BER_simulation,'r*');
grid on
legend('BER Theory', 'BER Simulation');
xlabel('Eb/No, dB');ylabel('BER');title('BER for BPSK using OFDM')
