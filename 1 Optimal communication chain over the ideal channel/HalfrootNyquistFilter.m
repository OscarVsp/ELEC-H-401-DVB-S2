function [filter] = HalfrootNyquistFilter(fs,T_symb)

% INPUTS:
% - input_signal : vector of input signal 
%
% OUTPUTS:
% - output_signal : vector of ouput signal
%
%This function creates the Nyquist filter that will be used to maximize the
%SNR. The sampling frequency has to be the same as the symbols to convolve
%correctly with this filter.

%% filter parameters

beta = 0.3; %Makes the window smoother as beta increases // roll-off factor given in the specifications
T = T_symb; %Sampling period to avoid ISI given by the f_cut (slide 29 p211) not sure about this relation !
N = 201; % Number of filter samples may be given by the fs

f_step = fs/N;
f_max = f_step*(N-1)/2;
f = linspace(-f_max,f_max,N);

t_step = 1/fs;
t = (-(N-1)/2:(N-1)/2)*t_step;

%% filter in the frequency domain
%go with fft and calculated function

%via fft
% H = fft(h);
% tailledeH = length(H);
% figure(3);
% title("Root raised cosine window by fft");
% plot(f,abs(H));grid on;

%By the course calculation

H=zeros(1,N);
%creating the fucntion by the cases given in slide 30 p212
lower_bound= (1-beta)/(2*T);
upper_bound= (1+beta)/(2*T);

for i = 1:N
    fi = abs(f(i));
    if fi <= lower_bound
        H(i)=T;
    elseif fi <= upper_bound
        H(i) = T*( 1 + cos( pi*T*(fi - lower_bound)/beta))/2;
    end
end


%figure(4);plot(f,H,'*'); grid on;title("RC filter window");
%% filter in time domain
%Passing from freq to time domain with IFFT


G = sqrt(H); %Root of the filter to implement it at transmiter and receiver
%figure(6);plot(f,G,'*'); grid on;title("RRC filter window");

G=ifftshift(G);
%hold on; plot(f,G,'O')
g_norm= ifft(ifftshift(H));
g = ifft(G);
%figure(11);plot(t,g); grid on; title("Impusle response g NOT normalised nor shifted") ; 
g= fftshift(g/sqrt(g_norm(1))); %shift to center the sinc
%figure(15);plot(t,g); grid on; title("Impusle response g NOT normalised") ; 

% h = ifft(sqrt(H));
% figure(4); plot(t,real(h)); grid on; title("Impulse response of the raised cosine h(t)");
% h = fftshift(h/h(1));
% figure(12); plot(t,real(h)); grid on; title("Impulse response h normalised and shifted"); 
% figure(5); plot(t,real(g));grid on; title("Impulse response of the raised cosine g(t) normalized and fftshifted"); 


%% Output
filter = g; %Output the filter, the convolution is done on the main function so that we only compute de filter once at the start

end

