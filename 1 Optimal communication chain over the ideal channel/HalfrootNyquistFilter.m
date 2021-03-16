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

H=zeros(1,N);
%creating the fucntion by the cases given in slide 30 p212
lower_bound = (1-beta)/(2*T);
upper_bound = (1+beta)/(2*T);

for i = 1:N
    fi = abs(f(i));
    if fi < lower_bound
        H(i)=T;
    elseif fi <= upper_bound
        H(i) = T*( 1 + cos( pi*T*(fi - lower_bound)/beta))/2;
    end
end



%% filter in time domain
%Passing from freq to time domain with IFFT

G = sqrt(H); %Root of the filter to implement it at transmiter and receiver

h = ifft(ifftshift(H));
g = ifft(ifftshift(G));
norm = h(1);
h = fftshift(h);
g = fftshift(g);
h_normed = (h/norm);
g_normed = (g/sqrt(norm));

% figure(1);plot(f,H,'*'); grid on;title("Frequency domaine RC filter");
% figure(2);plot(f,G,'*'); grid on;title("Frequency domaine RRC filter");hold on;
% figure(3);plot(t,h); grid on; title("Time domain RC filter");
% figure(4);plot(t,h_normed); grid on; title("Time domaine normed RC filter")
% figure(5);plot(t,g); grid on; title("Time domain RRC filter");
% figure(6);plot(t,g_normed); grid on; title("Time domaine normed RRC filter")

%% Output
filter = g_normed; %Output the filter, the convolution is done on the main function so that we only compute de filter once at the start

end

