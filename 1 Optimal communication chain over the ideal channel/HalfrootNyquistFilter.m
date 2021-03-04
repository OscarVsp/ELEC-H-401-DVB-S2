function [output_signal] = HalfrootNyquistFilter(input_signal)

% INPUTS:
% - input_signal : vector of input signal 
%
% OUTPUTS:
% - output_signal : vector of ouput signal


%% filter parameters

beta = 0.3; %Makes the window smoother as beta increases // roll-off factor
f_cut = 1e6; %Hz Cutoff frequency 
T = 1/(2*f_cut); %Sampling period to avoid ISI given by the f_cut (slide 29 p211) not sure about this relation !
N = 10001; % Number of filter samples may be given by the fs
fs = 10*f_cut; % Sampling frequency (rule of thumb for the 10 25 times f_cut)

f_step = fs/N;
f_max = f_step*(N-1)/2;
f = linspace(-f_max,f_max,N);

t_step = 1/fs;
t = (-(N-1)/2:(N-1)/2)*t_step;
%% filter in time domain
%Formula given at slide 30
% h = ( sinc(pi*t/T).*(cos(pi*beta*t/T)) )./( 1 - 4*(beta^2)*(t.^2)/(T^2) ); %Watch out beacause alpha at denom but is it a "coquille 
% figure(2);
% title("Impulse response of the raised cosine");
% 
% plot(t,h);grid on;
% hold on;
% plot(t,h,'.');
%% filter in the frequency domain
%go with fft and calculated function

%via fft
% H = fft(h);
% tailledeH = length(H);
% figure(3);
% title("Root raised cosine window by fft");
% plot(f,abs(H));grid on;

%By the course calculation

H=zeros(length(f));
%creating the fucntion by the cases given in slide 30 p212
lower_bound= (1-beta)/(2*T);
upper_bound= (1+beta)/(2*T);
tic
for i = 1:length(f)
    fi = abs(f(i));
    if fi < lower_bound
        H(i)=T;
    elseif fi <= upper_bound
        H(i) = (T/2)*( 1 + cos( (pi*T/beta)*fi - lower_bound));
    end
end
toc
figure(4);

plot(f,H); grid on;title("RRC filter window");
%% Output
output_signal = input_signal; %Do the filter

end

