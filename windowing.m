ts = 1/50; %Sampling time
t = -10 : ts : 10; 
x =  sin(2*pi*15*t); %A sine wave with f = 15 Hz

hamming_w = 0.54 - 0.46*cos(2*pi*t/(length(x)-1)); %Hamming Window

%Element-wise multiplication of signal x with Hamming's window
x_window_ham = x.*hamming_w; 

hanning_w = 0.5 + 0.5*cos(pi*t/length(x)-1); %Hanning window

%Element-wise multiplication of signal x with Hanning's window
x_window_han = x.*hanning_w;

fs = 1/ts; %Frequency samples


%32 point DFT of x
y_DFT_32 = DFT(x, 32);
f_dft_32 = (0:length(y_DFT_32)-1)*fs/length(y_DFT_32);

%64 point DFT of x
y_DFT_64 = DFT(x, 64);
f_dft_64 = (0:length(y_DFT_64)-1)*fs/length(y_DFT_64);

%FFT of x
y_FFT = FFT(x);
f_fft = (0:length(y_FFT)-1)*fs/length(y_FFT);

%FFT of Hamming windowed function
y_window_ham = FFT(x_window_ham);
f_window_ham = (0:length(y_window_ham)-1)*fs/length(y_window_ham);

%FFT of Hanning windowed function
y_window_han = FFT(x_window_han);
f_window_han = (0:length(y_window_han)-1)*fs/length(y_window_han);

%Plotting 32 point DFT of x
figure(1);
subplot(3,1,1);
plot(f_dft_32, y_DFT_32);
title('32-point DFT of x');

%Plotting 64 point DFT of x
subplot(3,1,2);
plot(f_dft_64, y_DFT_64);
title('64-point DFT of x');

%Plotting FFT of x
subplot(3,1,3);
plot(f_fft, y_FFT);
title('FFT of x (1024-point)');

figure(2);
%subplot(2,2,1);
plot(f_fft, y_FFT, f_window_ham, y_window_ham, f_window_han, y_window_han ...
    , 'linewidth', 1.5);
legend('FFT of x', 'FFT of windowed x (Hamming)', 'FFT of windowed x (Hanning)'); 
title('FFT of x, FFT of windowed x (Hamming), FFT of windowed x (Hanning)'); 

figure(3);
%subplot(2,2,2);
plot(f_fft, y_FFT, f_window_ham, y_window_ham);
legend('FFT of x', 'FFT of windowed x (Hamming)'); 
title('FFT of x, FFT of windowed x (Hamming)'); 

figure(4);
%subplot(2,2,3);
plot(f_fft, y_FFT, f_window_han, y_window_han);
legend('FFT of x', 'FFT of windowed x (Hanning)'); 
title('FFT of x, FFT of windowed x (Hanning)');

figure(5);
%subplot(2,2,4);
plot(f_window_ham, y_window_ham, f_window_han, y_window_han);
legend('FFT of windowed x (Hamming)', 'FFT of windowed x (Hanning)'); 
title('FFT of windowed x (Hamming), FFT of windowed x (Hanning)');



function [y1] = DFT(x,num)
x_m = zeros(1, num);
    for m = 1:num
        for n = 1:num
            x_m(m) = x_m(m) + x(n) .* exp(-1i*2*pi*(m-1)*(n-1)/num);
        end
     end
    
y1 = (x_m);
end

function [y] = FFT(x) 
p=nextpow2(length(x));            % checking the size of the input array
x=[x zeros(1,(2^p)-length(x))];   % making a zero array 
N=length(x);                      % computing the array size
S=log2(N);                        % computing the #conversion stages
Half=N/2;                         % half the length of the array
for stage=1:S                     % stages of transformation
    % series of "butterflies" for each stage
    for index=0:(N/(2^(stage-1))):(N-1) 
        for n=0:(Half-1)    % creating "butterfly" and saving the results
            pos=n+index+1;  % index of the data sample
            pow=(2^(stage-1))*n;  % part of power of the complex multiplier
            w=exp((-1i)*(2*pi)*pow/N); % complex multiplier
            a=x(pos)+x(pos+Half); % 1st part of the "butterfly" 
            b=(x(pos)-x(pos+Half)).*w; % 2nd part of the "butterfly" 
            x(pos)=a;               % saving computation of the 1-st part
            x(pos+Half)=b;          % saving computation of the 2-nd part
        end
    end
Half=Half/2;        % computing the next "Half" value
end
y=bitrevorder(x);   % performing bit-reverse operation and return 
end


