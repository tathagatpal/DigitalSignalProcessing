%CODE FOR DFT, IDFT and Radix-2 FFT 

n = 0:1:7; %FFT at 8 points
ts = 1/8000; %Sampling freq = 1/ts = 8000 Hz
n1 = 0:ts:7;

%x_n = sin(2*pi*1000*n*ts) + 0.5*sin(2*pi*2000*n*ts + 3*pi/4); %Ex in Lyon
x_n = sin(2*pi*1000*n*ts);



for i=1:7
    x_1 = sin(2*pi*i*1000*n1*ts);
    x_2 = 0.5*sin(2*pi*i*2000*n1*ts + 3*pi/4);
    subplot(2,4,i);
    plot(n1, x_1, n1, x_2, n, x_n, 'linewidth', 2);
    grid on;
    title("Value of m: ", i);
    legend('sin(2.pi.1000.n1.ts)', '0.5*sin(2.pi.2000.n1.ts + 3*pi/4)', ...
            'x_n');
end

%FFT
[y] = FFT(x_n);  
y

%DFT
[y1] = DFT(x_n, 8);
y1

%IDFT
[y2] = IDFT(y1, 8);
figure(2);
subplot(2,1,1);
plot(n, real(y2));
title('IDFT(y1), where y1 = FT(x_n)');
subplot(2,1,2);
plot(n, x_n);
title('Original signal x_n');

%Function for DFT
function [y1] = DFT(x,num)
x_m = zeros(1, num);
    for m = 1:num
        for n = 1:num
            x_m(m) = x_m(m) + x(n) .* exp(-1i*2*pi*(m-1)*(n-1)/num);
        end
     end
    
y1 = (x_m);
end

%Function for IDFT
function [y] = IDFT(x, num)
x_n = zeros(1, num);
    for m = 1:num
        for n = 1:num
            x_n(m) = x_n(m) + x(n) .* exp(1i*2*pi*(m-1)*(n-1)/num);
        end
        x_n(m) = x_n(m) * 1/num;
    end
y = x_n;
end

%Function for Radix-2 FFT
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
