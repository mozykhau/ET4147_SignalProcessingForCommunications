%%  Exam Assignment
%   ET4147: Signal Processing for Communications

%   Author: M.Kyryliuk
%   Email:  m.s.kyryliuk@student.tudelft.nl

%% Pre-work
clc;
clear all
close;



%% Instanteneous Model

%% Initialisation
M = 5;  %   number of elements in array
N = 20; %   number of samples
SNR = 100; % SNR [dB]
delta = 0.5;    %   spacing of elements in wavelengths
theta = [-20;30]; %	source direction in degrees
f = [0.1;0.3]; % normalised frequency
d = length(theta);
theta_rad = deg2rad(theta); %   source direction in radians

%% Array response
a_resp = gen_a(M,delta,theta_rad);


%% Beamformer
%w2 = gen_a(M,delta,deg2rad(theta));
w = ones(M,1);
theta_range = -90:90;
theta_range = deg2rad(theta_range);
%a_vector(:) = gen_a(M,delta,theta_range(:))
% l = length(theta_range);
% for i=1:l
%     a(:,i) = gen_a(M,delta,theta_range(i));
% end
%     y = abs(w'*a);

%check = spat_response(w,delta,theta_range,M);

%% Make data


X = gendata(M,N,delta,theta_rad,f,SNR);

%% Plot singular values
figure(1)
% H0
M = 5;  %   number of elements in array
N = 20; %   number of samples
SNR = 20; % SNR [dB]
delta = 0.5;    %   spacing of elements in wavelengths
theta = [-20;30]; %	source direction in degrees
f = [0.1;0.3]; % normalised frequency
theta_rad = deg2rad(theta); %   source direction in radians

X_H0 = gendata(M,N,delta,theta_rad,f,SNR);

% different number of antennas
M2 = M*2; %number of antennas
X_Antennas = gendata(M2,N,delta,theta_rad,f,SNR);

% different frequency
f2 = [0.1;0.3]; % normalised frequency

X_freq = gendata(M,N,delta,theta_rad,f2,SNR);

% different DOA
theta2 = [20;30]; %different DOA separation
theta_rad2 = deg2rad(theta2);
X_DOA = gendata(M,N,delta,theta_rad,f,SNR);

% number of samples
N2 = N*2;
X_Samples = gendata(M,N,delta,theta_rad,f,SNR);

% Singular value decomposition
singular_H0 = svd(X_H0);
singular_Antennas = svd(X_Antennas);
singular_DOA = svd(X_DOA);
singular_Samples = svd(X_Samples);
singular_freq = svd(X_freq);

%% plot task 1.3
subplot(1,2,1)
plot(singular_H0,'ko','MarkerSize',15);
hold on
plot(singular_Antennas,'b+','MarkerSize',15)
plot(singular_Samples,'rx','MarkerSize',15)
title({'Singular values of X.', 'Comparison of different number of antennas and samples'});
legend('M = 5; SNR = 20; N = 20; alpha = [-20 30]; f = [0.1 0.3]',...
    'M = 10',...
    'N = 40');
ylim([0 150]);
subplot(1,2,2)
plot(singular_H0,'ko','MarkerSize',15);
hold on
plot(singular_freq,'.','MarkerSize',15)
plot(singular_DOA,'b*','MarkerSize',15)
title({'Singular values of X.','Comparison of different normalised frequency and DOA'});
legend('M = 5; SNR = 20; N = 20; alpha = [-20 30]; f = [0.1 0.3]',...
    'f = [0.25, 0.3]',...
    'theta = [20, 30]');
ylim([0 150]);

%% MUSIC DOA
close all
DOA = music(X,d,delta)

%% MUSIC Frequency
Freqs = musicfreq(X,d)
close all
clear all
%% Comparison Performance (mean and std)

M = 3;  %   number of elements in array
N = 20; %   number of samples
SNR = [0,4,8,12,16,20]; % SNR [dB]
delta = 0.5;    %   spacing of elements in wavelengths
theta = [-20;30]; %	source direction in degrees
f = [0.1;0.12]; % normalised frequency
theta_rad = deg2rad(theta); %   source direction in radians
d = length(theta);

len = length(SNR);
num_tests = 1;

DOA = zeros(len,num_tests,d);
freqs = zeros(len,num_tests,d);
statDOA = struct('snr',SNR,'mean1',zeros(1,len),'std1',zeros(1,len),'mean2',zeros(1,len),'std2',zeros(1,len));
statFreq = struct('snr',SNR,'mean1',zeros(1,len),'std1',zeros(1,len),'mean2',zeros(1,len),'std2',zeros(1,len));
for i = 1:len
    for k = 1:num_tests
        X = gendata(M,N,delta,theta_rad,f,SNR(i));
        DOA(i,k,:) = music(X,d,delta);
        freqs(i,k,:) = musicfreq(X,d);
    end
    statDOA.mean1(i) = mean(DOA(i,:,1));
    statDOA.mean2(i) = mean(DOA(i,:,2));
    statDOA.std1(i) = std(DOA(i,:,1),1);
    statDOA.std2(i) = std(DOA(i,:,2),1);
    statFreq.mean1(i) = mean(freqs(i,:,1));
    statFreq.mean2(i) = mean(freqs(i,:,2));
    statFreq.std1(i) = std(freqs(i,:,1),1);
    statFreq.std2(i) = std(freqs(i,:,2),1);
end
name = fieldnames(statDOA);
figure()
subplot(1,2,1)
plot(SNR,statDOA.mean1)
hold on
plot(SNR, statDOA.std1)
plot(SNR, statDOA.mean2)
plot(SNR, statDOA.std2)
title('Performance of Music DOA estimation vs SNR')
xlabel('SNR [dB]')
ylabel('DOA')
legend(name(2:5))

subplot(1,2,2)
plot(SNR,statFreq.mean1)
hold on
plot(SNR, statFreq.std1)
plot(SNR, statFreq.mean2)
plot(SNR, statFreq.std2)
title('Performance of Music Spectrum Estimation vs SNR')
xlabel('SNR [dB]')
ylabel('Normalised frequency x 2\pi')
legend(name(2:5))

close all
clear all

%% Zero Forcing Beamformer Angle

% Inititalise
M = 3;  %   number of elements in array
N = 20; %   number of samples
SNR = 10;   % SNR [dB]
delta = 0.5;    %   spacing of elements in wavelengths
theta = [-20;30]; %	source direction in degrees
f = [0.1;0.12]; % normalised frequency
theta_rad = deg2rad(theta); %   source direction in radians
d = length(theta);
theta_range = -90:90;
theta_range = deg2rad(theta_range);
% Generate data
[X,A,S] = gendata(M,N,delta,theta_rad,f,SNR);
DOA = music(X,d,delta);
doa_rad = deg2rad(DOA);
% Generate beamformer
for i=1:length(doa_rad)
    beamformer(:,i) = gen_a(M,delta,doa_rad(i));
end
beamformer = sum(beamformer,2);
beamformer = beamformer/norm(beamformer);
newS = beamformer'*X;
figure()
subplot(1,2,1)
spat_response(beamformer,delta,theta_range,M);
title('Spatial response of ZF built on MUSIC DOA');
%% Zero Forcing beamformer based on frequency
Freqs = musicfreq(X,d);
% Generate data signal based on frequency
signal_freq = complex_sin_signal(Freqs,0,N);
A_freq = X*signal_freq';


% Try
beamformer_freq = 1./A_freq;
beamformer_freq = sum(beamformer_freq,2);
beamformer_freq = beamformer_freq/norm(beamformer_freq);

%Working version
% beamformer_freq = pinv(A_freq).';
% beamformer_freq = sum(beamformer_freq,2);
% beamformer_freq = beamformer_freq/norm(beamformer_freq);

% Plot spatial response
subplot(1,2,2)
spat_response(beamformer_freq,delta,theta_range,M);
title('Spatial response of ZF built on MUSIC Spectrum');

%% CMA algorithm
clear all
close all

% Inititalise
M = 3;  %   number of elements in array
N = 5000; %   number of samples
SNR = 10;   % SNR [dB]
delta = 0.5;    %   spacing of elements in wavelengths
theta = [-20;30]; %	source direction in degrees
f = [0.1;0.12]; % normalised frequency
theta_rad = deg2rad(theta); %   source direction in radians
d = length(theta);
theta_range = -90:90;
theta_range = deg2rad(theta_range);

X = gendata(M,N,delta,theta_rad,f,SNR);

mu = 0.001;
figure()
for i = 1:3
    w_init = rand(M,1)
    [w,y] = cma_try(X,mu,w_init);
    beamformer = w(:,N)
    spat_response(beamformer,delta,theta_range,M);
    hold on
    title('Spatial response of beamformer CMA');
end


clear all
close all
%% Channel equalization

%% Signal model
N = 500;
P = 4;
sigma = 0.5;

% Generate signal
source = sourceqpsk(N);

% Generate data
x = gendata_conv(source,P,N,sigma);

% Data matrix

for pl = 1:2*P
    for nl = 1:N-1
        X(pl,nl) = x(pl+P*(nl-1));
    end
end

rank = rank(X);
%% Zero-forcing and Wiener Receiver
% create matrix H
L = 1;
m = 2;
Hm = zeros(P,m+1);
center_tap = (L+m+1)/2;
delay = center_tap;
for i = 1:m
    for k = 1:L
        H_temp(:,k) = channel_exam(L,P,0);
    end
    Hm((i-1)*P+1:i*P,i*(1:L)) = H_temp((1:P),:)
end

% Zero Forcing
H_inverse = pinv(Hm);
H_inverse2 = H_inverse(delay,:);
newS_ZF = H_inverse2*X;

% Wiener
H = Hm;
I = eye(length(H));
wiener = inv(H*H'+I*sigma)*H;

Wiener_f = wiener(:,delay);
%Wiener_f = wiener;
newS_w = Wiener_f'*X;
%Plot
figure()
subplot(1,2,1)
scatter(real(newS_w),imag(newS_w))
title('Wiener')
xlabel('Real')
ylabel('Imaginary')
xlim([-1.5,1.5]);
ylim([-1.5,1.5]);
subplot(1,2,2)
scatter(real(newS_ZF.'),imag(newS_ZF.'))
title('Zero Forcing')
xlabel('Real')
ylabel('Imaginary')
xlim([-1.5,1.5]);
ylim([-1.5,1.5]);
hold off
close all

%% Channel estimation
h_est = channel_estimator(X,source,L);
h_est = h_est(:,1);
Hm = Hm(:,2);
figure()
subplot(1,2,1)
plot(real(Hm))
hold on
plot(real(h_est))
legend('True','Blind')
title('Real \delta = 1')
xlabel('Samples')
ylabel('Magnitude')
ylim([-1,1])
subplot(1,2,2)
plot(imag(Hm))
hold on
plot(imag(h_est))
title('Imag \delta = 1')
legend('True','Blind')
xlabel('Samples')
ylabel('Magnitude')
ylim([-0.1,0.1])
hold off

%% Blind Spatial filtering (Channel)
h = blind_channel(X);
%Hm = Hm(:,1:2)

% figure()
% subplot(1,2,1)
% plot(real(Hm))
% hold on
% plot(real(h))
% legend('true','blind')
% title('Real')
% subplot(1,2,2)
% plot(imag(Hm))
% hold on
% plot(imag(h))
% title('Imag')
% legend('true','blind')

%% Blind Spatial filtering (Signal)
s = blind_symbol(X);
% figure()
% scatter(real(s),imag(s))
% xlabel('Real')
% ylabel('Imaginary')
% title('Blind Signal Estimation')

