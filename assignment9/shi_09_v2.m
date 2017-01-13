function shi_09_v2
M = 4;
% SNR in db
SNRdb = -30 : 2 : 10;
% SNR
SNR = 10 .^ (SNRdb/10);
% shifted EBN0
Eb0 = SNR /log(M) *64;
% Sigema of noise
sigma = sqrt(2./SNR/2);

% Number of QAM symbols
N1 = 10^5;
% Number of Spreading sequence symbol
N2 = 64*N1;
% Generate the 2bits/symbols signel
a = randi(4,1,N1);
t = [1+1i -1+1i -1-1i 1-1i];
% 4 QAM signal
A = t(a);
% Spreading sequence
C = randi(2,1,N2)*2-3;
% Spread the symbol sequence
t = 1:N2;
t = ceil(t ./ 64);
S = A(t) .* C;

%-------------2 non-fading channel with impulse response---------
% Initialize SER
SER = zeros(1,length(SNRdb));
% Impulse response
H = [1,0,0,0,0,1];
% Pass signal through response
S1 = filter(H,1,S);
% Unit noise
n01 = randn(1,N2) + 1i .* randn(1,N2);
for l = 1: length(SNRdb)
    % Noise 
    n1 = sigma(l).* n01;
    % Received signal
    S2 = S1 + n1;
    % Path search
    % Since perfect knowledge of thechannel.
    % I assume receiver already know which channel carried signal
    % we have 2 branch
    S2 = [S2.*C;S2.*filter([0,0,0,0,0,1],1,C)];
    % MRC - Apply attenuation
    % Since the gain is one in this response
    % We pass this stage 
    % Integgrate
    S2 = [mean(reshape(S2(1,:),64,[]),1);mean(reshape(S2(2,:),64,[]),1)];
    % MRC - combaining
    S2 = sum(S2,1);
    S2 = S2(1,:);
    S2 = 1*(real(S2)>0&imag(S2)>=0) + 2*(real(S2)<=0&imag(S2)>0) + 3*(real(S2)<0&imag(S2)<=0) + 4*(real(S2)>=0&imag(S2)<0);
    SER(l) = (sum(S2 ~= a)) / length(S2);
end
figure('name','assignment 91');
semilogy(SNRdb,SER,'k') ;
hold on
SER = 2 .* qfunR (sqrt(2*Eb0));
semilogy(SNRdb,SER,'k*') ;
xlabel('SNR(dB)');
ylabel('symbol Error rate');
title('SER for RAKE receiver with channel h =[1, 0, 0, 0, 0, 1]');
legend('Sim','The');
xlim([-30 10])
ylim([0.00001 1])
hold 


%-------------2-fading channel ---------
% Initialize SER
SER = zeros(1,length(SNRdb));
% Impulse response
H = [1,0,0,0,0,1];
% Pass signal through response
g1 = abs(sqrt(0.5).*(randn(1,N1) + 1i * randn(1,N1)));
g1 = g1(t);
g2 = abs(sqrt(0.5).*(randn(1,N1) + 1i * randn(1,N1)));
g2 = g2(t);
S1 = filter([0,0,0,0,0,1],1,S).*g1 + S.*g2;
% Unit noise
n01 = randn(1,N2) + 1i .* randn(1,N2);
for l = 1: length(SNRdb)
    % Noise 
    n1 = sigma(l).* n01;
    % Received signal
    S2 = S1 + n1;
    % Path search
    % Since perfect knowledge of the channel.
    % I assume receiver already know which channel carried signal
    % we have 2 branch
    S2 = [S2.*C;S2.*filter([0,0,0,0,0,1],1,C)];
    % MRC - Apply attenuation
    % Since the gain is one in this response
    S2 = [S2(1,:).*g2;S2(2,:).*g1];
    % We pass this stage 
    % Integration
    S2 = [mean(reshape(S2(1,:),64,[]),1);mean(reshape(S2(2,:),64,[]),1)];
    % MRC - combaining
    S2 = sum(S2.*1,1);
    S2 = 1*(real(S2)>0&imag(S2)>=0) + 2*(real(S2)<=0&imag(S2)>0) + 3*(real(S2)<0&imag(S2)<=0) + 4*(real(S2)>=0&imag(S2)<0);
    SER(l) = (sum(S2 ~= a)) / length(S2);
end
figure('name','assignment 92');
semilogy(SNRdb,SER,'b') ;
hold on
%Theoretical BER
%Uncoded symbol error rate
D = sqrt(Eb0 ./(1 + Eb0));
SER = 2 * ((1-D)./2).^2.*(1+(2*((1+D)./2).^1));
semilogy(SNRdb,SER,'b*') ;
xlabel('SNR(dB)');
ylabel('symbol Error rate');
title('SER for RAKE receiver with channel h =[g1, 0, 0, 0, 0, g2]');
legend('Sim','The');
xlim([-30 10])
ylim([0.00001 1])
hold off

%-------------4-fading channel ---------
% Initialize SER
SER = zeros(1,length(SNRdb));
% Impulse response
H = [1,0,0,1,0,0,1,0,0,1];
% Pass signal through response
g1 = abs(sqrt(0.5).*(randn(1,N2) + 1i * randn(1,N2)));
g1 = g1(t);
g2 = abs(sqrt(0.5).*(randn(1,N2) + 1i * randn(1,N2)));
g2 = g2(t);
g3 = abs(sqrt(0.5).*(randn(1,N2) + 1i * randn(1,N2)));
g3 = g3(t);
g4 = abs(sqrt(0.5).*(randn(1,N2) + 1i * randn(1,N2)));
g4 = g4(t);
S1 = S.*g1 + filter([0,0,0,1,0,0,0,0,0,0],1,S).*g2 + filter([0,0,0,0,0,0,1,0,0,0],1,S).*g3 + filter([0,0,0,0,0,0,0,0,0,1],1,S).*g4;
% Unit noise
n01 = randn(1,N2) + 1i .* randn(1,N2);
for l = 1: length(SNRdb)
    % Noise 
    n1 = sigma(l).* n01;
    % Received signal
    S2 = S1 + n1;
    % Received signal
    S2 = S1 + n1;
    % Path search
    % Since perfect knowledge of thechannel.
    % I assume receiver already know which channel carried signal
    % we have 2 branch
    S2 = [S2.*C;S2.*filter([0,0,0,1,0,0,0,0,0,0],1,C);S2.*filter([0,0,0,0,0,0,1,0,0,0],1,C);S2.*filter([0,0,0,0,0,0,0,0,0,1],1,C)];
    % MRC - Apply attenuation
    % Since the gain is one in this response
    S2 = [S2(1,:).*g1;S2(2,:).*g2;S2(3,:).*g3;S2(4,:).*g4];
    % We pass this stage 
    % Integgrate
    S2 = [mean(reshape(S2(1,:),64,[]),1);mean(reshape(S2(2,:),64,[]),1);mean(reshape(S2(3,:),64,[]),1);mean(reshape(S2(4,:),64,[]),1)];
    % MRC - combaining
    S2 = sum(S2.*1,1);
    S2 = 1*(real(S2)>0&imag(S2)>=0) + 2*(real(S2)<=0&imag(S2)>0) + 3*(real(S2)<0&imag(S2)<=0) + 4*(real(S2)>=0&imag(S2)<0);
    SER(l) = (sum(S2 ~= a)) / length(S2);
end
figure('name','assignment 93');
semilogy(SNRdb,SER,'r') ;
hold on
D = sqrt(Eb0 ./(1 + Eb0));
SER = 2 .* ((1-D)./2).^4.*(1+(4*((1+D)./2)+(10*((1+D)./2).^2)+(20*((1+D)./2).^3)));
semilogy(SNRdb,SER,'r*') ;
hold off
xlabel('SNR(dB)');
ylabel('symbol Error rate');
title('SER for RAKE receiver for channel [g1, 0, 0, g2, 0, 0, g3, 0, 0, g4]');
xlim([-30 10])
ylim([0.00001 1])
legend('Sim','The');
end

function q = qfunR( x)
q = 0.5 * erfc (x/sqrt(2));
end