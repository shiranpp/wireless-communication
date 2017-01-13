function shi_07_v2


%Number of symbols
N = 10^6;
%Generate 2bits/symbols signal
b = randi(4,1,N/5);
%Coded signal to 5 continuous signal
t = 1:N;
t = ceil(t./5);
b1 = b (t);
%Transform ditital signal to 4 Qam signal
t = [1+1i -1+1i -1-1i 1-1i];
a = t(b1);

%For 4-Qam Signal
M=4;
%Eb/N0 in db
Eb0db = -10:2:40;
%Eb/N0
Eb0 = 10.^(Eb0db/10);
%SNR
SNR = Eb0 * log2(M);
%Signal power
Es = sum(abs(a).^2)/N;
%Sigema of nise
sigma = sqrt(Es./SNR/2);

%Unit complex AWGN
v = randn(1,N) + 1i * randn(1,N);

%-----------AWGN------------
%Initialize SER
SER1 = zeros(1,length(Eb0));
SER2 = zeros(1,length(Eb0));
figure('name','AWGN');
%Simulation BER
for l=1:length(Eb0);
    %Apply unit noise with sigema
    v1 = sigma (l) .* v;
    %Add noise and agin to 
    r = a + v1;
    %Hard decision decoding
    s2 = 1*(real(r)>0&imag(r)>=0) + 2*(real(r)<=0&imag(r)>0) + 3*(real(r)<0&imag(r)<=0) + 4*(real(r)>=0&imag(r)<0);
    %Decode the error control coding
    s3 = decode(s2);
    %Uncoded symbol error rate
    SER1(l) = (sum(s2 ~= b1)) / length(s2);
    %Coded symbol error rate
    SER2(l) = (sum(s3 ~= b)) / length(s3);
end
BER1 = SER1/2;
semilogy(Eb0db,BER1,'k') ;
hold on
BER2 = SER2/2;
semilogy(Eb0db,BER2,'r') ;

%Theoretical BER
%Uncoded symbol error rate
SER1 = 2 .* qfunR (sqrt(2*Eb0));
BER1 = SER1/2;
semilogy(Eb0db,BER1,'k*') ;
%Coded symbol error rate
SER2 = 2 .* qfunR (sqrt(6*Eb0));
BER2 = SER2/2;
semilogy(Eb0db,BER2,'r*') ;

xlabel('SNR(dB)');
ylabel('Error rate');
title('BER for AWGN');
xlim([-10 40])
ylim([0.00001 1])
legend('Simulation uncoded','Simulation coded','Theoretical uncoded','Theoretical uncoded');
hold off
%----------------------------------------------


%--------Rayleigh fading with MRC L=1--------
figure('name','Ray MRC1');
%Simulation BER
%Chi-Square gain with level 2
n1 = sqrt(0.5).*(randn(1,N) + 1i * randn(1,N));
g1 =  sqrt((abs(n1) .^ 2) .* (1));
for l=1:length(Eb0);
    %Apply unit noise with sigema
    v1 = sigma (l) .* v;
    %Add noise and agin to 
    r = a.*g1 + v1;
    %Hard decision decoding
    s2 = 1*(real(r)>0&imag(r)>=0) + 2*(real(r)<=0&imag(r)>0) + 3*(real(r)<0&imag(r)<=0) + 4*(real(r)>=0&imag(r)<0);
    %Decode the error control coding
    s3 = decode(s2);
    %Uncoded symbol error rate
    SER1(l) = (sum(s2 ~= b1)) / length(s2);
    %Coded symbol error rate
    SER2(l) = (sum(s3 ~= b)) / length(s3);
end
BER1 = SER1/2;
semilogy(Eb0db,BER1,'k') ;
hold on
BER2 = SER2/2;
semilogy(Eb0db,BER2,'r') ;

%Theoretical BER
%Uncoded symbol error rate
D = sqrt(Eb0 ./(1 + Eb0));
SER1 = 2 .* ((1-D) ./ 2);
BER1 = SER1/2;
semilogy(Eb0db,BER1,'k*') ;


xlabel('SNR(dB)');
ylabel('Error rate');
title('BER for Rayleigh fading with MRC L=1');
xlim([-10 40])
ylim([0.00001 1])
legend('Simulation uncoded','Simulation coded','Theoretical uncoded');
hold off
%-----------------------------------------------------

%--------Rayleigh fading with MRC L=2--------
figure('name','Ray MRC2');
%Simulation BER
%Chi-Square gain with level 4
n1 = sqrt(0.5).*(randn(1,N) + 1i * randn(1,N));
n2 = sqrt(0.5).*(randn(1,N) + 1i * randn(1,N));
g1 =  sqrt((abs(n1) .^ 2) + (abs(n2) .^ 2));
for l=1:length(Eb0);
    %Apply unit noise with sigema
    v1 = sigma (l) .* v;
    %Add noise and agin to 
    r = a.*g1 + v1;
    %Hard decision decoding
    s2 = 1*(real(r)>0&imag(r)>=0) + 2*(real(r)<=0&imag(r)>0) + 3*(real(r)<0&imag(r)<=0) + 4*(real(r)>=0&imag(r)<0);
    %Decode the error control coding
    s3 = decode(s2);
    %Uncoded symbol error rate
    SER1(l) = (sum(s2 ~= b1)) / length(s2);
    %Coded symbol error rate
    SER2(l) = (sum(s3 ~= b)) / length(s3);
end
BER1 = SER1/2;
semilogy(Eb0db,BER1,'k') ;
hold on
BER2 = SER2/2;
semilogy(Eb0db,BER2,'r') ;

%Theoretical BER
%Uncoded symbol error rate
D = sqrt(Eb0 ./(1 + Eb0));
SER1 = 2 .* ((1-D)./2).^2.*(1+(2*((1+D)./2)));
BER1 = SER1/2;
semilogy(Eb0db,BER1,'k*') ;


xlabel('SNR(dB)');
ylabel('Error rate');
title('BER for Rayleigh fading with MRC L=2');
xlim([-10 40])
ylim([0.00001 1])
legend('Simulation uncoded','Simulation coded','Theoretical uncoded');
hold off
%-----------------------------------------------------

%--------Rayleigh fading with MRC L=3--------
figure('name','Ray MRC3');
%Simulation BER
%Chi-Square gain with level 6
n1 = sqrt(0.5).*(randn(1,N) + 1i * randn(1,N));
n2 = sqrt(0.5).*(randn(1,N) + 1i * randn(1,N));
n3 = sqrt(0.5).*(randn(1,N) + 1i * randn(1,N));
g1 =  sqrt((abs(n1) .^ 2) + (abs(n2) .^ 2) + (abs(n3) .^ 2));
for l=1:length(Eb0);
    %Apply unit noise with sigema
    v1 = sigma (l) .* v;
    %Add noise and agin to 
    r = a.*g1 + v1;
    %Hard decision decoding
    s2 = 1*(real(r)>0&imag(r)>=0) + 2*(real(r)<=0&imag(r)>0) + 3*(real(r)<0&imag(r)<=0) + 4*(real(r)>=0&imag(r)<0);
    %Decode the error control coding
    s3 = decode(s2);
    %Uncoded symbol error rate
    SER1(l) = (sum(s2 ~= b1)) / length(s2);
    %Coded symbol error rate
    SER2(l) = (sum(s3 ~= b)) / length(s3);
end
BER1 = SER1/2;
semilogy(Eb0db,BER1,'k') ;
hold on
BER2 = SER2/2;
semilogy(Eb0db,BER2,'r') ;

%Theoretical BER
%Uncoded symbol error rate
D = sqrt(Eb0 ./(1 + Eb0));
SER1 = 2 .* ((1-D)./2).^3.*(1+(3*((1+D)./2)+(6*((1+D)./2).^2)));
BER1 = SER1/2;
semilogy(Eb0db,BER1,'k*') ;


xlabel('SNR(dB)');
ylabel('Error rate');
title('SER for Rayleigh fading with MRC L=3');
xlim([-10 40])
ylim([0.00001 1])
legend('Simulation uncoded','Simulation coded','Theoretical uncoded');
hold off
%-----------------------------------------------------

%--------Slow fading rayleigh fading with MRC L=1--------
figure('name','Slow-Ray MRC1');
%Simulation BER
%Chi-Square gain with level 2
n1 = sqrt(0.5).*(randn(1,N/5) + 1i * randn(1,N/5));
g =  sqrt((abs(n1) .^ 2) .* (1));
t = 1:N;
t = ceil(t./5);
g1 = g(t);
for l=1:length(Eb0);
    %Apply unit noise with sigema
    v1 = sigma (l) .* v;
    %Add noise and agin to 
    r = a.*g1 + v1;
    %Hard decision decoding
    s2 = 1*(real(r)>0&imag(r)>=0) + 2*(real(r)<=0&imag(r)>0) + 3*(real(r)<0&imag(r)<=0) + 4*(real(r)>=0&imag(r)<0);
    %Decode the error control coding
    s3 = decode(s2);
    %Uncoded symbol error rate
    SER1(l) = (sum(s2 ~= b1)) / length(s2);
    %Coded symbol error rate
    SER2(l) = (sum(s3 ~= b)) / length(s3);
end
BER1 = SER1/2;
semilogy(Eb0db,BER1,'k') ;
hold on
BER2 = SER2/2;
semilogy(Eb0db,BER2,'r') ;

%Theoretical BER
%Uncoded symbol error rate
D = sqrt(Eb0 ./(1 + Eb0));
SER1 = 2 .* ((1-D) ./ 2);
BER1 = SER1/2;
semilogy(Eb0db,BER1,'k*') ;


xlabel('SNR(dB)');
ylabel('Error rate');
title('SER for Slow fading rayleigh fading with MRC L=1');
xlim([-10 40])
ylim([0.00001 1])
legend('Simulation uncoded','Simulation coded','Theoretical uncoded');
hold off
%-----------------------------------------------------


end

%Qfunction
function q = qfunR( x)
q = 0.5 * erfc (x/sqrt(2));
end

%Decoder for error control coding
function s1=decode(s2)
s1 = zeros(1,ceil(length(s2)/5));
temp = 0;
temp1 = zeros(1,4);
for l=1:length(s2)
    temp1(s2(l)) =  temp1(s2(l)) + 1;
    temp = temp + 1;
    if (temp ==5) || (l==length(s2))
        [~,temp3] = max(temp1);
        s1(ceil(l/5)) = temp3;
        temp = 0;
        temp1 = zeros(1,4);
    end
end 
end