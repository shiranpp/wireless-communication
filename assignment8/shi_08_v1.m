function shi_08_v1
%I compute all possible signal in likehood
% and return the signal with minimum value
%What i seen in the figure is that MIMO can reach the same snr vs BER level
%of MRC and SIMO with same output branch.
%Which means the MIMO can reach the get same BER in same power of MRC while
%transmitting more information.
%Number of symbols
N1 = 10^2;
N2 = 10^4;
N = N1*N2;
%EbN0 in db
Eb0db = -10 : 2 : 30;
%EbN0
Eb0 = 10 .^ (Eb0db/10);
%S/N
SNR = Eb0;
%Signal power
Es = 1;
%Sigema of nise
sigma = sqrt(Es./SNR./2);
%unit signal
o1 = ones(1,N1);
o2 = -1 .* ones (1,N1);

%--------------MIMO 2*2-------------------
%Generate signal of 1,-1 bits
A = randi([0,1],2,N1) * 2 - 1;
%All signal state
s = zeros (2,2);
for d = 0:3
    s(d+1,:)=[-1 + 2*mod(d,2);-1 + 2*mod(floor(d/2),2)];
end
%Initialize BER
BER = zeros(1,length(Eb0db));
for d = 1:N2
    %Channel matrix 
    H = sqrt(0.5).*(randn(2,2) + randn(2,2) * 1i);
    %Unit complex gussian noise
    V =  randn(2,N1) + randn(2,N1) * 1i;
    for j = 1: length(Eb0db)
        %Reciveed signal
        Y = H * A + sigma(j) .* V;
        %Maximum likehood
        L = [sum(abs(Y-H*(s(1,:)'*o1)),1);sum(abs(Y-H*(s(2,:)'*o1)),1);sum(abs(Y-H*(s(3,:)'*o1)),1);sum(abs(Y-H*(s(4,:)'*o1)),1)];
        [~,l] = min(L);
        %Rebuild signal
        B = s(l,:)';
        BER(j) = BER(j)+ (sum(sum((A ~= B))))/2/N;
    end
end
semilogy(Eb0db,BER,'k-.') ;
hold on

%--------------MIMO 4*4-------------------
%Generate signal of 1,-1 bits
A = randi([0,1],4,N1) * 2 - 1;
%All signal state
s = zeros (16,4);
for d = 0:15
    s(d+1,:)=[-1 + 2*mod(d,2);-1 + 2*mod(floor(d/2),2);-1 + 2*mod(floor(d/4),2);-1 + 2*mod(floor(d/8),2)];
end

%Initialize BER
BER = zeros(1,length(Eb0db));
for d = 1:N2
    %Channel matrix 
    H = sqrt(0.5).*(randn(4,4) + randn(4,4) * 1i);
    %Unit complex gussian noise
    V =  randn(4,N1) + randn(4,N1) * 1i;
    for j = 1: length(Eb0db)
        %Reciveed signal
        Y = H * A + sigma(j) .* V;
        %Maximum likehood
        L = [sum(abs(Y-H*(s(1,:)'*o1)),1)];
        for z = 2:16
            L = [L;sum(abs(Y-H*(s(z,:)'*o1)),1)];
        end
        [~,l] = min(L);
        %Rebuild signal
        B = s(l,:)';
        BER(j) = BER(j)+ (sum(sum((A ~= B))))/4/N;
    end
end
semilogy(Eb0db,BER,'m-.') ;
hold on


%--------------SIMO 1*2-------------------
%Generate signal of 1,-1 bits
A = randi([0,1],1,N1) * 2 - 1;
%Initialize BER
BER = zeros(1,length(Eb0db));
for d = 1:N2
    %Channel matrix 
    H = sqrt(0.5).*(randn(2,1) + randn(2,1) * 1i);
    %Unit complex gussian noise
    V =  randn(2,N1) + randn(2,N1) * 1i;
    for j = 1: length(Eb0db)
        %Reciveed signal
        Y = H * A + sigma(j) .* V;
        %Maximum likehood
        L = [sum(abs(Y-H*o1),1);sum(abs(Y-H*o2),1)];
        [~,l] = min(L);
        %Rebuild signal
        B = 1* (l==1) + -1* (l==2);
        BER(j) = BER(j)+ (sum(sum((A ~= B))))/N;
    end
end
semilogy(Eb0db,BER,'r') ;
hold on
%Theoretical BER
%Uncoded symbol error rate
D = sqrt(Eb0 ./(1 + Eb0));
BER =   ((1-D)./2).^2.*(1+(2*((1+D)./2)));
semilogy(Eb0db,BER,'k*') ;

%--------------SIMO 1*4-------------------
%Generate signal of 1,-1 bits
A = randi([0,1],1,N1) * 2 - 1;
%Initialize BER
BER = zeros(1,length(Eb0db));
for d = 1:N2
    %Channel matrix 
    H = sqrt(0.5).*(randn(4,1) + randn(4,1) * 1i);
    %Unit complex gussian noise
    V =  randn(4,N1) + randn(4,N1) * 1i;
    for j = 1: length(Eb0db)
        %Reciveed signal
        Y = H * A + sigma(j) .* V;
        %Maximum likehood
        L = [sum(abs(Y-H*o1),1);sum(abs(Y-H*o2),1)];
        [~,l] = min(L);
        %Rebuild signal
        B = 1* (l==1) + -1* (l==2);
        BER(j) = BER(j)+ (sum(sum((A ~= B))))/N;
    end
end
semilogy(Eb0db,BER,'g') ;
hold on

%--------------SISO 1*1-------------------
%Generate signal of 1,-1 bits
A = randi([0,1],1,N1) * 2 - 1;
%Initialize BER
BER = zeros(1,length(Eb0db));
for d = 1:N2
    %Channel matrix 
    H = sqrt(0.5).*(randn(1,1) + randn(1,1) * 1i);
    %Unit complex gussian noise
    V =  randn(1,N1) + randn(1,N1) * 1i;
    for j = 1: length(Eb0db)
        %Reciveed signal
        Y = H * A + sigma(j) .* V;
        %Maximum likehood
        L = [sum(abs(Y-H*o1),1);sum(abs(Y-H*o2),1)];
        [~,l] = min(L);
        %Rebuild signal
        B = 1* (l==1) + -1* (l==2);
        BER(j) = BER(j)+ (sum(sum((A ~= B))))/N;
    end
end
semilogy(Eb0db,BER,'b') ;
hold on
%Theoretical BER
%Uncoded symbol error rate
D = sqrt(Eb0 ./(1 + Eb0));
BER = ((1-D) ./ 2);
semilogy(Eb0db,BER,'b*') ;

legend('MIMO 2X2','MIMO 4X4','Sim SIMO 1X2','The SIMO 1X2','SIMO 1X4','Sim SISO 1X1', 'The  SISO 1X1');
xlabel('Eb/No(dB)');
ylabel('Bit Error rate');
title({'BER for Different ?I?O System with Bit Pulse Signal ';' By using maximum liklihood for all possible signal'});
xlim([-10 30])
ylim([0.00001 1])
end

