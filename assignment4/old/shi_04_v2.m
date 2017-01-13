function ser = shi_04_v2
T=10;%total Times
a = 1/sqrt(2)+1/sqrt(2)*1i;
fs=2000;
n0 = randi(4,1,fs*T)-1;%generate 2 bits symbols with 2000 symbols/second and 4seconds
n1 = exp(-1i*(n0/4*2*pi-0.25*pi));%generate the QAM of signal
Eb0db = 0:5:60;%Eb/N0 in db
Eb = sum(abs(n1).^2)/T/fs;
N0 = Eb./10.^(Eb0db/10);
sigma = sqrt( N0 /2);% Sigma of noise
no = randn(length(N0),length(n1)) + randn(length(N0),length(n1)) * 1i;
no=no.*repmat(sigma,length(n1),1)';%Produce the noise for 
n3 = repmat(n1,length(sigma),1)+no;
n4 = Dec(n3);
ser = n4 -  repmat(n0,length(sigma),1);
ser = (0~=ser);
ser = (10^-5+sum(ser'))./fs/T;%Symbol error rate
figure('name','assignment4')
semilogy(Eb0db,ser) 
hold on
%semilogy(Eb0db,2*erfc(sqrt(Eb0*2)),'b--o') 
n2=zeros(1,length(n1));
n2(1)=a;
for l=2:length(n1)+1
    n2(l) = n2(l-1)*n1(l-1);
end
Eb = sum(abs(n2).^2)/T/fs;
N0 = Eb./10.^(Eb0db/10);
sigma = sqrt( N0./2 );% Sigma of noise
no = randn(length(N0),length(n2)) + randn(length(N0),length(n2)) * i;
no=no.*repmat(sigma,length(n2),1)';%Produce the noise for 
n3 = repmat(n2,length(sigma),1)+no;
n3 = (n3(:,2:length(n3))./n3(:,1:length(n3)-1));
n4 = Dec(n3);
ser = n4 -  repmat(n0,length(sigma),1);
ser = (0~=ser);
ser = (10^-5+sum(ser'))./fs/T;%Symbol error rate
semilogy(Eb0db,ser) 

for fd=[20 60 200]
    g = genG( fd , length(n2))
    n3 = repmat(n2,length(sigma),1).*repmat(g,length(sigma),1)+no;
    n31 = (n3(:,2:length(n3))./n3(:,1:length(n3)-1));
    n4 = Dec(n31);
    ser = n4 -  repmat(n0,length(sigma),1);
    ser = (0~=ser);
    ser = (10^-5+sum(ser'))./fs/T;%Symbol error rate
    semilogy(Eb0db,ser,'--') 
end

end

function n4 = Dec(n3)
n4 = atan2(imag(n3),real(n3));
n4 = 0.*(0<=n4 & n4 < pi/2) + 3.*(pi/2<=n4 & n4 <= pi) + 1.*(-pi/2<=n4 & n4<0) + 2.*(-pi<=n4 & n4<-pi/2) ;
end

function g = genG( fd , n1) %generate the gain with doppler frequency and sampling frequency and time T
fs=2000;
fx = - fs / 2 : 1 : fs / 2;
S = zeros( 1 , length( fx ) );
n = randn( 1 , n1 ) + randn( 1 , n1 ) * 1i;%white Noise with sigama 1 mu 0
n =n*.10;
S=zeros(1,length(fx));
for l=-fs/2:1:fs/2
    if abs(l)<fd
        S(l+fs/2+1)=1/4/pi/fd/sqrt(1-(l/fd)^2);%Produce the PSD of Doppler function
    else
        S(l+fs/2+1)=0;
    end
end
S1 = circshift(S,ceil(size(S)/2));
h = ifft(sqrt(S1));%Doppler filter
h = circshift(h,floor(size(h)/2));
g = filter ( h , 1 , n );
end