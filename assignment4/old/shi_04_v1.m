function n2=shi_04_v1
T=10;
a=0.25*pi;%initial phase of DIF encode
fd=30;
fs=2000;
n0 = randi(4,1,fs*T)-1;%generate 2 bits symbols with 2000 symbols/second and 4seconds
n1 = exp(-i*(n0/4*2*pi));%generate the QAM of signal
n2 = zeros(1,length(n1)+1);%initial the encoded signal
n2(1)=exp(-i*a);
for l=2:1:length(n1)+1
    n2(l) = n2(l-1)*n1(l-1);%encoding signal
end
Eb = sum(abs(n2).^2)/T/fs;
x=0:5:60;
x1=zeros(1,length(x));
for fd = [20 60 200]
    g = genG( fd , fs ,length(n2));
    for Eb0 = 0:5:60
        N0 = Eb/(10^(Eb0/10));
        sigma =  sqrt(N0/2);
        no = rand(1,length(n2)) * sigma + rand(1,length(n2))* sigma * i;
        y = n2 + no;
        n3 = difR(y);
        n4 = dis(n3);
        x1(Eb0/5+1) = 1 - ser(n0,n4)
    end
end
end

function x = ser(n0,n4)
x=0;
for l = 1:1:length(n4)
    if n0(l)==n4(l)
        x=x+1;
    end
end
x=x/length(n4);
end

function n3 = difR(n2)
n3=zeros(1,length(n2)-1);;
for l = 1:1:length(n2)-1
    n3(l) = n2(l+1)/n2(l);
end
end

function n4 = dis(n3)
for l=1:1:length(n3)
    if real(n3(l))>=0
        if imag(n3(l))>=0
            n4(l) = 0;
        else
            n4(l) = 1;
        end
    else
        if imag(n3(l))>=0
            n4(l) = 3;
        else
            n4(l) = 2;
        end  
    end
end
end

function g = genG( fd , fs , n1) %generate the gain with doppler frequency and sampling frequency and time T
fx = - fs / 2 : 1 : fs / 2;
S = zeros( 1 , length( fx ) );
n = randn(1,n1) + randn(1,n1) * i;%white Noise with sigama 1 mu 0
for j = -fs / 2 : 1 : fs / 2
    if abs( j ) < fd
        S( j + fs / 2 + 1 ) = 1 / sqrt( 1 -( j / fd ) ^ 2);%Produce the PSD of Doppler function
    else
        S( j + fs / 2 + 1 ) = 0;
    end
end
S1 = circshift(S,ceil(size(S)/2));
S2 = sqrt(S);%frequency respose of Doppler filter
h = ifft(sqrt(S1))%Doppler filter
h = circshift(h,floor(size(h)/2));
g = filter ( h , 1 , n );
end