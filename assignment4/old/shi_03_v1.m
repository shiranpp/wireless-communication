function shi_03_v1
T=4;
fd=30;
fs=2000;
n0 = genqam ( 4 , T*fs );
n1 = difR( n0 , 1 );
Eb = sum ( n1 .^2 )/T;
for fd = [20 60 100]
    g = genG( fd , fs ,T);
    for Eb0 = 0:5:60
        N0 = Eb/10^(Eb0/10);
        sigma =  sqrt(N0/2);
        no = randn(1,fs*T) * sigma + randn(1,fs*T)* sigma * i;
        y = n1 .* g + no;
        
    end
end
end

function n0 = genqam( m , n )
%generate an n symbols and m bits persymbols QAM sequence.%
k=log2(m);
n0=randi( 2 ^ ( k-1 ), 1 , n )*2 - 2 * ( k-1 ) - 1 + (randi( 2 ^ ( k - 1 ) , 1 , n )*2 - 2 * ( k-1 )-1) * i;
end

function n2 = difR( n1 , a )
%generate the Differentially encode of n0 and initial value a%
n2 = [ [a] n1 ];
n2 = diff( n2 );
end

function n1 = idifR( n0 , a )
%generate the Differentially encode of n0 and initial value a%
n1=zeros(1,length(n0));
n1(1)=a+n0(1)
for l=2:1:length(n0)
    n1(l)=n1(l-1)+n0(l);
end  
end

function g = genG( fd , fs , T) %generate the gain with doppler frequency and sampling frequency
fx = - fs / 2 : 1 : fs / 2;
S = zeros( 1 , length( fx ) );
n = randn(1,fs*T) + randn(1,fs*T) * i;%white Noise with sigama 1 mu 0
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

