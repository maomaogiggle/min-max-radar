%clairvoyant SINR loss
clc;clear;

radar = radar_init;
%inter-element space
d = radar.d;
%antenna index
nvs = radar.pos';
%pulse index
nvt = radar.pulse;
%antenna number
N = radar.N;
%pulse number
M = radar.M;
%wavelength
lambda = radar.lambda;
MN = M*N;
%selection vector
K = 16;
%number of range gates for true covariance matrix
L = 1000;

%read the file for optimum sub-configuration for each Doppler frequency
load('Sx.mat');

%non-optimum sub-configuration
x2_t = [ones(4,1);zeros(4,1)];
x2_s = [ones(4,1);zeros(4,1)];
x2 = reshape(x2_s*x2_t',MN,1);
It2 = find(x2);
Itt2 = It2*ones(1,L) + MN*ones(length(It2),1)*(0:L-1);

fs = 0;
f = [-0.1:0.01:0.1,0.1:0.02:0.5];
s_s = exp(1i*2*pi*nvs*fs);

%statistic matrix
SINR1 = zeros(1,length(f));
SINR = zeros(1,length(f));
SINR2 = zeros(1,length(f));

%generate the true covariance matrix for the whole array
CMR = clutter_gen(radar,0,L);
C_f = (CMR*CMR')/L;
Pc = trace(C_f)/MN;
Pn = Pc*(10^(-radar.CNR/10));

%covariance matrix for the testing data
CMRn = CMR + sqrt(Pn/2)*(randn(MN,L)+1i*randn(MN,L));
C = (CMRn*CMRn')/L;

C2 = (CMRn(Itt2)*CMRn(Itt2)')/L;

for i = 1:length(f) 
    s_t = exp(1i*2*pi*nvt*f(i));
    S = s_s*s_t;
    s = S(:); 

    %searching for the optimum configuration
    x1 = Sx(:,i);
    It1 = find(x1);
    Itt1 = It1*ones(1,L) + MN*ones(length(It1),1)*(0:L-1);
    %selected covariance matrix
    C1 = (CMRn(Itt1)*CMRn(Itt1)')/L;
     
    %clairvoyant SINR loss for the whole array
    SINR(i) = abs(s'*(C\s))/(MN/Pn);
        
    %clairvoyant SINR loss for the subarray 1
    sN = s(logical(x1),:);         
    SINR1(i) = abs(sN'*(C1\sN))/(K/Pn);
        
    %clairvoyant SINR loss for the subarray 2
    sN = s(logical(x2),:);          
    SINR2(i) = abs(sN'*(C2\sN))/(K/Pn);
       
end   

figure;
plot(f,10*log10(SINR));
hold on;
plot(f,10*log10(SINR1),'r');
hold on;
plot(f,10*log10(SINR2),'k');










