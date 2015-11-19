function CMR = clutter_gen(radar,nR,K)
lambda = radar.lambda;
Nphi = 720;
dphi = 2*pi/Nphi;
d = radar.d;
vp = radar.vp;
PRF = radar.PRF;
mv = radar.pulse;
nv = radar.pos';
M = radar.M;
N = radar.N;

%----------------------------------------------------
%generate the clutter for the required range
H = radar.H;
dR = radar.dR;
Rs = radar.Rs;
Amp = sqrt(radar.Pt/2); %No free space loss
CMR = zeros(M*N,K);
theta = asin(H/(Rs+nR*dR));

SM = zeros(M*N,Nphi);

for nphi = 1:Nphi
        phi = (nphi-1)*dphi;
        fd = (2*vp/lambda)*cos(phi)*cos(theta);
        fs = (d/lambda)*cos(phi)*cos(theta);
        s_t = exp(1i*2*pi*mv*fd/PRF);
        s_s = exp(1i*2*pi*nv*fs);
        S = s_s*s_t;
        SM(:,nphi) = S(:);
end

for k = 1:K
    Rf = Amp*(randn(Nphi,1)+1i*randn(Nphi,1));
    CMR(:,k) = SM*Rf;
end
