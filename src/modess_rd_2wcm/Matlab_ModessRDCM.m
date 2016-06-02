addpath(genpath('/home/doru/Doru_utilities'))
% fpth = '/home/doru/infra/Cpp_code/DoruV/cpp_tests/Modal_NBtest1/samples/2012-12-18_17_26_51' %0.01 Hz
% fpth = '/home/doru/infra/Cpp_code/DoruV/cpp_tests/Modal_NBtest1/samples/2012-12-18_12_55_00' % 0.5 Hz
% fpth = '/home/doru/infra/Cpp_code/DoruV/cpp_tests/Modal_NBtest1/samples/2012-12-19_11_54_39' % 0.01 Hz
% fpth = '/home/doru/infra/Cpp_code/DoruV/cpp_tests/Modal_NBtest1/samples/2012-12-19_13_30_58'
fpth = '/home/doru/infra/Cpp_code/DoruV/cpp_tests/Modal_NBtest1/samples/2012-12-19_16_29_01'
%% read eigen vals vecs
a = textread(fullfile(fpth, 'eigvalvecs_1.dat'));

%%
freq = a(1,1)
Nm = a(1,2)
Nz = a(1,3)
dz = a(1,4)
%%
kk = complex(a(2:2+Nm-1,2), a(2:2+Nm-1,3))
rho  = a(2+Nm:2+Nm+Nz-1, 2);
vs = reshape(a(2+Nm+Nz:end,2), Nz, Nm);
%% view nth mode
moden = 7;
I = (1+Nm)+Nz+(moden-1)*Nz+1: ((1+Nm)+Nz + moden*Nz);
figure
plot(a(I,2), a(I,1), 'm')
title('complex')
hold on; grid on
% plot(Re(I,2), Re(I,1), 'b')
title(sprintf('Mode %d', moden))


%% read Rmat
tmp = textread(fullfile(fpth,'Rmat_1.dat'));
N1 = sqrt(size(tmp,1))
%%
RR = reshape((complex(tmp(:,1), tmp(:,2))),N1,N1);
RR=RR.';
clear tmp
%% view R matrix
figure
% imagesc(abs(a))
imagesc(log10(abs(RR)))
hold on; zoom on;
colorbar

%%  zero very small elements
I = find(RR<1e-14);
RR(I) = 0;
%%
Nmin = 7;
%% initialize the S matrix
S = eye(2*Nmin);
%% list of Rmat files
% listf = list_files_in_folder(fpth, 'Rmat_')
flist = get_filename_by_patterns(fpth, 'Rmat_')
flist = sort(flist)
for j=1:length(flist)
    fprintf('%s\n', flist{j})
end

nregions = length(flist) + 1
Rv = [100 100 200 250 500]*1000;

%% update the S matrix
for j=1:length(flist)
    j
    tmp = textread(flist{j});
    N1 = sqrt(size(tmp,1));
    RR = reshape((complex(tmp(:,1), tmp(:,2))),N1,N1);
    RR=RR.';
%     I=find(RR<1e-14);
%     RR(I) = 0;
    
    
    S = RR*S;
end
%% view S matrix
figure
% imagesc(abs(a))
imagesc(log10(abs(S)))
hold on; zoom on;
colorbar
title('S')

%% read S outputted from cpp version
tmp = textread(fullfile('/home/doru/infra/Cpp_code/DoruV/cpp_tests/Modal_NBtest1/samples/','S_curr.txt'));
N1 = sqrt(size(tmp,1))
Scpp = reshape((complex(tmp(:,1), tmp(:,2))),N1,N1);
Scpp=Scpp.';
%% view Scpp matrix
figure
% imagesc(abs(a))
imagesc(log10(abs(Scpp)))
hold on; zoom on;
colorbar
title('Scpp')

%% use real part of k?
use_real = 1;



%% read eigenvecs, vals
a = textread(fullfile(fpth, 'eigvalvecs_1.dat'));
freq = a(1,1)
Nm = a(1,2)
Nz = a(1,3)
dz_km = a(1,4)
kk_curr = complex(a(2:2+Nm-1,2), a(2:2+Nm-1,3)); % Nm elements
rho_curr = a(2+Nm:2+Nm+Nz-1, 2);    % Nz elements
vs_curr = reshape(a(2+Nm+Nz:end,2), Nz, Nm);

if use_real
    kk_curr = real(kk_curr);
end

%% trim to Nmin
kk_curr = kk_curr(1:Nmin);
vs_curr = vs_curr(:, 1:Nmin);

%% form D
r1 = Rv(2);
dii = zeros(1,Nmin);
dii = exp(2*1i*(kk_curr*r1-pi/4));
D = diag(dii);

%%
nz_src = 1;
H01 = sqrt(2./(pi*kk_curr*r1)).*exp(1i*(kk_curr.*r1-pi/4));
ss = 1i/(4*rho_curr(nz_src))*sqrt(rho_curr(nz_src)).*(vs_curr(nz_src,:).').*H01;

% ss2 = 1i*exp(-1i*pi/4.0)./(rho_curr(nz_src)*sqrt(8*pi*kk_curr*r1)).*sqrt(rho_curr(nz_src)).*(vs_curr(nz_src,:).').*exp(1i*kk_curr*r1);
%       ss[m] = I*exp(-I*Pi/4.0)/(rho_curr[n_zsrc]*sqrt(8*Pi*k_curr[m]*Rv[0]))*sqrtrho_zs*v_curr[n_zsrc][m]*exp(I*k_curr[m]*Rv[0])
%% find b1: from (S4+S3*D)b1 = -S3*ss
% Nmin = Nm;
b1 = (S(Nmin+1:end, Nmin+1:end) + S(Nmin+1:end, 1:Nmin)*D)\(-S(Nmin+1:end, 1:Nmin)*ss);
a1 = ss + D*b1;
ab_curr = [a1; b1];

%% pressure in region 1
N = 500; % number of points in range
rngv = (1:N)*1000;
P = zeros(1,N);

j=1;
rng = rngv(1);
while(rng<=Rv(1))
    H11 = sqrt(Rv(1)/rng)*exp(1i*kk_curr*(rng-Rv(1)));
    H12 = sqrt(Rv(1)/rng)*exp(-1i*kk_curr*(rng-Rv(1)));

%     H11 = sqrt(Rv(1)/rng)*exp(1i*real(kk_curr)*(rng-Rv(1)));
%     H12 = sqrt(Rv(1)/rng)*exp(-1i*real(kk_curr)*(rng-Rv(1)));
    P(j) = sum((ab_curr(1:Nmin).*H11 + ab_curr(Nmin+1:end).*H12).*(vs_curr(nz_src, :).')*sqrt(rho_curr(nz_src)));
    j = j + 1
    rng = rngv(j);
end

%% plot pressure in first region
figure
plot(20*log10(abs(4*pi*P)))
hold on; zoom on

%%
for reg_ix = 2:nregions
    tmp = textread(flist{reg_ix-1});
    N1 = sqrt(size(tmp,1));
    RR = reshape((complex(tmp(:,1), tmp(:,2))),N1,N1);
    RR=RR.';  
    ab_next = RR*ab_curr;
    ab_curr = ab_next;
    % read from eigenvalvec_j.dat file
    fn = sprintf('eigvalvecs_%d.dat', reg_ix);
    % read eigenvecs, vals
    a = textread(fullfile(fpth, fn));
    freq = a(1,1);
    Nm = a(1,2);
    Nz = a(1,3);
    dz_km = a(1,4);
    kk_curr = complex(a(2:2+Nm-1,2), a(2:2+Nm-1,3)); % Nm elements
    rho_curr = a(2+Nm:2+Nm+Nz-1, 2);    % Nz elements
    vs_curr = reshape(a(2+Nm+Nz:end,2), Nz, Nm);
%     trim to Nmin
    kk_curr = kk_curr(1:Nmin);
    vs_curr = vs_curr(:, 1:Nmin);
    while rng<Rv(reg_ix+1)
        H11 = sqrt(Rv(reg_ix)/rng)*exp(1i*kk_curr*(rng-Rv(reg_ix)));
        H12 = sqrt(Rv(reg_ix)/rng)*exp(-1i*kk_curr*(rng-Rv(reg_ix)));
%         H11 = sqrt(Rv(reg_ix)/rng)*exp(1i*real(kk_curr)*(rng-Rv(reg_ix)));
%         H12 = sqrt(Rv(reg_ix)/rng)*exp(-1i*real(kk_curr)*(rng-Rv(reg_ix)));
        P(j) = sum((ab_curr(1:Nmin).*H11 + ab_curr(Nmin+1:end).*H12).*(vs_curr(nz_src, :).')*sqrt(rho_curr(nz_src)));
        j = j + 1
        rng = rngv(j);
    end
end

%% plot pressure
figure
%
plot(20*log10(abs(4*pi*P)))
hold on; zoom on


%% -----------------------------------------------------

%% read ab1_curr
samples_dir = '/home/doru/infra/Cpp_code/DoruV/cpp_tests/Modal_NBtest1/samples/'
tmp = textread(fullfile(samples_dir,'ab_curr.txt'))







%% eigen vecs are normalized
j = 1
p = sum(v_curr(:,j).*(v_curr(:,j)))*dz_km*1000



%% load tloss_rdcm_1d.nm
tl = textread(fullfile('/home/doru/infra/Cpp_code/DoruV/cpp_tests/Modal_NBtest1/samples/','tloss_rdcm_1d.nm'));
figure
plot(tl(:,1), 20*log10(abs(complex(tl(:,2),tl(:,3)))));
hold on


%% read Rmat
tmp = textread(fullfile(fpth,'Rmat_2.dat'));
N1 = sqrt(size(tmp,1))
%%
RR = reshape((complex(tmp(:,1), tmp(:,2))),N1,N1);
RR=RR.';
clear tmp
%% view R matrix
figure
% imagesc(abs(a))
imagesc(log10(abs(RRmlab)))
hold on; zoom on;
colorbar

%%  zero very small elements
I = find(RR<1e-14);
RR(I) = 0;

%% diagonal of RR
dRR = zeros(1,2*Nm);
for j=1:2*Nm
    dRR(j) = RR(j,j);
end

figure
plot(log10(dRR))
hold on; zoom on; grid on;


%% calculate RR

%% read eigen vals vecs
a = textread(fullfile(fpth, 'eigvalvecs_1.dat'));

%%
freq = a(1,1)
Nm = a(1,2)
Nz = a(1,3)
dz_km = a(1,4)
%%
k_curr   = complex(a(2:2+Nm-1,2), a(2:2+Nm-1,3))
rho_curr = a(2+Nm:2+Nm+Nz-1, 2);
v_curr   = reshape(a(2+Nm+Nz:end,2), Nz, Nm);
%% trim to Nmin
kk_curr = kk_curr(1:Nmin);
vs_curr = vs_curr(:, 1:Nmin);

%% read eigen vals vecs
a = textread(fullfile(fpth, 'eigvalvecs_2.dat'));
%%
freq = a(1,1);
Nm = a(1,2);
Nz = a(1,3);
dz_km = a(1,4);
%%
k_next = complex(a(2:2+Nm-1,2), a(2:2+Nm-1,3))
rho_next  = a(2+Nm:2+Nm+Nz-1, 2);
v_next = reshape(a(2+Nm+Nz:end,2), Nz, Nm);

%% trim to Nmin
kk_next = kk_next(1:Nmin);
vs_next = vs_next(:, 1:Nmin);

%%
if use_real
    kk_curr = real(kk_curr);
	kk_next = real(kk_next);
end

%% the H matrices
%  r1 = 100000; r0 = r1; %first region
r1 = 200000; r0=100000; % for 2nd region
sr1ovr2 = sqrt(r0/r1);
H1 = sr1ovr2*exp(1i*(r1-r0)*(k_curr));
H2 = sr1ovr2*exp(-1i*(r1-r0)*(k_curr));

H1mx = diag(H1);
H2mx = diag(H2);

%% C_LR
dz = dz_km*1000;
C_LR = zeros(Nm,Nm);
C_RL = zeros(Nm,Nm);
for l=1:Nm
    for m=1:Nm
        C_LR(l,m) = k_curr(m)/k_next(l).*sum(sqrt(rho_next./rho_curr).*v_next(:,l).*v_curr(:,m)*dz); % eq 63 in NMRD-CM notes
        C_RL(l,m) = sum(sqrt(rho_curr./rho_next).*v_next(:,l).*v_curr(:,m)*dz); %eq. 62
        
    end
end
R1 = 1/2*(C_LR + C_RL)*H1mx;
R2 = 1/2*(-C_LR + C_RL)*H2mx;
R3 = 1/2*(-C_LR + C_RL)*H1mx;
R4 = 1/2*(C_LR + C_RL)*H2mx;
RRmlab = [R1 R2;R3 R4];

%%
figure
imagesc(10*log10(abs(C_RL)))
hold on
colorbar

%% diagonal of RR and RRmlab
dRR = zeros(1,2*Nm);
dRRmlab = zeros(1,2*Nm);
for j=1:2*Nm
    dRR(j) = RR(j,j);
    dRRmlab(j) = RRmlab(j,j);
end

figure
plot((dRR))
hold on; zoom on; grid on;
plot(abs(dRRmlab), 'r')
legend('cpp', 'mlab')

%% load tloss_rdcm_1d.nm
% fpth2 = '/home/doru/infra/Cpp_code/DoruV/cpp_tests/Modal_NBtest1/samples/compare_1way_vs_2way_coupled_modes'
fpth2 = '/home/doru/infra/Cpp_code/DoruV/cpp_tests/Modal_NBtest1/samples/'
tl = textread(fullfile(fpth2,'tloss_rdcm_1d.nm'));
figure
plot(tl(:,1), 20*log10(abs(complex(tl(:,2),tl(:,3)))));
hold on

%% load tloss_rd_1d.nm
% tl1way = textread(fullfile(fpth2,'tloss_rd_1d.lossless.nm'));
tl1way = textread(fullfile(fpth2,'tloss_rd_1d.nm'));
figure
plot(tl1way(:,1), 20*log10(abs(complex(tl1way(:,2),tl1way(:,3)))));
hold on
%%
figure
plot(tl1way(:,1), 20*log10(abs(complex(tl1way(:,2),tl1way(:,3)))));
hold on; grid on; zoom on;
plot(tl(:,1), 20*log10(abs(complex(tl(:,2),tl(:,3)))), 'r');
legend('1-way coupled modes', '2-way coupled modes')

title(sprintf('Freq=0.5 Hz; atmosfile g2sgcp2011012606L.jordan.env'))

%% plot 2D TL
fpth1 = '/home/doru/infra/Cpp_code/DoruV/wip_Modal_Code_20121106/samples'
a = textread(fullfile(fpth1, 'tloss_rdcm_2d.nm'));
%%
figure
plot(



