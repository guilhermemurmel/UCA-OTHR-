clear; close all; clc;
%--------------------------------------------------------------------------
% [References] MURMEL, Guilherme, F. DoA Estimation Performance of UCAs
% with Reduced Number of Sensors using Phase-Mode Transformation and
% Small Sample Support. Dissertação de Mestrado, Instituto Militar de
% Engenharia, 2022.
%--------------------------------------------------------------------------
montecarlo=3;                    % Number of trials
%--------------------------------
mu=0.8;                          % Factor for leankage minimization (0<mu<1)
                                 % Reduce mu factor when reducing sample support
%--------------------------------
SIRdB=-3;                         % 10*log10(sigmaS1/sigmaS2)
%--------------------------------
SNRdB=-30;                       % 10*log10(sigmaS1/sigman2)
%--------------------------------
SampleSupport=100;               % Ensemble for batch processing
%--------------------------------
% signal parameters (for every trial)
%--------------------------------
fs = 16e7;                        % Sampling frequency
fo= 12.5e6;                      % Array operating frequency and SoI
f1= 15e6;                        % Interferer central frequency
dt = 1/fs;
duration = 1e-4;
t = (0:dt:duration-dt)';
%--------------------------------
% Bandpass Narrowband filters for signal sources
%--------------------------------
wpass1=[0.95*fo  1.05*fo];
wpass2=[0.90*f1  1.10*f1];
%--------------------------------
% Array parameters
%--------------------------------
M=25;                           % Number of sensors in the array
N=7;                            % Phase Mode index
MODES=2*N+1;                     % Number os Phase Modes (spatial harmonics)
%--------------------------------
% Propagation velocity
%--------------------------------
vpropag=299792458;
%--------------------------------
% UCA Radius and Wavenumbers used in the steering vector
%--------------------------------
lambda=vpropag/fo;
%--------------------------------
K=2*pi/lambda;                    % Wavenumber of SoI (center frequency)
R=lambda*M/(4*pi);                % Radius to assure a inter-element spacing of lambda/2
KR=K*R;
%--------------------------------
D=2;                              % Number of sources
%--------------------------------
phideg=0:0.1:360;                 % Discrete azimuth parameter scan
phi=deg2rad(phideg);
tetadeg=0:0.1:90;                 % Discrete zenith parameter scane
teta=deg2rad(tetadeg);
%--------------------------------
% Signal received by each sensor
%--------------------------------
% SoI source DoA
%--------------------------------
phisim1=37.9;                    % AZIMUTH 1
phisim1=deg2rad(phisim1);
% tetasim1=90;                   % Used only for acessing anisotropic behavior
tetasim1=19.7;                   % ZENITH 1
tetasim1=deg2rad(tetasim1);
%--------------------------------
% Interference source DoA
%--------------------------------
phisim2=120.5;                   % AZIMUTH 2
phisim2=deg2rad(phisim2);
tetasim2=45.7;                   % ZENITH 2
tetasim2=deg2rad(tetasim2);
%--------------------------------
% Initial presumed zenith angle (Initially, sources are in the same plane)
% step 1 of decoupled 2-D DoA
%--------------------------------
zenini=pi/2;
%--------------------------------
plotaz1=zeros(1,montecarlo);
plotaz2=zeros(1,montecarlo);
plotze1=zeros(1,montecarlo);
plotze2=zeros(1,montecarlo);
%--------------------------------
% Received signals (for each trial)
%--------------------------------
for y=1:montecarlo
    %--------------------------------
    S1 = cos(2*pi*fo*t)+1i*sin(2*pi*fo*t);
    % S1 = rand(duration*fs,1);
    % S1 = rand(duration*fs,1)+1i*rand(duration*fs,1);
    % S1 = bandpass(S1,wpass1,fs);
    S1=S1-mean(S1); S1=S1/std(S1);
    % S2 = cos(2*pi*f1*t+(pi/2)*(rand))+1i*sin(2*pi*f1*t+(pi/2)*(rand));
    % S2 = rand(duration*fs,1);
    S2 = rand(duration*fs,1)+1i*rand(duration*fs,1);
    S2 = bandpass(S2,wpass2,fs);
    sigmaS1=var(S1);
    sigmaS2=sigmaS1*10^(-SIRdB/10);
    S2=S2-mean(S2); S2=sqrt(sigmaS2)*S2/std(S2);
    %--------------------------------
    % Background AWGN
    %--------------------------------
    sigmax2=var(S1);
    sigman2=sigmax2*10^(-SNRdB/10);
    %--------------------------------
    % Reference signals in the center of the array
    %--------------------------------
    % In case of analitical (real) signals
    %--------------------------------
    % xcentro1=hilbert(S1(:,1));         % SoI (centro do arranjo)
    % xcentro2=hilbert(S2(:,1));         % Interferidor (centro do arranjo)
    %--------------------------------
    xcentro1=S1;
    xcentro2=S2;
    %--------------------------------
    % Steering Vectors
    %--------------------------------
    A1=zeros(M,1);
    A2=zeros(M,1);
    for m=1:M
        A1(m,1)=exp(1i*(KR)*sin(tetasim1)*cos(-(2*pi*(m-1)/M)+phisim1));
        A2(m,1)=exp(1i*(KR)*sin(tetasim2)*cos(-(2*pi*(m-1)/M)+phisim2));
    end
    %--------------------------------
    % Beamforming Matrix F (Phase Modes Matrix) (IDFT) generated based on Nr of
    % sensors and phase modes
    %--------------------------------
    F=zeros(MODES,M);
    w=exp(1i*2*pi/M);
    j=zeros(MODES,1);
    for n=1:(MODES)
        for m=1:M
            F(n,m)= (1/sqrt(M))*w^((m-1)*(-N-1+n));
        end
        j(n,1)= 1/(1i^(-N-1+n)*besselj((-N-1+n),KR*sin(zenini)));
    end
    %-------------------------------
    % Matrix J isolates the azimuth-dependent portion of the array manifold
    % vectors
    %-------------------------------
    J=diag(j);
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % TOEPLITZ CONDITIONING AND FBA
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % Sample Covariance Matrices
    %--------------------------------
    Rcorr=zeros(M,M);
    RcorrPHASE=zeros(MODES,MODES);
    for k=1:SampleSupport
        %--------------------------------
        % Input AWGN dos sensores
        %--------------------------------
        noise=randn(M,1)+1i*randn(M,1);
        noise=noise-mean(noise); noise=sqrt(sigman2)*noise/std(noise);
        %noise=hilbert(noise(:,1));
        %--------------------------------
        % Element-Space Sensor Output
        %--------------------------------
        snapshot=A1*xcentro1(k+1000)+A2*xcentro2(k+1000)+noise;
        %--------------------------------
        % Beamspace Sensor Output
        %--------------------------------
        snapshotPHASE=F*snapshot;
        %--------------------------------
        % TOEPLITZ CONDIDIONING - Sensor Output in BEAMSPACE with
        % azimuth parameter only (Vandermonde structure)
        %--------------------------------
        snapshotPHASE=J*snapshotPHASE;
        %--------------------------------
        Rcorr=Rcorr+snapshot*snapshot';
        RcorrPHASE=RcorrPHASE+snapshotPHASE*snapshotPHASE';
    end
    %----------------------------------
    % Element Space and Beamspace Spatial Covariance Matrices
    %----------------------------------
    Rcorr=Rcorr/SampleSupport;
    RcorrPHASE=RcorrPHASE/SampleSupport;
    %--------------------------------
    % Diagonal Loading
    %--------------------------------
    Rcorr=Rcorr+1*(eye(M));
    RcorrPHASE=RcorrPHASE+1*(eye(MODES));
    %--------------------------------
    % Forward-Backward Avereging (FBA)
    %--------------------------------
    Rcorr=(Rcorr+(flip(eye(M))*(conj(Rcorr))*flip(eye(M))))/2;
    RcorrPHASE=(RcorrPHASE+(flip(eye(MODES))*(conj(RcorrPHASE))*flip(eye(MODES))))/2;
    %--------------------------------
    % PREWHITENING
    %--------------------------------

    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % 2-D DoA in 2 Phases: Root-MUSIC for Azimuth and Angular Scan Method
    % (Element-space Spectral MUSIC) for Zenith
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    [EPHASE,LPHASE]=eig(RcorrPHASE);
    [val,pos] = sort(diag(LPHASE),'descend');
    LPHASE=diag(val);
    EPHASE(:,:)=EPHASE(:,pos);
    % ----------------------------------------------------------------
    % BEAMSPACE Signal and Noise Sub-Spaces
    % ----------------------------------------------------------------
    EsPHASE=EPHASE(:,1:D);
    EnPHASE=EPHASE(:,D+1:MODES);
    EnEnHPHASE=EnPHASE*EnPHASE';
    %-------------------------------------------------------------------------
    % Azimuth DoA with ROOT-MUSIC (PHI)
    %-------------------------------------------------------------------------
    % Signal Subspace Covariance Matrix
    %----------------------------------
    RTS=EnEnHPHASE;
    %----------------------------------
    Coefsup=zeros(1,MODES-1);
    Coefinf=zeros(1,MODES-1);
    for i=1:MODES-1
        Coefsup(i) = sum(diag(RTS,i));
    end

    for i=1:MODES-1
        Coefinf(i) = sum(diag(RTS,-i));
    end

    CoefTotal = [fliplr(Coefinf) sum(diag(RTS)) Coefsup];
    CoefTotal = fliplr(CoefTotal);
    %--------------------------------------------------------------------------
    % ROOT-MUSIC Polynomial Rooting
    %--------------------------------------------------------------------------
    rootMUSIC= roots(CoefTotal);
    DOA=rad2deg(angle(rootMUSIC));
    MAG=abs(rootMUSIC);
    [val,pos] = sort(MAG,'descend');
    DOA(:)=DOA(pos);
    ROOTS=[MAG DOA];
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % SMALL SAMPLE SUPPORT - LEAKAGE MINIMIZATION
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    DOAini1=abs(DOA(MODES,1));       % Initial DoAs are the first solutions of
    DOAini2=abs(DOA(MODES+1,1));     % Root-MUSIC without this procedure
    %--------------------------------
    % Array Manifold with Initial DoAs
    %--------------------------------
    aini1=zeros(M,1);
    aini2=zeros(M,1);
    for m=1:M
        aini1(m,1)=exp(1i*(KR)*sin(zenini)*cos(-(2*pi*(m-1)/M)+deg2rad(DOAini1)));
        aini2(m,1)=exp(1i*(KR)*sin(zenini)*cos(-(2*pi*(m-1)/M)+deg2rad(DOAini2)));
    end
    Aini = [aini1 aini2];    % Array Manifold with Initial DoAs
    plotazini1(y)=DOAini1;
    plotazini2(y)=DOAini2;
    %--------------------------------
    P = Aini*pinv(Aini'*Aini)*(Aini'); % Projector matrix P element-space
    W=F.';
    PPHASE = W'*Aini*pinv(Aini'*(W)*(W')*Aini)*(Aini'*W); % Projector matrix P beamspace
    %--------------------------------
    Tx = P*Rcorr*(eye(M)-P);       % Leakage Minimization in element-space
    TxPHASE = PPHASE*RcorrPHASE*(eye(MODES)-PPHASE);
    %--------------------------------
    % Updated Sample Covariance Matrix with LEAKAGE MINIMIZATION
    %--------------------------------
    Rcorrup = Rcorr-mu*(Tx + Tx');                      % Element Space
    RcorrPHASEup = RcorrPHASE-mu*(TxPHASE + TxPHASE');  % BEAMSPACE
    %--------------------------------
    % Forward-Backward Avereging (FBA)
    %--------------------------------
    Rcorrup=(Rcorrup+(flip(eye(M))*(conj(Rcorrup))*flip(eye(M))))/2;
    RcorrPHASEup=(RcorrPHASEup+(flip(eye(MODES))*(conj(RcorrPHASEup))*flip(eye(MODES))))/2;
    %--------------------------------
    % Subspaces in BEAMSPACE after LEAKAGE MINIMIZATION (up)
    %----------------------------------------------------------------
    [EPHASEup,LPHASEup]=eig(RcorrPHASEup);
    [val,pos] = sort(diag(LPHASEup),'descend');
    LPHASEup=diag(val);
    EPHASEup(:,:)=EPHASEup(:,pos);
    EsPHASEup=EPHASEup(:,1:D);
    EnPHASEup=EPHASEup(:,D+1:MODES);
    EnEnHPHASEup=EnPHASEup*EnPHASEup';
    %-------------------------------------------------------------------------
    %  ROOT-MUSIC para AZIMUTE (PHI) após LEAKAGE MINIMIZATION (up)
    %-------------------------------------------------------------------------
    RTSup=EnEnHPHASEup;
    % --------------------------------------------
    Coefsupup=zeros(1,MODES-1);
    Coefinfup=zeros(1,MODES-1);
    for i=1:MODES-1
        Coefsupup(i) = sum(diag(RTSup,i));
    end
    for i=1:MODES-1
        Coefinfup(i) = sum(diag(RTSup,-i));
    end
    CoefTotalup = [fliplr(Coefinfup) sum(diag(RTSup)) Coefsupup];
    CoefTotalup = fliplr(CoefTotalup);
    %--------------------------------------------------------------------------
    % Calculando as raízes do polinômio ROOT-MUSIC após LEAKAGE MINIMIZATION (up)
    %--------------------------------------------------------------------------
    rootMUSICup= roots(CoefTotalup);
    DOAup=rad2deg(angle(rootMUSICup));
    MAGup=abs(rootMUSICup);
    [val,pos] = sort(MAGup,'descend');
    DOAup(:)=DOAup(pos);
    ROOTSup=[MAGup DOAup];

    % DOAaz1=abs(DOAup(MODES,1));
    % DOAaz2=abs(DOAup(MODES+1,1));
    % ------------------------------------
    DOAaz1A=abs(DOAup(MODES,1));
    DOAaz2A=abs(DOAup(MODES+1,1));
    % ------------------------------------
    E1=abs(DOAaz1A-rad2deg(phisim1));
    E1A=abs(DOAaz2A-rad2deg(phisim1));
    if E1<E1A
        DOAaz1=DOAaz1A;
        DOAaz2=DOAaz2A;
    else
        DOAaz1=DOAaz2A;
        DOAaz2=DOAaz1A;
    end
    %------------------------------------
    plotaz1(y)=DOAaz1;
    plotaz2(y)=DOAaz2;
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % RESOLUÇÃO EM ZÊNITE: MUSIC Espectral em Element Space and Beamspace
    % em 1-D no azimute da solução Root-MUSIC
    %--------------------------------------------------------------------------
    %-------------------------------------------------------------------------
    % Subspaces in ELEMENT SPACE
    %----------------------------------------------------------------
    [E,L]=eig(Rcorrup);
    [val,pos] = sort(diag(L),'descend');
    L=diag(val);
    E(:,:)=E(:,pos);
    Esup=E(:,1:D);
    Enup=E(:,D+1:M);
    EnEnHup=Enup*Enup';
    %----------------------------------------------------------------
    PMUSICSoIze=zeros(length(teta));
    aSoI=zeros(M,1);
    PMUSICIntze=zeros(length(teta));
    aInt=zeros(M,1);
    %------------------------------------
    for contteta=1:length(teta)
        for m=1:M
            aSoI(m,1)= exp(1i*(KR)*sin(teta(contteta))*cos(-(2*pi*(m-1)/M)+deg2rad(DOAaz1)));
            aInt(m,1)= exp(1i*(KR)*sin(teta(contteta))*cos(-(2*pi*(m-1)/M)+deg2rad(DOAaz2)));
        end
        %--------------------------------
        PMUSICSoIze(contteta)=1/(aSoI'*EnEnHup*aSoI);
        PMUSICIntze(contteta)=1/(aInt'*EnEnHup*aInt);
    end
    %--------------------------------
    PMUSICSoIze=abs(PMUSICSoIze)/max(abs(PMUSICSoIze));
    PMUSICIntze=abs(PMUSICIntze)/max(abs(PMUSICIntze));
    %--------------------------------
    [M1,DOAze1]=max(PMUSICSoIze);    % Zenith DoA obtained from the maximum of
    [M2,DOAze2]=max(PMUSICIntze);    % MUSIC Pseudo-spectrum
    %--------------------------------
    DOAze1A=DOAze1/10;
    DOAze2A=DOAze2/10;
    %------------------------------------
    E1=abs(DOAze1A-rad2deg(tetasim1));
    E1A=abs(DOAze2A-rad2deg(tetasim1));
    if E1<E1A
        DOAze1=DOAze1A;
        DOAze2=DOAze2A;
    else
        DOAze1=DOAze2A;
        DOAze2=DOAze1A;
    end
    %------------------------------------
    plotze1(y)=DOAze1;
    plotze2(y)=DOAze2;
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % BEAMSPACE MVDR Beamformer for 2-D DoA
    % 
    % Use BEAMSPACE Root-MUSIC and ELEMENT SPACE Spectral MUSIC to construct the MVDR Beamformer
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % PMVDRPHASE=zeros(length(phi),length(teta));
    % aMVDRPHASE=zeros(MODES,1);
    % for contphi=1:length(phi)
    %     for contteta=1:length(teta)
    %         for n=1:MODES
    % 
    %             aMVDRPHASE(n,1)= 1i^(-N-1+n)*besselj((-N-1+n),KR*sin(teta(contteta)))*exp(1i*((-N-1+n))*(phi(contphi)));
    % 
    %         end
    % 
    %         PMVDRPHASE(contphi,contteta)=1/(aMVDRPHASE'*pinv(RcorrPHASEup)*aMVDRPHASE);
    % 
    %     end
    % end
    % PMVDRPHASE(:,:)=abs(PMVDRPHASE(:,:))/max(max(abs(PMVDRPHASE)));
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % Element Space MVDR Beamformer for 2-D DoA
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % PMVDR=zeros(length(phi),length(teta));
    % aMVDR=zeros(M,1);
    % for contphi=1:length(phi)
    %     for contteta=1:length(teta)
    %         for m=1:M
    %             aMVDR(m,1)=[exp(1i*(KR)*sin(teta(contteta))*cos(-(2*pi*(m-1)/M)+phi(contphi)))];
    %         end
    %         PMVDR(contphi,contteta)=1/(aMVDR'*pinv(Rcorrup)*aMVDR);
    % 
    %     end
    % end
    % PMVDR(:,:)=abs(PMVDR(:,:))/max(max(abs(PMVDR)));
end
%--------------------------------
% Results metrics
%--------------------------------
meanazini1=mean(plotazini1);
meanazini2=mean(plotazini2);
meanaz1=mean(plotaz1);
%stdaz1=std(plotaz1);
meanze1=mean(plotze1);
%stdze1=std(plotze1);
meanaz2=mean(plotaz2);
%stdaz2=std(plotaz2);
meanze2=mean(plotze2);
%stdze2=std(plotze2);
%--------------------------------
% Consolidated RMSE from all montecarlo trials 
%--------------------------------
E1AZini=abs(rad2deg(phisim1)-meanazini1)  % Initial Error-Azimuth-SoI
E2AZini=abs(rad2deg(phisim2)-meanazini2)  % Initial Error-Azimuth-Interferer
E1AZ=abs(rad2deg(phisim1)-meanaz1)  % Error-Azimuth-SoI
E2AZ=abs(rad2deg(phisim2)-meanaz2)  % Error-Azimuth-Interferer
E1ZE=abs(rad2deg(tetasim1)-meanze1) % Error-Zenith-SoI
E2ZE=abs(rad2deg(tetasim2)-meanze2) % Error-Zenith-Interferer
%--------------------------------
% Vectors of RMSE from each montecarlo trials 
%--------------------------------
plotaz1;
plotaz2;
plotze1;
plotze2;
