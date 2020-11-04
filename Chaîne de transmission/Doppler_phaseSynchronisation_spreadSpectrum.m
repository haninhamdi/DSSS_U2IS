clc;
close all;
clear all;

% Définir les paramètres de simulation 
trellis = poly2trellis(3,{'1 + x^2','1 + x + x^2'}); 
RSB=10;
N0 = 10^(-RSB/10);
sigma_bruit = sqrt(N0);
Nb = 100000 ;
Nb_p = 500 ;
ksi = 0.01 ;  
facteur_etalement = 8 ;
Ts = 32 ;
fs = 1/Ts ;
Tc = Ts/facteur_etalement;
fc = 1/Tc;
f0 = 0.4 ;
T0 = 1/f0 ;
alpha = 0.22 ;
Fe = f0*3 ;
Seuil=0 ;
iterations = 3;
mu_DA = (ksi/N0)^(3/2) ;
doppler = 0.003 ;
sigma_phase = sqrt(0.005);  % A changer

for N0_bruit = -5:25
    % Générer aléatroirement les bits 
    data = randi([0,1], Nb, 1);
    data_p = randi([0,1], Nb_p, 1);
    
    % Codage canal 
    dataEnc = convenc(data,trellis);
    dataEnc_p = convenc(data_p,trellis);
    
    % Modulation BPSK de dataEnc : 1 et -1 
    for j = 1 : length(dataEnc)
        if dataEnc(j) == 0
            dataEnc(j) = -1;
        end
    end
    
    for j = 1 : length(dataEnc_p)
        if dataEnc_p(j) == 0
            dataEnc_p(j) = -1;
        end
    end
    
    % Etalement de spectre 
    [spreadData,spreadSequence] = EtalementDeSpectre(facteur_etalement,dataEnc,length(dataEnc));
    [spreadData_p,spreadSequence_p] = EtalementDeSpectre(facteur_etalement,dataEnc_p,length(dataEnc_p));
    
    % conversion en un signal physique
    N = length(spreadData);
    peigne = zeros(1,N*Tc) ;
    for i= 0: N-1
        peigne(1+i*Tc) = spreadData(1+i);
    end
    N_p = length(spreadData_p);
    peigne_p = zeros(1,N_p*Tc) ;
    for i= 0: N_p-1
        peigne_p(1+i*Tc) = spreadData_p(1+i);
    end

    % Construction du filtre d'émission en racine de cosinus suréléve 
    axe_temps = -10*Tc:10*Tc;
    for i=1:length(axe_temps)
        filtre_emission(i) = cosSurreleve(alpha,axe_temps(i),Tc);
    end

    % effectuer la convolution entre le signal (peigne) et la fonction du
    % filtre
    signal_filtre = conv(filtre_emission,peigne);
    signal_filtre_p = conv(filtre_emission,peigne_p);

    % effectuer le produit entre le signal réel et la porteuse et ajouter l'effet Doppler 
    N2 = length(signal_filtre);
    axe_temps1 = -N2/2:(N2-1)/2;
    porteuse = cos(2*pi*(f0+doppler)*axe_temps1);
    signal_transmis = signal_filtre.*porteuse;
    signal_transmis_sans_doppler = signal_filtre.*cos(2*pi*f0*axe_temps1);
    
    N2_p = length(signal_filtre_p);
    axe_temps1_p = -N2_p/2:(N2_p-1)/2;
    porteuse_p = cos(2*pi*(f0+doppler)*axe_temps1_p);
    signal_transmis_p = signal_filtre_p.*porteuse_p;
    signal_transmis_p_sans_doppler = signal_filtre_p.*cos(2*pi*f0*axe_temps1_p);
    
    % ajout d'un déphasage 
    
    bruit_phase = sigma_phase * randn(1,length(signal_transmis)+length(signal_transmis_p)-1)*pi/180;
    dephasage = zeros(1,length(signal_transmis)+length(signal_transmis_p));
    for k1 = 2 : length(dephasage)
        c = dephasage(k1-1) +ksi + bruit_phase(k1-1);
        if (c>=-pi/4)&&(c<=pi/4) dephasage(k1)=c;
        elseif c>pi/4 dephasage(k1) = c - dephasage(k1-1)/4;
        elseif c<-pi/4 dephasage(k1) = c + dephasage(k1-1)/4;
        end
        %elseif c>pi/4 dephasage(k1) = pi/4;
        %elseif c<-pi/4 dephasage(k1) = -pi/4;
        %end 
    end
    erreurDePhase = exp(1i*dephasage);

    % ajout d'un BBAG  
    bruit = sigma_bruit * randn(1,length(signal_transmis));
    signal_canal_non_deph = signal_transmis + bruit;
    signal_canal_non_deph_sans_doppler = signal_transmis_sans_doppler + bruit;
    signal_canal_deph = signal_transmis.*erreurDePhase(length(signal_transmis_p)+1:end) + bruit;
    signal_canal_deph_sans_doppler = signal_transmis_sans_doppler.*erreurDePhase(length(signal_transmis_p)+1:end) + bruit;
    
    bruit_p = sigma_bruit * randn(1,length(signal_transmis_p));
    signal_canal_p_non_deph = signal_transmis_p + bruit_p;
    signal_canal_p_deph = signal_transmis_p.*erreurDePhase(1:length(signal_transmis_p)) + bruit_p;
    signal_canal_p_non_deph_sans_doppler = signal_transmis_p_sans_doppler + bruit_p;
    signal_canal_p_deph_sans_doppler = signal_transmis_p_sans_doppler.*erreurDePhase(1:length(signal_transmis_p)) + bruit_p;
    
    % Bloc pour le brouillage 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Pb = 10^(N0_bruit/10);
    sigma_brouilleur = sqrt(Pb);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Brouillage gaussien
    [b,a] = butter(8,fc) ; 
    
    brouilleur_gaussien = 1.2*sigma_bruit*randn(1,length(signal_canal_deph));
    brouilleur_gaussien = sigma_brouilleur*filter(b,a,brouilleur_gaussien); 
    signal_recu_gaus_non_deph = signal_canal_non_deph + brouilleur_gaussien.*porteuse;
    signal_recu_gaus_deph = signal_canal_deph + brouilleur_gaussien.*porteuse;
    signal_recu_gaus_non_deph_sans_doppler = signal_canal_non_deph_sans_doppler + brouilleur_gaussien.*cos(2*pi*f0*axe_temps1);
    signal_recu_gaus_deph_sans_doppler = signal_canal_deph_sans_doppler + brouilleur_gaussien.*cos(2*pi*f0*axe_temps1);
    
    brouilleur_gaussien_p = 1.2*sigma_bruit*randn(1,length(signal_canal_p_deph));
    brouilleur_gaussien_p = sigma_brouilleur*filter(b,a,brouilleur_gaussien_p);
    signal_recu_gaus_p_non_deph = signal_canal_p_non_deph + brouilleur_gaussien_p.*porteuse_p;
    signal_recu_gaus_p_deph = signal_canal_p_deph + brouilleur_gaussien_p.*porteuse_p;
    signal_recu_gaus_p_non_deph_sans_doppler = signal_canal_p_non_deph_sans_doppler + brouilleur_gaussien_p.*cos(2*pi*f0*axe_temps1_p);
    signal_recu_gaus_p_deph_sans_doppler = signal_canal_p_deph_sans_doppler + brouilleur_gaussien_p.*cos(2*pi*f0*axe_temps1_p);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % estimation de l'effet Doppler à partir du préambule
    epsilon = 0.001 ;
    t = -10:10 ; 
    test = 0;
    correlation = epsilon ;
    
    %figure
    %hold on
    for j=1:11
        testPorteuse = cos(2*pi*(f0+(j-1)*epsilon)*axe_temps1_p);
        signal_gaus = conv(fliplr(signal_filtre_p.*testPorteuse),signal_recu_gaus_p_non_deph);
        
        if test < max(signal_gaus)
            test = max(signal_gaus);
            correlation = (j-1)*epsilon ;
        end
        %axe_freq1= -Fe/2 : Fe/length(Y): Fe/2-Fe/length(Y);
        %chr = ['ksi = ',num2str(epsilon*(j-1))];
        %plot(signal_gaus,'DisplayName',chr);
        %grid on 
    end
     %legend show
     %hold off
    
    porteuse_correle = cos(2*pi*(f0+doppler)*axe_temps1);
    porteuse_p_correle = cos(2*pi*(f0+doppler)*axe_temps1_p);
    
    signal_sans_porteuse_gaus_non_deph = signal_recu_gaus_non_deph.*porteuse_correle;
    signal_sans_porteuse_gaus_p_non_deph = signal_recu_gaus_p_non_deph.*porteuse_p_correle;
    signal_sans_porteuse_gaus_non_deph_sans_doppler = signal_recu_gaus_non_deph_sans_doppler.*cos(2*pi*f0*axe_temps1);
    signal_sans_porteuse_gaus_p_non_deph_sans_doppler = signal_recu_gaus_p_non_deph_sans_doppler.*cos(2*pi*f0*axe_temps1_p);
    
    porteuseQ = -sin(2*pi*(f0+correlation)*axe_temps1);
    porteuseQ_p = -sin(2*pi*(f0+correlation)*axe_temps1_p);
    
    signal_sans_porteuse_gaus_deph_I = 2*signal_recu_gaus_deph.*porteuse_correle;
    signal_sans_porteuse_gaus_deph_Q = 2*signal_recu_gaus_deph.*porteuseQ;
    signal_sans_porteuse_gaus_deph_sans_synch = signal_sans_porteuse_gaus_deph_I+1i*signal_sans_porteuse_gaus_deph_Q;
    beta_estime_NDA = NDA(signal_sans_porteuse_gaus_deph_sans_synch,mu_DA,iterations,N0);
    signal_sans_porteuse_gaus_deph = signal_sans_porteuse_gaus_deph_sans_synch.*exp(-1i*beta_estime_NDA*1.1);
    
    signal_sans_porteuse_gaus_deph_I_sans_doppler = 2*signal_recu_gaus_deph_sans_doppler.*cos(2*pi*f0*axe_temps1);
    signal_sans_porteuse_gaus_deph_Q_sans_doppler = -2*signal_recu_gaus_deph_sans_doppler.*sin(2*pi*f0*axe_temps1);
    signal_sans_porteuse_gaus_deph_sans_synch_sans_doppler = signal_sans_porteuse_gaus_deph_I_sans_doppler+1i*signal_sans_porteuse_gaus_deph_Q_sans_doppler;
    beta_estime_NDA_sans_doppler = NDA(signal_sans_porteuse_gaus_deph_sans_synch_sans_doppler,mu_DA,iterations,N0);
    signal_sans_porteuse_gaus_deph_sans_doppler = signal_sans_porteuse_gaus_deph_sans_synch_sans_doppler.*exp(-1i*beta_estime_NDA_sans_doppler*1.1);
    
    signal_sans_porteuse_gaus_p_deph_I = 2*signal_recu_gaus_p_deph.*porteuse_p_correle;
    signal_sans_porteuse_gaus_p_deph_Q = 2*signal_recu_gaus_p_deph.*porteuseQ_p;
    signal_sans_porteuse_gaus_p_deph_sans_synch = signal_sans_porteuse_gaus_p_deph_I+1i*signal_sans_porteuse_gaus_p_deph_Q;
    beta_estime_DA = DA(signal_sans_porteuse_gaus_p_deph_sans_synch,signal_transmis_p,mu_DA,iterations);
    signal_sans_porteuse_gaus_p_deph = signal_sans_porteuse_gaus_p_deph_sans_synch.*exp(-1i*beta_estime_DA);
    
    signal_sans_porteuse_gaus_p_deph_I_sans_doppler = 2*signal_recu_gaus_p_deph_sans_doppler.*cos(2*pi*f0*axe_temps1_p);
    signal_sans_porteuse_gaus_p_deph_Q_sans_doppler = -2*signal_recu_gaus_p_deph_sans_doppler.*sin(2*pi*f0*axe_temps1_p);
    signal_sans_porteuse_gaus_p_deph_sans_synch_sans_doppler = signal_sans_porteuse_gaus_p_deph_I_sans_doppler+1i*signal_sans_porteuse_gaus_p_deph_Q_sans_doppler;
    beta_estime_DA_sans_doppler = DA(signal_sans_porteuse_gaus_p_deph_sans_synch_sans_doppler,signal_transmis_p_sans_doppler,mu_DA,iterations);
    signal_sans_porteuse_gaus_p_deph_sans_doppler = signal_sans_porteuse_gaus_p_deph_sans_synch_sans_doppler.*exp(-1i*beta_estime_DA_sans_doppler);
    
    
    % filtre à la réception adapté (en racine de cosinus surréléve)
    for j = 1:length(filtre_emission)
        filtre_reception(j)= cosSurreleve(alpha,Tc-axe_temps(j),Tc);
    end
    
    % filtrer le signal recu 
    signal_filtre_adapte_gaus_non_deph = conv(filtre_reception,signal_sans_porteuse_gaus_non_deph);
    signal_filtre_adapte_gaus_p_non_deph = conv(filtre_reception,signal_sans_porteuse_gaus_p_non_deph);
    
    signal_filtre_adapte_gaus_non_deph_sans_doppler = conv(filtre_reception,signal_sans_porteuse_gaus_non_deph_sans_doppler);
    signal_filtre_adapte_gaus_p_non_deph_sans_doppler = conv(filtre_reception,signal_sans_porteuse_gaus_p_non_deph_sans_doppler);

    signal_filtre_adapte_gaus_deph = conv(filtre_reception,signal_sans_porteuse_gaus_deph);
    signal_filtre_adapte_gaus_p_deph = conv(filtre_reception,signal_sans_porteuse_gaus_p_deph);
    
    signal_filtre_adapte_gaus_deph_sans_doppler = conv(filtre_reception,signal_sans_porteuse_gaus_deph_sans_doppler);
    signal_filtre_adapte_gaus_p_deph_sans_doppler = conv(filtre_reception,signal_sans_porteuse_gaus_p_deph_sans_doppler);
    
    signal_filtre_adapte_gaus_deph_sans_synch = conv(filtre_reception,signal_sans_porteuse_gaus_deph_sans_synch);
    signal_filtre_adapte_gaus_p_deph_sans_synch = conv(filtre_reception,signal_sans_porteuse_gaus_p_deph_sans_synch);
    
    signal_filtre_adapte_gaus_deph_sans_synch_sans_doppler = conv(filtre_reception,signal_sans_porteuse_gaus_deph_sans_synch_sans_doppler);
    signal_filtre_adapte_gaus_p_deph_sans_synch_sans_doppler = conv(filtre_reception,signal_sans_porteuse_gaus_p_deph_sans_synch_sans_doppler);
    
    % Décoder le signal reçu filtré
    k=0 ;
    for i=1 : length(spreadData)
        k=k+1 ;
        m_gaus_non_deph=real(signal_filtre_adapte_gaus_non_deph(i*Tc+length(filtre_reception))) ;
        m_gaus_non_deph_sans_doppler=real(signal_filtre_adapte_gaus_non_deph_sans_doppler(i*Tc+length(filtre_reception))) ;
        
        if m_gaus_non_deph>= Seuil
            bit_estime_gaus_non_deph(k)=1 ;
        else
            bit_estime_gaus_non_deph(k)=-1 ;
        end
        
        if m_gaus_non_deph_sans_doppler>= Seuil
            bit_estime_gaus_non_deph_sans_doppler(k)=1 ;
        else
            bit_estime_gaus_non_deph_sans_doppler(k)=-1 ;
        end
    end
    
    k=0 ;
    for i=1 : length(spreadData)
        k=k+1 ;
        m_gaus_deph=signal_filtre_adapte_gaus_deph(i*Tc+length(filtre_reception)) ;
        m_gaus_deph_sans_synch = signal_filtre_adapte_gaus_deph_sans_synch(i*Tc+length(filtre_reception)) ;
        
        m_gaus_deph_sans_doppler=signal_filtre_adapte_gaus_deph_sans_doppler(i*Tc+length(filtre_reception)) ;
        m_gaus_deph_sans_synch_sans_doppler = signal_filtre_adapte_gaus_deph_sans_synch_sans_doppler(i*Tc+length(filtre_reception)) ;
        
        if m_gaus_deph>= Seuil
            bit_estime_gaus_deph(k)=1 ;
        else
            bit_estime_gaus_deph(k)=-1 ;
        end 
        
        if m_gaus_deph_sans_doppler>= Seuil
            bit_estime_gaus_deph_sans_doppler(k)=1 ;
        else
            bit_estime_gaus_deph_sans_doppler(k)=-1 ;
        end 
        
        if m_gaus_deph_sans_synch>= Seuil
            bit_estime_gaus_deph_sans_synch(k)=1 ;
        else
            bit_estime_gaus_deph_sans_synch(k)=-1 ;
        end 
        
        if m_gaus_deph_sans_synch_sans_doppler>= Seuil
            bit_estime_gaus_deph_sans_synch_sans_doppler(k)=1 ;
        else
            bit_estime_gaus_deph_sans_synch_sans_doppler(k)=-1 ;
        end 
    end
    
    k=0 ;
    for i=1 : length(spreadData_p)
        k=k+1 ;
        m_gaus_p_non_deph =signal_filtre_adapte_gaus_p_non_deph(i*Tc+length(filtre_reception)) ;
        m_gaus_p_non_deph_sans_doppler =signal_filtre_adapte_gaus_p_non_deph_sans_doppler(i*Tc+length(filtre_reception)) ;
        
        if m_gaus_p_non_deph>= Seuil
            bit_estime_gaus_p_non_deph(k)=1 ;
        else
            bit_estime_gaus_p_non_deph(k)=-1 ;
        end
        
        if m_gaus_p_non_deph_sans_doppler>= Seuil
            bit_estime_gaus_p_non_deph_sans_doppler(k)=1 ;
        else
            bit_estime_gaus_p_non_deph_sans_doppler(k)=-1 ;
        end
    end
    
    k=0 ;
    for i=1 : length(spreadData_p)
        k=k+1 ;
        m_gaus_p_deph=signal_filtre_adapte_gaus_p_deph(i*Tc+length(filtre_reception)) ;
        m_gaus_p_deph_sans_synch=signal_filtre_adapte_gaus_p_deph_sans_synch(i*Tc+length(filtre_reception)) ;
        
        m_gaus_p_deph_sans_doppler=signal_filtre_adapte_gaus_p_deph_sans_doppler(i*Tc+length(filtre_reception)) ;
        m_gaus_p_deph_sans_synch_sans_doppler=signal_filtre_adapte_gaus_p_deph_sans_synch_sans_doppler(i*Tc+length(filtre_reception)) ;
        
        if m_gaus_p_deph>= Seuil
            bit_estime_gaus_p_deph(k)=1 ;
        else
            bit_estime_gaus_p_deph(k)=-1 ;
        end 
        if m_gaus_p_deph_sans_synch>= Seuil
            bit_estime_gaus_p_deph_sans_synch(k)=1 ;
        else
            bit_estime_gaus_p_deph_sans_synch(k)=-1 ;
        end 
        
        if m_gaus_p_deph_sans_doppler>= Seuil
            bit_estime_gaus_p_deph_sans_doppler(k)=1 ;
        else
            bit_estime_gaus_p_deph_sans_doppler(k)=-1 ;
        end 
        if m_gaus_p_deph_sans_synch_sans_doppler>= Seuil
            bit_estime_gaus_p_deph_sans_synch_sans_doppler(k)=1 ;
        else
            bit_estime_gaus_p_deph_sans_synch_sans_doppler(k)=-1 ;
        end
    end
    
    
    bit_estime_gaus_deph = real(bit_estime_gaus_deph);
    bit_estime_gaus_p_deph = real(bit_estime_gaus_p_deph);
    
    bit_estime_gaus_deph_sans_doppler = real(bit_estime_gaus_deph_sans_doppler);
    bit_estime_gaus_p_deph_sans_doppler = real(bit_estime_gaus_p_deph_sans_doppler);
    
    bit_estime_gaus_sans_synch = real(bit_estime_gaus_deph_sans_synch); 
    bit_estime_gaus_p_sans_synch = real(bit_estime_gaus_p_deph_sans_synch);
    
    bit_estime_gaus_sans_synch_sans_doppler = real(bit_estime_gaus_deph_sans_synch_sans_doppler); 
    bit_estime_gaus_p_sans_synch_sans_doppler = real(bit_estime_gaus_p_deph_sans_synch_sans_doppler);

    % Inverser l'étalement de spectre 
    Data_estime_gaus_non_deph = InverseEtalementDeSpectre(bit_estime_gaus_non_deph,spreadSequence,facteur_etalement,length(dataEnc));
    Data_estime_gaus_p_non_deph = InverseEtalementDeSpectre(bit_estime_gaus_p_non_deph,spreadSequence_p,facteur_etalement,length(dataEnc_p));
    
    Data_estime_gaus_non_deph_sans_doppler = InverseEtalementDeSpectre(bit_estime_gaus_non_deph_sans_doppler,spreadSequence,facteur_etalement,length(dataEnc));
    Data_estime_gaus_p_non_deph_sans_doppler = InverseEtalementDeSpectre(bit_estime_gaus_p_non_deph_sans_doppler,spreadSequence_p,facteur_etalement,length(dataEnc_p));
    
    Data_estime_gaus_deph_sans_doppler = InverseEtalementDeSpectre(bit_estime_gaus_deph_sans_doppler,spreadSequence,facteur_etalement,length(dataEnc));
    Data_estime_gaus_p_deph_sans_doppler = InverseEtalementDeSpectre(bit_estime_gaus_p_deph_sans_doppler,spreadSequence_p,facteur_etalement,length(dataEnc_p));
    
    Data_estime_gaus_deph = InverseEtalementDeSpectre(bit_estime_gaus_deph,spreadSequence,facteur_etalement,length(dataEnc));
    Data_estime_gaus_p_deph = InverseEtalementDeSpectre(bit_estime_gaus_p_deph,spreadSequence_p,facteur_etalement,length(dataEnc_p));
    
    Data_estime_gaus_sans_synch = InverseEtalementDeSpectre(bit_estime_gaus_sans_synch,spreadSequence,facteur_etalement,length(dataEnc));
    Data_estime_gaus_p_sans_synch = InverseEtalementDeSpectre(bit_estime_gaus_p_sans_synch,spreadSequence_p,facteur_etalement,length(dataEnc_p));
    
    Data_estime_gaus_sans_synch_sans_doppler = InverseEtalementDeSpectre(bit_estime_gaus_sans_synch_sans_doppler,spreadSequence,facteur_etalement,length(dataEnc));
    Data_estime_gaus_p_sans_synch_sans_doppler = InverseEtalementDeSpectre(bit_estime_gaus_p_sans_synch_sans_doppler,spreadSequence_p,facteur_etalement,length(dataEnc_p));
    % Décodage canal 
    tbdepth = 34;
    
    Data_coded_gaus_non_deph = vitdec(Data_estime_gaus_non_deph,trellis,tbdepth,'trunc','hard');
    Data_coded_gaus_p_non_deph = vitdec(Data_estime_gaus_p_non_deph,trellis,tbdepth,'trunc','hard');
    Ne_gaus_non_deph=sum(abs(data-Data_coded_gaus_non_deph))+sum(abs(data_p-Data_coded_gaus_p_non_deph));
    TEB_gaus_non_deph(N0_bruit+6) = Ne_gaus_non_deph/(Nb+Nb_p);
    
    Data_coded_gaus_non_deph_sans_doppler = vitdec(Data_estime_gaus_non_deph_sans_doppler,trellis,tbdepth,'trunc','hard');
    Data_coded_gaus_p_non_deph_sans_doppler = vitdec(Data_estime_gaus_p_non_deph_sans_doppler,trellis,tbdepth,'trunc','hard');
    Ne_gaus_non_deph_sans_doppler=sum(abs(data-Data_coded_gaus_non_deph_sans_doppler))+sum(abs(data_p-Data_coded_gaus_p_non_deph_sans_doppler));
    TEB_gaus_non_deph_sans_doppler(N0_bruit+6) = Ne_gaus_non_deph_sans_doppler/(Nb+Nb_p);
    
    
    
    Data_coded_gaus_deph = vitdec(Data_estime_gaus_deph,trellis,tbdepth,'trunc','hard');
    Data_coded_gaus_p_deph = vitdec(Data_estime_gaus_p_deph,trellis,tbdepth,'trunc','hard');
    Ne_gaus_deph=sum(abs(data-Data_coded_gaus_deph))+sum(abs(data_p-Data_coded_gaus_p_deph));
    TEB_gaus_deph(N0_bruit+6) = Ne_gaus_deph/(Nb+Nb_p);
    
    Data_coded_gaus_deph_sans_doppler = vitdec(Data_estime_gaus_deph_sans_doppler,trellis,tbdepth,'trunc','hard');
    Data_coded_gaus_p_deph_sans_doppler = vitdec(Data_estime_gaus_p_deph_sans_doppler,trellis,tbdepth,'trunc','hard');
    Ne_gaus_deph_sans_doppler=sum(abs(data-Data_coded_gaus_deph_sans_doppler))+sum(abs(data_p-Data_coded_gaus_p_deph_sans_doppler));
    TEB_gaus_deph_sans_doppler(N0_bruit+6) = Ne_gaus_deph_sans_doppler/(Nb+Nb_p);
    
    
    
    Data_coded_gaus_sans_synch = vitdec(Data_estime_gaus_sans_synch,trellis,tbdepth,'trunc','hard');
    Data_coded_gaus_p_sans_synch = vitdec(Data_estime_gaus_p_sans_synch,trellis,tbdepth,'trunc','hard');
    Ne_gaus_sans_synch=sum(abs(data-Data_coded_gaus_sans_synch))+sum(abs(data_p-Data_coded_gaus_p_sans_synch));
    TEB_gaus_sans_synch(N0_bruit+6) = Ne_gaus_sans_synch/(Nb+Nb_p);
    
    Data_coded_gaus_sans_synch_sans_doppler = vitdec(Data_estime_gaus_sans_synch_sans_doppler,trellis,tbdepth,'trunc','hard');
    Data_coded_gaus_p_sans_synch_sans_doppler = vitdec(Data_estime_gaus_p_sans_synch_sans_doppler,trellis,tbdepth,'trunc','hard');
    Ne_gaus_sans_synch_sans_doppler=sum(abs(data-Data_coded_gaus_sans_synch_sans_doppler))+sum(abs(data_p-Data_coded_gaus_p_sans_synch_sans_doppler));
    TEB_gaus_sans_synch_sans_doppler(N0_bruit+6) = Ne_gaus_sans_synch_sans_doppler/(Nb+Nb_p);
end

figure 
semilogy([-5:25],TEB_gaus_deph,'-o r')
hold on
semilogy([-5:25],TEB_gaus_deph_sans_doppler,'-* k')
hold on
semilogy([-5:25],TEB_gaus_non_deph,'-p r')
hold on
semilogy([-5:25],TEB_gaus_non_deph_sans_doppler,'- k')
hold on
semilogy([-5:25],TEB_gaus_sans_synch,'-< r')
hold on
semilogy([-5:25],TEB_gaus_sans_synch_sans_doppler,'-s k')

%title("TEB aprés synchronisation de phase pour ksi = 0.1 (brouilleur gaussien)(Pu=1)") 
xlabel('Pb (dB)')
ylabel('TEB')
%legend('avec dephasage et synch sans doppler','sans dephasage sans doppler','avec déphasage et sans synch sans doppler')
legend('avec dephasage et synch','avec dephasage et synch sans doppler','sans dephasage','sans dephasage sans doppler','avec déphasage et sans synch','avec déphasage et sans synch sans doppler')
grid on