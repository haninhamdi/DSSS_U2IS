%%%%%%%%%%%%%%%%%%%% PRE_U2IS : Quatrième semaine %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Objectif : Etalement de spectre + brouillage %%%%%%%%%%%%%%%%
%*************************************************************************%
% Nb : Le nombre des bits générés
% data : vecteur représentant les bits transmis 
% Ts : fréquence symbole, binaire
% fs : fréquence binaire
% fc : fréquence chip ( de l'étalement de spectre )
% f0 : fréquence porteuse 
% trellis : utilisé pour le codage canal, 2/3 rate
% T0 : période porteuse
%*************************************************************************%
%****************************  Les vecteurs ******************************%
% data : les bits générés initialement 
% dataEnc : les bits obtenus après un codage canal de data 
% S : représente la modulation BPSK du vecteur dataEnc ( en 1 et -1)
% spreadData : le résultat de l'étalement de spectre appliqué à S 
% spreadSequence : la séquence détalement utilisée dans l'algorithme
%                  d'étalement de spectre 
% peigne : peigne de dirac représentant le signal continu construit à
%          partir des bits de spreadData
% filtre_emission : le filtre appliqué à l'émission
% signal_filtre : le résultat du filtrage de peigne par le filtre_emission
% porteuse : le signal de la porteuse
% signal_transmis : le signal transmis par la porteuse : résultat du
%                   produit entre signal_filtre et porteuse 
% bruit : le bruit (son enveloppe complexe) ajouté par le canal 
% signal_canal : le résultat de l'ajout de bruit au signal_transmis
% signal_recu : le signal reçu après bruitage ( pour le moment égal au
%               signal_canal)
% signal_sans_porteuse : le signal réel obtenu après élimination de la
%                         porteuse 
% filtre_reception : filtre adapté à filtre_émission 
% signal_filtre_adapte : signal reçu après filtrage par filtre_reception
% bit_estime : le signal_filtre_adapte décodé en 1 et -1
% Data_estime : résultat de l'inversement de l'étalement de spectre du
%               vecteur bit_estime
% Data_coded : résultat du décodage canal : le résultat final cherché 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
close all;
clear all;

% Définir les paramètres de simulation 
trellis = poly2trellis(3,{'1 + x^2','1 + x + x^2'}); 
RSB=10;
N0 = 10^(-RSB/10);
sigma_bruit = sqrt(N0);
Nb = 50000 ;
Nb_p = 100 ;
ksi = 0.004 ;  
facteur_etalement = 4 ;
Ts = 32 ;
fs = 1/Ts ;
Tc = Ts/facteur_etalement;
fc = 1/Tc;
f0 = 0.4 ;
T0 = 1/f0 ;
alpha = 0.22 ;
Fe = 50*(10^6) ;
Seuil=0 ;
for N0_bruit = 0:0
    % Générer aléatroirement les bits 
    data = randi([0,1], Nb+Nb_p, 1);
    data_p = data(1:Nb_p);
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
    N2_p = length(signal_filtre_p);
    axe_temps1_p = -N2_p/2:(N2_p-1)/2;
    
    N2 = length(signal_filtre);
    axe_temps1 = -N2/2:(N2-1)/2;
    porteuse = cos(2*pi*(f0+ksi)*axe_temps1);
    signal_transmis = signal_filtre.*porteuse;
    
    porteuse_non_correlee = cos(2*pi*f0*axe_temps1);
    signal_transmis_non_correle = signal_filtre.*porteuse_non_correlee;

    % ajout d'un BBAG  
    bruit = sigma_bruit * randn(1,length(signal_transmis));
    signal_canal = signal_transmis + bruit;
    signal_canal_non_correle = signal_transmis_non_correle + bruit;
    
    % Bloc pour le brouillage 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Pb = 10^(N0_bruit/10);
    sigma_brouilleur = sqrt(Pb);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Brouillage gaussien
    [b,a] = butter(8,fc) ; 
    
    brouilleur_gaussien = 1.1*sigma_bruit*randn(1,length(signal_canal));
    brouilleur_gaussien = sigma_brouilleur*filter(b,a,brouilleur_gaussien); 
    signal_recu_gaus = signal_canal + brouilleur_gaussien.*porteuse;
    signal_recu_gaus_non_correle = signal_canal_non_correle + brouilleur_gaussien.*porteuse_non_correlee;
    
    % Brouillage par balayage 
    
    % delta_f = fc/20 ;
    % f = -fc+6*delta_f:delta_f:fc-6*delta_f ; 
    % for i = 1:length(axe_temps1)
    %    somme = 0 ;
    %    for j=1:length(f)
    %        somme = somme + cos(2*pi*f(j)*axe_temps1(i));
    %    end
    %    brouilleur_balayage(i)= sqrt(2*Pb/length(f))*somme;
    % end
   % brouilleur_balayage =  brouilleur_balayage.*porteuse;
    brouilleur_balayage = 1.7*sigma_bruit*randn(1,length(signal_canal));
    brouilleur_balayage = sigma_brouilleur*filter(b,a,brouilleur_balayage); 
    signal_recu_balay = brouilleur_balayage.*porteuse + signal_canal;
    signal_recu_balay_non_correle = brouilleur_balayage.*porteuse_non_correlee + signal_canal_non_correle;
    
    % Brouilleur répéteur
    % générer un retard Tc< tau < Ts aléatoire 
    
     %tau = randi([Tc,Ts]);
     %brouilleur_repeteur = zeros(1,length(signal_canal));
     %for j = tau:length(signal_canal)-1
     %    brouilleur_repeteur(j) = signal_canal(j-tau+1);
     %end
     %brouilleur_repeteur = sigma_brouilleur*brouilleur_repeteur;
     brouilleur_repeteur = 2.3*sigma_bruit*randn(1,length(signal_canal));
     brouilleur_repeteur = sigma_brouilleur*filter(b,a,brouilleur_repeteur); 
     signal_recu_rep = brouilleur_repeteur.*porteuse+ signal_canal;
     signal_recu_rep_non_correle = brouilleur_repeteur.*porteuse_non_correlee+ signal_canal_non_correle;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % estimation de l'effet Doppler à partir du préambule
    epsilon = 0.001 ;
    t = -10:10 ; 
    test = 0;
    correlation = epsilon ;
    
    figure
    hold on
    for j=1:10
        testPorteuse = cos(2*pi*(f0+(j-1)*epsilon)*axe_temps1_p);
        signal_gaus = conv(fliplr(signal_filtre_p.*testPorteuse),signal_recu_gaus(1:length(peigne_p)+length(filtre_emission)-1));
        
        if test < max(signal_gaus)
            test = max(signal_gaus);
            correlation = (j-1)*epsilon ;
        end
        %axe_freq1= -Fe/2 : Fe/length(Y): Fe/2-Fe/length(Y);
        chr = ['ksi = ',num2str(epsilon*(j-1))];
        plot(signal_gaus,'DisplayName',chr);
        grid on 
    end
     %legend show
     hold off
     
    % élimination de la porteuse 
    porteuse_correlee = cos(2*pi*(f0+correlation)*axe_temps1);
    
    signal_sans_porteuse_gaus = signal_recu_gaus.*porteuse_correlee;
    signal_sans_porteuse_balay = signal_recu_balay.*porteuse_correlee;
    signal_sans_porteuse_rep = signal_recu_rep.*porteuse_correlee;
    
    
    
    signal_sans_porteuse_gaus_non_correle = signal_recu_gaus_non_correle.*porteuse_non_correlee;
    signal_sans_porteuse_balay_non_correle = signal_recu_balay_non_correle.*porteuse_non_correlee;
    signal_sans_porteuse_rep_non_correle = signal_recu_rep_non_correle.*porteuse_non_correlee;
    
    % Représentation fréquencielle du signal transmis et du brouilleur
    %X =fftshift(fft(brouilleur_balayage)); 
    %X=X.*conj(X)/length(X); 
    %axe_freq= -Fe/2 : Fe/length(X): Fe/2-Fe/length(X);
    %Y =fftshift(fft(signal_transmis)); 
    %Y=Y.*conj(Y)/length(Y); 
    %axe_freq1= -Fe/2 : Fe/length(Y): Fe/2-Fe/length(Y);
    %figure
    %semilogy(axe_freq, X);
    %hold on
    %semilogy(axe_freq1, Y);
    %legend('bruit','signal_transmis')
    
    
    % filtre à la réception adapté (en racine de cosinus surréléve)
    for j = 1:length(filtre_emission)
        filtre_reception(j)= cosSurreleve(alpha,Tc-axe_temps(j),Tc);
    end
    
    % filtrer le signal recu 
    signal_filtre_adapte_gaus = conv(filtre_reception,signal_sans_porteuse_gaus);
    signal_filtre_adapte_balay = conv(filtre_reception,signal_sans_porteuse_balay);
    signal_filtre_adapte_rep = conv(filtre_reception,signal_sans_porteuse_rep);
    
    
    signal_filtre_adapte_gaus_non_correle = conv(filtre_reception,signal_sans_porteuse_gaus_non_correle);
    signal_filtre_adapte_balay_non_correle = conv(filtre_reception,signal_sans_porteuse_balay_non_correle);
    signal_filtre_adapte_rep_non_correle = conv(filtre_reception,signal_sans_porteuse_rep_non_correle);
    
    % Décoder le signal reçu filtré
    k=0 ;
    for i=1 : length(spreadData)
        k=k+1 ;
        m_gaus=real(signal_filtre_adapte_gaus(i*Tc+length(filtre_reception))) ;
        m_balay=signal_filtre_adapte_balay(i*Tc+length(filtre_reception)) ;
        m_rep=signal_filtre_adapte_rep(i*Tc+length(filtre_reception)) ;
        
        if m_gaus>= Seuil
            bit_estime_gaus(k)=1 ;
        else
            bit_estime_gaus(k)=-1 ;
        end 
        
        if m_balay>= Seuil
            bit_estime_balay(k)=1 ;
        else
            bit_estime_balay(k)=-1 ;
        end
        
        if m_rep>= Seuil
            bit_estime_rep(k)=1 ;
        else
            bit_estime_rep(k)=-1 ;
        end
        
        m_gaus_non_correle=signal_filtre_adapte_gaus_non_correle(i*Tc+length(filtre_reception)) ;
        m_balay_non_correle=signal_filtre_adapte_balay_non_correle(i*Tc+length(filtre_reception)) ;
        m_rep_non_correle=signal_filtre_adapte_rep_non_correle(i*Tc+length(filtre_reception)) ;
        
        if m_gaus_non_correle>= Seuil
            bit_estime_gaus_non_correle(k)=1 ;
        else
            bit_estime_gaus_non_correle(k)=-1 ;
        end 
        
        if m_balay_non_correle>= Seuil
            bit_estime_balay_non_correle(k)=1 ;
        else
            bit_estime_balay_non_correle(k)=-1 ;
        end
        
        if m_rep_non_correle>= Seuil
            bit_estime_rep_non_correle(k)=1 ;
        else
            bit_estime_rep_non_correle(k)=-1 ;
        end
    end
    

    % Inverser l'étalement de spectre 
    Data_estime_gaus = InverseEtalementDeSpectre(bit_estime_gaus,spreadSequence,facteur_etalement,length(dataEnc));
    Data_estime_balay = InverseEtalementDeSpectre(bit_estime_balay,spreadSequence,facteur_etalement,length(dataEnc));
    Data_estime_rep = InverseEtalementDeSpectre(bit_estime_rep,spreadSequence,facteur_etalement,length(dataEnc));
    
    Data_estime_gaus_non_correle = InverseEtalementDeSpectre(bit_estime_gaus_non_correle,spreadSequence,facteur_etalement,length(dataEnc));
    Data_estime_balay_non_correle = InverseEtalementDeSpectre(bit_estime_balay_non_correle,spreadSequence,facteur_etalement,length(dataEnc));
    Data_estime_rep_non_correle = InverseEtalementDeSpectre(bit_estime_rep_non_correle,spreadSequence,facteur_etalement,length(dataEnc));
    
    
    % Décodage canal 
    tbdepth = 34;
    
    Data_coded_gaus = vitdec(Data_estime_gaus,trellis,tbdepth,'trunc','hard');
    Ne_gaus=sum(abs(data-Data_coded_gaus));
    TEB_gaus(N0_bruit+6) = Ne_gaus/(Nb+Nb_p);
    
    Data_coded_balay = vitdec(Data_estime_balay,trellis,tbdepth,'trunc','hard');
    Ne_balay=sum(abs(data-Data_coded_balay));
    TEB_balay(N0_bruit+6) = Ne_balay/(Nb+Nb_p);
    
    Data_coded_rep = vitdec(Data_estime_rep,trellis,tbdepth,'trunc','hard');
    Ne_rep=sum(abs(data-Data_coded_rep));
    TEB_rep(N0_bruit+6) = Ne_rep/(Nb+Nb_p);
    
    Data_coded_gaus_non_correle = vitdec(Data_estime_gaus_non_correle,trellis,tbdepth,'trunc','hard');
    Ne_gaus_non_correle=sum(abs(data-Data_coded_gaus_non_correle));
    TEB_gaus_non_correle(N0_bruit+6) = Ne_gaus_non_correle/(Nb+Nb_p);
    
    Data_coded_balay_non_correle = vitdec(Data_estime_balay_non_correle,trellis,tbdepth,'trunc','hard');
    Ne_balay_non_correle=sum(abs(data-Data_coded_balay_non_correle));
    TEB_balay_non_correle(N0_bruit+6) = Ne_balay_non_correle/(Nb+Nb_p);
    
    Data_coded_rep_non_correle = vitdec(Data_estime_rep_non_correle,trellis,tbdepth,'trunc','hard');
    Ne_rep_non_correle=sum(abs(data-Data_coded_rep_non_correle));
    TEB_rep_non_correle(N0_bruit+6) = Ne_rep_non_correle/(Nb+Nb_p);
end

%figure 
%semilogy([-5:25],TEB_gaus,'-o k')
%hold on
%semilogy([-5:25],TEB_gaus_non_correle,'-o r')
%hold on
%semilogy([-5:25],TEB_balay,'-* k')
%hold on
%semilogy([-5:25],TEB_balay_non_correle,'-* r')
%hold on
%semilogy([-5:25],TEB_rep,'-< k')
%hold on
%semilogy([-5:25],TEB_rep_non_correle,'-< r') 
%xlabel('Pb (dB)')
%ylabel('TEB')
%legend('gaussien avec Doppler','gaussien sans Doppler','peigne de raies avec Doppler','peigne de raies sans Doppler','répéteur avec Doppler','répéteur sans Doppler')
%legend('gaussien','peigne de raies','répéteur')
%grid on
