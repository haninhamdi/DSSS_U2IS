%%%%%%%%%%%%%%%%%%%% PRE_U2IS : Quatrième semaine %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Objectif : Etalement de spectre + brouillage %%%%%%%%%%%%%%%%
%*************************************************************************%
% NbPk : Le nombre des packets
% Nb : nombre des symboles (=bits dans le cas BPSK) par packet 
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
NbPk = 2500 ; 
Nb = 1000 ; 
facteur_etalement = 16 ; 
Ts = 32 ;
fs = 1/Ts ;
Tc = Ts/facteur_etalement;
fc = 1/Tc ;
f0 = 0.4 ;
T0 = 1/f0 ;
alpha = 0.22 ;
Fe = 50*(10^6) ;
Seuil=0 ;
for N0_bruit = 0:25
    
    ErTotale_gaus = 0 ; 
    ErTotale_balay = 0 ; 
    ErTotale_rep = 0 ; 
    for pk =1:NbPk
        % Générer aléatroirement les bits 
    data = randi([0,1], Nb, 1);
    
    % Codage canal 
    dataEnc = convenc(data,trellis);
    % Modulation BPSK de dataEnc : 1 et -1 
    for j = 1 : length(dataEnc)
        if dataEnc(j) == 0
            dataEnc(j) = -1;
        end
    end
    
    % Etalement de spectre 
    [spreadData,spreadSequence] = EtalementDeSpectre(facteur_etalement,dataEnc,length(dataEnc));
  
    % conversion en un signal physique
    N = length(spreadData);
    peigne = zeros(1,N*Tc) ;
    for i= 0: N-1
        peigne(1+i*Tc) = spreadData(1+i);
    end

    % Construction du filtre d'émission en racine de cosinus suréléve 
    axe_temps = -10*Tc:10*Tc;
    for i=1:length(axe_temps)
        filtre_emission(i) = cosSurreleve(alpha,axe_temps(i),Tc);
    end

    % effectuer la convolution entre le signal (peigne) et la fonction du
    % filtre
    signal_filtre = conv(filtre_emission,peigne);

    % effectuer le produit entre le signal réel et la porteuse 
    N2 = length(signal_filtre);
    axe_temps1 = -N2*Tc/2:Tc:(N2-1)*Tc/2;
    porteuse = cos(2*pi*f0*axe_temps1);
    signal_transmis = signal_filtre.*porteuse;

    % ajout d'un BBAG  
    RSB=10;
    N0 = 10^(-RSB/10);
    sigma_bruit = sqrt(N0);
    bruit = sigma_bruit * randn(1,length(signal_transmis));
    signal_canal = signal_transmis + bruit;

    % Bloc pour le brouillage 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Pb = 10^(N0_bruit/10);
    sigma_brouilleur = sqrt(Pb);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Brouillage gaussien
    brouilleur_gaussien = 1.2*sigma_bruit*randn(1,length(signal_canal));
    [b,a] = butter(8,fc) ; 
    brouilleur_gaussien = sigma_brouilleur*filter(b,a,brouilleur_gaussien).*porteuse; 
    signal_recu_gaus = signal_canal + brouilleur_gaussien;
    
     % Brouillage par balayage 
    %delta_f = fc/20 ;
    %f = -fc+4*delta_f:delta_f:fc-4*delta_f;
    %somme = 0;
    %for i = 1:length(axe_temps1)
    %    somme = 0;
    %    for j=1:length(f)
    %        somme = somme + cos(2*pi*f(j)*axe_temps1(i));
    %    end
    %    brouilleur_balayage(i)= sqrt(Pb/length(f))*somme;
    %end
    brouilleur_balayage = 1.7*sigma_bruit*randn(1,length(signal_canal));
    [b,a] = butter(8,fc) ; 
    brouilleur_balayage = sigma_brouilleur*filter(b,a,brouilleur_balayage).*porteuse; 
    signal_recu_balay = brouilleur_balayage + signal_canal;
    
    % Brouilleur répéteur
    % générer un retard Tc< tau < Ts aléatoire 
     %tau = randi([Tc,Ts]);
     %brouilleur_repeteur = zeros(1,length(signal_canal));
     %for j=tau:length(signal_canal)-1
     %    brouilleur_repeteur(j) = sigma_brouilleur*signal_canal(j-tau+1);
     %end
     brouilleur_repeteur = 2.3*sigma_bruit*randn(1,length(signal_canal));
     [b,a] = butter(8,fc) ; 
     brouilleur_repeteur = sigma_brouilleur*filter(b,a,brouilleur_repeteur).*porteuse; 
     signal_recu_rep = brouilleur_repeteur + signal_canal;
     
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % élimination de la porteuse 
    signal_sans_porteuse_gaus = signal_recu_gaus.*porteuse;
    signal_sans_porteuse_balay = signal_recu_balay.*porteuse;
    signal_sans_porteuse_rep = signal_recu_rep.*porteuse;
    
    % Représentation fréquencielle du signal transmis et du brouilleur
    %X =fftshift(fft(brouilleur_balayage)); 
    %X=X.*conj(X)/length(X); 
    %axe_freq= -Fe/2 : Fe/length(X): Fe/2-Fe/length(X);
    %Y =fftshift(fft(signal_filtre)); 
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
    
    % Décoder le signal reçu filtré
    k=0 ;
    for i=1 : length(spreadData)
        k=k+1 ;
        m_gaus= signal_filtre_adapte_gaus(i*Tc+length(filtre_reception)) ;
        m_balay= signal_filtre_adapte_balay(i*Tc+length(filtre_reception)) ;
        m_rep= signal_filtre_adapte_rep(i*Tc+length(filtre_reception)) ;
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
    end 

    % Inverser l'étalement de spectre 
    Data_estime_gaus = InverseEtalementDeSpectre(bit_estime_gaus,spreadSequence,facteur_etalement,length(dataEnc));
    Data_estime_balay = InverseEtalementDeSpectre(bit_estime_balay,spreadSequence,facteur_etalement,length(dataEnc));
    Data_estime_rep = InverseEtalementDeSpectre(bit_estime_rep,spreadSequence,facteur_etalement,length(dataEnc));
    
    
    % Décodage canal 
    tbdepth = 34;
        Data_coded_gaus = vitdec(Data_estime_gaus,trellis,tbdepth,'trunc','hard');
        erPkt_gaus=sum(abs(data-Data_coded_gaus));
        ErTotale_gaus = ErTotale_gaus + (erPkt_gaus > 0) ;
        
        Data_coded_balay = vitdec(Data_estime_balay,trellis,tbdepth,'trunc','hard');
        erPkt_balay=sum(abs(data-Data_coded_balay));
        ErTotale_balay = ErTotale_balay + (erPkt_balay > 0) ;
    
        Data_coded_rep = vitdec(Data_estime_rep,trellis,tbdepth,'trunc','hard');
        erPkt_rep=sum(abs(data-Data_coded_rep));
        ErTotale_rep = ErTotale_rep + (erPkt_rep > 0) ;
    end
    TEP_gaus(N0_bruit+10) = ErTotale_gaus/NbPk ;
    TEP_balay(N0_bruit+1) = ErTotale_balay/NbPk ;
    TEP_rep(N0_bruit+1) = ErTotale_rep/NbPk ;
end

figure 
semilogy([0:25],TEP_gaus)
hold on
semilogy([0:25],TEP_balay)
hold on
semilogy([0:25],TEP_rep)
%title('Le taux d erreur par paquet en fonction de Pb (Pu = 1) (RSB=10)') 
xlabel('Pb (dB)')
ylabel('TEP')
legend('gaussien','peigne de raies','répéteur')
grid on
