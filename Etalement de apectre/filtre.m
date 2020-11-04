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
Nb = 100 ;
facteur_etalement = 8 ;
Ts = 128 ;
fs = 1/Ts ;
Tc = Ts/facteur_etalement;
fc = 1/Tc;
f0 = 0.4 ;
T0 = 1/f0 ;
alpha = 0.22 ;
Fe = 50 ;
Seuil=0 ;
for N0_bruit = 0:0
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
    Nt = length(dataEnc);
    peignet = zeros(1,Nt*Ts) ;
    for i= 0: Nt-1
        peignet(1+i*Ts) = dataEnc(1+i);
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
    axe_temps1 = -N2/2:(N2-1)/2;
    porteuse = cos(2*pi*f0*axe_temps1);
    signal_transmis = signal_filtre;
   
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
    brouilleur_gaussien = 1.1*sigma_bruit*randn(1,length(signal_canal));
    [b,a] = butter(12,fc,'low') ; 
    brouilleur_gaussien = sigma_brouilleur*filter(b,a,brouilleur_gaussien); 
    signal_recu_gaus = signal_canal + brouilleur_gaussien;
    
    % Brouillage par balayage 
    
     delta_f = fc/40 ;
     f = -fc+12*delta_f:delta_f:fc-12*delta_f ; 
     for i = 1:length(axe_temps1)
        somme = 0 ;
        for j=1:length(f)
            somme = somme + cos(2*pi*f(j)*axe_temps1(i));
        end
        brouilleur_balayage(i)= somme;
     end
    brouilleur_balayage =  brouilleur_balayage;
    
    %brouilleur_balayage = 1.7*sigma_bruit*randn(1,length(signal_canal));
    %[b,a] = butter(8,fc) ; 
    %brouilleur_balayage = sigma_brouilleur*filter(b,a,brouilleur_balayage).*porteuse; 
    brouilleur_balayage= sqrt(Pb/length(f))*brouilleur_balayage ;
    %brouilleur_balayage= brouilleur_balayage.*porteuse ;
    
    signal_recu_balay = brouilleur_balayage + signal_canal;
    
    % Brouilleur répéteur
    % générer un retard Tc< tau < Ts aléatoire 
    
     %tau = randi([Tc,Ts]);
     %brouilleur_repeteur = zeros(1,length(signal_canal));
     %for j = tau:length(signal_canal)-1
     %    brouilleur_repeteur(j) = signal_canal(j-tau+1);
     %end
     %brouilleur_repeteur = sigma_brouilleur*brouilleur_repeteur;
     brouilleur_repeteur = 2.3*sigma_bruit*randn(1,length(signal_canal));
     %[b,a] = butter(8,fc) ; 
     brouilleur_repeteur = sigma_brouilleur*filter(b,a,brouilleur_repeteur).*porteuse; 
     signal_recu_rep = brouilleur_repeteur+ signal_canal;
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    % Représentation fréquencielle du signal transmis et du brouilleur
    % signal en bande de base avant étalement 

    
    figure
    plot(axe_temps,filtre_emission);
    xlabel('Temps (\mus)')
    ylabel('Amplitude')
    grid on 
    
    
    
end
