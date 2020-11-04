% Cettte fonction effectue l'étalement de spectre
function [spreadData,PN_sequence] = EtalementDeSpectre(facteur,data,Nb)
PN_sequence=randsrc(facteur,1,[0:1]); 
for j=1:facteur
    if PN_sequence(j)==0
        PN_sequence(j)=-1;
    end
end
j=1;
for i = 1: Nb
    for k = j:j+facteur
        spreadData(k) = data(i);
    end
    spreadData(j:j+facteur-1) = spreadData(j:(j+facteur-1))'.*PN_sequence(1:facteur);
    j = facteur*i+1;    
end
end