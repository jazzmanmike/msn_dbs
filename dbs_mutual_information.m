function [ MI ] = dbs_mutual_information( grotI, grotJ )

%DBS_MUTUL_INFORMATION calculates the mutual information between two variables
%   Based on Simas and Salvador

n = size(grotI,1);

[Pxy,~]=cpsd(grotI, grotJ, [], [], n); 
[Pxx,~]=cpsd(grotI, grotI, [], [], n);
[Pyy,~]=cpsd(grotJ, grotJ, [], [], n);

Rxy=Pxy./(sqrt(Pxx.*Pyy)); % estimate from Salvador et al 2005

PCohxy=(abs(Rxy)).^2;  % Coherence [0,1]

%PCohxy(find(F<lmin | F>lmax))=[]; % frequency band filter

dxy=-(1/(2*pi))*mean(log(1-PCohxy)); %Mutual Information

MI=sqrt(1-exp(-2*dxy)); %Mutual Information normalized [0,1]

end

