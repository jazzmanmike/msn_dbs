function [ RDN, DGN, CMP ] = dbs_complexity( network )
%DBS_COMPLEXITY    computes redundancy & degeneracy
%
%   [RDN, DGN, CMP] = dbs_complexity(network); 
%
%   Redundancy: a parcel with similar connections in-situ 
%   i.e. its function is redundant compared to the original parcel
%   
%   Degeneracy: a parcel with similar connections to another after removal
%   i.e. it reconstitutes function
%
%   Input:      network,    connectivity matrix
%
%   Output:     RDN,        redundancy
%               DGN,        degeneracy
%               CMP,        complexity
%
%   Michael Hart, University of Cambridge, July 2018

%% Baseline definitions

A = network; %original connectivity matrix re-defined
N = size(network, 2); %number of parcels
RDN = zeros(size(network)); % initialise output - matrix 
DGN = zeros(size(network)); % initialise output 

%% Redundancy

for i = 1:N; % for each parcel in turn
    grotI = A(i,:); % parcel of interest
    for j = 1:N; % work through all other parcels in turn
        grotJ = A(j,:); % specific parcel for comparison
        RDN(i,j) = dbs_mutual_information(grotI, grotJ); %gives a scalar value
    end
end

%% Degeneracy

for i = 1:N; % for each parcel in turn
    grotI = A(i,:); % parcel of interest
    grotA = A; grotA(i,:) = 0; grotA(:,i) = 0; % new *lesioned* network
    for j = 1:N; % work through all other parcels in turn
        grotJ = grotA(j,:); % specific parcel for comparison
        DGN(i,j) = dbs_mutual_information(grotI, grotJ); %gives a scalar value
    end
end

%% Complexity

CMP = RDN - DGN;

end

