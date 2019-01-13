function [ CostMeasures ] = dbs_network_cost( CIJ, maxCost )
%DBS_NETWORK_COST Cost function plots using a MST approach
%   Measures calculated in steps of 100 edges
%
%   CostMeasures = dbs_network_cost(CIJ);
%
%   Inputs:     CIJ,            weighted connectivity matrix
%               maxCost,        maximumCost of runs (optional)
%
%   Outputs:    CostMeasures,    metrics that vary with cost
%
% Michael Hart, University of Cambridge, August 2017

%% Define & initialise

if nargin <2
    maxCost = 0.3;
end

% Check matrix
nNodes = size(CIJ,1);           %Number of node
CIJ(CIJ<0) = 0;                 %Take out negative correlations
CIJ(1:nNodes+1:end) = 0;        %Set diagonal to 0
CIJ = max(CIJ, CIJ');           %trick to make symmetric

% Outputs
%nodalMeasures = zeros();
%globalMeasures = zeros();
%hubScore = zeros();

%% Create MST matrix
% Uses brainGL (which works with sparse matrices)

MST = kruskal_mst(sparse((1-CIJ))); %120 edges (67 * 67) ?why sqrt(2*(1-CIJ) ?trying to make mean0/std1?

%Store Initial MST in the adjacency matrix A that defines the network
A = full(MST);                  %Make unsparse
[iCol, jRow] = find(MST);       %Matrix co-ordinates of MST
for iMST = 1:length(jRow)
    A(iCol(iMST), jRow(iMST)) = CIJ(iCol(iMST), jRow(iMST));  
    A(jRow(iMST), iCol(iMST)) = CIJ(jRow(iMST), iCol(iMST));  
end

%% Order weights in matrix (for adding later)
%Order C according to decreasing wieghts in the correlation matrix
Co = triu(CIJ, 1);                                      %Only upper triangle
ind = find(Co);                                         %To get rid of 1s in matrix
Clist = Co(ind);                                        %Values of matrix   
[~, I] = sort(Clist,'descend');                         %Finds indices of values in order
[row, col] = ind2sub([nNodes,nNodes],ind(I));           %Puts indices into row/column co-ordinates   

%% Run loop to grow edges
% Records measures after every 100 edges
% Output will be measures at certain cost

%Initially, with just the MST: set counters and calculate cost and all measures
t = 1;          %edge index
eNum = nnz(A);  %number of edges
g = 1;          %counter
cost = zeros(); %cost

%Now add edges in correct order until all possible edges exist
while (eNum < (maxCost*nNodes*(nNodes-1)/2)) %maximum cost
    % if edge wasn't initially included in MST
    if A(row(t),col(t)) == 0
        %add edge weights according to co-ordinates
        A(row(t),col(t)) = Co(row(t),col(t)); 
        A(col(t),row(t)) = Co(row(t),col(t)); 
        eNum = eNum + 1; %increment edge counter
        if mod(eNum, 10) == 0 %after every 100 edges
            %Increment counter
            g = g + 1; %for plotting
            %calculate cost
            cost(g,1) = 2 * eNum / (nNodes * (nNodes - 1)); %calculates cost of measures
            %Call function that calculates all measures
            [nodalMeasures(:,:,g), globalMeasures(g,:), hubScore(:,g)] = dbs_quick_measures(A);
        end
    end
    t = t+1; %add next edge
end

%% Parse outputs

CostMeasures.global = globalMeasures;
CostMeasures.nodes = nodalMeasures;
CostMeasures.hubs = hubScore;
CostMeasures.cost = cost;

dbs_cost_plots %plot the cost - note called within a function cf. variable scope

end

