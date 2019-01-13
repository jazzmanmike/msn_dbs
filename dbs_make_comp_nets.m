function [ graphsArray, graphsCode ] = dbs_make_comp_nets( CIJ )
%DBS_MAKE_COMP_NETS Makes a selection of graphs for comparisons
%
%   Binary only
%   
%   [graphsArray, graphsCode] = dbs_make_comp_nets(CIJ, threshold);
%
%   Makes a selection of comparison graphs for a set matrix CIJ
%
%   Inputs:     CIJ,            weighted matrix
%
%   Outputs:    graphsArray,    array of 8 comparison graphs
%               graphsCode,     string of comparison graphs
%
%   Version 2:  thresholding / cost removed
%
% Michael Hart, University of Cambridge, May 2018

%% Basic checks on CIJ

if max(rem(CIJ,1)) == 0
    disp('input is binary')
else
    disp('input is weighted - binarising network')
    CIJ = double(CIJ>0);
end

%% Threshold CIJ and initialise parameters

nNodes = size(CIJ,1); %number of nodes
nEdges = nnz(CIJ); %number of edges
degree_dist = sum(CIJ); %degree distribution
maxDegree = max(sum(CIJ)); %maximum degree
degreeRange = 1:maxDegree; %range of degrees, ignoring isolated nodes

graphsArray = zeros(nNodes,nNodes,8); % define array (116,116,8);
graphsCode = {'empirical'; 'lattice'; 'smallworld'; 'random'; ...
    'modular'; 'exponential'; 'scalefree'; ...
    'preferential_attachment'};

%% 2. Generate test graphs, including CIJ as #1

%% i. Empirical 

graphsArray(:,:,1) = CIJ;

%% ii. Lattice

graphsArray(:,:,2) = makelatticeCIJ(nNodes, nEdges);

%% iii. Small world

SW = smallw(nNodes,8,10); %re-wiring empirically adjusted: contest toolbox
graphsArray(:,:,3) = full(SW); %NB: edges will not match

%% iv. Random weights & topology

randomNetwork = zeros(nNodes,nNodes,10); %100 iterations
for j = 1:10
    randomNetwork(:,:,j) = randmio_und(CIJ, 100);
end
graphsArray(:,:,4) = mean(randomNetwork, 3);

%% v. Hierarchical modular

%Make network
HM_network = makefractalCIJ(7, 3, 4); %heuristical - give 2862 edges but 128 nodes

nRemovals = length(HM_network) - nNodes;

select_removals = randperm(size(HM_network,1), nRemovals);

HM_network(:,select_removals) = []; HM_network(select_removals,:) = [];

%Initialise parameters
%HM_nodes = 0;
%tmp = 2;
%while HM_nodes < nNodes
%    HM_nodes = 2^tmp;
%    tmp = tmp+1;
%end
%mx_lvl = tmp-1; %number of nodes in original
%modules = 6; %corresponds to roughly 12 modules at nodes=1024
%now threshold number of nodes
%select_nodes = randperm(size(HM_network,1),nNodes); %randomly N nodes of HM_node
%HM_network = HM_network(1:4:end, 1:4:end);

%now threshold number of edges
%if nnz(HM_network)>nEdges %therefore HM_network must be big to start
%    ratio = nnz(HM_network) / nEdges;
%    HM_network = threshold_proportional(HM_network, ratio);
%end

graphsArray(:,:,5) = HM_network;

%% vi. Exponential

mu = 1/(sum(degree_dist)/nNodes);
probx = exp(-mu*(1:maxDegree)); %generates probability of node degrees

expDegrees = randsample(degreeRange,nNodes,true,probx); %2161 edges

[exponentialNetwork,b] = makerandCIJdegreesfixed(expDegrees, expDegrees);
if b==1
    disp 'exponential network made successfully'
end

%now check edges are the same
%if nnz(exponentialNetwork)>nEdges %usually larger than CIJ
%    ratio = nnz(exponentialNetwork) / nEdges;
%    exponentialNetwork = threshold_proportional(exponentialNetwork, ratio);
%end

graphsArray(:,:,6) = exponentialNetwork;

%% Amendments

%use histogram to get number of bins and number of counts in each bin
%divide number of counts in each bin by total counts to get PDF summing to 1

%next, create ft of function (SF, exp, exp-trunc) with fttype(degree range, degree probabilities, function)
%get specific values with fit (also gives new function f)
%use new function f with degrees f(degree range)

%now, the probability of each degree with this function is given
%% vii. Scale-free

nb = nnz(degree_dist); 

gamma = 1 + nb/(sum(log(degree_dist(degree_dist>0)))); %not needed
probx = degreeRange.^(-gamma); %exponential probability

%s = zeros(nNodes);
%stop = 0;
%while(stop==0)
     s = randsample(degreeRange,nNodes,true,probx);
%    if sum(s)==nEdges 
%        stop=1;
%    end
%end

[sfNetwork,b] = makerandCIJdegreesfixed(s,s);
if b==1
    disp 'scale-free network made successfully'
end

graphsArray(:,:,7) = sfNetwork; %NB: edges will not match

%% viii. Exponentially truncated
%tmp = threshold_proportional(CIJ,0.1);
%degree_dist = sum(tmp);
%x = degree_dist;
%x = x(x>0); %only keep if not 0
%n = length(x); 
%n = degreeRange;
%p = 1;

%fn = -(-n * p * log(sum(x)/(n * p)) - n * log(gamma(p)) + (p - 1) * sum(log(x)) - n * p);
   
%out = nlm(fn, p = 1, hessian = TRUE) %non-linear minimisation
%alpha <- out$estimate %0.6149
%beta <- sum(degree.dist)/(n.regions*alpha) %6.6899
%gamma.trace <- 1-pgamma((0:nmax), shape = d$alpha, scale = d$beta)

%need to make graphs out of given degree distributions - use

%% x. Preferential attachment

paNetwork = pref(nNodes); %contest toolbox
d=3; %baseline is degree 2

while nnz(paNetwork)<nEdges
    paNetwork = pref(nNodes,d); %will only be approximate
    d = d+1;
end

message = sprintf('mimimum degree in preferential attachment network is %d', d);
disp(message);

graphsArray(:,:,8) = paNetwork;

%% Time to plot

figure1 = figure('Name','Binary comparison matrices');

subplot1 = subplot(2,4,1,'Parent',figure1);
hold(subplot1,'on');
image(graphsArray(:,:,1),'Parent',subplot1,'CDataMapping', 'scaled');
title({'CIJ binary thresholded'});
xlim(subplot1,[0 nNodes]);
ylim(subplot1,[0 nNodes]);

subplot2 = subplot(2,4,2,'Parent',figure1);
hold(subplot2,'on');
image(graphsArray(:,:,2),'Parent',subplot2,'CDataMapping', 'scaled');
title({'lattice'});
xlim(subplot2,[0 nNodes]);
ylim(subplot2,[0 nNodes]);

subplot3 = subplot(2,4,3,'Parent',figure1);
hold(subplot3,'on');
image(graphsArray(:,:,3),'Parent',subplot3,'CDataMapping', 'scaled');
title({'small world'});
xlim(subplot3,[0 nNodes]);
ylim(subplot3,[0 nNodes]);

subplot4 = subplot(2,4,4,'Parent',figure1);
hold(subplot4,'on');
imagesc(graphsArray(:,:,4),'Parent',subplot4, 'CDataMapping', 'scaled');
title({'random'});
xlim(subplot4,[0 nNodes]);
ylim(subplot4,[0 nNodes]);

subplot5 = subplot(2,4,5,'Parent',figure1);
hold(subplot5,'on');
imagesc(graphsArray(:,:,5),'Parent',subplot5, 'CDataMapping', 'scaled');
title({'hierarchical modular'});
xlim(subplot5,[0 nNodes]);
ylim(subplot5,[0 nNodes]);

subplot6 = subplot(2,4,6,'Parent',figure1);
hold(subplot6,'on');
imagesc(graphsArray(:,:,6),'Parent',subplot6, 'CDataMapping', 'scaled');
title({'exponential'});
xlim(subplot6,[0 nNodes]);
ylim(subplot6,[0 nNodes]);

subplot7 = subplot(2,4,7,'Parent',figure1);
hold(subplot7,'on');
imagesc(graphsArray(:,:,7),'Parent',subplot7, 'CDataMapping', 'scaled');
title({'scalefree'});
xlim(subplot7,[0 nNodes]);
ylim(subplot7,[0 nNodes]);

subplot8 = subplot(2,4,8,'Parent',figure1);
hold(subplot8,'on');
imagesc(graphsArray(:,:,8),'Parent',subplot8, 'CDataMapping', 'scaled');
title({'preferential attachment'});
xlim(subplot8,[0 nNodes]);
ylim(subplot8,[0 nNodes]);
end

