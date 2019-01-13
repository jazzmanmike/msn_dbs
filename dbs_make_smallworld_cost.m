function dbs_make_smallworld_cost( CIJ )
%DBS_MAKE_SMALLWORLD_COST small world depending on cost
%
%   dbs_make_smallworld_cost(CIJ);
%
%   Inputs:     CIJ,            weighted adjacency matrix (diagonal = zero)
%
%   Outputs:    a series of figures plotting smallworldness with cost
%
% Michael Hart, University of Cambridge, July 2018
%% Initialise

nNodes = size(CIJ, 1);
density = density_und(CIJ);
cost_range = 0.01:0.01:density;

%% Generate cost data

Humphries = zeros(length(cost_range), 1); 
Latora = zeros(length(cost_range), 1);
Telesford = zeros(length(cost_range), 1);
counter = 1;
for iCost = cost_range
    CIJ_iCost = threshold_proportional(CIJ, iCost);
    [Humphries(counter, 1), Latora(counter, 1), Telesford(counter, 1)] = dbs_make_SmallWorlds(CIJ_iCost);
    counter = counter + 1;
end

%% Plot data
figure1 = figure('Name','Measure cost plots');
hold on

subplot1 = subplot(3,1,1,'Parent',figure1);
hold(subplot1,'on');
plot(cost_range, Humphries,'Parent',subplot1);
xlim([0 density]);
title({'Humphries'});

subplot2 = subplot(3,1,2,'Parent',figure1);
hold(subplot2,'on');
plot(cost_range, Latora,'Parent',subplot2);
xlim([0 density]);
title({'Latora'});

subplot3 = subplot(3,1,3,'Parent',figure1);
hold(subplot3,'on');
plot(cost_range, Telesford,'Parent',subplot3);
xlim([0 density]);
title({'Telesford'});

end
