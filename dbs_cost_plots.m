%costplots script

figure1 = figure('Name','Measure cost plots');

nPoints = length(CostMeasures.global); 
nNodes = length(CostMeasures.nodes);
x = ((1:nPoints) * 10) ./ ((nNodes*(nNodes-1))/2);

subplot1 = subplot(5,2,1,'Parent',figure1);
hold(subplot1,'on');
plot(x, squeeze(mean(nodalMeasures(:,1,:))),'Parent',subplot1);
xlim([0 maxCost]);
title({'Degree'});

subplot2 = subplot(5,2,2,'Parent',figure1);
hold(subplot2,'on');
plot(x, globalMeasures(:,1),'Parent',subplot2);
xlim([0 maxCost]);
title({'Strength'});

subplot3 = subplot(5,2,3,'Parent',figure1);
hold(subplot3,'on');
plot(x, globalMeasures(:,4),'Parent',subplot3);
xlim([0 maxCost]);
title({'Clustering'});

subplot4 = subplot(5,2,4,'Parent',figure1);
hold(subplot4,'on');
plot(x, squeeze(mean(nodalMeasures(:,4,:))),'Parent',subplot4);
xlim([0 maxCost]);
title({'Closeness'});

subplot5 = subplot(5,2,5,'Parent',figure1);
hold(subplot5,'on');
plot(x, globalMeasures(:,5),'Parent',subplot5);
xlim([0 maxCost]);
title({'Global efficiency'});

subplot6 = subplot(5,2,6,'Parent',figure1);
hold(subplot6,'on');
plot(x, globalMeasures(:,6),'Parent',subplot6);
xlim([0 maxCost]);
title({'Betweenness'});

subplot7 = subplot(5,2,7,'Parent',figure1);
hold(subplot7,'on');
plot(x, squeeze(mean(nodalMeasures(:,6,:))),'Parent',subplot7);
xlim([0 maxCost]);
title({'Z score'});

subplot8 = subplot(5,2,8,'Parent',figure1);
hold(subplot8,'on');
plot(x, squeeze(mean(nodalMeasures(:,7,:))),'Parent',subplot8);
xlim([0 maxCost]);
title({'Participation'});

subplot9 = subplot(5,2,9,'Parent',figure1);
hold(subplot9,'on');
plot(x, squeeze(mean(nodalMeasures(:,8,:))),'Parent',subplot9);
xlim([0 maxCost]);
title({'Eigenvector'});

subplot10 = subplot(5,2,10,'Parent',figure1);
hold(subplot10,'on');
plot(x, globalMeasures(:,7),'Parent',subplot10);
xlim([0 maxCost]);
title({'Semi-metricity'});
