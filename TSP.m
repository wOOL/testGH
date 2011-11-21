D = load('tspadata2.txt');
numAnts = 48;
numIterations = 200;
pheromoneFactor = 1;
heuristicFactor = 5;
evaporationCoef = 0.1;
pheromoneBoost = 100;

dis = zeros(50,1);
tic;
for i=1:50
    [Route,Distance] = TSP_AC(D,numIterations,numAnts,pheromoneFactor,heuristicFactor,evaporationCoef,pheromoneBoost);
    dis(i) = Distance;
end;
toc;
hist(dis);