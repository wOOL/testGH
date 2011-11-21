D = load('tspadata2.txt');
numAnts = 48;
numIterations = 200;
pheromoneFactor = 1;
heuristicFactor = 5;
evaporationCoef = 0.1;
pheromoneBoost = 100;

tic;
[Route,Distance] = TSP_AC(D,numIterations,numAnts,pheromoneFactor,heuristicFactor,evaporationCoef,pheromoneBoost);
toc;