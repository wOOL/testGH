% Ant Colony algorithm for TSP
% Yongxin Yang <iamwool@gmail.com>
% Last edited: 16/11/2011
% I never hurt an ant when I was coding this :-)

function [resultRoute,resultDistance] = ...
    TSP_AC(D, numIterations, numAnts, pheromoneFactor,...
           heuristicFactor, evaporationCoef, pheromoneBoost)

% Explanation of parameters
% -------------------------------------------------------------------------
% D                  -    Distance matrix, i.e. tspdataN
% numIterations      -    number of iteration times
% numAnts            -    number of ants
% pheromoneFactor    -    the weight of pheromone factor
% heuristicFactor    -    the weight of heuristic factor
% evaporationCoef    -    pheromone evaporation coefficient 
% pheromoneBoost     -    pheromone boost coefficient
% -------------------------------------------------------------------------

% ---------------------- Step 1: Initialize variables --------------------- 

numCities = size(D,1);
% Get the number of cities
D(D==0) = 10e-10;
% Let the zero-value-element be an extremely small value such that 
% 1 divided by 0 cases will be avoided in following procedure
Alpha = 1 ./ D;         
% Alpha is heuristic factor, in TSP case, it is reciprocal for 
% the distance between two cities
pheromoneMatrix = ones(numCities,numCities);     
% pheromoneMatrix is the matrix recording pheromone left by ants
routeMatrix = zeros(numAnts,numCities);  
% routeMatrix is the matrix recording cities visited by each ant, it's also
% called tabu list which means cities on it can not be visited again
bestRoute = zeros(numIterations,numCities);       
% The best route of each iteration
shortestDistance = inf.*ones(numIterations,1);   
% The shortest distance corresponding to best route

% ---------------------------- End of Step 1 ------------------------------ 

for Iteration = 1:numIterations
    
    % --------- Step 2: Generate a random start city for each ant --------- 
    % --------- Number of ants isn't larger than number of cities ---------
    
    randStartPoint = randperm(numCities);
    routeMatrix(:,1) = randStartPoint(1:numAnts);
    
    % ------------------------ End of Step 2 ------------------------------ 

    % ----------- Step 3: Each ant chooses its next city based   ---------- 
    % ----------- on probability until its travel is finished    ----------
    
    for j = 2:numCities     
    % Start from the second city
        for i = 1:numAnts    
            visitedCities = routeMatrix(i,1:(j-1)); 
            % Record the cities which have been visited
            J = zeros(1,(numCities-j+1));       
            % Initialize a vector to store the cities which will be visited
            P = J;                      
            % Initialize a vector to store probability of corresponding cities
            Jc = 1;
            % Counter - number of visited cities
            for k = 1:numCities
                if isempty(find(visitedCities==k, 1))
                    J(Jc) = k;
                    Jc = Jc + 1;
                end
            end
            % If next city is not visited, add it into list,
            % then counter + 1
            
            for k = 1:length(J)
                P(k) = (pheromoneMatrix(visitedCities(end),J(k))^pheromoneFactor)...
                      *(Alpha(visitedCities(end),J(k))^heuristicFactor);
            end
            % Calculate the probability distribution P using pheromone
            % factor and heuristic factor
            
            P = P / sum(P);
            % Normalization
            
            Pcum=cumsum(P);
            Selected = find(Pcum>=rand);
            toVisit = J(Selected(1));
            routeMatrix(i,j) = toVisit;
            % Select the route using roulette wheel selection
        end
    end
    
    if Iteration >= 2
        routeMatrix(1,:) = bestRoute(Iteration-1,:);
    end
    
    % ------------------------ End of Step 3 ------------------------------

    % -------- Step 4: Calculate the best route for this iteration --------
    
    L = zeros(numAnts,1);
    % Initialize a vector to store the distance of route chosen by each ant
    for i = 1:numAnts
        R = routeMatrix(i,:);
            for j = 1:(numCities-1)
                L(i) = L(i) + D(R(j),R(j+1));
                % Every time I wrote the code like this, I wanna ask
                % Why? there is no 'i += i_new' stuff in Matlab. FML
            end
        L(i) = L(i)+D(R(1),R(numCities));
        % Add the distance from terminal point to starting point
    end
    shortestDistance(Iteration) = min(L);
    % Shortest Distance
    pos = find(L==shortestDistance(Iteration));
    bestRoute(Iteration,:) = routeMatrix(pos(1),:);
    % Best Route
    
    % ------------------------ End of Step 4 ------------------------------
    
    % --------- Step 5: Update pheromone using pheromone boost    ---------
    % --------- coefficient and pheromone evaporation coefficient ---------
    
    Delta = zeros(numCities,numCities);
    % Initialize a matrix to store the changes of pheromone matrix
    for i=1:numAnts
        for j=1:(numCities-1)
            Delta(routeMatrix(i,j),routeMatrix(i,j+1)) = ...
            Delta(routeMatrix(i,j),routeMatrix(i,j+1)) + pheromoneBoost/L(i);          
            % The addition of pheromone between city i and city j
            % caught by one ant.
        end
        Delta(routeMatrix(i,numCities),routeMatrix(i,1)) = ...
        Delta(routeMatrix(i,numCities),routeMatrix(i,1)) + pheromoneBoost/L(i);
        % The addition of pheromone caught by all ants
    end
    pheromoneMatrix = (1-evaporationCoef).*pheromoneMatrix+Delta; 
    % Some of pheromone will evaporate
    
    % ------------------------ End of Step 5 ------------------------------
    
    routeMatrix = zeros(numAnts,numCities);
    % Clear the route matrix preparing for next iteration
end

% Output
pos = find(shortestDistance==min(shortestDistance));
resultRoute = bestRoute(pos(1),:);
resultDistance = shortestDistance(pos(1)); 