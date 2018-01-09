 classdef GMOTSS < handle
    
    properties
        name_parameters = {'numIt','nPobl','k','pCross','pMut','seed','minSeg','maxSeg','polyDegree','sizeChromosome','iterClust'}
        dataFile
        data
        idealSegFile
        idealseg
        parameters
    end
    
    methods
        
        %% Constructor
        function obj = GMOTSS()
            obj.defaultParameters();
        end

        %% Default parameters of the class
        function obj = defaultParameters(obj)
            % Number of generations
            obj.parameters.numIt = 200;
            % Number of individuals of the population
            obj.parameters.nPobl = 80;
            % Number of clusters for k-means
            obj.parameters.k = 5;
            % Crossover probability
            obj.parameters.pCross = 0.8;
            % Mutation probability
            obj.parameters.pMut = 0.2;
            % Percentage to remove/add points
            obj.parameters.mutedPoints = 0.1;
            % Random number generation seed
            obj.parameters.seed = 1;
            % Minimum length of each segment
            obj.parameters.minSeg = 1;
            % Maximum length of each segment
            obj.parameters.maxSeg = 10;
            % Number of iterations of k-means algorithm
            obj.parameters.iterClust = 20;
            % Degree polynomial
            obj.parameters.polyDegree = 2;
            % Fitness for clustering 1-10
            obj.parameters.typeFitness = 1;
            % Fitness for error: 1-RMSE 2-RMSEp 3-MAXe
            obj.parameters.typeFitness2 = 1;
        end
        
        %% Parameters of the algorithm
        function [parameters_as_str] = getParameters(obj)
            parameters = obj.parameters;
            
            fields = fieldnames(parameters);
            parameters_as_str = '';
            
            for i = 1:numel(fields)
                parameters_as_str = [parameters_as_str sprintf('%s;%f\n', fields{i}, parameters.(fields{i}))];
            end
            
        end
        
        %% Main Evolutionary algorithm
        function [information,informationBestClustering,informationBestError] = runEvolutive(obj, serie)
            addpath kmeans/
            obj.data = serie;
            nOfData = length(serie);
            if strcmp(version('-release'),'2013a')
                s = RandStream('mt19937ar','Seed',obj.parameters.seed);
                RandStream.setGlobalStream(s);
            else
                s = RandStream.create('mt19937ar','seed',obj.parameters.seed);
                RandStream.setDefaultStream(s);
            end
            obj.parameters.sizeChromosome = nOfData;
            initialPopulation = obj.initialisePopulation();
            initialFitness = zeros(1,obj.parameters.nPobl)*NaN;
            [initialFitness1, initialFitness2] = obj.evaluate(initialPopulation,initialFitness);
            [f,currentPopulation,currentFitness1,currentFitness2,crowding_distances] = obj.non_domination_sort(initialPopulation, initialFitness1, initialFitness2);
                        
            for i=1:obj.parameters.numIt,
                fprintf('It: %d\n',i);
                % To use 'Tournament Selection'
                % [parentPopulation,parentFitness1,parentFitness2] = GMOTSS.binary_tournament(currentPopulation,currentFitness1,currentFitness2);
                parentPopulation = currentPopulation;
                parentFitness1 = currentFitness1;
                parentFitness2 = currentFitness2;
                %'Crossover'
                [newPopulation, changedFitness1] = obj.crossover(parentPopulation,parentFitness1);
                %'Mutation'
                [newPopulation, changedFitness1] = obj.mutation(newPopulation,changedFitness1);

                %'Evaluation'
                [newFitness1, newFitness2] = obj.evaluate(newPopulation,changedFitness1);
                [f,resultantPopulation,resultantFitness1,resultantFitness2,crowding_distances] = obj.non_domination_sort([parentPopulation; newPopulation], [parentFitness1 newFitness1], [parentFitness2 newFitness2]);
                
                 %'Selection'
                [currentPopulation, currentFitness1, currentFitness2] = obj.selection(resultantPopulation,resultantFitness1,resultantFitness2);
                
                % To see 'Repeated individuals' in each generation
                % [r, c] = size(unique(currentPopulation, 'rows'));
                % repetidos = obj.parameters.nPobl - r;
                % fprintf('Repeated Individuals: %d ----> Percentage: %f\n', repetidos, repetidos/obj.parameters.nPobl);
            end
            
            % Best global solution: Fitness closer to 1,1
            minimo = min(currentFitness1(1,:));
            maximo = max(currentFitness1(1,:));
            normFit1(1,:) = (currentFitness1(1,:)-minimo) / (maximo -minimo);
            minimo = min(currentFitness2(1,:));
            maximo = max(currentFitness2(1,:));
            normFit2(1,:) = (currentFitness2(1,:)-minimo) / (maximo -minimo);
            normFit1(1,:)= 1 - normFit1(1,:);
            normFit1(1,:)= normFit1(1,:).*normFit1(1,:);
            normFit2(1,:)= 1 - normFit2(1,:);
            normFit2(1,:)= normFit2(1,:).*normFit2(1,:);
            distances_to_ideal(1,:) = sqrt(normFit1(1,:)+normFit2(1,:));
            [mini fbestidx] = min(distances_to_ideal(1,:));
            fprintf('*******BEST SEGMENTATION IN GLOBAL TERMS*******\n');
            fprintf('Indx: %d\n', fbestidx);
            fprintf('Fitness Clustering: %.15f\n', currentFitness1(1,fbestidx));
            fprintf('Fitness Error: %.15f\n', currentFitness2(1,fbestidx));
            % Information for the reporter
            information.fbestClustering = currentFitness1(1,fbestidx);
            information.fbestError = currentFitness2(1,fbestidx);
            information.segmentation = currentPopulation(fbestidx,:);
            [information.RMSE, information.RMSEp, information.MAXe] = obj.computeErrors(information.segmentation);
            [information.features, errors] = obj.computeMetrics(information.segmentation);
            information.yEstimada = obj.estimationSerie(information.segmentation);
            information.cuts = find(information.segmentation==1);
            [normCharac] = GMOTSS.normalizeFunction(information.features);
            C = initCentroids(normCharac,obj.parameters.k);
            [information.L, information.C] = dcKMeans(normCharac,obj.parameters.k,C,obj.parameters.iterClust);
            information.parameters = obj.parameters;
            information.fitness = [currentFitness1; currentFitness2];
            information.number_of_firstFront = numel(find(f==1));
            if information.number_of_firstFront > obj.parameters.nPobl,
                information.number_of_firstFront = obj.parameters.nPobl;
            end
            [spacing,m3,hv,HRS]=obj.evaluateFronts(information.fitness(:,1:information.number_of_firstFront));
            information.spacing = spacing;
            information.m3 = m3;
            information.hv = hv;
            information.HRS = HRS;
            information.fbest = fbestidx;
            information.degree = obj.parameters.polyDegree;
            
            % Best solution in Clustering terms
            [maxClustering fbestidx] = max(currentFitness1(1,:));
            fprintf('*******BEST SEGMENTATION IN CLUSTERING TERMS*******\n');
            fprintf('Indx: %d\n', fbestidx);
            fprintf('Fitness Clustering: %.15f\n', currentFitness1(1,fbestidx));
            fprintf('Fitness Error: %.15f\n', currentFitness2(1,fbestidx));
            % Information for the reporter
            informationBestClustering.fbestClustering = currentFitness1(1,fbestidx);
            informationBestClustering.fbestError = currentFitness2(1,fbestidx);
            informationBestClustering.segmentation = currentPopulation(fbestidx,:);
            [informationBestClustering.RMSE, informationBestClustering.RMSEp, informationBestClustering.MAXe] = obj.computeErrors(informationBestClustering.segmentation);
            [informationBestClustering.features, errors] = obj.computeMetrics(informationBestClustering.segmentation);
            informationBestClustering.yEstimada = obj.estimationSerie(informationBestClustering.segmentation);
            informationBestClustering.cuts = find(informationBestClustering.segmentation==1);
            [normCharac] = GMOTSS.normalizeFunction(informationBestClustering.features);
            C = initCentroids(normCharac,obj.parameters.k);
            [informationBestClustering.L, informationBestClustering.C] = dcKMeans(normCharac,obj.parameters.k,C,obj.parameters.iterClust);
            informationBestClustering.parameters = obj.parameters;
            informationBestClustering.fitness = [currentFitness1; currentFitness2];
            informationBestClustering.number_of_firstFront = numel(find(f==1));
            if informationBestClustering.number_of_firstFront > obj.parameters.nPobl,
                informationBestClustering.number_of_firstFront = obj.parameters.nPobl;
            end
            informationBestClustering.spacing = spacing;
            informationBestClustering.m3 = m3;
            informationBestClustering.hv = hv;
            informationBestClustering.HRS = HRS;
            informationBestClustering.fbest = fbestidx;
            informationBestClustering.degree = obj.parameters.polyDegree;
            
           % Best solution in Error terms
            [maxError fbestidx] = max(currentFitness2(1,:));
            fprintf('*******BEST SEGMENTATION IN ERROR TERMS*******\n');
            fprintf('Indx: %d\n', fbestidx);
            fprintf('Fitness Clustering: %.15f\n', currentFitness1(1,fbestidx));
            fprintf('Fitness Error: %.15f\n', currentFitness2(1,fbestidx));
            % Information for the reporter
            informationBestError.fbestClustering = currentFitness1(1,fbestidx);
            informationBestError.fbestError = currentFitness2(1,fbestidx);
            informationBestError.segmentation = currentPopulation(fbestidx,:);
            [informationBestError.RMSE, informationBestError.RMSEp, informationBestError.MAXe] = obj.computeErrors(informationBestError.segmentation);
            [informationBestError.features, errors] = obj.computeMetrics(informationBestError.segmentation);
            informationBestError.yEstimada = obj.estimationSerie(informationBestError.segmentation);
            informationBestError.cuts = find(informationBestError.segmentation==1);
            [normCharac] = GMOTSS.normalizeFunction(informationBestError.features);
            C = initCentroids(normCharac,obj.parameters.k);
            [informationBestError.L, informationBestError.C] = dcKMeans(normCharac,obj.parameters.k,C,obj.parameters.iterClust);
            informationBestError.parameters = obj.parameters;
            informationBestError.fitness = [currentFitness1; currentFitness2];
            informationBestError.number_of_firstFront = numel(find(f==1));
            if informationBestError.number_of_firstFront > obj.parameters.nPobl,
                informationBestError.number_of_firstFront = obj.parameters.nPobl;
            end
            informationBestError.spacing = spacing;
            informationBestError.m3 = m3;
            informationBestError.hv = hv;
            informationBestError.HRS = HRS;
            informationBestError.fbest = fbestidx;
            informationBestError.degree = obj.parameters.polyDegree;
        end
        
        %% Initialize the population
        function [newPopulation] = initialisePopulation(obj)
            newPopulation = false(obj.parameters.nPobl,obj.parameters.sizeChromosome);
            cont = 0;
            
            for i=1:obj.parameters.nPobl
                cont = 0;
                while (obj.parameters.sizeChromosome - cont) > obj.parameters.maxSeg
                    cont = cont + randi([obj.parameters.minSeg,obj.parameters.maxSeg]);
                    newPopulation(i,cont)=1;
                end
            end     
            
        end
        
        %% Non domination sort
        function [f, sorted_population, changedFitness1, changedFitness2, crowding_distances] = non_domination_sort(obj, population, fitness1, fitness2)
            N = numel(fitness1(1,:));
            crowding_distances = zeros(N,1);
            f = zeros(N,1);
            front = 1;
            F(front).f = [];
            individual = [];
            for i = 1 : N,
                individual(i).n = 0;
                individual(i).p = [];
                for j = 1 : N,
                    dom_less = 0;
                    dom_equal = 0;
                    dom_more = 0;
                    %We compare fitness 1 (clustering)
                    if fitness1(1,i) < fitness1(1,j),
                        dom_less = dom_less + 1;
                    elseif fitness1(1,i) == fitness1(1,j),
                        dom_equal = dom_equal+1;
                    else
                        dom_more = dom_more +1;
                    end
                    % We compare fitness 2 (error aproximation)
                    if fitness2(1,i) < fitness2(1,j),
                        dom_less = dom_less + 1;
                    elseif fitness2(1,i) == fitness2(1,j),
                        dom_equal = dom_equal+1;
                    else
                        dom_more = dom_more +1;
                    end
                    
                    if dom_more == 0 && dom_equal ~= 2,
                        individual(i).n = individual(i).n + 1;
                    elseif dom_less == 0 && dom_equal ~=2,
                        individual(i).p = [individual(i).p j];
                    end
                end
                if individual(i).n == 0,
                    f(i,1)=front;
                    F(front).f = [F(front).f i];
                end
            end
            
            while ~isempty(F(front).f)
                Q = [];
                for i = 1 : length(F(front).f)
                    if ~isempty(individual(F(front).f(i)).p)
                        for j = 1 : length(individual(F(front).f(i)).p)
                            individual(individual(F(front).f(i)).p(j)).n = individual(individual(F(front).f(i)).p(j)).n - 1;
                            if individual(individual(F(front).f(i)).p(j)).n == 0
                                f(individual(F(front).f(i)).p(j),1) = front + 1;
                                Q = [Q individual(F(front).f(i)).p(j)];
                            end
                        end
                    end
                end
                front = front + 1;
                F(front).f = Q;
            end
            sorted_population = zeros(numel(population(:,1)),numel(population(1,:)));
            changedFitness1 = zeros(1,numel(fitness1(1,:)));
            changedFitness2 = zeros(1,numel(fitness2(1,:)));
            [temp, index_of_fronts] = sort(f(:,1));
            for i = 1 : length(index_of_fronts),
                sorted_population(i,:)=population(index_of_fronts(i),:);
                changedFitness1(1,i)=fitness1(1,index_of_fronts(i));
                changedFitness2(1,i)=fitness2(1,index_of_fronts(i));
                f(i,1)=temp(i,1);
            end
            
            % Apply crowding
            current_index = 0;
            
            for front = 1 : (length(F) - 1),
                distance = 0;
                y = [];
                yFitness1 = [];
                yFitness2 = [];
                yf = [];
                previous_index = current_index + 1;
                for i = 1: length(F(front).f)
                    y(i,:) = sorted_population(current_index + i, :);
                    yFitness1(1,i) = changedFitness1(1,current_index + i);
                    yFitness2(1,i) = changedFitness2(1,current_index + i);
                    yf(i,1) = f(current_index + i,1);                    
                end
                current_index = current_index + i;
                
                %%FITNESS 1
                [sorted_based_on_objective, index_of_objectives] = sort(yFitness1(1,:));
                sorted_based_on_objective = [];
                sorted_based_on_objectiveFitness1 = [];
                sorted_based_on_objectiveFitness2 = [];
                for j = 1 : length(index_of_objectives)
                    sorted_based_on_objective(j,:) = y(index_of_objectives(j),:);
                    sorted_based_on_objectiveFitness1(1,j)=yFitness1(1,index_of_objectives(j));
                    sorted_based_on_objectiveFitness2(1,j)=yFitness2(1,index_of_objectives(j));
                end
                f_max = sorted_based_on_objectiveFitness1(1,length(index_of_objectives));
                f_min = sorted_based_on_objectiveFitness1(1,1);
                y(index_of_objectives(length(index_of_objectives)),obj.parameters.sizeChromosome+1) = Inf;
                y(index_of_objectives(1), obj.parameters.sizeChromosome + 1) = Inf;
                for j = 2 : length(index_of_objectives) - 1
                   next_obj  = sorted_based_on_objectiveFitness1(1,j + 1);
                   previous_obj  = sorted_based_on_objectiveFitness1(1,j - 1);
                   if (f_max - f_min == 0)
                       y(index_of_objectives(j),obj.parameters.sizeChromosome+1) = Inf;
                   else
                       y(index_of_objectives(j),obj.parameters.sizeChromosome+1) = (next_obj - previous_obj)/(f_max - f_min);
                   end
                end
                %%FITNESS 2
                [sorted_based_on_objective, index_of_objectives] = sort(yFitness2(1,:));
                sorted_based_on_objective = [];
                sorted_based_on_objectiveFitness1 = [];
                sorted_based_on_objectiveFitness2 = [];
                for j = 1 : length(index_of_objectives)
                    sorted_based_on_objective(j,:) = y(index_of_objectives(j),:);
                    sorted_based_on_objectiveFitness1(1,j)=yFitness1(1,index_of_objectives(j));
                    sorted_based_on_objectiveFitness2(1,j)=yFitness2(1,index_of_objectives(j));
                end
                f_max = sorted_based_on_objectiveFitness2(1,length(index_of_objectives));
                f_min = sorted_based_on_objectiveFitness2(1,1);
                y(index_of_objectives(length(index_of_objectives)),obj.parameters.sizeChromosome+2) = Inf;
                y(index_of_objectives(1), obj.parameters.sizeChromosome + 2) = Inf;
                for j = 2 : length(index_of_objectives) - 1
                   next_obj  = sorted_based_on_objectiveFitness2(1,j + 1);
                   previous_obj  = sorted_based_on_objectiveFitness2(1,j - 1);
                   if (f_max - f_min == 0)
                       y(index_of_objectives(j),obj.parameters.sizeChromosome+2) = Inf;
                   else
                       y(index_of_objectives(j),obj.parameters.sizeChromosome+2) = (next_obj - previous_obj)/(f_max - f_min);
                   end
                end
                distance = [];
                distance(:,1)= zeros(length(F(front).f),1);
                distance(:,1) = distance(:,1) + y(:,obj.parameters.sizeChromosome+1) + y(:,obj.parameters.sizeChromosome+1);

                y(:,obj.parameters.sizeChromosome+1)=[];
                y(:,obj.parameters.sizeChromosome+1)=[];
                
                sorted_populationAux=[];
                changedFitness1Aux=[];
                changedFitness2Aux=[];
                fAux=[];

                [sorted_based_on_distances, index_of_distances] = sort(distance(:,1),'descend');
                for j = 1: length(index_of_distances),
                    sorted_populationAux(j,:) = y(index_of_distances(j),:);
                    changedFitness1Aux(1,j) = yFitness1(1,index_of_distances(j));
                    changedFitness2Aux(1,j) = yFitness2(1,index_of_distances(j));
                    fAux(j,1) = yf(index_of_distances(j),1);
                end
                sorted_population(previous_index:current_index,:)=sorted_populationAux(:,:);
                changedFitness1(1,previous_index:current_index)=changedFitness1Aux(1,:);
                changedFitness2(1,previous_index:current_index)=changedFitness2Aux(1,:);
                f(previous_index:current_index,:)=fAux(:,1);
                crowding_distances(previous_index:current_index,:)=sorted_based_on_distances(:,1);
            end
            sorted_population; changedFitness1; changedFitness2; f; crowding_distances;
            % When the function is finished we have everything in:
            % sorted_population, changedFitness1, changedFitness2, f, crowding_distances
        end

        %% Single point cross over operator
        function [crossedPopulation,fitnessChanged] = crossover(obj,population,fitness)
            crossedPopulation = population;
            fitnessChanged = fitness;
            [nPop, nCutPoints] = size(population);
            for i=1:nPop
                %Crossover
                if rand()<obj.parameters.pCross,
                    % Find individuals to apply crossover
                    ind1 = i;
                    attempt = 1;
                    while attempt == 1,
                        ind2 = randi(nPop,1,1);
                        while ind2==i,
                            ind2 = randi(nPop,1,1);
                        end
                        attempt2=0;
                        while attempt2<=2,
                            crossPoint = randi(nCutPoints-3,1,1) + 1;
                            crossedPopulation(ind1,:) = [population(ind1,1:crossPoint) population(ind2,crossPoint+1:end)];
                            crossedPopulation(ind2,:) = [population(ind2,1:crossPoint) population(ind1,crossPoint+1:end)];
                            point1= find(crossedPopulation(ind1,1:crossPoint)==1,1,'last');
                            if isempty(point1)
                                point1=1;
                            end
                            point2= find(crossedPopulation(ind1,crossPoint+1:end)==1,1,'first');
                            point2= point2 + crossPoint;
                            if isempty(point2)
                                point2= obj.parameters.sizeChromosome;
                            end
                            point3= find(crossedPopulation(ind2,1:crossPoint)==1,1,'last');
                            if isempty(point3)
                                point3=1;
                            end
                            point4= find(crossedPopulation(ind2,crossPoint+1:end)==1,1,'first');
                            point4= point4 + crossPoint;
                            if isempty(point4)
                                point4= obj.parameters.sizeChromosome;
                            end
                           
                            if((point2 - point1)< obj.parameters.minSeg) || ((point4 - point3)< obj.parameters.minSeg),
                               attempt2=attempt2+1;
                               crossedPopulation(ind1,:)=population(ind1,:);
                               crossedPopulation(ind2,:)=population(ind2,:);
                            elseif ((point2 - point1)> obj.parameters.maxSeg) || ((point4 - point3)> obj.parameters.maxSeg),
                               attempt2=attempt2+1;
                               crossedPopulation(ind1,:)=population(ind1,:);
                               crossedPopulation(ind2,:)=population(ind2,:);
                             else
                               attempt2=3;
                               attempt=0;
                            end
                        end
                    end
                  
                    population(ind1,:) = crossedPopulation(ind1,:);
                    population(ind2,:) = crossedPopulation(ind2,:);
                    fitnessChanged(ind1) = NaN;
                    fitnessChanged(ind2) = NaN;                    
                end
            end
        end
        
        %% Mutation: Two types Add/Remove segment
        function [mutatedPopulation,changedFitness] = mutation(obj,population,fitness)
            mutatedPopulation = population;
            changedFitness = fitness;
            
            %t= cputime
            for i=1:size(population,1),
                %Mutate?
                if rand()<obj.parameters.pMut,
                    type = randi(2,1,1);
                    cLength=sum(population(i,:));
                    if type == 1,
                        for j=1:(obj.parameters.mutedPoints*cLength),
                            % add point
                            if rand()>0.5,
                                [ind] = find(mutatedPopulation(i,:)==0);
                                attempts=0;
                                flag = 1;
                                while attempts<3 && flag==1,
                                    point = ind(randi(numel(ind),1,1));
                                    mutatedPopulation(i,point) = 1;
                                    previousPoint = find(mutatedPopulation(i,1:point-1)==1,1,'last');
                                    nextPoint = find(mutatedPopulation(i,point+1:end)==1,1,'first');
                                    nextPoint = nextPoint + point;
                                    if isempty(previousPoint),
                                        previousPoint=1;
                                    end
                                    if isempty(nextPoint),
                                        nextPoint = obj.parameters.sizeChromosome;
                                    end
                                    if ((nextPoint - point)< obj.parameters.minSeg || (point - previousPoint)< obj.parameters.minSeg),
                                        mutatedPopulation(i,point)=0;
                                        attempts=attempts+1;
                                    else
                                        flag=0;
                                    end
                                end
                                if attempts < 3,
                                    changedFitness(i) = NaN;
                                end
                            % remove point
                            else
                                [ind] = find(mutatedPopulation(i,:)==1);
                                attempts=0;
                                flag=1;
                                
                                while attempts < 3 && flag==1,
                                    choice = randi(numel(ind),1,1);
                                    point = ind(choice);
                                    mutatedPopulation(i,point) = 0;
                                    if choice == 1,
                                        previousPoint = 1;
                                    else
                                        previousPoint = ind(choice-1);
                                    end
                                    if choice == numel(ind),
                                        nextPoint = obj.parameters.sizeChromosome;
                                    else
                                        nextPoint = ind(choice+1);
                                    end
                                    if nextPoint - previousPoint > obj.parameters.maxSeg,
                                        mutatedPopulation(i,point)=1;
                                        attempts=attempts+1;
                                    else
                                        flag=0;
                                    end
                                   
                                end
                                if attempts < 3,
                                    changedFitness(i) = NaN;
                                end
                            end
                        end
                    else
                        for j=1:(obj.parameters.mutedPoints*cLength),
                            [ind] = find(mutatedPopulation(i,:)==1);
                            [choice] = randi(numel(ind),1);
                            % Desplacement to the right
                            if rand()>0.5,
                                if choice == numel(ind),
                                    difference = numel(mutatedPopulation(i,:)) - ind(choice);
                                else
                                    difference = ind(choice+1) - ind(choice);
                                end
                                if (difference > obj.parameters.minSeg),
                                    attempts=0;
                                    flag=1;
                                    while flag==1 && attempts<3,
                                        mutatedPopulation(i,ind(choice)) = 0;
                                        displacement = randi(difference-obj.parameters.minSeg,1);
                                        mutatedPopulation(i,ind(choice)+displacement) = 1;
                                        if choice==1,
                                            previousPoint = 1;
                                        else
                                            previousPoint = ind(choice-1);
                                        end
                                        if ind(choice) - previousPoint + displacement > obj.parameters.maxSeg,
                                            mutatedPopulation(i,ind(choice)) = 1;
                                            mutatedPopulation(i,ind(choice)+displacement) = 0;
                                            attempts=attempts+1;
                                        else
                                            flag=0;
                                        end
                                    end
                                    if attempts < 3,
                                    changedFitness(i) = NaN;
                                    end
                                end 
                            else %Desplacement to the left
                                if choice == 1,
                                    difference = ind(choice);
                                else
                                    difference = ind(choice) - ind(choice-1);
                                end
                                if (difference > obj.parameters.minSeg),
                                    attempts=0;
                                    flag=1;
                                    while flag==1 && attempts<3,
                                        mutatedPopulation(i,ind(choice)) = 0;
                                        displacement = randi(difference-obj.parameters.minSeg,1);
                                        mutatedPopulation(i,ind(choice)-displacement) = 1;
                                        if choice == numel(ind),
                                            nextPoint = obj.parameters.sizeChromosome;
                                        else
                                            nextPoint = ind(choice+1);
                                        end
                                        if nextPoint - ind(choice) - displacement > obj.parameters.maxSeg,
                                            mutatedPopulation(i,ind(choice)) = 1;
                                            mutatedPopulation(i,ind(choice)-displacement) = 0;
                                            attempts=attempts+1;
                                        else
                                            flag=0;
                                        end
                                    end
                                    if attempts < 3,
                                    changedFitness(i) = NaN;
                                    end
                                end
                            end
                            
                        end
                        
                        
                    end
                end
                
            end
        end
        
        %% Selection of individuals
        function [newPopulation,newFitness1,newFitness2] = selection(obj,population,fitness1,fitness2)
              newFitness1 = zeros(1,obj.parameters.nPobl);
              newFitness2 = zeros(1,obj.parameters.nPobl);
              newPopulation = zeros(obj.parameters.nPobl,obj.parameters.sizeChromosome);
              
              newPopulation(1:obj.parameters.nPobl,:)=population(1:obj.parameters.nPobl,:);
              newFitness1(1,1:obj.parameters.nPobl)=fitness1(1,1:obj.parameters.nPobl);
              newFitness2(1,1:obj.parameters.nPobl)=fitness2(1,1:obj.parameters.nPobl);
        end 
        
        %% Binary tournament
        function [newPopulation, newFitness1, newFitness2] = binary_tournament(obj,population,fitness1,fitness2)
            nRows = obj.parameters.nPobl;
            nCols = obj.parameters.sizeChromosome;
            newPopulation = zeros(nRows,nCols);
            newFitness1 = zeros(1, nRows);
            newFitness2 = zeros(1, nRows);
            
            for i=1:numel(population(:,1)),
                ind1 = randi(nRows);
                ind2 = randi(nRows);
                while ind1 == ind2,
                    ind2 = randi(nRows);
                end
                if ind1 < ind2,
                    newPopulation(i,:)=population(ind1,:);
                    newFitness1(1,i)=fitness1(1,ind1);
                    newFitness2(1,i)=fitness2(1,ind1);
                else
                    newPopulation(i,:)=population(ind2,:);
                    newFitness1(1,i)=fitness1(1,ind2);
                    newFitness2(1,i)=fitness2(1,ind2);
                end
            end
        end
        
        %% Evaluation method
        function [fitness, fitness2] = evaluate(obj,population,oldFitness)
            fitness = zeros(1,size(population,1));
            fitness2 = zeros(1,size(population,1));
            for i=1:size(population,1),
                if isnan(oldFitness(i)),
                    [charac,errores] = obj.computeMetrics(population(i,:));
                    [normCharac] = GMOTSS.normalizeFunction(charac);
                    C = initCentroids(normCharac,obj.parameters.k);
                    [assignation, centroids] = dcKMeans(normCharac,obj.parameters.k,C,obj.parameters.iterClust);
                    fitness(i) = obj.fitnessF(normCharac,assignation,centroids);
                    fitness2(i) = obj.fitnessError(errores);
                end
            end
            
        end
        
        %% Value of the fitness function for a clustering assignment
        function [fitness] = fitnessF(obj, charac, assignation, centroids)
            
            realTargets = unique(assignation);
            newk = numel(realTargets);
            distanceToCentroid = zeros(1,newk);
            fitness = 0;
            for i=1:newk
                distanceToCentroid(i) = (mean(pdist2(charac(assignation==realTargets(i),:),centroids(realTargets(i),:))));
            end
            
            if newk < obj.parameters.k,
                fitness = 0;
                
            elseif obj.parameters.typeFitness == 1,
                %%  typeFitness == 1 => Davies-Bouldin Criterion => Minimize
                for i=1:newk
                    selected = 0;
                    for j=(i+1):newk
                        aux = (distanceToCentroid(i) + distanceToCentroid(j))/pdist2(centroids(realTargets(i),:),centroids(realTargets(j),:));
                        if aux>selected,
                            selected = aux;
                        end
                    end
                    fitness = fitness + selected;
                end
                fitness = fitness / newk;
                % Change the fitness so as we have to maximise it
                fitness = 1/(1+fitness);
           
            elseif obj.parameters.typeFitness == 2,
                %%  typeFitness == 2 => Dunn index => Maximize
                bad_solution=0;
                for i=1:newk,
                    num_elem = numel(find(assignation==i));
                    if num_elem < 2,
                        bad_solution=1;
                    end
                end
                if bad_solution==0,
                    selected = inf;
                    diam = zeros(1,newk);
                    indexesI = false(newk,size(assignation,1));
                    allDistances = squareform(pdist(charac));
                    for i=1:newk
                        indexesI(i,:) = assignation==realTargets(i);
                        numElI = sum(indexesI(i,:));
                        diam(i) = (1/(2*numElI*(numElI-1)))*sum(sum(allDistances(indexesI(i,:),indexesI(i,:))));
                    end
                    maxintraclusterdistance = max(diam);
                    for i=1:newk
                        for j=(i+1):newk
                            aux = min(min(allDistances(indexesI(i,:),indexesI(j,:))))/maxintraclusterdistance;
                            if aux < selected
                                selected = aux;
                            end

                        end
                    end
                    fitness = selected;
                else
                    fitness = 0;
                end
            
            elseif obj.parameters.typeFitness == 3,
                %% typeFitness == 3 => CALINSKI AND HARABASZ INDEX 1974 => Maximize
                Sw=zeros(4+obj.parameters.polyDegree,4+obj.parameters.polyDegree); %Degree
                Sb=zeros(4+obj.parameters.polyDegree,4+obj.parameters.polyDegree);
                m = mean(charac);
                totalNum = numel(assignation);
                for currentClass = 1:newk,
                    indexesI = assignation == realTargets(currentClass);
                    mk = mean(charac(indexesI,:));
                    Sw = Sw + (1/totalNum)*cov( (charac(indexesI,:)),1);
                    Sb = Sb + sum(indexesI)/totalNum*(mk-m)'*(mk-m);
                end

                
                fitness = (trace(Sb)/(obj.parameters.k-1))/(trace(Sw)/(totalNum-obj.parameters.k));
            
            elseif obj.parameters.typeFitness == 4,
                %% typeFitness == 4 => SSE => Minimize
                for i=1:newk
                    fitness = fitness + sum(pdist2(charac(assignation==realTargets(i),:),centroids(realTargets(i),:)));
                end
                fitness = fitness/numel(assignation);
                fitness = 1/(1+fitness);
            
            elseif obj.parameters.typeFitness == 5,
                %% typeFitness == 5 => SSE/DentreK => Minimize
                for i=1:newk
                    fitness = fitness + sum(pdist2(charac(assignation==realTargets(i),:),centroids(realTargets(i),:)));
                end
                dist_entre=0;
                for i=1:newk,
                    for j=(i+1):newk,
                        dist_entre = dist_entre + sum(pdist2(centroids(i,:),centroids(j,:)));  
                    end
                end
                fitness = fitness/numel(assignation);
                dist_entre = dist_entre / factorial(newk-1);
                
                fitness = fitness / dist_entre;
                
                fitness = 1/(1+fitness);
            
            elseif obj.parameters.typeFitness == 6,
                %% typeFitness == 6 => Silhouette index => Maximize
                fitness = 0;
                allDistances = squareform(pdist(charac));
                for k=1:newk,
                   indK = find(assignation==realTargets(k));
                   for i=1:numel(indK),
                       b=zeros(1,newk);
                       for l=1:newk,
                           b(1,l)=sum(allDistances(indK(i),assignation==realTargets(l)));
                           b(1,l)=b(1,l)/(sum(assignation==realTargets(l)));
                       end
                       a=b(1,k);
                       b(k)=[];
                       b=min(b);
                       fitness=fitness + ( (b-a)/(max(a,b)));
                   end                                 
                end
                fitness = fitness / numel(assignation);
                
            elseif obj.parameters.typeFitness == 7,
                %% typeFitness == 7 => gD33 index => Maximize
                fitness=obj.generalizedDunn(charac, assignation, centroids, 3, 3);
                
            elseif obj.parameters.typeFitness == 8,
                %% typeFitness == 8 => gD43 index => Maximize
                fitness=obj.generalizedDunn(charac, assignation, centroids, 4, 3);
                
            elseif obj.parameters.typeFitness == 9,
                %% typeFitness == 9 => gD53 index => Maximize
                fitness=obj.generalizedDunn(charac, assignation, centroids, 5, 3);
                
            else
                %% typeFitness == 10 => COP index => Minimize
                allDistances = squareform(pdist(charac));
                for k=1:newk,
                   num = sum(pdist2(charac(assignation==k,:),centroids(k,:)));
                   num = num / sum(assignation==k);
                   
                   maximos=zeros(1,newk);
                   
                   for i=1:newk,
                       maximos(1,i)=max(max(allDistances(assignation==k,assignation==i)));
                   end
                   maximos(k)=[];
                   den=min(maximos);
                   fitness=fitness + num/den;
                end
                fitness=fitness/numel(assignation);
                fitness = 1/(1+fitness);
            
            end
           
        end
        
        %% Value of the fitness function for error
        function [fitness2] = fitnessError(obj, errores)
            if obj.parameters.typeFitness2 == 1,
               fitness2 = mean(errores);
               fitness2 = sqrt(fitness2);
               fitness2 = 1/(1+fitness2);
            elseif obj.parameters.typeFitness2 == 2,
               fitness2 = sum(errores)/obj.parameters.sizeChromosome;
               fitness2 = sqrt(fitness2);
               fitness2 = 1/(1+fitness2); 
            else
               fitness2 = max(errores);
               fitness2 = sqrt(fitness2);
               fitness2 = 1/(1+fitness2);
            end
        end
        
        %% Value for \delta
        function [value] = deltaValue(obj,charac,assignation,centroids,type,k,l)
            value=0;
            if type == 3,
                value=pdist2(charac(assignation==k,:),charac(assignation==l,:));
                value=sum(sum(value))/(sum(assignation==k)*sum(assignation==l));
            elseif type == 4,
                value=pdist2(centroids(k,:),centroids(l,:));    
            else
                value=sum(pdist2(charac(assignation==k,:),centroids(k,:)));
                value2=sum(pdist2(charac(assignation==l,:),centroids(l,:)));
                value=value+value2;
                value=value/(sum(assignation==k)+sum(assignation==l));
            end
                    
            
        end
        
        %% Value for \Delta
        function [value] = deltaValueG(obj, charac, assignation, centroids, type, k)
            value=0;
            % Types 1 and 3, is correct for generalizedDunn
            if type==1,
                allDistances=squareform(pdist(charac(assignation==k)));
                value=max(max(allDistances));
            else
                value=sum(pdist2(charac(assignation==k,:),centroids(k,:)));
                value=2*value/(sum(assignation==k));
            end
            
        end
        
        %% Generalized Dunn
        function [gDunn] = generalizedDunn(obj, charac, assignation, centroids, type_delta, type_Delta)
           newk=numel(unique(assignation));
           min = obj.deltaValue(charac,assignation,centroids,type_delta,1,2);
           max = obj.deltaValueG(charac,assignation,centroids,type_Delta,1);
           for i=1:newk,
               for j=1:newk,
                  if i~=j,
                      newMin=obj.deltaValue(charac,assignation,centroids,type_delta,i,j);
                      if newMin < min,
                          min=newMin;
                      end
                  end
               end
               newMax=obj.deltaValueG(charac,assignation,centroids,type_Delta,i);
               if newMax > max,
                   max=newMax;
               end
           end
           gDunn=min/max;
        end
        
        %% Obtaining the metrics for all the segments
        function [characteristics, errores] = computeMetrics(obj,individual)
            ind = find(individual==1);
            nOfSegments = size(ind,2);
            characteristics = zeros(nOfSegments+1,4+obj.parameters.polyDegree);
            errores = zeros(nOfSegments+1,1);
            [characteristics(1,:), errores(1,1)] = obj.metrics(obj.data(1:ind(1)));
            
            % Number of Segments - 1
            for j=1:nOfSegments-1,
                [characteristics(j+1,:), errores(j+1,1)] = obj.metrics(obj.data(ind(j):ind(j+1)));
            end
            [characteristics(end,:), errores(end,1)] = obj.metrics(obj.data(ind(end):end));
        end 
        
        %% Values of the metrics for a given segment
        function [values, error] = metrics(obj,segment)
            % To extract statistic metrics
            values = zeros(1,4+obj.parameters.polyDegree);
            X = 1:numel(segment);
            X = transpose(X);
            c = (1/numel(segment));
            m = sum(segment)*c;
            a = (segment-m);
            varsegment = c * sum(a.^2);
            s = sqrt(varsegment);
            if s==0,
                % Variance and Skewness
                values([1,2]) = 0;
                % Kurtosis
                values(3) = -3;
                % Autocorrelation
                values(4) = 0;
                for i=5:(4+obj.parameters.polyDegree),
                    values(i) = 0;
                end
                error=0;
            else
                % Variance
                values(1) = varsegment;
                % Skewness
                values(2) = c*(sum(a.^3)/(s.^3));
                % Kurtosis
                values(3) = c*(sum(a.^4)/varsegment.^2) - 3;
                % Autocorrelation
                values(4) = sum((segment(1:end-1) - m) .* (segment(2:end) - m))/varsegment; %cambio el 6 por 5 para quitar el MSE

                p = polyfit(X,segment,obj.parameters.polyDegree);
                for i=5:(4+obj.parameters.polyDegree),
                    values(i) = p(i-4);
                end
                % Error
                % estimated(:,1)= p(1)*X(:,1).*X(:,1) + p(2)*X(:,1) + p(3);
                estimated(:,1) = polyval(p,X(:,1));
                error = estimated(:,1) - segment(:,1);
                error = error.*error;
                if obj.parameters.typeFitness2 == 1,
                    N=numel(error);
                    error=sum(error);
                    error=error/N;
                elseif obj.parameters.typeFitness2 == 2,
                    error=sum(error);
                else
                    error=max(error);
                end
            end    
            
        end
        
        %% Error values for reporter
        function [RMSE, RMSEp, MAXe] = computeErrors(obj, individual)
            ind = find(individual==1);
            nOfSegments = size(ind,2);
            errorsRMSE = zeros(nOfSegments+1,1);
            errorsRMSEp = zeros(nOfSegments+1,1);
            errorsMAXe = zeros(nOfSegments+1,1);
            [errorsRMSE(1,1), errorsRMSEp(1,1), errorsMAXe(1,1)] = obj.errors(obj.data(1:ind(1)));
            
            % Number of Segments - 1
            for j=1:nOfSegments-1,
                [errorsRMSE(j+1,1), errorsRMSEp(j+1,1), errorsMAXe(j+1,1) ] = obj.errors(obj.data(ind(j):ind(j+1)));
            end
            [errorsRMSE(end,1), errorsRMSEp(end,1), errorsMAXe(end,1)] = obj.errors(obj.data(ind(end):end));
            
            MSE=mean(errorsRMSE);
            RMSE=sqrt(MSE);
            
            MSEp=sum(errorsRMSEp)/obj.parameters.sizeChromosome;
            RMSEp=sqrt(MSEp);
            
            MAXe=max(errorsMAXe);
            MAXe=sqrt(MAXe);
            
        end
        
        %% Error values for a given segment
        function [errorMSE, errorSSE, errorMAXe] = errors(obj, segment)
            X=1:numel(segment);
            X=transpose(X);
			p = polyfit(X,segment,obj.parameters.polyDegree);
            estimated(:,1) = polyval(p,X(:,1));
            error = estimated(:,1) - segment(:,1);
            error = error.*error;
            %RMSE
            N=numel(error);
            errorMSE=sum(error);
            errorMSE=errorMSE/N;
            %RMSEp
            errorSSE=sum(error);
            %MAXe
            errorMAXe=max(error);            
        end
        
        %% Estimating time series
        function [yEstimado] = estimationSerie(obj, individual)
            ind = find(individual==1);
            nOfSegments = size(ind,2);
            yEstimado = zeros(numel(obj.data),1);
            yEstimado(1:ind(1))= obj.estimationSegment(obj.data(1:ind(1)));
           
            % Number of Segments - 1
            for j=1:nOfSegments-1,
                yEstimado(ind(j):ind(j+1)) = obj.estimationSegment(obj.data(ind(j):ind(j+1)));
            end
            yEstimado(ind(end):end) = obj.estimationSegment(obj.data(ind(end):end));
        end
        
        %% Estimation of the segment
        function [yEstimado] = estimationSegment(obj, segment)
            yEstimado = zeros(numel(segment),1);
				X=1:numel(segment);
                X=transpose(X);
				p = polyfit(X,segment,obj.parameters.polyDegree);
                yEstimado(:,1) = polyval(p,X);
 				% yEstimado(:,1) = p(1)*X(:,1).*X(:,1) + p(2)*X(:,1) + p(3);
        end
        
        %% Metrics for the evaluation of the Pareto front
        function [spacing,m3,hv,HRS] = evaluateFronts(obj,fitness)
            spacing = GMOTSS.evaluateSpacing(fitness);
            m3 = GMOTSS.evaluateM3(fitness);
            hv = GMOTSS.evaluateHyperarea(fitness);
            HRS = GMOTSS.evaluateHRS(fitness);
        end
        
    end
    
    methods (Static = true)
        
        %% Standarize a matrix
        function [XN, XMeans, XStds] = standarizeFunction(X,XMeans,XStds)
            
            if (nargin<3)
                XStds = std(X);
            end
            if (nargin<2)
                XMeans = mean(X);
            end
            XN = zeros(size(X));
            for i=1:size(X,2)
                XN(:,i) = (X(:,i) - XMeans(i)) / XStds(i);
            end
        end
        
        %% Normalize a matrix
        function xnorm = normalizeFunction(x)
            mins = min(x);
            maxs = max(x);
            tama = size(x,1);
            xnorm = (x-repmat(mins,tama,1))./(repmat(maxs,tama,1)-repmat(mins,tama,1));
            
        end
        
        %% Spacing Metric
        function [spacing] = evaluateSpacing(fitness)
            % Minimize
            num = numel(fitness(1,:));
            auxiliar_distances = zeros(1,num);
            d = zeros(1,num);

            for i=1:num,
                auxiliar_distances(1,:)=0;
                for j=1:num,
                    auxiliar_distances(1,j)=abs(fitness(1,i)-fitness(1,j))+abs(fitness(2,i)-fitness(2,j));
                    if i==j,
                       auxiliar_distances(1,j)=inf;
                    end
                end
                d(i)=min(auxiliar_distances);
            end

            spacing = sum((d - mean(d)).*(d - mean(d)));
            spacing = sqrt(spacing/(num-1));    
        end

        %% M3 Metric
        function [m3] = evaluateM3(fitness)
            % Maximize
            maximo1 = max(fitness(1,:));
            minimo1 = min(fitness(1,:));

            maximo2 = max(fitness(2,:));
            minimo2 = min(fitness(2,:));

            m3 = sqrt((maximo1-minimo1) + (maximo2-minimo2));

        end

        %% HyperArea Metric
        function [H] = evaluateHyperarea(fitness)
            % Maximize
            [fit1 ind]=sort(fitness(1,:));

            fit2 = fitness(2,ind);

            H = fit1(1)*fit2(1);
            for i=2:numel(fit1),
                H = H + ((fit1(i) - fit1(i-1))*fit2(i)); 
            end

        end

        %% HRS Metric
        function [HRS] = evaluateHRS(fitness)
            % Minimize
            num = numel(fitness(1,:));
            auxiliar_distances = zeros(1,num);
            d = zeros(1,num);

            for i=1:num,
                auxiliar_distances(1,:)=0;
                for j=1:num,
                    auxiliar_distances(1,j)=abs(fitness(1,i)-fitness(1,j))+abs(fitness(2,i)-fitness(2,j));
                    if i==j,
                       auxiliar_distances(1,j)=inf;
                    end
                end
                d(i)=min(auxiliar_distances);
            end

            HRS = max(d)/mean(d);   
        end
 
    end
    
    
end
