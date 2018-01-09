 classdef TRADTSS < handle
    
    properties
        name_parameters = {'iterClust','k','maxError', 'polyDegree', 'sizeChromosome'}
        dataFile
        data
        parameters
    end
    
    methods
        
        %% Constructor
        function obj = TRADTSS()
            obj.defaultParameters();
        end

        %% Default parameters of the class
        function obj = defaultParameters(obj)
            % Number of iterations for k-means
            obj.parameters.iterClust = 20;
            % k for k-means
            obj.parameters.k = 5;
            % Maximum error of each segment
            obj.parameters.maxError = 10;
            % Degree polynomial
            obj.parameters.polyDegree = 2;
            % Type of Algorithm
            % Type 1: Sliding Windows
            % Type 2: Top-Down
            % Type 3: Bottom-Up
            % Type 4: SWAB
            obj.parameters.typeAlgorithm = 1;
            % Fitness for clustering 1-10
            obj.parameters.typeFitness = 2;
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
        
        %% Value of the error for a given segment
        function [error] = calculate_error(obj,segment)
            X=1:numel(segment);
            X=transpose(X);
            p=polyfit(X,segment,obj.parameters.polyDegree);
            estimated(:,1) = polyval(p,X(:,1));
            error = estimated(:,1) - segment(:,1);
            error = error.*error;
            N=numel(error);
            error=sum(error);
            error=error/N;
        end
        
        %% Value of the error for a given chromosome
        function [error] = calculate_full_error(obj,chromosome)
            ind = find(chromosome==1);
            error=0;
            error = error + obj.calculate_error(obj.data(1:ind(1)));
            for i=1:numel(ind)-1,
                error = error + obj.calculate_error(obj.data(ind(i):ind(i+1)));
            end
            error = error + obj.calculate_error(obj.data(ind(end):end));
            error = error/(numel(ind)+1);
            error = sqrt(error);
        end
        
        %% Sliding Window algorithm
        function [chromosome] = Sliding_Window(obj,segment,max_error)
            chromosome = zeros(1,numel(segment));
            left=1;
            right=2;
            while right < numel(segment),
                while (right <= numel(segment)) && (obj.calculate_error(segment(left:right)) < max_error),
                    right=right+1;
                end
                left=right-1;
                chromosome(left)=1;
            end
            chromosome(numel(segment))=0;
        end
        
        %% Top-Down algorithm
        function [chromosome] = Top_Down(obj,segment,max_error)
            chromosome = zeros(1,numel(segment));
            best_so_far = inf;
            for i=2:numel(segment)-1,
                error_partition = obj.calculate_error(segment(1:i)) + obj.calculate_error(segment(i:end));
                if error_partition < best_so_far,
                    breakpoint = i;
                    best_so_far = error_partition;
                end
            end
            
            if obj.calculate_error(segment(1:breakpoint)) > max_error,
                chromosome(1:breakpoint) = obj.Top_Down(segment(1:breakpoint),max_error);
            end
            
            if obj.calculate_error(segment(breakpoint:end)) > max_error,
                chromosome(breakpoint:end) = obj.Top_Down(segment(breakpoint:end),max_error);
            end
            chromosome(1,breakpoint)=1;            
        end
        
        %% Bottom_Up algorithm
        function [chromosome] = Bottom_Up(obj,segment,max_error)
            chromosome = ones(1,numel(segment));
            merge_cost = zeros(1,numel(segment));
            merge_cost(1,1)=inf;
            merge_cost(1,end)=inf;
            indexes = find(chromosome==1);
            for i=2:numel(segment)-1,
                merge_cost(i)=obj.calculate_error(segment(indexes(i-1):indexes(i+1))); 
            end
            
            while min(merge_cost) < max_error,
                ind = find(merge_cost==min(merge_cost),1);
                if ind > 2,
                    merge_cost(ind-1)=obj.calculate_error(segment(indexes(ind-2):indexes(ind+1)));
                end
                if ind < (numel(merge_cost)-1),
                    merge_cost(ind+1)=obj.calculate_error(segment(indexes(ind-1):indexes(ind+2)));
                end
                chromosome(1,indexes(ind))=0;
                merge_cost(ind)=[];
                indexes(ind)=[];                
            end
            chromosome(1,1)=0;
            chromosome(1,end)=0;           
        end
        
        %% SWAB algorithm
        function [chromosome] = SWAB(obj,segment,max_error,buffer_size)
            chromosome = zeros(1,numel(segment));
            left = 1;
            right = buffer_size;
            salida = 0;
            while salida==0 && left~=right,
                initial_left = left;
                initial_right = right;
                %chromosome(left:right) = 0;
                chromosome(left:right) = obj.Bottom_Up(segment(left:right),max_error);
                chromosome(left)=1;
                chromosome(right)=1;
                ind = find(chromosome(left+1:right)==1,1);
                chromosome(left+ind)=1;
                left=left+ind;
                %chromosome(right:end) = 0;
                chromosome(right:end) = obj.Sliding_Window(segment(right:end),max_error);
                chromosome(end) = 1;
                ind = find(chromosome(right:end)==1,1);
                right = right + ind - 1; 
                if(initial_left == left) && (initial_right == right),
                    salida=1;
                end
            end
            chromosome(1,1)=0;
            chromosome(1,end)=0;
        end
        
        %% Main Traditional algorithm
        function [information] = runTraditionals(obj,serie)
            addpath kmeans/
            obj.data = serie;
            obj.parameters.sizeChromosome = length(serie);
            
            switch obj.parameters.typeAlgorithm
                case 1
                    chromosome = obj.Sliding_Window(obj.data,obj.parameters.maxError);
                case 2
                    chromosome = obj.Top_Down(obj.data,obj.parameters.maxError);
                case 3
                    chromosome = obj.Bottom_Up(obj.data,obj.parameters.maxError);
                case 4
                    chromosome = obj.SWAB(obj.data,obj.parameters.maxError,295);
                otherwise
                    disp('This method does not exist');
            end
            
            error = obj.calculate_full_error(chromosome);
            information.fBestError=1/(1+error);
            information.segmentation=chromosome;
            [information.RMSE, information.RMSEp, information.MAXe] = obj.computeErrors(information.segmentation);
            information.features = obj.computeMetrics(information.segmentation);
            information.yEstimada = obj.estimationSerie(information.segmentation);
            information.cuts = find(information.segmentation==1);
            [normCharac] = TRADTSS.normalizeFunction(information.features);
            C = initCentroids(normCharac,obj.parameters.k);
            [information.L, information.C] = dcKMeans(normCharac,obj.parameters.k,C,obj.parameters.iterClust);
            information.fBestClustering = obj.fitnessF(normCharac,information.L,information.C);
            information.parameters = obj.parameters;
            information.degree = obj.parameters.polyDegree;
            
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
            [characteristics(j+2,:), errores(j+2,1)] = obj.metrics(obj.data(ind(end):end));
        end 
        
        %% Values of the metrics for a given segment
        function [values, error] = metrics(obj,segment)
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
                    N=numel(error);
                    error=sum(error);
                    error=error/N;
                    

                end
            
        end
        
        %% Errors values for reporter
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
        
        %% Errors values for a given segment
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
        
        function saveAll(obj,model,dataset,repsuffix, paramstr)
            outputPDF = true;
            outputSVG = true;

            addpath external_tools/export_fig/
            addpath external_tools/plot2svg/
            
            % X is the original time series
            X = load(dataset);
            % X2 is the estimated time series
            X2 = zeros(numel(X(:,1)),2);
            X2(:,1)=X(:,1);
            X2(:,2)=model.yEstimada(:,1);
            
            cuts = model.cuts;
            cutscero = [1 cuts length(X)];
            
            nOfClusters = size(model.C,1);
            nOfCuts = size(cuts,2);
            
            fe=model.features;
            class = model.L;
            centroids = model.C;            
            degree = model.degree;
            outputFile = [repsuffix filesep dataset];
            
            if numel(unique(model.L))==nOfClusters
                
                % Unnormalized clusters
                centroidsNorm = centroids;
                for j=1:(4+degree),
                   maximo=max(fe(:,j));
                   minimo=min(fe(:,j));
                   for i=1:nOfClusters,
                       centroids(i,j)=(centroids(i,j)*(maximo-minimo))+minimo;
                   end
                    
                end
              
                % Segments
                f = fopen([outputFile 'segments.csv'], 'wt');
                for i=1:nOfCuts+1,
                    fprintf(f, '%d;%d;', i, cutscero(i));
                    for j=1:(4+degree),
                        fprintf(f, '%f;', fe(i,j));
                    end
                    fprintf(f, '%d\n', class(i));
                end
                fclose(f);

                % Clusters
                f = fopen([outputFile '_centroids.csv'], 'wt');
                for i=1:nOfClusters
                    fprintf(f, '%d;', i);
                    for j=1:(4+degree-1),
                        fprintf(f, '%f;', centroids(i,j));
                    end
                    fprintf(f, '%f\n', centroids(i,(4+degree)));
                end
                fclose(f);
                
                % Clusters Norm
                f = fopen([outputFile '_centroidsNorm.csv'], 'wt');
                for i=1:nOfClusters
                    fprintf(f, '%d;', i);
                    for j=1:(4+degree-1),
                        fprintf(f, '%f;', centroidsNorm(i,j));
                    end
                    fprintf(f, '%f\n', centroidsNorm(i,(4+degree)));
                end
                fclose(f);
                
                % Estimated Time Series
                f = fopen([outputFile '_estimated.txt'], 'wt');
                for i=1:numel(X(:,1)),
                   fprintf(f, '%d\t%f\n', X2(i,1), X2(i,2)); 
                end
                fclose(f);

                % Other information
                f = fopen([outputFile '_info.csv'], 'wt');
                fprintf(f, 'Number of Cuts;%d\n',nOfCuts);
                fprintf(f, 'Number of Segments;%d\n',nOfCuts+1);
                fprintf(f, 'Clustering Fitness value;%f\n',model.fbestClustering);
                fprintf(f, 'Error Fitness value;%f\n', model.fbestError);
                fprintf(f, 'RMSE;%f\n', model.RMSE);
                fprintf(f, 'RMSEp;%f\n', model.RMSEp);
                fprintf(f, 'MAXe;%f\n', model.MAXe);
                fprintf(f, 'TRAD parameters\n');
                fprintf(f, '%s\n',paramstr);
                fclose(f);
                
                

               %% GRAPHICS %%
                markers = {'r-','g-','b-','m-','c-','y-','k-'};
                % Set lim of all graphs
                xmin=0;
                xmax=max(max(X(:,1)),max(X2(:,1)));
                ymin=min(min(X(:,2)),min(X2(:,2)));
                ymax=max(max(X(:,2)),max(X2(:,2)));
                
                
                %% ORIGINAL TIME SERIES WITHOUT CUT POINTS
                initialPoint = 1;
                f=figure;
                set(f, 'Position', [50 50 1200 400])
                hold on;
                set(gca,'fontsize',14,'LineWidth',1) ;
                ylabel('Y Axis','fontsize',14')
                xlabel('X Axis','fontsize',14)

                xlim([xmin xmax]);
                ylim([ymin ymax]);

                h = zeros((nOfClusters),1);
                cuts = [cuts length(X)];
                nOfCuts = size(cuts,2);

                for c=1:(nOfCuts)
                    ht = plot(X(initialPoint:cuts(c),1),...
                         X(initialPoint:cuts(c),2),markers{model.L(c)},'linewidth',1);
                    initialPoint = cuts(c);

                    % To build the correct legend
                    h(model.L(c))=ht;
                end

                % Build legend string
                legendStr = cell(nOfClusters,1);
                for j=1:nOfClusters
                    legendStr{j,1} = sprintf('Cluster %d',j); 
                end

                legend(h,legendStr,'Location','NorthWest')

                hold off;

                if outputPDF
                    export_fig([outputFile '.pdf'],'-pdf','-transparent');
                end
                if outputSVG
                    plot2svg([outputFile '.svg']);
                end
                close all;

                %% ORIGINAL TIME SERIES WITH CUT POINTS
                f2=figure;
                initialPoint = 1;
                set(f2, 'Position', [50 50 1200 400])
                hold on;
                set(gca,'fontsize',14,'LineWidth',1) ;
                ylabel('Y Axis','fontsize',14)
                xlabel('X Axis','fontsize',14)
                % Set lim of graph
                xlim([xmin xmax]);
                ylim([ymin ymax]);

                h2 = zeros((nOfClusters),1);
                for c=1:(nOfCuts)
                    ht2 = plot(X(initialPoint:cuts(c),1),...
                         X(initialPoint:cuts(c),2),markers{model.L(c)},'linewidth',1);
                    initialPoint = cuts(c);

                    % To build the correct legend
                    h2(model.L(c))=ht2;
                end
                legend(h2,legendStr,'Location','NorthWest')


                timecuts = X(cuts,1);
                timecuts = timecuts/1.^3;
                for i=1:numel(timecuts),
                    plot([timecuts(i) timecuts(i)], [min((X(:,2))) max((X(:,2)))], 'k--');
                end

                if outputPDF
                    export_fig([outputFile 'cuts.pdf'],'-pdf','-transparent');
                end
                if outputSVG
                    plot2svg([outputFile 'cuts.svg']);
                end

                hold off;
                close all;
                
                %% ESTIMATED TIME SERIES WITHOUT CUT POINTS
                initialPoint = 1;
                f3=figure;
                set(f3, 'Position', [50 50 1200 400])
                hold on;
                set(gca,'fontsize',14,'LineWidth',1) ;
                ylabel('Y Axis','fontsize',14')
                xlabel('X Axis','fontsize',14)
                % Set lim of graph
                xlim([xmin xmax]);
                ylim([ymin ymax]);

                h3 = zeros((nOfClusters),1);
                nOfCuts = size(cuts,2);

                for c=1:(nOfCuts)
                    ht3 = plot(X2(initialPoint:cuts(c),1),...
                         X2(initialPoint:cuts(c),2),markers{model.L(c)},'linewidth',1);
                    initialPoint = cuts(c);

                    % To build the correct legend
                    h3(model.L(c))=ht3;
                end
                legend(h3,legendStr,'Location','NorthWest')

                hold off;

                if outputPDF
                    export_fig([outputFile '_estimated.pdf'],'-pdf','-transparent');
                end
                if outputSVG
                    plot2svg([outputFile '_estimated.svg']);
                end
                close all;

                %% ESTIMATED TIME SERIES WITH CUT POINTS
                f4=figure;
                initialPoint = 1;
                set(f4, 'Position', [50 50 1200 400])
                hold on;
                set(gca,'fontsize',14,'LineWidth',1) ;
                ylabel('Y Axis','fontsize',14)
                xlabel('X Axis','fontsize',14)
                % Set lim of graph
                xlim([xmin xmax]);
                ylim([ymin ymax]);

                h4 = zeros((nOfClusters),1);
                for c=1:(nOfCuts)
                    ht4 = plot(X2(initialPoint:cuts(c),1),...
                         X2(initialPoint:cuts(c),2),markers{model.L(c)},'linewidth',1);
                    initialPoint = cuts(c);

                    % To build the correct legend
                    h4(model.L(c))=ht4;
                end
                legend(h4,legendStr,'Location','NorthWest')

                timecuts = X2(cuts,1);
                timecuts = timecuts/1.^3;
                for i=1:numel(timecuts),
                    plot([timecuts(i) timecuts(i)], [min((X2(:,2))) max((X2(:,2)))], 'k--');
                end

                if outputPDF
                    export_fig([outputFile '_estimatedCuts.pdf'],'-pdf','-transparent');
                end
                if outputSVG
                    plot2svg([outputFile '_estimatedCuts.svg']);
                end

                hold off;
                close all;               
            end

            %% Limpiar classpath
            rmpath external_tools/plot2svg/ 
            rmpath external_tools/export_fig/ 
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
 
    end
    
    
end
