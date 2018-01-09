function saveAll(model,dataset,repsuffix, paramstr)
            outputPDF = true;
            outputSVG = true;

            addpath external_tools/export_fig/
            addpath external_tools/plot2svg/
            
            % X is the original time series
            X = load(['time_series' filesep dataset]);
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
                fprintf(f, 'Front Evaluation Metrics\n');
                fprintf(f, 'Spacing;%f\n', model.spacing);
                fprintf(f, 'M3;%f\n', model.m3);
                fprintf(f, 'HV;%f\n', model.hv);
                fprintf(f, 'HRS;%f\n', model.HRS);
                fprintf(f, 'GMO parameters\n');
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
                
                %% FRONT
                f5=figure;
                xmin=min(model.fitness(1,:));
                xmax=max(model.fitness(1,:));
                ymin=min(model.fitness(2,:));
                ymax=max(model.fitness(2,:));
                xlabel('Clustering Fitness','fontsize',13);
                ylabel('Error Fitness','fontsize',13);
                hold on;
                xlim([xmin xmax]);
                ylim([ymin ymax]);
                
                plot(model.fitness(1,1:model.number_of_firstFront),model.fitness(2,1:model.number_of_firstFront),'bo');
                plot(model.fitness(1,model.fbest:model.fbest),model.fitness(2,model.fbest:model.fbest),'b*');
                if model.number_of_firstFront < numel(model.fitness(1,:)),
                    plot(model.fitness(1,model.number_of_firstFront+1:end),model.fitness(2,model.number_of_firstFront+1:end),'co');
                end
                
                if outputPDF
                    export_fig([outputFile '_FRONT.pdf'],'-pdf','-transparent');
                end
                if outputSVG
                    plot2svg([outputFile '_FRONT.svg']);
                end

                hold off;
                close all;
            end

            %% Limpiar classpath
            rmpath external_tools/plot2svg/ 
            rmpath external_tools/export_fig/ 
        end