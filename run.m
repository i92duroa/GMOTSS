files = {'IBEX.txt', 'Donoho-Johnstone.txt', 'b46001.txt', 'MIT_BIH_Arrhythmia_108.txt'};
warning('off')
nOfruns = 1;

for file=1:numel(files),    
    for typeFitness1 = 1:10,
        for typeFitness2 = 1:3,
            serie=load(['time_series' filesep char(files(file))]);
            c = clock;
            folder = ['reports' filesep num2str(c(1)) '-' num2str(c(2)) '-'  num2str(c(3)) '-' num2str(c(4)) '-' num2str(c(5)) '-' num2str(uint8(c(6)))];
            mkdir('.',folder)
            for i=1:nOfruns,
                alg(i) = GMOTSS;
                alg(i).dataFile = char(files(file));
                alg(i).parameters.maxSeg=150;
                alg(i).parameters.polyDegree=2;
                alg(i).parameters.minSeg=alg(i).parameters.polyDegree+2;
                alg(i).parameters.numIt = 3;
                alg(i).parameters.k = 5;
                alg(i).parameters.nPobl = 100;
                % Fitness Clustering
                alg(i).parameters.typeFitness = typeFitness1;
                % Fitness Error
                alg(i).parameters.typeFitness2 = typeFitness2;
                alg(i).parameters.seed = i*10;
                alg(i).parameters.mutedPoints = 0.2;
                %0.8%
                alg(i).parameters.pCross = 0.8;
                %0.2%
                alg(i).parameters.pMut = 0.2;
                [information(i),informationClustering(i),informationError(i)] = alg(i).runEvolutive(serie(:,2));
                mkdir('.',[folder filesep num2str(i)])
                mkdir('.',[folder filesep num2str(i) filesep 'bestGlobal'])
                mkdir('.',[folder filesep num2str(i) filesep 'bestClustering'])
                mkdir('.',[folder filesep num2str(i) filesep 'bestError'])
                saveAll(information(i),char(files(file)),[folder filesep num2str(i) filesep 'bestGlobal'], alg(i).getParameters());
                saveAll(informationClustering(i),char(files(file)),[folder filesep num2str(i) filesep 'bestClustering'], alg(i).getParameters());
                saveAll(informationError(i),char(files(file)),[folder filesep num2str(i) filesep 'bestError'], alg(i).getParameters());
            end

            masterSaveAll(folder,'resultsBestGlobal',information);
            masterSaveAll(folder,'resultsBestClustering',informationClustering);
            masterSaveAll(folder,'resultsBestError',informationError);
        end
    end
end


