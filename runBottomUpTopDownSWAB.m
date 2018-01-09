file = ['time_series' filesep 'synthetic.txt'];
warning('off')
serie = load(file);
i=1;
for i=1:35,
    alg(i) = TRADTSS;
    alg(i).dataFile = file;
    alg(i).parameters.k = 4;
    alg(i).parameters.maxError = 0.000038+(i*0.000002);
%    alg(i).parameters.maxError = 0.000057;
    alg(i).parameters.typeFitness = 6;
    alg(i).parameters.typeAlgorithm = 4;
    [information] = alg(i).runTraditionals(serie(:,2));
    fprintf('*****************\n');
    fprintf('Resultados con max error: %f\n',alg(i).parameters.maxError);
    NSEG = size(information.cuts,2)+1;
    fprintf('NSEG: %f\n', NSEG);
    fprintf('SI: %f\n', information.fBestClustering);
    fprintf('RMSE: %f\n', information.RMSE);
    fprintf('RMSEp: %f\n', information.RMSEp);
    fprintf('MAXe: %f\n', information.MAXe);
    fprintf('*****************\n');
end