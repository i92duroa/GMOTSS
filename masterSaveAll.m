%% Summary of the seeds
function masterSaveAll(folder,filename,model)
    if numel(model)>1,        
        for i=1:numel(model),
            resultsRMSE(i) = model(i).RMSE;
            resultsRMSEp(i) = model(i).RMSEp;
            resultsMAXe(i) = model(i).MAXe;
            resultsClustering(i) = model(i).fbestClustering;
            numSegments(i) = size(model(i).cuts,2) + 1;
            spacing(i) = model(i).spacing;
            m3(i) = model(i).m3;
            hv(i) = model(i).hv;
            HRS(i) = model(i).HRS;
        end

        fid = fopen([folder filesep filename '.txt'],'wt');
        fprintf(fid,'MeanNumberSegments, StdNumberSegments, MeanFitnessClustering, StdFitnessClustering, MeanRMSE, StdRMSE, MeanRMSEp, StdRMSEp, MeanMAXe, StdMAXe\n');
        fprintf(fid,'%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n',mean(numSegments),std(numSegments),mean(resultsClustering),std(resultsClustering),mean(resultsRMSE),std(resultsRMSE),mean(resultsRMSEp),std(resultsRMSEp),mean(resultsMAXe),std(resultsMAXe));
        fprintf(fid,'MeanSpacing, StdSpacing, MeanM3, StdM3, MeanHV, StdHV, MeanHRS, StdHRS\n');
        fprintf(fid,'%f,%f,%f,%f,%f,%f,%f,%f\n', mean(spacing),std(spacing),mean(m3),std(m3),mean(hv),std(hv),mean(HRS),std(HRS));
        fclose(fid);
        save([folder filesep filename '.mat'], 'model');
    else
        fid = fopen([folder filesep filename '.txt'],'wt');
        fprintf(fid,'MeanNumberSegments, StdNumberSegments, MeanFitnessClustering, StdFitnessClustering, MeanRMSE, StdRMSE, MeanRMSEp, StdRMSEp, MeanMAXe, StdMAXe\n');
        numSegments = size(model.cuts,2) + 1;
        fprintf(fid,'%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n',numSegments,0,model.fbestClustering,0,model.RMSE,0,model.RMSEp,0,model.MAXe,0);
        fprintf(fid,'MeanSpacing, StdSpacing, MeanM3, StdM3, MeanHV, StdHV, MeanHRS, StdHRS\n');
        fprintf(fid,'%f,%f,%f,%f,%f,%f,%f,%f\n', model.spacing,0,model.m3,0,model.hv,0,model.HRS,0);
        fclose(fid);
    end
end            