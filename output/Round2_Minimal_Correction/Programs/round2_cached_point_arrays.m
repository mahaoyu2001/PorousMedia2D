function [meanMean,meanSD,varMean,varSD,skeMean,skeSD,kurMean,kurSD,sampleCount] = ...
    round2_cached_point_arrays(vol,methods,componentNames)
t=round2_read_statistics('point');
fields={'Mean','Mean_SD','Variance','Variance_SD','Skewness','Skewness_SD','Kurtosis','Kurtosis_SD','NumRealizations'};
out=cell(1,9);
for k=1:9, out{k}=nan(numel(methods),numel(componentNames),numel(vol)); end
for im=1:numel(methods)
    for i=1:numel(componentNames)
        for j=1:numel(vol)
            mask=string(t.Method)==string(methods{im}) & ...
                string(t.StressComponent)==string(componentNames(i)) & abs(t.Porosity-vol(j))<1e-12;
            assert(nnz(mask)==1,'Missing or duplicate point-process summary row.');
            for k=1:9, out{k}(im,i,j)=t.(fields{k})(mask); end
        end
    end
end
[meanMean,meanSD,varMean,varSD,skeMean,skeSD,kurMean,kurSD,sampleCount]=deal(out{:});
end
