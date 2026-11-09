function [meanMean,meanSD,varMean,varSD,skeMean,skeSD,kurMean,kurSD] = ...
    round2_cached_arrays(kind,values,componentNames)
t=round2_read_statistics(kind);
key='Porosity'; if strcmpi(kind,'load'), key='LoadRatio'; end
fields={'Mean','Mean_SD','Variance','Variance_SD','Skewness','Skewness_SD','Kurtosis','Kurtosis_SD'};
out=cell(1,8);
for k=1:8, out{k}=nan(numel(componentNames),numel(values)); end
for i=1:numel(componentNames)
    for j=1:numel(values)
        mask=abs(t.(key)-values(j))<1e-12 & string(t.StressComponent)==string(componentNames(i));
        assert(nnz(mask)==1,'Missing or duplicate cached %s row for %s, %g.',kind,componentNames(i),values(j));
        for k=1:8, out{k}(i,j)=t.(fields{k})(mask); end
    end
end
[meanMean,meanSD,varMean,varSD,skeMean,skeSD,kurMean,kurSD]=deal(out{:});
end
