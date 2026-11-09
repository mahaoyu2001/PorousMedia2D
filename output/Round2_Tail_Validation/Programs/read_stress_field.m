function data=read_stress_field(fileName)
% Four columns: S11 S22 S12 EVOL, optionally preceded by text headers.
fid=fopen(fileName,'rt');assert(fid>=0,'Cannot open field file: %s',fileName);
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
headerLines=0;
while true
    line=fgetl(fid);
    assert(ischar(line),'No numeric field data in %s',fileName);
    text=strtrim(line);
    if ~isempty(text) && ~isempty(regexp(text,'^[+\-0-9.]','once')),break,end
    headerLines=headerLines+1;
end
frewind(fid);
parsed=textscan(fid,'%f%f%f%f','HeaderLines',headerLines, ...
    'CollectOutput',true,'ReturnOnError',false);
data=parsed{1};
assert(~isempty(data) && size(data,2)==4 && all(isfinite(data(:))) ...
    && all(data(:,4)>0),'Expected S11 S22 S12 and positive EVOL in %s',fileName);
end
