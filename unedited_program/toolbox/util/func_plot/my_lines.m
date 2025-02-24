function cmap = my_lines(nline, pat)

if ~exist('pat','var') || isempty(pat)
    pat = 1;
end

% Register new pattern here.
switch pat
    case 1
        pattern = [
            0    0.4470    0.7410;
            0.8500    0.3250    0.0980;
            0.9290    0.6940    0.1250;
            0.4940    0.1840    0.5560;
            0.4660    0.6740    0.1880;
            0.3010    0.7450    0.9330;
            0.6350    0.0780    0.1840;
            ];
        npat = size(pattern,1);
end

if ~exist('nline','var') || isempty(nline)
    nline = npat;
end

nrep = ceil(nline/npat);
cmap = repmat(pattern, [nrep,1]);
cmap = cmap(1:nline,:);

