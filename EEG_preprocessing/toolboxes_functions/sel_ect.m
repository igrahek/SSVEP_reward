function elec = sel_ect(cluster,channs,num_channs)

% selects a predefined number of electrodes with the maximum magnitude
% chosen from a given cluster of electrodes, if the third input is not specified
% 4 electrodes per subject are chosen

if nargin < 3
    num_channs = 4;
end

coll_elec = numel(cluster);
elec = zeros(size(channs,1),numel(cluster));
for isub = 1:size(channs,1)
    for ichann = 1:numel(cluster)
        coll_elec(:,ichann) = find(channs(isub,:,:) == cluster(ichann));
    end
    elec(isub,:) = channs(isub,sort(coll_elec),:);
end
elec = elec(:,1:num_channs);
