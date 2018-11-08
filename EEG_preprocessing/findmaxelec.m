function elec = findmaxelec(chanmat,numelec,plotrank)

% returns the ranking of the electrodes with maximum magnitude across subjects
% chanmat (subjects x electrodes x conditions) = matrix with all the channels sorted as a function of magnitude
% numeelec = scalar; number of electrodes to insert in the ranking of the strongest electrodes
% plotrank = 0 or not specified, it does not plot the ranking

if nargin < 3
    plotrank = 0;
elseif nargin < 2
    numelec = 4;
end

numsubjs = size(chanmat,1);
chanindex = 1:prod(numel(chanmat));
channels = 1:64;
pos = zeros(numsubjs,size(channels,2));
for ielec = 1:size(channels,2)
    where = chanindex(chanmat == channels(ielec));
    pos(:,ielec) = ceil(where/numsubjs);
end
elecweight = sum(1./pos)/numsubjs; % normalized across subjects
[weights elecrank] = sort(elecweight,'descend');
elec = elecrank(1:numelec);

if plotrank
    figure
    h = bar(weights(1:numelec),'r');
    set(gca,'fontsize',15)
    set(h,'LineWidth',2.5)
    title('Electrodes rank','FontSize',20)
    ylabel('electrode weight','FontSize',20)
    set(gca, 'XTickLabel',{elecrank(1:numelec)},'FontSize',15)
end


