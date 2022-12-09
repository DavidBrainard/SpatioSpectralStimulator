% t_WatsonCSF

sfs = linspace(0,20,100);
tf = 1;
fakeContrastdB = 1;


% Model params
minThreshDb = -1000;
c0s = linspace(-10,10,10);
cfs = linspace(-10,10,10);
c2s = linspace(-10,10,10);

for ii = 1:length(sfs)
    stimParams(ii,1) = sfs(ii);
    stimParams(ii,2) = tf;
    stimParams(ii,3) = fakeContrastdB;
end


figure; clf; hold on;
for kk = 1:length(c0s)
    for ll = 1:length(cfs)
        for mm = 1:length(c2)
            psiParams = [minThreshDb c0s(kk) cfs(ll) c2s(mm)];

            [~,predictedContrastThreshold] = qpPFSTCSF(stimParams,psiParams);
            sensitivity = -predictedContrastThreshold;

            plot(sfs,sensitivity,'r');
        end
    end
end