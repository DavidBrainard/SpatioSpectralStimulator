% ContrastMeasureDataForLens1
%
% Contrast measurement data to decide the Lens1 position on the cage rod
% system for SACC project.
 
% History:
%    12/06/21 smo   Measured and noted.

%% Measured data.
% Raw data
distanceInches = [0.0000, 0.5000, 1.0000, 1+1/16, 1+2/16, 1+3/16, 1+4/16, 1+5/16, 1+6/16, 1+7/16, 1+8/16, 1+9/16, 1+10/16, 1+11/16, 1+12/16, 1+13/16, 1+14/16, 1+15/16, 2.0];  
contrast =       [0.1635, 0.2952, 0.9262, 0.9535, 0.9537, 0.9547, 0.9539, 0.9545, 0.9450, 0.9444, 0.9070, 0.8461, 0.7788, 0.7037, 0.5392, 0.3674, 0.2796, 0.2093, 0.1268];
idxMaxContrast = find(contrast == max(contrast));

% Plot it.
figure; clf; hold on;
plot(distanceInches, contrast, 'k.-', 'markersize', 12);
plot(distanceInches(idxMaxContrast), contrast(idxMaxContrast), 'ro', 'markersize', 14);
xlabel('Distance between Lens 1 and the nominal position (inch)');
ylabel('Contrast');
legend('All data','Max Contrast','location','northwest');

% Add text to the max point.
text(distanceInches(idxMaxContrast)+0.05, contrast(idxMaxContrast)+0.02, ...
    append('Contrast: ',num2str(round(contrast(idxMaxContrast),2))),'Color','red','FontSize',14);
text(distanceInches(idxMaxContrast)+0.05, contrast(idxMaxContrast)+0.06, ...
    append('Distance: ',num2str(round(distanceInches(idxMaxContrast),2))),'Color','black','FontSize',14);