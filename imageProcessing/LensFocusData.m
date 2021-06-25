%% Width 20
distance = [0.5,0.75,1.0,1+1/16,1+2/16,1+3/16,1+7/32,1+8/32,1+9/32,1+9.5/32,1+10/32,1+11/32,1+12/32,1+7/16,1+8/16,1.75,2.00];
contrast = [0.19799,0.66611,0.92,0.93624,0.9405,0.94702,0.95017,0.95033,0.94702,0.9438,0.94711,0.94719,0.9403,0.92799,0.90728,0.52174,0.087894];
fftarea = [1.4796,2.695,3.7637,5.2641,5.013,3.3249,3.7601,4.2912,5.3609,5.6031,5.4896,5.5493,5.121,3.2796,2.8132,1.719,1.2619];
fftarea_normalized = fftarea./max(fftarea);

figure(1); hold on;
plot(distance,contrast,'k--.','MarkerSize',15);
plot(distance(find(contrast==max(contrast))),max(contrast),'ko','MarkerSize',10); % Best point
plot(distance,fftarea_normalized,'b--.','MarkerSize',15);
plot(distance(find(fftarea_normalized==max(fftarea_normalized))),max(fftarea_normalized),'bo','MarkerSize',10); % Best point

title('Contrast / FFT area over distance')
xlabel('Distance from the nominal (inch)')
ylabel('Contrast or Normalized FFT area')
legend('Contrast','Best contrast','FFT area','Best FFT area')

% figure(3); hold on;
% plot(contrast,fftarea_normalized,'k--.','MarkerSize',15);
% title('Contrast vs. FFT area')
% xlabel('Contrast')
% ylabel('Normalized FFT area')
% legend('All data')

%% Raw Data (as of 0616) - Width 50
% Nominal=2.6582
% 0.5inch=3.0357
% 0.75inch=4.272
% 1-inch =3.6342
% 1.25-inch = 6.05430
% 1.5-inch=3.6289 (0.97701 contrast)
% 2-inch=4.8405 (0.95025)
% 
% %% Raw Data (as of 0616) - Width 20
% 
% 0.5-inch=1.4796 (0.19799)
% 0.75-inch = 2.695 (0.66611)
% 
% 1.0-inch = 3.7637 (0.92)
% % 1.0-inch_2nd = 3.6743 (0.92053)
% 1 1/16 = 5.2641 (0.93624) 
% 1 2/16 = 5.013 (0.9405)
% 1 3/16 = 3.3249 (0.94702)
% % 1 3/16 = 3.4485 (0.94745)
% 1 7/32 = 3.7601 (0.95017)
% 1 8/32 = 4.2912 (0.95033)
% % 1 4/16 = 4.77 (0.94719)
% % 1.25-inch = 5.0002 (0.94745)
% 
% ***Final position*** (as of June 16)
% 1 9/32 = 5.3609 (0.94702)
% % 1 9/32 = 5.1131 (0.9509)
% % 1 9/32 = 5.4034 (0.94702)
% % 1 9/32 = 5.4013 (0.94719)
% 
% 1 9.5/32 = 5.6031 (0.9438)/ 5.5986 (0.94702)
% 1 10/32 = 5.4896 (0.94711) 
% % 1 5/16 = 5.6455 (0.94684)
% 1 11/32 = 5.5493 (0.94719)
% 1 12/32 = 5.121 (0.9403)
% % 1 6/16 = 5.0404 (0.94)
% 1 7/16 = 3.2796 (0.92799)
% 1 8/16 = 2.8132 (0.90728)
% % 1.5-incn = 3.1738 (0.90713)
% 1.75-incn=1.719 (0.52174)
% 2-inch=1.2619 (0.087894)
% 
% 
