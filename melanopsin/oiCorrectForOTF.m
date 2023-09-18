function oi = oiCorrectForOTF(oi, otfM)


% Get the current data set. It has the right size. We over-write it below.
p = oiGet(oi, 'photons');
wave = oiGet(oi,'wave');

for ii = 1:length(wave)
    img = p(:, :, ii);
    % figure(1);
    % imagesc(img);
    % colormap(gray);

    otf = otfM(:, :, ii);
    % vcNewGraphWin;
    % mesh(fftshift(otf));
    % otf(1, 1)

    % Put the image center in (1, 1) and take the transform.
    imgFFT = fft2(fftshift(img));
    % figure(1);
    % imagesc(abs(imgFFT));
    % figure(2);
    % imagesc(abs(otf));
    % colormap(gray)

    % Divide the transformed image by the otf
    filteredIMG = abs(ifftshift(ifft2(imgFFT ./ otf)));

    % Put it back
    p(:, :, ii) = filteredIMG;
end

% Store correced image in oi for return
oi = oiSet(oi, 'photons', p);

