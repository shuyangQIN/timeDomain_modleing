phantom_image = phantom('Modified Shepp-Logan',128);
subplot(1,2,1)
imagesc(phantom_image)
colormap jet
colorbar
dim = size(phantom_image);
win = getWin(dim, 'Hanning', 'Rotation', true , 'Symmetric', true);
image_sm = real(ifftn(fftn(phantom_image) .* ifftshift(win)));
subplot(1,2,2)
imagesc(image_sm)
colormap jet
colorbar
