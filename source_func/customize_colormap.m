function cmap = customize_colormap(index)
switch index
    case 1
        % colormap -- yellow to blue
        d1 = [255 255-128 0]/(255*63);
        d2 = [0 0 255]/(255*63);

        c1 = zeros(64,3);
        c1(1,:) = [0 128 255]/255;
        c2 = zeros(64,3);
        c2(1,:) = [255 255 0]/255;
        for j = 2:64
            c1(j,:) = c1(1,:)+(j-1)*d1;
            c2(j,:) = c2(1,:)+(j-1)*d2;
        end

        c2 = flipud(c2);
        cmap = [c1;c2(2:64,:)];
    
    case 2
        % colormap -- same as kwave
        c = hot;
        c = flipud(c);
        cmap = [bone;c];
        
    otherwise
        error('Unsupported colormap')
end