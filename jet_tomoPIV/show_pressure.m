FileName = 'JetAR1_03_conAMIC_7';
load([FileName, '.mat']);

iFrame = 6;
U = reshape(U_array(:,iFrame), size(X));
V = reshape(V_array(:,iFrame), size(X));
W = reshape(W_array(:,iFrame), size(X));
P = reshape(P_array(:,iFrame), size(X));
E = reshape(Enabled(:,iFrame), size(X));


iSlice = 1;        figure;
subplot(2,2,1); pcolor(X(:,:,iSlice), Y(:,:,iSlice), U(:,:,iSlice));
shading flat; colorbar('northoutside'); colormap(jet); axis equal;
subplot(2,2,2); pcolor(X(:,:,iSlice), Y(:,:,iSlice), V(:,:,iSlice));
shading flat; colorbar('northoutside'); colormap(jet); axis equal;
subplot(2,2,3); pcolor(X(:,:,iSlice), Y(:,:,iSlice), W(:,:,iSlice));
shading flat; colorbar('northoutside'); colormap(jet); axis equal;
subplot(2,2,4); pcolor(X(:,:,iSlice), Y(:,:,iSlice), P(:,:,iSlice));
shading flat; colorbar('northoutside'); colormap(jet); axis equal;

% set(gca, 'fontsize', 36);
set(gcf,'outerposition',get(0,'screensize'));



%%
iSlice = 18;
% figure;
close all;

% iFrame = 3;
for iFrame = 100:10:1100
    figure;
    U = reshape(U_array(:,iFrame), size(X));
    V = reshape(V_array(:,iFrame), size(X));
    W = reshape(W_array(:,iFrame), size(X));
    P = reshape(P_array(:,iFrame), size(X));
    E = reshape(Enabled(:,iFrame), size(X));

    subplot(2,2,1); pcolor(X(:,:,iSlice), Y(:,:,iSlice), U(:,:,iSlice));
    shading flat; colorbar('northoutside'); colormap(jet); axis equal;
    clim([-0.02 0.16]);
    subplot(2,2,2); pcolor(X(:,:,iSlice), Y(:,:,iSlice), V(:,:,iSlice));
    shading flat; colorbar('northoutside'); colormap(jet); axis equal;
    clim([-0.06 0.06]);
    subplot(2,2,3); pcolor(X(:,:,iSlice), Y(:,:,iSlice), W(:,:,iSlice));
    shading flat; colorbar('northoutside'); colormap(jet); axis equal;
    clim([-0.06 0.06]);
    subplot(2,2,4); pcolor(X(:,:,iSlice), Y(:,:,iSlice), P(:,:,iSlice));
    shading flat; colorbar('northoutside'); colormap(jet); axis equal;
    clim([-5 5]);

    set(gcf,'outerposition',get(0,'screensize'));
    drawnow
    pause(1);
    img{iFrame} = frame2im(getframe(gcf));
    close(gcf);
    % clf;
end

%%

filename = [FileName, '.gif'];
ffra = 0;
for idx = 1:size(img,2)
    if ~isempty(img{idx})
        [A,map] = rgb2ind(img{idx},256);
        if ffra == 0
            imwrite(A,map,filename,"gif","LoopCount",0,"DelayTime",1);
            ffra = 1;
        else
            imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",1);
        end
    end
end