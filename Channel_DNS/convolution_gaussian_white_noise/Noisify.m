% Woii 231208
% 221212 v6 including POD spectrum
SourceFolder = '../../Channel/volumn_edge/Fields_testing/';
OutputFolder = 'Field/';
AFrame = 1:501;
SourceFolder2= '../../Channel/volumn_edge/Fields/';
AFrame2= 1:1899;
OutputFile = 'V_NOISE.mat';

ROI = [1 88 1 88 1 12]; % [begin end begin end begin end] of the 1/2/3 dim
% parameters of adding noise and outliers
NoiseLevel = 0.25*6;
FreqOutlier = 0;
SizeOutlier = 240;
% outlier detecting parameters
U_range = [0 1.4];   V_range = [-0.25 0.25];   W_range = [-0.25 0.25];
Size_Median = 9;   Eps_Median = 0.1;   Th_Median = 3;
tic

% processing
n_nodes = (ROI(2)-ROI(1)+1)*(ROI(4)-ROI(3)+1)*(ROI(6)-ROI(5)+1);
U_noise = zeros(n_nodes,length(AFrame));
V_noise = U_noise;   W_noise = U_noise;   Enabled = U_noise;
U_array = U_noise;   V_array = U_noise;   W_array = U_noise;
U_ext   = zeros(n_nodes,length(AFrame2));
V_ext   = U_ext;     W_ext   = U_ext;
disp('processing flow field data...');
iCount = 0;
for iFrame = AFrame
    sub_name = ['Field_', num2str(iFrame,'%06u'), '.mat'];
    % loading field
    load([SourceFolder, sub_name], 'u', 'v', 'w');
    if iFrame == AFrame(1),    FieldStd = std(v, 0, 'all');    end
    % croping
    u = u(ROI(1):ROI(2),ROI(3):ROI(4),ROI(5):ROI(6));
    v = v(ROI(1):ROI(2),ROI(3):ROI(4),ROI(5):ROI(6));
    w = w(ROI(1):ROI(2),ROI(3):ROI(4),ROI(5):ROI(6));

    % adding noise
    u = u + FieldStd*NoiseLevel*imgaussfilt3(randn(size(u)),1);
    v = v + FieldStd*NoiseLevel*imgaussfilt3(randn(size(v)),1);
    w = w + FieldStd*NoiseLevel*imgaussfilt3(randn(size(w)),1);
    % adding outliers
    % bOutlier = binornd(1,FreqOutlier, size(u));
    % tOutlier = find(bOutlier == 1);
    % for iSub = 1:sum(bOutlier, 'all')
    %     mapOutlier = GenOutlierMap(SizeOutlier,u,tOutlier(iSub));
    %     shift = @()(poissrnd(4)+randn)*((randn > 0)*2-1);
    %     u(mapOutlier) = u(mapOutlier) + shift();
    %     v(mapOutlier) = v(mapOutlier) + shift();
    %     w(mapOutlier) = w(mapOutlier) + shift();
    % end

    % detecting outliers by velocity range
    e1= u > U_range(1) & u < U_range(2) & v > V_range(1) &...
        v < V_range(2) & w > W_range(1) & w < W_range(2);
    % detecting outliers by median filter
    e2= BMedian(u, Size_Median, Eps_Median, Th_Median) &...
        BMedian(v, Size_Median, Eps_Median, Th_Median) &...
        BMedian(w, Size_Median, Eps_Median, Th_Median);
    e = e1 & e2;
    % 
    u_mf = FMedian(u, e, Size_Median);
    v_mf = FMedian(v, e, Size_Median);
    w_mf = FMedian(w, e, Size_Median);

    iCount = iCount + 1;
    U_noise(:,iCount) = u(:);
    V_noise(:,iCount) = v(:);
    W_noise(:,iCount) = w(:);
    Enabled(:,iCount) = e(:);
    U_array(:,iCount) = u_mf(:);
    V_array(:,iCount) = v_mf(:);
    W_array(:,iCount) = w_mf(:);
end
toc

iCount = 0;
disp('processing more flow field data for POD...');
for iFrame = AFrame2
    sub_name = ['Field_', num2str(iFrame,'%06u'), '.mat'];
    % loading field
    load([SourceFolder2, sub_name], 'u', 'v', 'w');
    % croping
    u = u(ROI(1):ROI(2),ROI(3):ROI(4),ROI(5):ROI(6));
    v = v(ROI(1):ROI(2),ROI(3):ROI(4),ROI(5):ROI(6));
    w = w(ROI(1):ROI(2),ROI(3):ROI(4),ROI(5):ROI(6));

    % adding noise
    u = u + FieldStd*NoiseLevel*imgaussfilt3(randn(size(u)),1);
    v = v + FieldStd*NoiseLevel*imgaussfilt3(randn(size(v)),1);
    w = w + FieldStd*NoiseLevel*imgaussfilt3(randn(size(w)),1);
    % adding outliers
    % bOutlier = binornd(1,FreqOutlier, size(u));
    % tOutlier = find(bOutlier == 1);
    % for iSub = 1:sum(bOutlier, 'all')
    %     mapOutlier = GenOutlierMap(SizeOutlier,u,tOutlier(iSub));
    %     shift = @()(poissrnd(4)+randn)*((randn > 0)*2-1);
    %     u(mapOutlier) = u(mapOutlier) + shift();
    %     v(mapOutlier) = v(mapOutlier) + shift();
    %     w(mapOutlier) = w(mapOutlier) + shift();
    % end

    % detecting outliers by velocity range
    e1= u > U_range(1) & u < U_range(2) & v > V_range(1) &...
        v < V_range(2) & w > W_range(1) & w < W_range(2);
    % detecting outliers by median filter
    e2= BMedian(u, Size_Median, Eps_Median, Th_Median) &...
        BMedian(v, Size_Median, Eps_Median, Th_Median) &...
        BMedian(w, Size_Median, Eps_Median, Th_Median);
    e = e1 & e2;
    % 
    u_mf = FMedian(u, e, Size_Median);
    v_mf = FMedian(v, e, Size_Median);
    w_mf = FMedian(w, e, Size_Median);

    iCount = iCount + 1;
    U_ext(:,iCount) = u_mf(:);
    V_ext(:,iCount) = v_mf(:);
    W_ext(:,iCount) = w_mf(:);
end
toc

% POD processing
disp('generating POD spectrum...');
Um = mean([U_array, U_ext], 2);
Vm = mean([V_array, V_ext], 2);
Wm = mean([W_array, W_ext], 2);
[PsiU, SigmaU, PhiU] = svd([[U_array, U_ext]-Um; [V_array, V_ext]-Vm;...
    [W_array, W_ext]-Wm]', 'econ');

% saving
disp('saving data...');
save(OutputFile, 'U_noise','V_noise','W_noise','Enabled',...
    'U_array','V_array','W_array','U_ext','V_ext','W_ext',...
    'Um','Vm','Wm','SigmaU','PhiU', '-v7.3');
figure; plot(PsiU(1:100, 1:2));
s_sigma = sum(diag(SigmaU));
figure; plot(diag(SigmaU)/s_sigma);
yyaxis right; plot(cumsum(diag(SigmaU))/s_sigma);
set(gca, 'xscale', 'log')
toc

%% functions
function mapOutlier = GenOutlierMap(SizeOutlier,u,position)
% gen map of outliers
CubicTable = (1:2:9).^3;
mapOutlier = false(size(u));
sOutlier = poissrnd(SizeOutlier);
lOutlier = sum(sOutlier > CubicTable)*2+1;
kOutlier = zeros(lOutlier, lOutlier, lOutlier);
kOutlier(2:end-1,2:end-1,2:end-1) = 1;
kOutlier2 = kOutlier;
kOutlier2(:,1,1) = 1;          kOutlier2(:,end,end) = 1;
kOutlier2(1,:,1) = 1;          kOutlier2(end,:,end) = 1;
kOutlier2(1,1,:) = 1;          kOutlier2(end,end,:) = 1;
tmp0 = find(kOutlier2 == 0);
tmpr = sOutlier-(lOutlier-2)^3;
tmpE = randperm(length(tmp0),min(tmpr,length(tmp0)));
kOutlier(tmp0(tmpE)) = 1;
if tmpr > length(tmp0)
    tmpr = tmpr - length(tmp0);
    tmp0 = find(kOutlier == 0);
    tmpE = randperm(length(tmp0),tmpr);
    kOutlier(tmp0(tmpE)) = 1;
end
suby = mod(position-1,size(u,1))+1;
subx = ceil((mod(position-1,size(u,1)*size(u,2))+1)/size(u,2));
subz = ceil(position/size(u,1)/size(u,2));
hl = sum(sOutlier > CubicTable);
mapOutlier2 = false(size(u)+hl*2);
mapOutlier2(1+hl:end-hl,1+hl:end-hl,1+hl:end-hl) = mapOutlier;
mapOutlier2(suby:suby+2*hl,subx:subx+2*hl,subz:subz+2*hl) = kOutlier |...
    mapOutlier2(suby:suby+2*hl,subx:subx+2*hl,subz:subz+2*hl);
mapOutlier = mapOutlier2(1+hl:end-hl,1+hl:end-hl,1+hl:end-hl);
end

function enabled = BMedian(u, FiltSize, E, Th)
u_med = medfilt3(u, [FiltSize, FiltSize, FiltSize]);
r     = abs(u - u_med);
r_med = medfilt3(r, [FiltSize, FiltSize, FiltSize]);
enabled = r./(r_med+E) < Th;
end

function u = FMedian(u, enabled, FiltSize)
u_med = medfilt3(u, [FiltSize, FiltSize, FiltSize]);
u(enabled == 0) = u_med(enabled == 0);
end