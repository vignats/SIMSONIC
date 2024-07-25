addpath(genpath('~/Documents'));

% Pathway to the bone image to process
bone = {'271G' '267G' '251G' '245D' '227G'};
bone_nb = 5;
dirname = ['~/Documents/BoneRugosity/RMS/', bone{bone_nb}, '/'];
num_slice = 1226;
file = ['SAMPLE_', bone{bone_nb}, '_SLICE_', int2str(num_slice), '.bmp']; 
filename = [dirname , file];

% Load image
bone_bmp = imread(filename); 

%% INITIAL COMPUTATION OF RMS HEIGTH
% THRESHOLD TO CONVERT TO BINARY IMAGE
threshold = graythresh(bone_bmp); % Find an automatic threshold
binary_image_init = imbinarize(bone_bmp, threshold);

% Plot to verifiy
print = false;
if print 
    figure,
    subplot(2,1,1);
    imshow(bone_bmp); 
    title('Original bone image');
    subplot(2,1,2);
    imshow(binary_image_init); 
    title(sprintf('Binary bone image with a %0.1f threshold ', threshold));
end 
% BOUNDARIES COMPUTATION
[boundary_endost_init, xy] = ExtractBoundary(binary_image_init);

% POLYNOMIAL FITTING IF THE ENDOST
surface = FitParabola(binary_image_init, boundary_endost_init, true);

% COMPUTATION OF RMS 
rms_0 = rms(surface(:, 2) - boundary_endost_init(:, 2))*9e-3;       % Rms (mm) - Pixel size = 9e-3 mm
%% COMPUTE RMS HEIGTH FOR ALL IMAGES
% Compute the rms for all the other images, using the parabola and endost
% limitation defined for the first image 
tic
disp('------- Computation of the mean rms -------')
% Initialise the rms_mean with the rms of the first image
rms_all = [num_slice, rms_0];

% To start the loop with the second image, change filename to point it to
% the second image
num_slice = num_slice + 100;
file = ['SAMPLE_', bone{bone_nb}, '_SLICE_', int2str(num_slice), '.bmp']; 
filename = [dirname , file];

while exist(filename, 'file')
    % Load image
    bone_bmp = imread(filename); 

    % Convert to binary image
    binary_image = imbinarize(bone_bmp, threshold);

    % Extract the boundary of the endost and periost
    boundaries = bwboundaries(binary_image,8,'noholes');
    [~, max_index] = max(cellfun('size', boundaries, 1));
    boundary = flip(boundaries{max_index}, 2);

    % Extract the endost boundary
    % The xy limits obtained with the first image need to be indicated for
    % the other images
    [boundary_endost, ~] = ExtractBoundary(binary_image, xy);
    surface = FitParabola(binary_image_init, boundary_endost, false);
    
    % Compute the rms
    [~, idx_surface, idx_boundary_endost] = intersect(surface(:, 1), boundary_endost(:, 1));
    rms_all(end +1,:) = [num_slice, rms(surface(idx_surface, 2) - boundary_endost(idx_boundary_endost, 2))*9e-3];       % Rms (mm) - Pixel size = 9e-3 mm

    num_slice = num_slice + 1;
    file = ['SAMPLE_', bone{bone_nb}, '_SLICE_', int2str(num_slice), '.bmp']; 
    filename = [dirname , file];
end
toc 

%% SAVE RMS VALUE IN .mat FILE
save(fullfile(dirname, 'rms.mat'), "rms_all");
figure 
plot(rms_all(:,1), rms_all(:,2));