function [ENV_bf, PeriParab,  EndoParab, Angles, Times]=DAS_boneimage_Salome( ...
    bone_image_output_data, fnumber, ...
    XMinSegmentationEndo,XMaxSegmentationEndo, plotMap)
%
apodization_start_end_YES = 1;

SIG = bone_image_output_data.rf_data;

bone_image_output_data.C_RADIAL;
%--------------------------------------------------%
%           P4-1 probe
%--------------------------------------------------%
% lens parameters found with semblance on single reflection on interface
% silicone/water (15 oct 2019)
FS                  = bone_image_output_data.PROBE_PARAM.Fs;%10e6; % [Hz]
[~,~,NSOURCES]           = size(SIG);
FREQ_Transducteur   = bone_image_output_data.PROBE_PARAM.FREQ_Transducteur;% 2.5e6; % [Hz]


% for synthetic aperture imaging, zero transmit delay added (only used for virtual sources)
add_to_delay_firing = zeros(1,NSOURCES);

% horizontal positions of the elements
XR =  bone_image_output_data.PROBE_PARAM.XS;%XR;% (0:NELEMENTS-1) *PITCH; XR = XR - mean(XR); % [m]
% define dimension and pixel size of reconstructed image

% create image axes/pixels coordinates
X       = bone_image_output_data.X;%xmin:resolution:xmax; % image width [m]
Z       = bone_image_output_data.Z;%zmin:resolution:zmax; % image depth [m]
%
% Sources to be used during reconstruction
Tx_to_be_used       = 1:NSOURCES;

% Remark: no acceptance angle HalfOpeningAng used in transmit

% Receive angle(s)
ReceiveAngle    = deg2rad(0); % receive angles, can be a vector for vector flow imaging [rad]

% vertical positions of the elements
ZR = zeros(1,length(XR));% [m]
% coordinates of sources
XS = XR;
ZS = ZR;
% backward shift of raw echo signals by time to peak of the round-trip waveform
SIG = SIG(bone_image_output_data.PROBE_PARAM.offset+1:end,:,:);

if apodization_start_end_YES
    muting_duration                     = 1.2e-6;
    SIG(1:ceil(muting_duration*FS),:,:) = 0;
end

% I/Q separation.
I_SIG = ones(size(SIG,1),size(SIG,2),size(SIG,3));
Q_SIG = ones(size(SIG,1),size(SIG,2),size(SIG,3));
for Tx = 1:size(SIG,3)
    tmp =  hilbert(squeeze(SIG(:,:,Tx)));
    I_SIG(:,:,Tx) = real(tmp);
    Q_SIG(:,:,Tx) = imag(tmp);
end



%
FnumberMin=fnumber;1.9;


Setup = bone_image_output_data.Setup;
% Parabolas (if needed the user can force the parabolas of periosteum and/or endosteum)
if isfield(bone_image_output_data, 'PeriParabIn')
    PeriParabIn = bone_image_output_data.PeriParabIn;
else
    PeriParabIn = [ 0 0 0 ];%0 0 ]; % if automatic segmentation desired, must be a vector with 5 zeros
end
EndoParabIn = [ 0 0 0 ];%0 0 ]; % if automatic segmentation desired, must be a vector with 5 zeros
Setup(8) = FnumberMin;

ReceiveAngAtPix_rad = ReceiveAngle;

SetupNEW = [Setup(2) Setup(4:7) Setup(7:15)...D_PROBE_PERIOS 
    ...MIN_CORTICAL MAX_CORTICAL 
    ...TxHalfOpeningAngInLensRad RxHalfOpeningAngInLensRad...
    ...FNumberMin Max_err_allowed_ReceiveAngle_rad FS 
    ...C_LENS C_TISSUE C_MARROW...
    ...C_AXIAL C_RADIAL 
    Setup(1) Setup(16:end)...ANISO_SHAPE_COEF NCORE ReConTo NeedTravelTime ReceiveSubApertureApodis...
    XMinSegmentationEndo XMaxSegmentationEndo];
% rad2deg(SetupNEW(5))
% rad2deg(SetupNEW(6))
[~, Time_T_Tissue, Time_R_Tissue, Angle_T_Tissue,...
    Angle_R_Tissue, Beam_I_Tissue, Beam_Q_Tissue, APix_T_Tissue, APix_R_Tissue,...
    ~, Time_T_Bone, Time_R_Bone, Angle_T_Bone,...
    Angle_R_Bone, Beam_I_Bone, Beam_Q_Bone,...
    APix_T_Bone, APix_R_Bone, PeriParab, EndoParab] =...
    PostProcessingFNumber_WeakAniso_SUPERFAST_RxApodAll_v4(...
        SetupNEW, XR, ZR, XS, ZS, X, Z,...
        I_SIG, Q_SIG, ...
        Tx_to_be_used, add_to_delay_firing,...
        ReceiveAngAtPix_rad, PeriParabIn, EndoParabIn);

fit_curve_Periosteum = polyval(PeriParab(1:3),X);
fit_curve_Endosteum = polyval(EndoParab(1:3),X);

Angles.Angle_T_Tissue = Angle_T_Tissue;
Angles.Angle_R_Tissue = Angle_R_Tissue;
Angles.APix_T_Tissue = APix_T_Tissue;
Angles.APix_R_Tissue = APix_R_Tissue;
Angles.Angle_T_Bone = Angle_T_Bone;
Angles.Angle_R_Bone = Angle_R_Bone;
Angles.APix_T_Bone = APix_T_Bone;
Angles.APix_R_Bone = APix_R_Bone;

Times.Time_T_Tissue = Time_T_Tissue;
Times.Time_R_Tissue = Time_R_Tissue;
Times.Time_T_Bone = Time_T_Bone;
Times.Time_R_Bone = Time_R_Bone;
% sum over transmissions for images with fixed receive angle and F-number
% if multiple receive angles, then receive angle must be chosen with
% dimension 3 in Beam_I and Beam_Q
C_TISSUE = Setup(12); C_RADIAL = Setup(14);
size(APix_R_Bone);
size(APix_T_Bone);
ImageTissue_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Tissue,4).^2+sum(Beam_Q_Tissue,4).^2));
ImageBone_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Bone,4).^2+sum(Beam_Q_Bone,4).^2));

%% Show
display_dynamic_range=40;
[stitching_mask_Tissue, stitching_mask_Bone]=get_stitching_masks...
    (X,Z, PeriParab(1:3), FREQ_Transducteur,C_TISSUE,C_RADIAL);
FullImage_FixedReceiveAng = stitching_mask_Tissue.*ImageTissue_FixedReceiveAng + ...
    stitching_mask_Bone.*ImageBone_FixedReceiveAng;

ENV_bf = FullImage_FixedReceiveAng;
ENV_bf = ENV_bf/max(ENV_bf(:));
show_image = 0;
if ~show_image
    return
else %show_image
    figure
    subplot 211 
    imagesc(X*1e3,Z*1e3,20*log10(ImageTissue_FixedReceiveAng))%/max(FullImage_FixedReceiveAng(:))))%))
    axis image
    title('Image tissus mous')
    colormap gray
    colorbar
    hold on
        plot(X*1e3,fit_curve_Periosteum*1e3,'k:', 'LineWidth',1.0)
    subplot 212
    imagesc(X*1e3,Z*1e3,20*log10(ImageBone_FixedReceiveAng))%/max(FullImage_FixedReceiveAng(:))))%))
    axis image
    title('Image Os')
    colormap gray
    colorbar
    hold on
        plot(X*1e3,fit_curve_Endosteum*1e3,'k:', 'LineWidth',1.5)

end    


end



