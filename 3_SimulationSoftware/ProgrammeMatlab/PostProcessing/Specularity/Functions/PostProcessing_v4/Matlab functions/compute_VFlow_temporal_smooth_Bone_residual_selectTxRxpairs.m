function [vz,vx,LS_sd_residuals,LS_sd_error_Vz,LS_sd_error_Vx]=compute_VFlow_temporal_smooth_Bone_residual_selectTxRxpairs(I,Q,transmitAngles,receiveAngles,C_RADIAL,C_AXIAL,ANISO_SHAPE_COEF,f0,PRF,Average_kernel_size_x,Average_kernel_size_z,kernel_size_temporal_smoothing,selected_iTx,selected_iRx)

%---------------------------------------------------------------------------------------------%
% kernel_size_temporal_smoothing = 1 does nothing, 0 average all frames, >1 temporal smoothing
%---------------------------------------------------------------------------------------------%

if length(selected_iTx)~=length(selected_iTx)
    disp('selected_iTx and selected_iRx must have the same size')
    return
end

sz = size(I);
NZ = sz(1);
NX = sz(2);
NT = sz(5);
N_TxRx_pairs = length(selected_iTx);

% for calculating phase angle in cortical bone from the group angle output by
% reconstruction code:
ANISO_COEF = (C_AXIAL - C_RADIAL)/C_AXIAL;
epsilon = -ANISO_COEF; % Thomsen weak aniso parameter 
delta_weak_aniso = -ANISO_COEF*ANISO_SHAPE_COEF; % Thomsen weak aniso parameter
A = 1+2*delta_weak_aniso;
B = 4*(epsilon-delta_weak_aniso);


lag = 1;

% kernel_size_temporal_smoothing = 1 does nothing, 0 average all frames, >1 temporal smoothing
if (kernel_size_temporal_smoothing==0) % temporal averaging on all frames
    Umat = zeros(NZ,NX,N_TxRx_pairs);
    vz = zeros(NZ,NX);
    vx = zeros(NZ,NX);
    LS_sd_residuals = zeros(NZ,NX);
    LS_sd_error_Vz = zeros(NZ,NX);
    LS_sd_error_Vx = zeros(NZ,NX);
else
    Umat = zeros(NZ,NX,N_TxRx_pairs,NT-1);
    vz = zeros(NZ,NX,NT-1);
    vx = zeros(NZ,NX,NT-1);
    LS_sd_residuals = zeros(NZ,NX,NT-1);
    LS_sd_error_Vz = zeros(NZ,NX,NT-1);
    LS_sd_error_Vx = zeros(NZ,NX,NT-1);
end


for iTxRx=1:N_TxRx_pairs
    fprintf(['Computing phase shifts for TxRx pair ' int2str(iTxRx) '/' int2str(N_TxRx_pairs) ' . . .  \n']);
        tx = selected_iTx(iTxRx);
        rx = selected_iRx(iTxRx);
        
        xmn = squeeze(I(:,:,rx,tx,:) + 1j*Q(:,:,rx,tx,:)); % IQ signal
        R=xmn(:,:,1:end-lag).*conj(xmn(:,:,lag+1:end))/lag;
        
        if (Average_kernel_size_z~=1 && Average_kernel_size_x~=1)
            %-----------------------------------------%
            %---------- Spatial smoothing ------------%
            %-----------------------------------------%
            sAvg = ones(Average_kernel_size_z,Average_kernel_size_x)./Average_kernel_size_x/Average_kernel_size_z;
            for iF=1:size(R,3)    
                R(:,:,iF) = filter2(sAvg,R(:,:,iF)); 
            end
        end
        
        if (kernel_size_temporal_smoothing~=1)&&(kernel_size_temporal_smoothing~=0)
            %-----------------------------------------%
            %--------- Temporal smoothing ------------%
            %-----------------------------------------%
            ht=ones(1,kernel_size_temporal_smoothing)/kernel_size_temporal_smoothing;
            Umat(:,:,iTxRx,:) = angle(filter(ht,1,R,[],3)); % phase shift is a 3D array
        elseif (kernel_size_temporal_smoothing==0)
            % temporal averaging on all frames
            Umat(:,:,iTxRx) = angle(mean(R,3)); % phase shift is a 2D array
        elseif (kernel_size_temporal_smoothing==1)
            Umat(:,:,iTxRx,:) = angle(R); % phase shift is a 3D array
        end

end



scaling_factor = (2*pi*f0)/PRF;

if (kernel_size_temporal_smoothing==0) % temporal averaging on all frames, vz and vx are 2D arrays
   
    % Inversion of Velocity components
    for z = 1:NZ
        for x = 1:NX
            AMat = nan(N_TxRx_pairs,2); % matrix A for this pixel
            for iTxRx=1:N_TxRx_pairs
                thetaM = transmitAngles(selected_iTx(iTxRx),z,x); % transmit group angle in cortical bone output by reconstruction code            
                X = tan(thetaM); % tan of transmit group angle
                tan_phase_ang = ((A^3*X^3)/27 + (A*X)/3 + (B*X)/2 + (((A^3*X^3)/27 + (A*X)/3 + (B*X)/2)^2 - ((A^2*X^2)/9 - 1/3)^3)^(1/2))^(1/3) + ((A^2*X^2)/9 - 1/3)/((A^3*X^3)/27 + (A*X)/3 + (B*X)/2 + (((A^3*X^3)/27 + (A*X)/3 + (B*X)/2)^2 - ((A^2*X^2)/9 - 1/3)^3)^(1/2))^(1/3) + (A*X)/3;
                S2 = sin(atan(tan_phase_ang)).^2;
                C2 = 1 - S2;
                V_bone_T = C_AXIAL - ( C_AXIAL - C_RADIAL ) * ( ANISO_SHAPE_COEF.*S2.*C2 + C2.^2 ); % approximation of weak transverse isotropy
                thetaM = atan(tan_phase_ang); % phase angle (transmit)
                
                phiN = receiveAngles(selected_iRx(iTxRx),z,x); % receive group angle in cortical bone output by reconstruction code
                X = tan(phiN); % tan of receive group angle
                tan_phase_ang = ((A^3*X^3)/27 + (A*X)/3 + (B*X)/2 + (((A^3*X^3)/27 + (A*X)/3 + (B*X)/2)^2 - ((A^2*X^2)/9 - 1/3)^3)^(1/2))^(1/3) + ((A^2*X^2)/9 - 1/3)/((A^3*X^3)/27 + (A*X)/3 + (B*X)/2 + (((A^3*X^3)/27 + (A*X)/3 + (B*X)/2)^2 - ((A^2*X^2)/9 - 1/3)^3)^(1/2))^(1/3) + (A*X)/3;
                S2 = sin(atan(tan_phase_ang)).^2;
                C2 = 1 - S2;
                V_bone_R = C_AXIAL - ( C_AXIAL - C_RADIAL ) * ( ANISO_SHAPE_COEF.*S2.*C2 + C2.^2 ); % approximation of weak transverse isotropy
                phiN = atan(tan_phase_ang); % phase angle (receive)
                        
                AMat(iTxRx,1) = (cos(thetaM)/V_bone_T + cos(phiN)/V_bone_R);
                AMat(iTxRx,2) = (sin(thetaM)/V_bone_T + sin(phiN)/V_bone_R);
            end
            
            % Pseudo inverse of over etsimated problem in Least Square (No
            % Damping nor Weighting so far)
            % pseudoInvA           = (AMat'*AMat)\(AMat');
            % VectorFlowMat_LS(z,x,:) = pseudoInvA*squeeze(UMat(z,x,:)); % should be [ vz and vx ]!
            % Supposed to be more accurate.
            idx_good    = find(~isnan(AMat(:,1)));
            AMat        = AMat(idx_good,:)*scaling_factor;
            
            pseudoInvA           = pinv(AMat); %(AMat'*AMat)\(AMat');
            y = squeeze(Umat(z,x,idx_good));
            beta = pseudoInvA*y;
            vz(z,x) = beta(1);
            vx(z,x) = beta(2);

            residuals =  y - AMat*beta;
            n = length(y); % number of residuals
            p = 2; % number of physical parameters
            var_res = ((residuals.')*residuals)/(n-p); % variance of the residuals
            X_t_X_inv = inv(AMat'*AMat);
            var_estimates = X_t_X_inv*var_res; % variance-covariance matrix of the estimates of the physical variables 
            % The square roots of the diagonal elements of this matrix are
            % the standard deviations of the estimates, that is the errors
            % to be quoted on the determinations of the physical variables.
            LS_sd_error_Vz(z,x) = sqrt(var_estimates(1,1));
            LS_sd_error_Vx(z,x) = sqrt(var_estimates(2,2));
            LS_sd_residuals(z,x) = sqrt(var_res);%sqrt(mean((AMat*squeeze(vector_flow_mat(z,x,:))-squeeze(Umat(z,x,idx_good,it))).^2));
            
        end
    end

    

else % vz and vx are 3D arrays (image depth,image width, slow time)

    
    % Inversion of Velocity components
    for z = 1:NZ
    disp(['Vector Doppler: pixel row ' int2str(z) '/' num2str(NZ)])
        for x = 1:NX
            AMat = nan(N_TxRx_pairs,2); % matrix A for this pixel
            for iTxRx=1:N_TxRx_pairs
                thetaM = transmitAngles(selected_iTx(iTxRx),z,x); % transmit group angle in cortical bone output by reconstruction code            
                X = tan(thetaM); % tan of transmit group angle
                tan_phase_ang = ((A^3*X^3)/27 + (A*X)/3 + (B*X)/2 + (((A^3*X^3)/27 + (A*X)/3 + (B*X)/2)^2 - ((A^2*X^2)/9 - 1/3)^3)^(1/2))^(1/3) + ((A^2*X^2)/9 - 1/3)/((A^3*X^3)/27 + (A*X)/3 + (B*X)/2 + (((A^3*X^3)/27 + (A*X)/3 + (B*X)/2)^2 - ((A^2*X^2)/9 - 1/3)^3)^(1/2))^(1/3) + (A*X)/3;
                S2 = sin(atan(tan_phase_ang)).^2;
                C2 = 1 - S2;
                V_bone_T = C_AXIAL - ( C_AXIAL - C_RADIAL ) * ( ANISO_SHAPE_COEF.*S2.*C2 + C2.^2 ); % approximation of weak transverse isotropy
                thetaM = atan(tan_phase_ang); % phase angle (transmit)
                
                phiN = receiveAngles(selected_iRx(iTxRx),z,x); % receive group angle in cortical bone output by reconstruction code
                X = tan(phiN); % tan of receive group angle
                tan_phase_ang = ((A^3*X^3)/27 + (A*X)/3 + (B*X)/2 + (((A^3*X^3)/27 + (A*X)/3 + (B*X)/2)^2 - ((A^2*X^2)/9 - 1/3)^3)^(1/2))^(1/3) + ((A^2*X^2)/9 - 1/3)/((A^3*X^3)/27 + (A*X)/3 + (B*X)/2 + (((A^3*X^3)/27 + (A*X)/3 + (B*X)/2)^2 - ((A^2*X^2)/9 - 1/3)^3)^(1/2))^(1/3) + (A*X)/3;
                S2 = sin(atan(tan_phase_ang)).^2;
                C2 = 1 - S2;
                V_bone_R = C_AXIAL - ( C_AXIAL - C_RADIAL ) * ( ANISO_SHAPE_COEF.*S2.*C2 + C2.^2 ); % approximation of weak transverse isotropy
                phiN = atan(tan_phase_ang); % phase angle (receive)
                        
                AMat(iTxRx,1) = (cos(thetaM)/V_bone_T + cos(phiN)/V_bone_R);
                AMat(iTxRx,2) = (sin(thetaM)/V_bone_T + sin(phiN)/V_bone_R);
            end
            
            % Pseudo inverse of over etsimated problem in Least Square (No
            % Damping nor Weighting so far)
            % pseudoInvA           = (AMat'*AMat)\(AMat');
            % VectorFlowMat_LS(z,x,:) = pseudoInvA*squeeze(UMat(z,x,:)); % should be [ vz and vx ]!
            % Supposed to be more accurate.
            idx_good    = find(~isnan(AMat(:,1)));
            AMat        = AMat(idx_good,:)*scaling_factor;            
            pseudoInvA  = pinv(AMat); %(AMat'*AMat)\(AMat');
            
            
            n = size(AMat,1); % length(y); % number of residuals
            p = 2; % number of physical parameters
            X_t_X_inv = inv(AMat'*AMat);
            
            Umat_at_pixel = squeeze(Umat(z,x,idx_good,:));

            for it=1:NT-1

                y = squeeze(Umat_at_pixel(:,it));
                beta = pseudoInvA*y;

                residuals =  y - AMat*beta;
                var_res = ((residuals.')*residuals)/(n-p); % variance of the residuals
                var_estimates = X_t_X_inv*var_res; % variance-covariance matrix of the estimates of the physical variables 
                % The square roots of the diagonal elements of this matrix are
                % the standard deviations of the estimates, that is the errors
                % to be quoted on the determinations of the physical variables.

                vz(z,x,it) = beta(1);
                vx(z,x,it) = beta(2); 
                LS_sd_residuals(z,x,it) = sqrt(var_res);
                LS_sd_error_Vz(z,x,it) = sqrt(var_estimates(1,1));
                LS_sd_error_Vx(z,x,it) = sqrt(var_estimates(2,2));

            end
        end
    end
    
end

