function CORR = xcorr_nested(SPECULAR_TRANSFORM,SPECULAR_MODEL)
    [NTilts,Nz,Nx] = size(SPECULAR_TRANSFORM);
    CORR = zeros([Nz Nx 2*NTilts-1]);
    for iz =1:Nz
        for ix = 1:Nx
            x = SPECULAR_TRANSFORM(:,iz,ix);
            y = SPECULAR_MODEL(:,iz,ix);
            CORR(iz,ix,:)  = xcorr(x,y,'normalized');
        end
    end
end