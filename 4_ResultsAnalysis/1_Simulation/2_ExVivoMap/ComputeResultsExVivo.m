clear; 
close all;
clc;
addpath(genpath('~/Documents'));
addpath(genpath('/calculSSD/salome'));
%%
pathSimuAll = '/calculSSD/salome/Simulation-10mai'; 
simuDirAll = dir(sprintf('%s/Bone*', pathSimuAll));

dirFlags = [simuDirAll.isdir]; 
simuDirAll = simuDirAll(dirFlags); 

bones = {'245D', '227G', '267G'};
slices = {'Z1', 'Z2', 'Z3', 'Z4', 'Z0'};
slicesNb = {{1134, 3195, 3852, 5511, 0000}, {2002, 3595, 5614, 6721, 0000}, {1700, 3410, 5478, 6716, 0000}}; 
%%
Results = struct();
%% Compute results
Results.Map = ComputeMap('Map', bones, slices, slicesNb, simuDirAll, Results);
Results.SpecuMap = ComputeMap('SpecuMap', bones, slices, slicesNb, simuDirAll, Results);

%% Plot maps
MapTitle = 'Interface profile for different bones';
PlotAll('Map', Results, bones, slices, MapTitle);

%%
MapTitle = 'Specular Index';
PlotAll('SpecuMap', Results, bones, slices, MapTitle);

%%
Param = ComputeMap('InterfaceParam', bones, slices, slicesNb, simuDirAll, Results);
% Results.InitParam = struct();
% Results.InitParam.Corr = table('Size', [numel(rmsAll), numel(corrAll)], ...
%                 'VariableType', repmat({'double'}, 1, numel(corrAll)), ...
%                 'VariableNames', cellstr(string(corrAll)), ...
%                 'RowNames', cellstr(string(rmsAll))); 
% Results.InitParam.Rms = table('Size', [numel(rmsAll), numel(corrAll)], ...
%                 'VariableType', repmat({'double'}, 1, numel(corrAll)), ...
%                 'VariableNames', cellstr(string(corrAll)), ...
%                 'RowNames', cellstr(string(rmsAll))); 
% 
% for i = 1:numel(bones)
%     for j = 1:numel(slices)
%         % Results.InitParam.Rms{i,j} = Param{i,j}{1}{1};
%         % Results.InitParam.Corr{i,j} = Param{i,j}{1}{2};
%         Results.InitParam.Rms{i,1} = Param{i,1}{1}{1};
%         Results.InitParam.Corr{i,1} = Param{i,1}{1}{2};
%     end
% end

%%
function[Table] = ComputeMap(MapType, bones, slices, slicesNb, simuDirAll, Results)
    % Compute parameters once 
    load('/calculSSD/salome/Simulation-10mai/boneSpeed.mat', 'boneSpeed')
    kcAll = [0.041 0.046 0.057 0.064; 0.051 0.046 0.033 0.037; 0.033 0.027 0.030 0.027];
    
    % Generate Map and probability regarding the case 
    Table = table('Size', [numel(bones), numel(slices)], ...
                'VariableType', repmat({'cell'}, 1, numel(slices)), ...
                'VariableNames', cellstr(string(slices)), ...
                'RowNames', cellstr(string(bones)));

    for i = 1:numel(bones)  
        bone.id = bones{i};         % Bone from ex-vivo files
        for j = 1:numel(slicesNb{i})
            bone.image = slicesNb{i}{j};        % Slice selected
            
            simuName = fullfile(simuDirAll(1).folder, sprintf('Bone%s-Image%04d/', bone.id, bone.image));
            try
            % Compute the required type of simulation 
                switch MapType
                    case 'Map'
                        [Map,~,~] = SimSonic2DReadMap2D(fullfile(simuName, 'Geometry.map2D'));
                        Table{i,j}{1} = Map;
                       
                    case 'SpecuMap'
                        if exist(fullfile(simuName, 'postProcess.mat'))
                            postProcess = load(fullfile(simuName, 'postProcess.mat'));
                            Table{i,j}{1} = postProcess.SpecularProbaMap;
                            disp(i)
                            disp(j)
                        end
                    % case 'Stat'
                    % 
                    %     Table{i,j}{1} = probability;
                    case 'InterfaceParam'
                        binaryImage = logical(Results.Map{i,j}{1});
                     
                        [~, roughness, ~, ~] = FilterProfil(binaryImage, kcAll(i,j));
                        % fprintf(['\nThe cut off frequency for the bone ' bone.id ' and the slices %04d is %.03f mm -1'], bone.image,kc)
                
                        % Compute the correlation length and rms of the roughness profile 
                        Rq = rms(roughness);
                        Corr = ComputeCorr(roughness, grid.step);
                
                        % Compute the porosity in a surface of one wavelength around the boundary
                        nbWavelength = 1/2;
                        lambda =  boneSpeed.(simuName(36:43)).(simuName(45:53)){1} / (probe.fc)*1e3; %(medium.cp(2) / (probe.fc))*1e3; 
                        width = lambda * nbWavelength;
                        [boundaryEndost, boundaryPores] = ExtractBoundary(binaryImage);
                        [porosity, poreSize, ~] = ComputePorosity(binaryImage, boundaryEndost, boundaryPores, width, false);
                        
                        Table{i,j}{1} = {Rq, Corr, porosity, poreSize};
                    otherwise 
                        warning('Unexpected Map type. No data computed.')
                end
            end
        end
    end
end

function[] = PlotAll(MapType, Struct, bones, slices, MapTitle)
    Table = getfield(Struct, MapType);
    
    figure
    t = tiledlayout(numel(bones),numel(slices),'TileSpacing','tight');
    % Plot Map
    for i = 1:numel(bones)
        for j = 1:numel(slices)
            % Plot each map
            nexttile(t)
            switch MapType
                case 'Map'
                    try
                        imagesc(Table{i,j}{1});
                    end
                case 'SpecuMap'
                    try
                        imagesc(Table{i,j}{1}.Bone);
                        axis image
                        colormap jet
                    end
                % case 'Stat'
                %     try
                %         recorded = LoadRfData(parameters.probe, simuDir1);
                %         [~, reconstruction] = GenerateParamRecon(recorded);
                %         plot(reconstruction.Xmm, Table.linearROI{i,j}{1});
                %         ylim([0, 1]);
                %     end
                otherwise
                     warning('Unexpected Map type. No data ploted.')
            end
            % Remove x and y axes
            set(gca, 'xtick', [], 'ytick', []);
        
        
            if j == 1
                ylabel(bones{i}, 'Interpreter', 'latex', 'FontSize', 16);
            end
            if i == 1
                title(slices{j}, 'Interpreter', 'latex', 'FontSize', 16);
            end
        end
    end
    % Set axis
    title(t, 'Zone', 'Interpreter', 'latex', 'FontSize', 22)
    ylabel(t, 'Bone', 'Interpreter', 'latex', 'FontSize', 22)
    xlabel(t, MapTitle, 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold')

end 