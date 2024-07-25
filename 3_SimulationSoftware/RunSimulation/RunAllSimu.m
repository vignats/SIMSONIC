pathSimsonic = '~/Documents/BoneRugosity/SIMSONIC/';    % Dossier contenant le script de simulation
pathSimu = '/calculSSD/salome/Simulation-10mai';        % Répertoire contenant tes dossiers de simulation 

% Récupérer la liste de tous les fichier et dossier dans ton répertoire de simul.
simuDir = dir(sprintf('%s/Bone*', pathSimu));

dirFlags = [simuDir.isdir]; % Flag qui te dit si c'et un répertoire ou un fichier
simuDir = simuDir(dirFlags); % Juste pour eter sur qu'on ne traite que les dossiers

for idx = length(simuDir)
    simu_dir = fullfile(pathSimu, simuDir(idx).name);
	terminal_string = sprintf('%s./run_simul.sh %s',pathSimsonic,simu_dir);
	system(terminal_string)
end

%%
path_to_run_simul = '~/Documents/BoneRugosity/SIMSONIC/'; % Ton dossier contenant le script de simulation
repertoire = '/calculSSD/salome/Simulation-10mai'; % Ton répertoire contenant tes dossiers de simulation 

% Récupérer la liste de tous les fichier et dossier dans ton répertoire de simul.
simuDir = dir(sprintf('%s/Bone*',repertoire));

dirFlags = [simuDir.isdir]; % Flag qui te dit si c'et un répertoire ou un fichier
simuDir = simuDir(dirFlags); % Juste pour eter sur qu'on ne traite que les dossiers

for simu_dir = simuDir
	terminal_string = sprintf('%s/run_simul.sh %s',path_to_run_simul, simu_dir.name);
	system(terminal_string)
end