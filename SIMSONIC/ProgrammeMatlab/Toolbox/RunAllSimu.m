path_to_run_simul = '~/Documents/BoneRugosity/SIMSONIC/'; % Ton dossier contenant le script de simulation
repertoire = fullfile(path_to_run_simul, 'Simulation'); % Ton répertoire contenant tes dossiers de simulation 

% Récupérer la liste de tous les fichier et dossier dans ton répertoire de simul.
dossiers_de_simul = dir(sprintf('%s/simulation_rms_*',repertoire));

dirFlags = [dossiers_de_simul.isdir]; % Flag qui te dit si c'et un répertoire ou un fichier

dossiers_de_simul = dossiers_de_simul(dirFlags); % Juste pour eter sur qu'on ne traite que les dossiers

for simu_dir=dossiers_de_simul
	parameters = load(fullfile(simu_dir, 'parameters.mat'));
    recorded = LoadRfData(parameters.probe, simu_dir);
end