#!/bin/bash

# Répertoire parent des simulations
parent_dir="/calculSSD/salome"

# Boucle pour parcourir les dossiers de simulation
for simu_dir in "$parent_dir"/"simulation_"*; do
	# Extraire le nom du dossier de simulation
	simu_name=$(basename "$simu_dir")
        
        # Vérifier si le fichier postProcess.mat existe dans le dossier de simulation
        if [ -f "$simu_dir/postProcess.mat" ]; then
            mkdir -p "$simu_name"
            
            # Copier le fichier .mat du sous-dossier source vers le sous-dossier correspondant dans le répertoire de destination
            cp "$simu_dir/postProcess.mat" "$simu_name/postProcess.mat"
    fi
done

