# Compactness Measure and Classification

Ce projet implémente des mesures de compacité et des algorithmes de classification pour des nuages de points 3D, en utilisant la bibliothèque CGAL.

## Fonctionnalités

- Calcul de mesures de compacité pour des nuages de points 3D
- Classification de points 3D en utilisant Random Forest
- Génération de caractéristiques géométriques pour la classification
- Évaluation des résultats de classification

## Prérequis

- CMake (version 3.1 ou supérieure)
- CGAL
- Boost (avec les composants serialization et iostreams)
- OpenMP (optionnel, pour le parallélisme)

## Compilation

1. Clonez le dépôt :
   ```
   git clone https://github.com/votre-nom-utilisateur/compactness-measure-classification.git
   cd compactness-measure-classification
   ```

2. Créez un répertoire de build et compilez :
   ```
   mkdir build
   cd build
   cmake ..
   make
   ```

## Utilisation

### Mesure de compacité

./compactness_measure input.ply scale alpha offset

### Classification 

./compactness_classifier input.ply

### Évaluation de la Classification

./compactness_classifier_eval input.ply

## Structure du projet

- `compactness_measure.cpp`: Calcul des mesures de compacité
- `compactness_classifier.cpp`: Classification des points 3D
- `compactness_classifier_eval.cpp`: Évaluation des résultats de classification
- `compactness_measure_planimetric.cpp`: Version planimétrique de la mesure de compacité
- `example_ethz_random_forest.cpp`: Exemple d'utilisation du classificateur Random Forest ETHZ
