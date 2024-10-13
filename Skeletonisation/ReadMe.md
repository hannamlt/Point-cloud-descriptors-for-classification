# Compactness Measure, Classification, and Skeletonization

Ce projet implémente des mesures de compacité, des algorithmes de classification pour des nuages de points 3D, ainsi que des techniques de squelettisation, en utilisant principalement la bibliothèque CGAL.

## Fonctionnalités

- Calcul de mesures de compacité pour des nuages de points 3D
- Classification de points 3D en utilisant Random Forest
- Génération de caractéristiques géométriques pour la classification
- Évaluation des résultats de classification
- Squelettisation de maillages 3D
- Extraction de caractéristiques basées sur le squelette

## Prérequis

- CMake (version 3.1 ou supérieure)
- CGAL
- Boost (avec les composants serialization et iostreams)
- OpenMP (optionnel, pour le parallélisme)

## Compilation

1. Clonez le dépôt :
   ```
   git clone https://github.com/votre-nom-utilisateur/compactness-measure-classification-skeletonization.git
   cd compactness-measure-classification-skeletonization
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

### Évaluation de la Classfication

./compactness_classifier input.ply

### Squelettisation et extraction de caractéristiques

./skeletonise_feature_extraction input.ply scale alpha offset

### Conversion de polylignes en nuage de points

./skeletonise_polylines_to_point_cloud input.ply scale alpha offset



## Structure du projet

- `compactness_measure.cpp`: Calcul des mesures de compacité
- `compactness_classifier.cpp`: Classification des points 3D
- `compactness_classifier_eval.cpp`: Évaluation des résultats de classification
- `compactness_measure_planimetric.cpp`: Version planimétrique de la mesure de compacité
- `example_ethz_random_forest.cpp`: Exemple d'utilisation du classificateur Random Forest ETHZ
- `skeletonise_classifier.cpp`: Classification basée sur la squelettisation
- `skeletonise_mean_curvature.cpp`: Squelettisation par flux de courbure moyenne
- `skeletonise_feature_extraction.cpp`: Extraction de caractéristiques basées sur le squelette
- `skeletonise_polylines_to_point_cloud.cpp`: Conversion de polylignes en nuage de points





