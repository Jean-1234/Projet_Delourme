# Projet Delourme — Méthode LOD (Localized Orthogonal Decomposition)

Implémentation MATLAB/Octave de la méthode des éléments finis standard et de la méthode **LOD (Localized Orthogonal Decomposition)** pour la résolution numérique de problèmes elliptiques à coefficients fortement oscillants (homogénéisation périodique en 1D).

Ce projet illustre, de manière progressive, comment la méthode LOD permet d'obtenir une précision d'ordre optimal sur un maillage grossier, là où les éléments finis standards échouent à capturer les oscillations à l'échelle ε sans raffiner massivement le maillage.

## Contexte

On considère le problème modèle 1D :

```
-d/dx( a_ε(x) du/dx ) = f(x)   sur (0,1)
u(0) = u(1) = 0
```

où `a_ε` est un coefficient oscillant à l'échelle `ε` (périodique ou discontinu). La méthode des éléments finis classiques nécessite `h << ε` pour converger correctement, ce qui devient prohibitif quand `ε → 0`. La méthode LOD construit un espace grossier enrichi (fonctions de base "LOD") qui retrouve un taux de convergence optimal en `H` indépendamment de `ε`.

## Structure du dépôt

```
.
├── FEM_Standard/        Éléments finis P1 classiques (référence / mise en évidence du problème)
├── LOD_Galerkin/        Méthode LOD "globale" (patchs = domaine entier)
└── LOD_Localisation/    Méthode LOD localisée (patchs de taille k, phases offline/online)
```

### `FEM_Standard/`

Solveur éléments finis P1 standard, utilisé comme référence et pour mettre en évidence la perte de convergence en présence d'un coefficient oscillant.

| Fichier | Rôle |
|---|---|
| `main.m` / `mainv2.m` | Scripts principaux : choix du cas test, étude de convergence |
| `parameters.m` | Définition des cas tests (coefficient `a_ε`, second membre `f`, valeurs de `ε`) |
| `mesh1D.m` | Génération du maillage 1D |
| `assembleFEM.m` / `assembleFEMv2.m` | Assemblage des matrices de rigidité/masse et du second membre |
| `Rigidite.m`, `Masse.m` | Matrices élémentaires de rigidité et de masse |
| `solveFEM.m` | Résolution du système linéaire |
| `exactSolution.m` | Solution exacte / solution de référence |
| `plotSolutions.m` | Visualisation des solutions et courbes de convergence |

### `LOD_Galerkin/`

Première implémentation de la méthode LOD, où le correcteur de chaque fonction de base est calculé en résolvant un problème sur le domaine entier (pas encore de localisation par patch).

| Fichier | Rôle |
|---|---|
| `mainLOD.m` | Script principal : choix du cas test, calcul LOD, étude de convergence |
| `parametersLOD.m` | Définition des cas tests |
| `assembleLOD.m` | Construction de l'espace LOD et assemblage du système grossier |
| `phiLOD.m` | Construction des fonctions de base LOD (correcteurs) |
| `reconstructLOD.m` | Reconstruction de la solution fine à partir de la solution grossière |
| `erreurLOD.m` | Calcul des erreurs L² et H¹ |
| `Plotcomparaison.m` | Comparaison graphique LOD vs solution exacte |

### `LOD_Localisation/`

Implémentation complète de la méthode LOD avec **localisation par patchs** (paramètre `k`), suivant la décomposition classique en phase *offline* (construction des correcteurs locaux) et phase *online* (résolution du problème grossier).

| Fichier | Rôle |
|---|---|
| `mainLOD1.m` | Script principal : choix du cas test, étude de convergence en fonction de `H` et `ε` |
| `construire_maillages.m` | Génération des maillages fin (`Xh`) et grossier (`XH`) |
| `assembler_systeme_fin.m` | Assemblage des matrices de rigidité/masse sur le maillage fin |
| `construire_patch.m` | Construction du patch de taille `k` autour d'un nœud grossier |
| `construire_IH.m` | Opérateur d'interpolation / quasi-interpolation grossier |
| `phase_offline.m` | Calcul des correcteurs LOD locaux par patch (résolution des problèmes locaux dans le noyau de l'interpolation) |
| `phase_online.m` | Assemblage et résolution du système LOD grossier, reconstruction de la solution fine |
| `solveLOD.m` | Résolution du système linéaire grossier |
| `plotsolution.m` | Visualisation des résultats |
| `test_*.m` | Scripts de test unitaires (maillage, assemblage, méthode LOD) |

## Cas tests disponibles

1. **Coefficient constant** (`a_ε = 1`) — vérification de l'ordre de convergence théorique (référence sans oscillation).
2. **Coefficient périodique oscillant** `a_ε(x) = 1 / (2 - cos(2πx/ε))` — cas d'homogénéisation classique, testé pour plusieurs valeurs de `ε` (`2⁻², 2⁻⁴, 2⁻⁶, 2⁻⁸`).
3. **Coefficient discontinu oscillant** (cas binaire par morceaux) — *(disponible dans certains scripts, ex. `parameters.m`)*.

Chaque script principal (`main.m`, `mainLOD.m`, `mainLOD1.m`) demande à l'utilisateur de choisir le cas test au lancement.

## Utilisation

Les scripts sont écrits en MATLAB/Octave et ne nécessitent aucune toolbox externe.

```bash
# Dans Octave ou MATLAB, depuis le dossier souhaité :
cd FEM_Standard
main          # ou mainv2

cd ../LOD_Galerkin
mainLOD

cd ../LOD_Localisation
mainLOD1
```

Chaque script demande de choisir un cas test, puis :
- calcule la solution numérique sur une série de maillages grossiers,
- la compare à une solution de référence (exacte ou fine),
- affiche les courbes d'erreur L²/H¹ en échelle log-log avec estimation de l'ordre de convergence,
- trace les solutions obtenues.

## Résultats attendus

- **FEM standard** : convergence en `O(h)` (H¹) / `O(h²)` (L²) uniquement quand `h << ε` ; dégradation nette sinon.
- **LOD (Galerkin global et localisé)** : convergence en `O(H)` (H¹) / `O(H²)` (L²) **indépendamment de ε**, dès que le patch est suffisamment grand (paramètre `k` dans `LOD_Localisation`).

Les figures générées (`*.png`) dans chaque dossier illustrent ces comportements pour les différents cas tests.

## Prérequis

- MATLAB (R2018b ou supérieur recommandé) **ou** GNU Octave ≥ 6.0
- Aucune toolbox supplémentaire requise

## Auteur

Projet réalisé dans le cadre d'un cours sur l'homogénéisation numérique et la méthode LOD (Målqvist–Peterseim), basé sur les supports de B. Verfürth et S. Fliss.

## Références

- A. Målqvist, D. Peterseim, *Numerical homogenization by localized orthogonal decomposition*, SIAM, 2020.
- P. Henning, A. Målqvist, *Localized Orthogonal Decomposition*, arXiv:1308.3379.
- B. Verfürth, *Analytical and Numerical Homogenization*, Winter Semester 2022/23.
- S. Fliss, *Homogénéisation périodique* (ENSTA / Polytechnique).
