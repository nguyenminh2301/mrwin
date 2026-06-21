# mrwin : Statistiques de Gain Causales pour Critères Composites Hiérarchisés

[English](README.md) | [Tiếng Việt](README.vi.md) | [中文](README.zh.md) | [日本語](README.ja.md) | [한국어](README.ko.md) | [Français](README.fr.md) | [Deutsch](README.de.md)

## Ce que fait ce package

`mrwin` est un package R qui répond à une question spécifique en épidémiologie :

**Une exposition génétiquement prédite déplace-t-elle les personnes vers des profils de résultats globaux meilleurs ou pires, lorsque les résultats sont classés cliniquement par gravité ?**

Dans de nombreuses maladies, les patients ne subissent pas un seul événement. Une personne atteinte d'une maladie cardio-rénale peut mourir, être hospitalisée pour insuffisance cardiaque ou développer un déclin rénal. Ces événements ne sont pas interchangeables : le décès est pire que l'hospitalisation, qui est pire qu'une anomalie biologique. Les analyses standard traitent tous les événements de manière égale ou les analysent un par un, perdant ainsi l'ordre clinique qui importe le plus aux patients et aux cliniciens.

`mrwin` combine deux cadres statistiques pour résoudre ce problème :

1. **Statistiques de gain (Win statistics)** — comparaisons par paires qui respectent les hiérarchies de gravité clinique (le décès l'emporte sur l'hospitalisation, qui l'emporte sur le déclin des biomarqueurs).
2. **Randomisation mendélienne** — variables instrumentales génétiques qui permettent des affirmations causales même lorsque les facteurs de confusion ne sont pas mesurés.

Le résultat est un **ratio de gain causal standardisé par dose (DS-CWR)** : un nombre unique qui résume si une exposition génétiquement prédite rend le profil de résultats priorisés meilleur ou pire.

---

## Randomisation mendélienne : une brève introduction

### Le problème : la confusion non mesurée

En épidémiologie observationnelle, nous voulons souvent savoir si une exposition (ex. : cholestérol LDL) cause un résultat (ex. : insuffisance cardiaque). La difficulté est que les personnes avec un LDL élevé diffèrent des personnes avec un LDL bas de bien d'autres manières — alimentation, exercice, statut socio-économique, autres médicaments. Ce sont des facteurs de confusion. Même après ajustement pour les facteurs de confusion mesurés, des facteurs de confusion non mesurés peuvent subsister.

### La solution : la génétique comme expérience naturelle

À la conception, chaque personne hérite aléatoirement d'un allèle de chaque parent à chaque locus génétique. C'est ce qu'on appelle la **ségrégation mendélienne**. Si une variante génétique affecte une exposition (ex. : une variante du gène LDLR augmente le cholestérol LDL), alors les personnes porteuses de cette variante sont, en moyenne, exposées à un LDL plus élevé tout au long de leur vie — non pas à cause de leur alimentation ou de leur mode de vie, mais à cause de leur génotype.

Cette assignation aléatoire est la base de la **randomisation mendélienne (MR)**. En utilisant des variantes génétiques comme **variables instrumentales**, la MR estime l'effet causal de l'exposition sur le résultat, en contournant la confusion non mesurée.

### Hypothèses clés de la MR

Pour que l'instrument génétique soit valide, trois conditions doivent être remplies :

1. **Pertinence** : La variante génétique doit être associée à l'exposition.
2. **Indépendance** : La variante génétique ne doit pas être associée aux facteurs de confusion (cela découle de la ségrégation mendélienne dans les études bien conçues).
3. **Restriction d'exclusion** : La variante génétique doit affecter le résultat uniquement via l'exposition, et non par d'autres voies.

### Types de données en MR

| Type de données | Description | Exemple |
|---|---|---|
| Niveau individuel (Individual-level) | Une ligne par personne avec génotype, exposition, résultat | Une cohorte de 10 000 personnes avec données SNP, mesures de LDL et dossiers hospitaliers |
| Niveau résumé (Summary-level) | Estimations d'effet par SNP issues de GWAS | « Le SNP rs12345 a un effet de 0,05 (SE 0,01) sur le LDL » |
| Score de risque polygénique (PRS) | Somme pondérée de nombreux SNPs | PRS pour le cholestérol LDL, calculé à partir de 100 SNPs |

`mrwin` utilise des **données de niveau individuel** avec un **instrument PRS**. Le PRS sert d'instrument continu qui capture la prédisposition génétique cumulative envers l'exposition.

---

## Statistiques de gain : l'estimande cible

### Qu'est-ce qu'une statistique de gain ?

Une statistique de gain compare chaque paire d'individus dans une étude. Pour chaque paire (A, B), la comparaison suit une hiérarchie clinique :

1. D'abord, comparer sur le résultat le plus important (ex. : décès). Si A est décédé et B non, B gagne.
2. Si égalité sur le premier résultat, passer au suivant (ex. : hospitalisation).
3. Continuer à descendre la hiérarchie jusqu'à ce qu'une personne gagne ou que tous les résultats soient à égalité.

Le **ratio de gain (win ratio)** est le nombre de gains pour le groupe le plus exposé divisé par le nombre de gains pour le groupe le moins exposé.

### Pourquoi les statistiques de gain sont importantes

| Méthode | Traite les événements également ? | Respecte l'ordre clinique ? | Gère les risques concurrents ? |
|---|---|---|---|
| Délai jusqu'au premier événement (Time-to-first-event) | Oui | Non | Non |
| Analyse de critère composite (Composite endpoint analysis) | Oui | Non | Partiellement |
| Ratio de gain / statistique de gain (Win ratio / win statistic) | Non | Oui | Oui |

Les statistiques de gain ont été introduites par Pocock et al. (2012) et ont été adoptées dans plus de 36 essais cliniques randomisés entre 2022 et 2024.

### Le ratio de gain causal (cCWR)

Les statistiques de gain standard nécessitent une randomisation ou une forte ignorabilité. En épidémiologie observationnelle, cette hypothèse est rarement défendable.

`mrwin` définit le **ratio de gain causal continu (cCWR)** comme un estimande de gradient à travers les quantiles de l'instrument génétique :

- Diviser la population en strates ordonnées selon le PRS (ex. : déciles).
- Dans chaque paire de strates adjacentes, calculer le ratio de gain.
- Regrouper à travers les strates en utilisant les gradients standardisés par instrument et la méta-analyse GLS.

Le résultat est un **effet causal marginal au niveau de la population** qui ne nécessite pas une forte ignorabilité.

---

## Formules fondamentales : explication en langage clair

### 1. Le noyau de comparaison hiérarchique

Pour deux individus *i* et *j*, comparer les résultats de la priorité la plus élevée à la plus basse :

- Si *j* a eu l'événement à la priorité *k* mais pas *i* (ou *i* a survécu plus longtemps) : *i* gagne (+1).
- Si *i* a eu l'événement à la priorité *k* mais pas *j* (ou *j* a survécu plus longtemps) : *i* perd (-1).
- Si égalité : passer à la priorité suivante.

**En clair** : « Cette personne avait-elle un meilleur profil de résultats que l'autre, en respectant l'ordre de gravité clinique ? »

### 2. Le ratio de gain spécifique à la strate

Au sein d'une strate (ex. : le 6e décile de PRS contre le 5e) :

```
theta_d = (nombre de gains pour la strate supérieure) / (nombre de pertes)
log_theta_d = log(theta_d)
```

**En clair** : « Parmi les personnes ayant une prédisposition génétique légèrement plus élevée à l'exposition, avaient-elles tendance à gagner plus souvent qu'à perdre ? »

### 3. Le gradient standardisé par instrument (ISG)

```
Delta_X_d = moyenne(exposition dans la strate d) - moyenne(exposition dans la strate d-1)
ISG_d = log_theta_d / Delta_X_d
```

**En clair** : « De combien le log du ratio de gain change-t-il par augmentation d'une unité de l'exposition, tel que prédit par la génétique ? »

### 4. Regroupement GLS

Les valeurs ISG à travers toutes les paires de strates sont regroupées à l'aide des moindres carrés généralisés (Generalized Least Squares) avec un estimateur de covariance à rétrécissement :

```
delta_GLS = ISG regroupé à travers toutes les strates
DS-CWR = exp(delta_GLS)
```

**En clair** : « Quel est l'effet causal global de l'exposition sur le profil de résultats priorisés ? »

- **DS-CWR > 1** : L'exposition est associée à un meilleur profil de résultats.
- **DS-CWR < 1** : L'exposition est associée à un moins bon profil de résultats.
- **DS-CWR = 1** : Aucune preuve d'un effet causal.

### 5. Le bootstrap multiplicatif

L'incertitude est estimée par :
1. Perturbation des poids GWAS (incertitude externe).
2. Tirage de poids multiplicateurs aléatoires (incertitude interne).
3. Recalcul du pipeline complet pour chaque itération bootstrap.

Cela donne un intervalle de confiance qui tient compte des deux sources d'incertitude.

### 6. Le diagnostic de pléiotropie par étapes (SDPD)

Le SDPD teste si l'instrument génétique affecte le résultat de plus haute priorité par des voies autres que l'exposition (pléiotropie directe). Cela se fait en ajustant une régression MR-Egger sur les effets par SNP de l'exposition et du résultat.

**En clair** : « L'instrument génétique est-il valide, ou affecte-t-il le résultat par d'autres voies ? »

---

## Exemples cliniques

### Exemple 1 : Maladie cardiovasculaire

**Question de recherche** : Le cholestérol LDL génétiquement prédit aggrave-t-il la trajectoire cardio-rénale globale ?

**Hiérarchie des critères** (de la priorité la plus élevée à la plus basse) :
1. Décès
2. Hospitalisation pour insuffisance cardiaque
3. Déclin rénal

**Interprétation** : Si DS-CWR = 0,82 (IC à 95 % : 0,70 à 0,96), cela signifie que par augmentation d'une unité du LDL génétiquement prédit, le profil de résultats cardio-rénaux priorisés s'aggrave de 18 %. L'IC ne franchit pas 1, suggérant un effet nocif statistiquement significatif.

**Signification clinique** : Un LDL plus élevé ne provoque pas seulement des crises cardiaques — il déplace toute la trajectoire de la maladie vers un décès plus précoce, plus d'hospitalisations et un déclin rénal plus rapide, dans cet ordre clinique.

### Exemple 2 : Démence

**Question de recherche** : Les troubles du sommeil génétiquement prédits déplacent-ils les personnes vers des résultats neurocognitifs sévères plus précoces ?

**Hiérarchie des critères** :
1. Décès
2. Diagnostic de démence
3. Admission en maison de retraite

**Interprétation** : Si DS-CWR = 0,91 (IC à 95 % : 0,78 à 1,06), l'estimation ponctuelle suggère un effet nocif, mais l'IC franchit 1. Il n'y a pas de preuve solide d'un effet causal.

**Pourquoi c'est important** : Dans la MR standard du sommeil sur la démence, les personnes qui décèdent avant de développer une démence sont perdues pour l'analyse. C'est une forme de biais du survivant : une variante de sommeil nocive peut sembler protectrice parce que ses porteurs ne vivent jamais assez longtemps pour être diagnostiqués. `mrwin` évite cela en traitant le décès comme l'événement de plus haute priorité, de sorte que la survie différentielle est absorbée dans l'estimande plutôt que conditionnée.

### Exemple 3 : Épidémiologie du cancer

**Question de recherche** : L'IMC génétiquement prédit aggrave-t-il les résultats de cancer priorisés ?

**Hiérarchie des critères** :
1. Décès par cancer
2. Progression ou métastase
3. Récidive

**Interprétation** : Si DS-CWR = 0,75 (IC à 95 % : 0,60 à 0,94), cela signifie qu'un IMC génétiquement prédit plus élevé déplace la trajectoire du cancer vers des résultats moins bons à chaque niveau de la hiérarchie.

**Signification clinique** : L'IMC n'augmente pas seulement le risque de cancer — il aggrave tout le cours de la maladie, de la récidive jusqu'au décès.

---

## Installation

```r
# Depuis GitHub (version de développement)
# install.packages("devtools")
devtools::install_github("nguyenminh2301/mrwin")
```

---

## Démarrage rapide

### Données simulées (aucune cohorte réelle nécessaire)

```r
library(mrwin)

# Simuler un jeu de données cardio-rénal
cfg <- mrwin_config(n_outcome = 500, m_snps = 20, seed = 1)
dat <- mrwin_simulate(cfg)

# Ajuster le modèle
fit <- mrwin(
  endpoint = mrwin_endpoint(dat$time, dat$status, c("death", "hf", "renal")),
  genotype = dat$G,
  exposure = dat$X,
  gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
  controls = mrwin_controls(n_strata = 5, bootstrap = 100, seed = 2)
)

# Examiner les résultats
print(fit)
summary(fit)
plot(fit)
```

### Lire les résultats

```
mrwin fit
  N: 500
  SNPs: 20
  Priorities: 3
  Strata: 5
  Bootstrap: 100 (100 valid)
  Backend: dense
  Adjustment: none
  delta_GLS: -0.2000
  DS-CWR: 0.8187
  95% CI: 0.7000 to 0.9600
  Q: 1.2345 ( df = 2 , p = 0.5394 )
```

| Champ | Signification | Ce qu'il faut rechercher |
|---|---|---|
| `delta_GLS` | Effet regroupé sur échelle logarithmique | Négatif = nocif, positif = bénéfique |
| `DS-CWR` | Effet exponentiel (ratio de gain) | < 1 = nocif, > 1 = bénéfique, 1 = pas d'effet |
| `95% CI` | Intervalle de confiance bootstrap | Franchit-il 1 ? |
| `Q` | Statistique d'hétérogénéité | Petite valeur p = l'effet varie selon les strates |
| `Warnings` | Mises en garde structurées | Toujours vérifier avant d'interpréter |

### Interprétation du résultat

Si DS-CWR = 0,82 (IC à 95 % : 0,70 à 0,96) :

- L'exposition est associée à un profil de résultats priorisés 18 % moins bon.
- L'IC ne franchit pas 1, l'effet est donc statistiquement significatif.
- Vérifier la valeur p de Q : si significative, l'effet peut varier selon la plage d'exposition.
- Vérifier les avertissements : si `sdpd_rejected`, la pléiotropie peut invalider l'estimation.

---

## Modèle pour cohorte réelle

```r
snp_cols <- paste0("rs", 1:20)

fit <- mrwin(
  data = dat,
  endpoint = mrwin_endpoint(
    time = c("t_death", "t_hf", "t_renal"),
    status = c("d_death", "d_hf", "d_renal"),
    priority = c("death", "heart_failure", "renal_decline")
  ),
  genotype = snp_cols,
  exposure = "ldl_cholesterol",
  gwas = mrwin_gwas(beta = beta_hat, se = se_beta, snp = snp_cols),
  covariates = c("age", "sex", "PC1", "PC2"),
  controls = mrwin_controls(
    n_strata = 10,
    bootstrap = 500,
    seed = 20260510,
    adjustment = "ordinal_iptw",
    sdpd_scale = "both"
  )
)

summary(fit)
plot(fit, type = "forest")
mrwin_report(fit, file = "analysis_report.md", format = "markdown")
tidy(fit)
```

---

## Formats de sortie

### Print

```r
print(fit)
```

Affiche le résultat principal : taille de l'échantillon, SNPs, strates, validité du bootstrap, DS-CWR, IC, statistique Q et avertissements.

### Summary

```r
summary(fit)
```

Produit un rapport formaté avec :
- Estimation principale DS-CWR, IC primaire (de Fieller), valeur p, et IC borné
  par pléiotropie optionnel
- Gradients des strates adjacentes (log-theta, Delta-X, ISG, CWR)
- Statistique d'hétérogénéité Q
- L'intervalle delta-method comme référence étiquetée (l'intervalle de Fieller
  est l'intervalle primaire ; voir *Intervalles de confiance* ci-dessous)
- Diagnostics SDPD (lorsqu'activés)
- Avertissements structurés

### Tidy (pour tableaux et analyses complémentaires)

```r
tidy(fit)
#   term estimate    delta se_delta statistic p_value  ci_low ci_high ...
# DS-CWR    0.82 -0.2000   0.075    -2.667  0.0077  0.7000  0.9600 ...
```

### Rapport

```r
mrwin_report(fit, format = "markdown")
mrwin_report(fit, file = "report.md", format = "markdown")
```

### Graphiques

```r
plot(fit, type = "isg")       # Gradients ISG adjacents (par défaut)
plot(fit, type = "forest")    # Graphique en forêt des CWR adjacents
plot(fit, type = "bootstrap") # Distribution bootstrap de log-theta
```

---

## Ajustement par covariables

Lorsque les facteurs de confusion sont mesurés (âge, sexe, composantes principales), utilisez l'ajustement IPTW ordinal :

```r
fit_adj <- mrwin(
  endpoint = mrwin_endpoint(dat$time, dat$status, c("death", "hf", "renal")),
  genotype = dat$G,
  exposure = dat$X,
  gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
  covariates = cbind(u_proxy = scale(dat$U)),
  controls = mrwin_controls(
    n_strata = 5,
    bootstrap = 100,
    seed = 2,
    adjustment = "ordinal_iptw",
    sdpd_scale = "both",
    pleiotropy_bias_radius = 0.05
  )
)
```

L'ajustement IPTW améliore la précision en équilibrant les covariables entre les strates de PRS. Le SDPD teste la pléiotropie directe.

---

## Forces et limites

### Forces

- **Inférence causale sous confusion non mesurée** : Utilise des instruments génétiques, pas seulement l'exposition observée.
- **Respecte la hiérarchie clinique** : Le décès est priorisé sur l'hospitalisation, qui est priorisée sur le déclin des biomarqueurs.
- **Gère naturellement les risques concurrents** : La survie différentielle est absorbée dans l'estimande, pas conditionnée.
- **Mesure synthétique unique** : Un seul DS-CWR remplace de multiples analyses par composante.
- **Cadre diagnostique** : Le SDPD teste la restriction d'exclusion ; la statistique Q teste l'hétérogénéité ; les avertissements structurés signalent les problèmes.
- **Reproductible** : Le bootstrap est contrôlé par graine (seed) ; toutes les sorties sont déterministes étant donné la graine.

### Limites

- **Données de niveau individuel requises** : Le package nécessite actuellement des données avec une ligne par participant contenant le génotype, l'exposition et les résultats. La MR sur données résumées uniquement n'est pas prise en charge.
- **Non-collapsibilité** : Le ratio de gain n'est pas collapsible. Le DS-CWR marginal diffère des ratios de gain conditionnels. C'est une propriété de l'estimande, pas un bug.
- **Instrument PRS linéaire** : Le package suppose un score de risque polygénique linéaire. Les interactions gène-gène non linéaires ne sont pas modélisées.
- **Covariance GWAS diagonale** : Le bootstrap actuel utilise les erreurs standard SNP par SNP, pas la matrice de covariance LD complète. Cela peut sous-estimer l'incertitude lorsque les SNPs sont en déséquilibre de liaison.
- **Vulnérabilité à la pléiotropie** : Le cCWR est sensible à la pléiotropie contaminant la hiérarchie. Une pléiotropie au niveau de la mortalité aussi faible que gamma = 0,05 peut réduire la couverture à 12 %.
- **Exigences de taille d'échantillon** : Des échantillons de taille biobanque (N > 100 000) sont un prérequis statistique strict pour une inférence fiable.
- **Nombre de strates `D`** : l'estimateur de type décile nécessite de choisir un nombre de strates. Un estimateur de gradient continu, contrôlé par largeur de bande (dont l'estimateur décile est le cas particulier boxcar exact), a été implémenté et validé, mais des vérifications externes ont montré qu'il est plus bruité que le décile et ne réduit pas la sensibilité à `D` (le gradient standardisé par instrument est un rapport, et un lissage plus fin rétrécit son dénominateur). La pratique recommandée est l'estimateur discret avec une analyse de sensibilité sur `D`, et non une reparamétrisation continue (voir `dev/wp15-continuous-isg.md`).

### Backends de performance et d'inférence (en option)

Le backend par défaut est inchangé, mais le travail de la Phase II (`dev/acceleration-roadmap.md`) a ajouté des alternatives validées et optionnelles via `mrwin_controls()` :

- `backend = "fast"` — un noyau win/loss sous-quadratique, compilé (Rcpp) (`Theta(N log^{K-1} N)` pour `K` niveaux de priorité) qui est identique bit à bit au backend dense et porte l'estimateur à l'échelle biobanque (ex. : `K = 3`, `N = 80 000` en moins d'une seconde).
- `inference = "analytic"` — une variance par fonction d'influence en forme close qui reproduit le bootstrap multiplicatif (avec un terme de Monte-Carlo exact pour l'incertitude des poids GWAS), supprimant la boucle bootstrap pour la composante d'échantillonnage.

> **Stratification doublement ordonnée — non recommandée (résultat négatif).**
> `mrwin_doubly_ranked_strata()` (Tian/Burgess) est implémentée, mais une
> calibration externe a montré que `stratification = "doubly_ranked"` est
> **incompatible** avec l'estimande DS-CWR inter-strates : elle équilibre
> l'instrument entre les strates, de sorte que le contraste entre strates
> adjacentes devient piloté par les facteurs de confusion (erreur de type I
> ~1,0 sous confusion). `mrwin()` l'exécute toujours mais émet un avertissement
> `doubly_ranked_invalid` ; la valeur par défaut est `"prs_rank"`. Voir
> `inst/spec/validation-findings.md`.

#### Intervalles de confiance

L'**intervalle rapporté en priorité est l'intervalle de Fieller**. Une simulation de calibration externe a montré que l'intervalle de ratio par delta-method d'origine est trop large (il n'atteint pas l'erreur de type I nominale), tandis que la construction de Fieller sur la même covariance est correctement calibrée ; l'intervalle delta-method n'est conservé que comme référence étiquetée. Sous un instrument faible, l'intervalle de Fieller est honnêtement rapporté comme non borné plutôt que comme un ensemble faussement étroit.

---

## Hypothèses

Le cadre `mrwin` repose sur les hypothèses suivantes :

1. **Pertinence** : Le PRS doit être associé à l'exposition. Les instruments faibles produisent des estimations instables (vérifier l'avertissement `weak_instrument`).

2. **Indépendance** : Le PRS ne doit pas être associé aux facteurs de confusion. Ceci est satisfait par la ségrégation mendélienne dans les populations homogènes. La stratification de population peut violer cette hypothèse ; utilisez les composantes principales comme covariables.

3. **Restriction d'exclusion** : Le PRS doit affecter le résultat uniquement via l'exposition. Le SDPD teste cela pour le critère de plus haute priorité. Si rejeté, l'estimation peut être invalide.

4. **Monotonicité** : L'instrument doit déplacer l'exposition dans la même direction pour tous les individus. Les violations peuvent biaiser l'estimation.

5. **Pas de modification d'effet par la survie** : L'effet causal ne doit pas varier entre ceux qui survivent et ceux qui ne survivent pas. Ceci est impossible à tester.

6. **Ordre de priorité correct** : La hiérarchie clinique doit être spécifiée avant l'analyse. Le package ne la déduit pas des données.

---

## Codes d'avertissement

| Code | Signification | Que faire |
|---|---|---|
| `weak_instrument` | Déplacement phénotypique proche de zéro ou IC de Fieller non borné | Traiter l'estimation comme instable |
| `sdpd_rejected` | L'intercept MR-Egger suggère une pléiotropie | Ajouter une discussion de sensibilité |
| `sdpd_underpowered` | Trop peu de SNPs pour le SDPD | Signaler la limitation de puissance |
| `positivity_failure` | ESS IPTW trop petit dans une strate | Inspecter le chevauchement des covariables |
| `bridged_strata` | Strates sautées en raison d'un échec de positivité | Signaler le schéma de pontage |
| `discordant_components` | Les directions MR des composantes sont conflictuelles | Ne pas revendiquer un effet global |

---

## Moteur de simulation

Exécutez des grilles de scénarios pour valider la méthode :

```r
scenarios <- mrwin_scenarios(
  mrwin_config(n_outcome = 200, m_snps = 10, seed = 10),
  scenarios = c("A_null", "B_valid_IV", "C_pleiotropy", "D_hierarchy_discordant")
)

grid <- mrwin_run_simulation_grid(
  scenarios = scenarios, n_iter = 2,
  controls = mrwin_controls(n_strata = 4, bootstrap = 20, seed = 11)
)

mrwin_simulation_summary(grid)
```

| Scénario | Signification | Comportement attendu |
|---|---|---|
| A_null | Aucun effet causal | DS-CWR proche de 1, taux de rejet proche de 5 % |
| B_valid_IV | Effet causal valide | DS-CWR dans la direction attendue |
| C_pleiotropy | Pléiotropie directe | Le rejet SDPD augmente |
| D_hierarchy_discordant | Effets opposés des composantes | Avertissement émis |

---

## Citation de ce package

Si vous utilisez `mrwin` dans vos recherches, veuillez citer :

> Nguyen Thien Minh, N. Ahmad Aziz. "Causal Win Statistics: Integrating Instrumental Variable Estimation with Hierarchical Composite Endpoints." *arXiv preprint*, 2026.

BibTeX :

```bibtex
@article{NguyenAziz2026,
  title={Causal Win Statistics: Integrating Instrumental Variable
         Estimation with Hierarchical Composite Endpoints},
  author={Nguyen Thien Minh and N. Ahmad Aziz},
  journal={arXiv preprint},
  year={2026}
}
```

**Auteurs :**
- Nguyen Thien Minh, University of Medicine and Pharmacy at Ho Chi Minh City, Vietnam (minhnt@ump.edu.vn)
- N. Ahmad Aziz (auteur correspondant), German Center for Neurodegenerative Diseases (DZNE) and University of Bonn, Germany (Ahmad.Aziz@dzne.de)

---

## Références associées

**Statistiques de gain :**
- Pocock SJ, et al. The win ratio: a new approach to the analysis of composite endpoints in clinical trials. *Eur Heart J*. 2012;33(14):1744-1749.
- Bebu I, Lachin JM. Large sample inference for a win ratio analysis of a composite endpoint based on prioritized components. *Biostatistics*. 2016;17(1):178-191.
- Even Z, Josse A. Causal win ratio. *arXiv preprint*. 2025.

**Randomisation mendélienne :**
- Lawlor DA, et al. Mendelian randomization: using genes as instruments for making causal inferences in epidemiology. *Stat Med*. 2008;27(8):1133-1163.
- Davey Smith G, Hemani G. Mendelian randomization: genetic anchors for causal inference in epidemiological studies. *Hum Mol Genet*. 2014;23(R1):R89-R98.
- Bowden J, et al. Mendelian randomization with invalid instruments: effect estimation and bias detection through Egger regression. *Int J Epidemiol*. 2015;44(2):512-525.

**Scores de risque polygénique :**
- Choi SW, et al. Tutorial: a guide to performing polygenic risk score analyses. *Nat Protoc*. 2020;15(9):2759-2772.

**Modèles multi-états :**
- Putter H, Fiocco M, Geskus RB. Tutorial in biostatistics: competing risks and multi-state models. *Stat Med*. 2007;26(11):2389-2430.

---

## Structure du package

```
mrwin/
  R/
    api.R              # Workflow haut niveau mrwin()
    kernel.R           # Noyau de comparaison hiérarchique
    estimate.R         # Estimateur DS-CWR, regroupement GLS
    bootstrap.R        # Inférence par bootstrap multiplicatif
    adjustment.R       # IPTW ordinal, ESS, pontage
    sdpd.R             # Diagnostic de pléiotropie par étapes (SDPD)
    simulate.R         # Processus générateur de données de simulation
    simulation_engine.R # Exécuteur de grille de scénarios
    validation.R       # Validation des entrées
    methods.R          # print, summary, plot, tidy, report
    benchmark.R        # Benchmarking de performance
    backend_sparse.R   # Backend sparse
    config.R           # Configuration de simulation
    strata.R           # Assignation des strates de PRS
  tests/
    testthat/          # 273 tests unitaires
  vignettes/           # 3 vignettes
  inst/spec/           # Spécifications d'implémentation
```

---

## Licence

MIT
