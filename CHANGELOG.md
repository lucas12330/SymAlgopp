# Changelog

Ce document répertorie tous les ajouts et correctifs majeurs de la bibliothèque SymAlgo++.

## [Unreleased] - Branche feature/equation-classique

### Ajouts majeurs (Major Additions)
* **Intégration Symbolique** : Implémentation d'un moteur de calcul formel d'intégrales par reconnaissance de motifs pour les expressions classiques (`ASTNode::integrer`). Prise en charge des polynômes, fonctions trigonométriques (sinus, cosinus), et de la linéarité.
* **Limites Symboliques (L'Hôpital)** : Ajout d'un système de calcul de limites formelles (`ASTNode::limite`). Pour le cas des fractions avec indétermination de type $\frac{0}{0}$, l'algorithme dérive symboliquement le numérateur et le dénominateur selon la célèbre règle de L'Hôpital.
* **Développements Limités (DL)** : Implémentation de la génération symbolique de développements de Taylor et Maclaurin (`ASTNode::DL`) via le calcul récursif de dérivées formelles (ex: approximation polynomiale).
* **Tracé et Tableaux de Points Adaptatifs** : Ajout de la méthode `genererPointsTrace` utilisant un algorithme récursif d'échantillonnage adaptatif. Cela permet de réduire radicalement la quantité de points calculés sur les sections linéaires des courbes tout en gardant une excellente précision dans les variations.
* **Dérivation Généralisée** : Le système dérive désormais les fonctions puissances de la forme $u(x)^{v(x)}$ en utilisant rigoureusement le logarithme népérien, tout en optimisant explicitement le comportement pour les exposants constants.
* **Nouvelle primitive mathématique** : Ajout du nœud `Logarithme` (`ast_ln`).

### Correctifs majeurs (Major Patches)
* *Aucun correctif majeur pour le moment.*
