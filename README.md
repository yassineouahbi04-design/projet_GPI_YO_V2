# **Projet_GPI_YO**

## **Présentation du projet**

Ce projet est un script *Python* établi sur de la **programmation orientée objet**. Il a pour objectif de lire un fichier de **structure 3D** d'ARN au format **PDB** et d'en déterminer la **structure secondaire 2D**. Le script réalise sa tache en extrayant spécifiquement les nucléotides de l'ARN, puis en détectant les appariements afin de générer leur structure secondaire sous format Dot-Bracket.

## **Paramètrage biologiques et physique**

Pour garantir la justesse scientifique, le script s'appuie sur des règles définies en début de fichier :

- **Les sites de liaison (BOND_SITES)** :

Le programme connaît la structure exacte de chaque nucléotide et sait précisément quels atomes de ces derniers peuvent agir comme donneurs ou accepteurs.

- **Les paires autorisées (REQUIRED_BONDS)** : 

Le script ne valide pas n'importe quel appariement. Il reconnaît uniquement les paires classiques de Watson-Crick (G-C nécessitant 3 liaisons hydrogènes, A-U nécessitant 2 liaisons hydrogènes) ainsi que la paire wobble (G-U nécessitant 2 liaisons hydrogènes).

- **Les seuils de distance (DIST_MIN, DIST_MAX)** : 

La longueur moyenne d'une liaison hydrogène est d'environ 3 Ångströms. On applique une tolérance de $\pm$ 0.4 Å. Pour qu'une liaison soit valide, la distance entre l'atome donneur et l'accepteur doit impérativement être comprise entre 2.6 et 3.4 Å.

## **Architecture du code**

Le code est structuré en **classes** qui vont s'imbriquer comme des poupées russes où chaque classe a un rôle précis :

- **Atom** (la base) : représente un **atome unique**. Elle stocke son nom et ses coordonnées cartésiennes ($X, Y, Z$). Elle intègre une méthode qui applique le théorème de *Pythagore* pour calculer la distance euclidienne entre deux atomes.

- **Nucleotide** (regroupement d'atomes) : représente une **base de l'ARN (A, U, G, C)**. Elle possède un dictionnaire pour stocker tous les atomes qui la composent, permetant de retrouver un atome instantanément par son nom.

- **RNA_Molecule** (regroupement de nucléotides) : représente la **molécule d'ARN complète**. Elle contient un dictionnaire de tous les nucléotides rangés par leur identifiant numérique (position dans la chaine). Elle gère la **lecture du fichier**, la **recherche des appariement** et la **génération du résultat**.

## **Fonctionnement du code**

Lorsque le script est exécuté, il déroule automatiquement les étapes suivantes :

- **Etape 1 : Lecture et filtrage (méthode load_pdb)** 

Le script ouvre le fichier PDB et lit les lignes une par une. Il applique un **double filtre** en ne gardant que les lignes commençant par ATOM et en vérifiant que le résidu est bien de l'ARN (A, U, G, C). Il extrait ensuite les coordonnées des atomes à des positions de texte fixes (dictées par le format PDB) et crée les objets en mémoire.

- **Étape 2 : Détection des liaisons (méthode check_hydrogen_bond)**

Le programme teste toutes les combinaisons possibles entre les atomes "Donneurs" du premier nucléotide et les atomes "Accepteurs" du second. Si la distance entre deux atomes compatibles entre dans la fourchette autorisée, la liaison est validée.

- **Étape 3 : Recherche des paires de bases (méthode find_base_pairs)**

Une **double boucle** teste toutes les combinaisons possibles de nucléotides (sans jamais tester deux fois la même paire ni comparer un nucléotide à lui-même). Si deux nucléotides forment un appariement valide (ex : $G-C$ ou $A-U$) et possèdent le nombre de liaisons hydrogène requises, ils sont considerés comme appariés.

- **Étape 4 : Traduction textuelle (méthode generate_output)**

Le script parcourt l'ARN du début à la fin pour écrire la séquence en lettres. En parallèle, il dessine la **structure en Dot-Bracket**. Si un nucléotide est seul, il écrit un point. S'il est apparié, il compare leurs positions. Le plus petit indice reçoit une parenthèse ouvrante "(" et le plus grand reçoit une parenthèse fermante ).

## **Optimisation**

Pour éliminer les **artefacts (faux couples)** que j'ai pu constater, deux filtre ont été appliqué avant de generer le Dot-Bracket:

- Un premier qui **interdit à deux nucléotides trop proches de s'apparier** (il faut au moins 3 nucléotides d'écart entre eux). La molécule n'a ni la **l'espace** ni la **fléxibilité** pour faire un virage aussi serré et apparié deux nucléotides trop proche. 

- Un deuxième qui supprime les **couples isolé**. Si une parenthèe (appariement) est toute seule, elle est remplacée par un point. **Il faut au moins un voisin apparié juste à coté pour avoir un structure (appariement) stable**.  

## **Prérequis et lancement du script**

- **Python 3** doit etre installé sur la machine.

- Le script utilise uniquement les **modules natifs de Python** (sys et math). Aucune installation supplémentaire n'est necessaire. 

- Un **fichier au format *PDB*** contenant une molécule d'ARN doit etre mis à disposition. Ce fichier doit se trouver dans le **meme dossier** que le script.  

- Poour lancé le script:
    - Ouvrir le terminal de l'appareil. 
    - Se déplacer dans le dossier ou se trouve le script et le/les fichier(s) *PDB*. 
    - Exécuter la commande suivant : ***python3 projet.py file_name.pdb***