# **Projet_GPI_YO**

## **Présentation du projet**

Ce projet est un script *Python* établi sur de la **programmation orientée objet**. Il a pour objectif de lire un fichier de structure moléculaire 3D d'ARN au format **PDB** et d'en déterminer la **structure secondaire 2D**. Le script réalise sa tache en extrayant spécifiquement les nucléotides de l'ARN, puis en détectant leurs liaisons hydrogène afin de générer leur structure secondaire sous format dot-bracket.

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

Le programme compare les atomes des nucléotides deux par deux. Il regarde dans un dictionnaire de référence quels atomes sont donneurs ou accepteurs. Si la distance entre deux atomes compatibles entre dans la fourchette autorisée (DIST_MIN et DIST_MAX), la liaison est validée.

- **Étape 3 : Recherche des paires de bases (méthode find_base_pairs)**

Une **double boucle** compare tous les nucléotides uniques de la molécule. Si deux nucléotides forment un appariement valide (ex : $G-C$ ou $A-U$) et possèdent le nombre minimal de liaisons hydrogène requises, ils sont déclarés "appariés".

- **Étape 4 : Traduction textuelle (méthode generate_output)**

Le script parcourt l'ARN du début à la fin pour écrire la séquence de lettres. En parallèle, il dessine la **structure en Dot-Bracket**. Si un nucléotide est seul, il écrit un point. S'il est apparié, il compare leurs positions. Le plus petit indice reçoit une parenthèse ouvrante "(" et le plus grand reçoit une parenthèse fermante ).