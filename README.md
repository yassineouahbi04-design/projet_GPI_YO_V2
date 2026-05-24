# **Projet_GPI_YO**

## **Présentation du projet**

Ce projet est un script *Python* établi sur de la **programmation orientée objet**. Il a pour objectif de lire un fichier de structure moléculaire 3D d'ARN au format **PDB** et d'en déterminer la **structure secondaire 2D**. Le script réalise sa tache en extrayant spécifiquement les nucléotides de l'ARN, puis en détectant leurs liaisons hydrogène afin de générer leur structure secondaire sous format dot-bracket.

## **Architecture du code**

Le code est structuré en **classes** qui vont s'imbriquer comme des poupées russes où chaque classe a un rôle précis :

- **Atom** (la base) : représente un **atome unique**. Elle stocke son nom et ses coordonnées cartésiennes ($X, Y, Z$). Elle intègre une méthode qui applique le théorème de *Pythagore* pour calculer la distance euclidienne entre deux atomes.

- **Nucleotide** (regroupement d'atomes) : représente une **base de l'ARN (A, U, G, C)**. Elle possède un dictionnaire pour stocker tous les atomes qui la composent, permetant de retrouver un atome instantanément par son nom.



