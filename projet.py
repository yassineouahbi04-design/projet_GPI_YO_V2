#!/usr/bin/env python3

""" On importe des modules nécessaires (sys pour la gestion des arguments, math pour les calculs de distance). """ 
import sys
import math

# ==========================================
# CONFIGURATION DES SITES DE LIAISON ET CRITÈRES DES LIAISONS HYDROGÈNE (A, U, G, C). 
# ==========================================

# On liste pour chaque nuclétide les atomes qui peuvent participer à des liaisons hydrogènes, classés en donneurs et accepteurs.
BOND_SITES = {
    'A': {'donors': ['N6'], 'acceptors': ['N1', 'N7', 'N3']},
    'G': {'donors': ['N1', 'N2'], 'acceptors': ['O6', 'N7', 'N3']},
    'C': {'donors': ['N4'], 'acceptors': ['N3', 'O2']},
    'U': {'donors': ['N3'], 'acceptors': ['O4', 'O2']}
}

# On définit le nombre de liaisons hydrogènes requises pour valider un appariement entre les nucléotides. 
REQUIRED_BONDS = {
    ('C', 'G'): 3, ('G', 'C'): 3,
    ('A', 'U'): 2, ('U', 'A'): 2,
    ('G', 'U'): 2, ('U', 'G'): 2
}

# On définit la distance entre donneurs et accepteurs pour considerer une liaison hydrogène (3 Å). 
# On s'accorde une marge de tolérance de 0.6 Å pour considérer qu'une liaison hydrogène est valide.
DIST_MIN = 2.4
DIST_MAX = 3.6

# ==========================================
# MODELISATION DES CLASSES (Atom, Nucleotide, RNAMolecule)
# ==========================================

class Atom:
    """ On définit ce qu'est un atome, extrait du fichier PDB, avec un nom (Carbone, Azote, etc.) et des coordonnées. """
    def __init__(self, name, x, y, z):
        self.name = name
        self.x = x
        self.y = y
        self.z = z

    def distance_to(self, other_atom):
        """ On calcule la distance euclidienne 3D entre l'atome actuel (self) et un autre atome (other_atom). On applique la formule qui est la suivante : d (ATOM 1 /ATOM 2)= sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2). L'ordre des atomes n'a pas d'importance pour le calcul de la distance. """ 
        return math.sqrt(
            (self.x - other_atom.x) ** 2 +
            (self.y - other_atom.y) ** 2 +
            (self.z - other_atom.z) ** 2
        )


class Nucleotide:
    """ On définit ce qu'est un nucléotide de la chaîne d'ARN. On lui donne un type (A, U, G, C), un identifiant numérique (position dans la chaîne) et un dictionnaire pour stocker les atomes qui lui sont associés. """
    def __init__(self, nt_name, nt_num):
        self.type = nt_name    # A, U, G, C
        self.id = nt_num       # Position numérique dans la chaîne.
        self.atoms = {}        # Dictionnaire stockant les objets atom qui constituent le nucléotide.

    def add_atom(self, atom):
        """ On ajoute un atome à ce nucléotide en le stockant dans le dictionnaire avec comme clé le nom de l'atome et comme valeur l'objet atom correspondant."""  
        self.atoms[atom.name] = atom  


class RNA_Molecule:
    """On définit un objet pour gérer l'ensemble de la structure de l'ARN et ses analyses."""
    def __init__(self):
        self.nucleotides = {}   # Dictionnaire pour stocker les nucléotides de l'ARN, avec comme clé leur identifiant numérique et comme valeur l'objet Nucleotide correspondant.

    def load_pdb(self, file_path):
        """ On lit le fichier PDB et extrait uniquement les atomes utiles pour l'ARN. """ 
        """ On se base sur les spécifications du format PDB pour extraire les informations de manière robuste, en utilisant les positions fixes des champs dans la ligne. On ne traite que les lignes commençant par "ATOM" et on ignore les autres types d'entrées (HETATM, etc.). """ 
        """ On ne considère que les nucléotides A, U, G, C et on ignore les autres résidus présentes dans le fichier. """
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith("ATOM"):
                        """Extraction stricte des données selon les spécifications des colonnes du format PDB officiel"""
                        atom_name = line[12:16].strip() # On extrait le nom de l'atome en utilisant les positions fixes (12-16) et on supprime les espaces superflus avec strip().
                        res_name = line[17:20].strip() # On extrait le nom du résidu (nucléotide) en utilisant les positions fixes (17-20) et on supprime les espaces superflus avec strip().
                        res_num = int(line[22:26].strip()) # On convertit le numéro de résidu en entier pour faciliter les comparaisons et le tri.
                        """On convertit les coordonnées en flottants pour les calculs de distance ultérieurs."""
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())

                        """On ne traite que les bases principales de l'ARN."""
                        if res_name in ['A', 'U', 'G', 'C']:
                            if res_num not in self.nucleotides:
                                self.nucleotides[res_num] = Nucleotide(res_name, res_num) # On crée un nouvel objet qui prend la classe Nucleotide pour ce résidu s'il n'existe pas déjà dans le dictionnaire.
                        
                            """On ajoute l'atome au nucléotide correspondant."""
                            self.nucleotides[res_num].add_atom(Atom(atom_name, x, y, z))
        except FileNotFoundError: # En cas d'erreur de lecture du fichier, on affiche un message d'erreur et on termine le programme avec un code de sortie pour indiquer une erreur.
            print(f"Error: The file '{file_path}' is not found.")
            sys.exit(1)

    def check_hydrogen_bond(self, nuc1, nuc2):
        """On compte le nombre de liaisons hydrogènes valides entre deux nucléotides."""
        bonds_count = 0
        
        """On récupère les rôles (donneurs / accepteurs) propres à chaque type de nucléotide."""
        sites1 = BOND_SITES[nuc1.type] # On accède au dictionnaire BOND_SITES pour obtenir les sites de liaison du nucléotide 1 en fonction de son type (A, U, G, C).
        sites2 = BOND_SITES[nuc2.type]

        """Possibilité 1 : Donneur du nuc1 vers accepteur du nuc2"""
        for d_name in sites1['donors']:
            for a_name in sites2['acceptors']:
                if d_name in nuc1.atoms and a_name in nuc2.atoms: # On vérifie que les atomes donneurs et accepteurs existent bien dans les nucléotides respectifs avant de calculer la distance.
                    dist = nuc1.atoms[d_name].distance_to(nuc2.atoms[a_name]) # On calcule la distance entre le donneur du nucléotide 1 et l'accepteur du nucléotide 2 en utilisant la méthode distance_to de la classe Atom.
                    if DIST_MIN <= dist <= DIST_MAX:
                        bonds_count += 1

        """Possibilité 2 : Donneur du nuc2 vers accepteur du nuc1"""
        for d_name in sites2['donors']:
            for a_name in sites1['acceptors']:
                if d_name in nuc2.atoms and a_name in nuc1.atoms:
                    dist = nuc2.atoms[d_name].distance_to(nuc1.atoms[a_name])
                    if DIST_MIN <= dist <= DIST_MAX:
                        bonds_count += 1

        return bonds_count
    

    def find_base_pairs(self):
        """On parcourt toutes la molécule pour trouver les nucléotides appariés."""
        pairs = {} # Dictionnaire pour stocker les appariements trouvés, avec comme clé l'identifiant numérique d'un nucléotide et comme valeur l'identifiant de son partenaire d'appariement.
        sorted_ids = sorted(self.nucleotides.keys()) # On trie les identifiants numériques des nucléotides pour garantir un ordre cohérent lors de la comparaison des paires.
        
        """Double boucle pour comparer chaque paire de nucléotides uniques (i < j)."""
        for i in range(len(sorted_ids)): # On parcourt les indices des nucléotides triés, en commençant par le premier (i) et en comparant avec tous les suivants (j > i) pour éviter les comparaisons redondantes et les auto-appariements.
            for j in range(i + 1, len(sorted_ids)):
                id1, id2 = sorted_ids[i], sorted_ids[j] # On récupère les identifiants numériques des deux nucléotides à comparer.
                nuc1, nuc2 = self.nucleotides[id1], self.nucleotides[id2]
                
                """On Vérifie si le couple est autorisé (CG, AU, GU)"""
                pair_type = (nuc1.type, nuc2.type)
                if pair_type in REQUIRED_BONDS:
                    nb_bonds = self.check_hydrogen_bond(nuc1, nuc2)

                    if nb_bonds >= REQUIRED_BONDS[pair_type]: # Si le nombre de liaisons trouvées correspond au critère requis pour ce type de paire, on les considère comme appariés et on les stocke dans le dictionnaire des paires.
                        pairs[id1] = id2
                        pairs[id2] = id1
        return pairs
    

def generate_output(self):
        """On génère la séquence d'ARN et la notation Dot-Bracket correspondante."""
        sorted_ids = sorted(self.nucleotides.keys()) # On trie les identifiants numériques des nucléotides pour garantir un ordre cohérent dans la séquence et la structure générées.
        pairs = self.find_base_pairs() # Retourne un dictionnaire des appariements trouvés.
        
        sequence = ""
        dot_bracket = []

        """Construction de la structure secondaire"""
        for res_num in sorted_ids:
            sequence += self.nucleotides[res_num].type # On ajoute le type de nucléotide (A, U, G, C) à la séquence finale.
            
            if res_num in pairs: # Si ce nucléotide est apparié, on vérifie son partenaire d'appariement pour déterminer s'il doit être représenté par une parenthèse ouvrante ou fermante dans la notation Dot-Bracket.
                partner_num = pairs[res_num] # On récupère l'identifiant numérique du partenaire d'appariement.
                if res_num < partner_num: # Si l'identifiant du nucléotide actuel est inférieur à celui de son partenaire, on le représente par une parenthèse ouvrante "(" dans la notation Dot-Bracket, sinon par une parenthèse fermante ")".
                    dot_bracket.append("(")
                else:
                    dot_bracket.append(")")
            else:
                dot_bracket.append(".") # Si le nucléotide n'est pas apparié, on le représente par un point ".". 

        return sequence, "".join(dot_bracket) # On retourne la séquence d'ARN complète et la notation Dot-Bracket correspondante en joignant les éléments de la liste dot_bracket en une chaîne de caractères.


