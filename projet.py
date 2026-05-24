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