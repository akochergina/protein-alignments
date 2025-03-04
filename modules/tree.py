""" Defines the tree class
"""

class Tree:
    def __init__(self, val = None, left=None, right=None):
        if val != None:
            self.val = val
        else:
            self.val = None

        if left != None:
            self.left = left
        else:
            self.left = None

        if right != None:
            self.right = right
        else:
            self.right = None

    def print_tree(self, level=0, prefix="Root: "):
        """
        Représente un arbre binaire de manière hiérarchique.
        :param node: Le nœud actuel (de type Tree)
        :param level: Niveau d'indentation actuel
        :param prefix: Préfixe à afficher avant chaque nœud
        """
        if self is not None:
            print(" " * (level * 4) + prefix + str(self.val))
            if self.left or self.right:  # Vérifie s'il y a des enfants avant d'afficher "None"
                if self.left:
                    self.left.print_tree(level + 1, "Left: ")
                else:
                    print(" " * ((level + 1) * 4) + "Left: None")
                if self.right:
                    self.right.print_tree(level + 1, "Right: ")
                else:
                    print(" " * ((level + 1) * 4) + "Right: None")
    

    def print_tree_sequences(self, liste_de_sequences, level=0, prefix="Root: "):
        """
        Affiche un arbre binaire en remplaçant les listes d'indices par les séquences correspondantes.
        
        :param liste_de_sequences: Liste des séquences utilisées pour construire l'arbre.
        :param level: Niveau d'indentation actuel.
        :param prefix: Préfixe à afficher avant chaque nœud.
        """
        if self is not None:
            # Vérifie si self.val est une liste d'indices
            if isinstance(self.val, list):
                sequences = [liste_de_sequences[i] for i in self.val if isinstance(i, int) and i < len(liste_de_sequences)]
                sequence_str = ", ".join(sequences)
            else:
                sequence_str = str(self.val)

            print(" " * (level * 4) + prefix + sequence_str)
            
            # Vérifie la présence d'enfants et affiche l'arbre récursivement
            if self.left or self.right:
                if self.left:
                    self.left.print_tree_sequences(liste_de_sequences, level + 1, "Left: ")
                else:
                    print(" " * ((level + 1) * 4) + "Left: None")
                
                if self.right:
                    self.right.print_tree_sequences(liste_de_sequences, level + 1, "Right: ")
                else:
                    print(" " * ((level + 1) * 4) + "Right: None")
