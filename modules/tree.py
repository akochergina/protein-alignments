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