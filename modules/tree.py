""" 
Defines the tree class
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
        Prints the binary tree in a hierarchical structure.

        Parameters:
        ----------
        level : int, optional
            The current indentation level (default is 0).
        prefix : str, optional
            The prefix to display before each node (default is "Root: ").

        Returns:
        -------
        None
            The function prints the tree structure to the console.

        Notes:
        ------
        - If a node has no left or right child, "None" is printed for clarity.
        - The function automatically calls `visualize_tree()` to display the tree graphically if used in a Jupyter Notebook.
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
        self.visualize_tree()
    

    def print_tree_sequences(self, liste_de_sequences, level=0, prefix="Root: "):
        """
        Prints a binary tree, replacing index lists with their corresponding sequences.

        Parameters:
        ----------
        sequence_list : list of str
            A list of sequences used to construct the tree.
        level : int, optional
            The current indentation level (default is 0).
        prefix : str, optional
            The prefix to display before each node (default is "Root: ").

        Returns:
        -------
        None
            The function prints the tree structure to the console.

        Notes:
        ------
        - If `self.val` is a list of indices, it is replaced with the corresponding sequences from `sequence_list`.
        - If a node has no left or right child, "None" is printed for clarity.
        - The function recursively traverses the tree and prints it in a hierarchical format.
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
