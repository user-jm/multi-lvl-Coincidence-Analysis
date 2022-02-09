#!/usr/bin/env python3

# simple tree with deep search function 

class Tree:

    def __init__(self, name):
        # initiate function
        self.name = name
        self.children = {}

    def __iter__(self):
        # how to iterate through the tree
        return iter(self.children)

    def __str__(self):
        # output of str(tree)
        return self.name

    def __repr__(self):
        # output of print(tree)
        return 'Tree("{}")'.format(self.name)

    def add_child(self, child):
        # add child below the current node
        self.children[child] = child
        return child
    
    def deep_search(self, include_self=True):
        # deep search through the tree
        if include_self:
            yield self
        for child in self.children:
            yield child
            yield from child.deep_search(False)
