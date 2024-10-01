import numpy as np
from random import randrange

class Tree:
    def __init__(self, n):
        self.edges = {i: [] for i in range(n)}
        self.labels = {}

    def add_edge(self, u, v):
        self.edges[u].append((v, None))
        self.edges[v].append((u, None))

    def is_leaf(self, v):
        return v in self.labels

    def create_topological_order(self):
        visited = set()
        order = []
        def dfs(v):
            if v in visited:
                return
            visited.add(v)
            for neighbor, _ in self.edges[v]:
                dfs(neighbor)
            order.append(v)
        dfs(max(self.edges.keys()))  # Assuming the last node is the root
        return order

class LabeledTree:
    def __init__(self, n):
        self.labels = {}

    def initialize_from(self, T):
        self.labels = {v: T.labels[v] for v in T.labels}

    def set_label(self, v, label):
        self.labels[v] = label

def SmallParsimony(T, alphabet='ATGC'):
    '''
    SmallParsimony

    Find the most parsimonious labeling of the internal nodes of a rooted tree.

    Given: An integer n followed by an adjacency list for a rooted binary tree with n leaves labeled by DNA strings.

    Return: The minimum parsimony score of this tree, followed by the adjacency list of the tree
            corresponding to labeling internal nodes by DNA strings in order to minimize the parsimony score of the tree.

    '''
    def calculate_score(v, n, scores):
        child_scores = np.full((2, n), np.nan)
        for i, (e, _) in enumerate(T.edges[v]):
            child_scores[i, :] = np.min(scores[e] + Delta, axis=1)
        return np.sum(child_scores, axis=0)

    def backtrack(v, current_assignment, scores):
        for v_next, _ in T.edges[v]:
            if T.is_leaf(v_next):
                continue
            if v_next not in assignments.labels:
                assignments.labels[v_next] = ''
            min_score = np.min(scores[v_next, :])
            indices = [i for i in range(n) if scores[v_next, i] == min_score]
            matched = False
            for i in indices:
                if alphabet[i] == current_assignment:
                    matched = True
                    assignments.set_label(v_next, assignments.labels[v_next] + current_assignment)
                    backtrack(v_next, current_assignment, scores)

            if not matched:
                next_assignment = alphabet[indices[randrange(0, len(indices))]]
                assignments.set_label(v_next, assignments.labels[v_next] + next_assignment)
                backtrack(v_next, next_assignment, scores)

    def update_assignments(v, s):
        c = alphabet[np.argmin(s)]
        assignments.set_label(v, assignments.labels[v] + c if v in assignments.labels else c)
        return c

    def SmallParsimonyC(Character):
        # Assign character indices only for leaf nodes
        leaf_nodes = {v: alphabet.index(c) for v, c in zip(T.labels.keys(), Character)}
        scores = np.full((len(T.edges), n), np.inf)

        for v in Nodes:
            if T.is_leaf(v):
                scores[v, leaf_nodes[v]] = 0  # Use the leaf node's correct index in character list
            else:
                scores[v, :] = calculate_score(v, n, scores)

        v = Nodes[-1]
        backtrack(v, update_assignments(v, scores[v, :]), scores)
        return np.min(scores[v, :])

    Nodes = T.create_topological_order()
    assignments = LabeledTree(len(T.edges))
    assignments.initialize_from(T)
    n = len(alphabet)
    Delta = np.ones((n, n))
    np.fill_diagonal(Delta, 0)

    # Iterate over each character in the sequence and calculate parsimony for all positions
    parsimony_score = sum([SmallParsimonyC([v[i] for _, v in T.labels.items()]) for i in range(len(next(iter(T.labels.values()))))])

    # Calculate the adjacency list with differences in labels
    edges = []
    for parent, children in T.edges.items():
        for child, _ in children:
            if parent in assignments.labels and child in assignments.labels:
                parent_label = assignments.labels[parent]
                child_label = assignments.labels[child]
                diff = sum([parent_label[i] != child_label[i] for i in range(len(parent_label))])
                edges.append(f'{parent_label}->{child_label}:{diff}')
    
    return parsimony_score, edges

# Test data

T = Tree(7)
T.add_edge(4, 6)
T.add_edge(5, 6)
T.labels[4] = "CAAATCCC"
T.labels[4] = "ATTGCGAC"
T.labels[5] = "CTGCGCTG"
T.labels[5] = "ATGGACGA"

parsimony_score, adjacency_list = SmallParsimony(T)

# Output the result
print(parsimony_score)
for edge in adjacency_list:
    print(edge)