from collections import defaultdict

# SECTION 1: Parsing the input into a tree structure
def parse_input(data):
    tree = defaultdict(list)
    dna_strings = {}
    n = int(data[0])
    
    # Fill the tree and DNA string dictionary
    for line in data[1:]:
        parent, child = line.split("->")
        parent = int(parent.strip())
        
        if child.strip().isdigit():
            child = int(child.strip())
            tree[parent].append(child)
        else:
            tree[parent].append(len(dna_strings))
            dna_strings[len(dna_strings)] = child.strip()

    # Checkpoint 1: Verify the parsed tree and DNA strings
    print("\nCheckpoint 1: Tree and DNA strings parsed")
    print("Tree:", tree)
    print("DNA strings:", dna_strings)
    
    return tree, dna_strings, n

# SECTION 2: Initialize scores for leaf nodes, character by character
def initialize_leaves(dna_strings, alphabet, pos):
    num_symbols = len(alphabet)
    s = defaultdict(lambda: [float('inf')] * num_symbols)  # Each node gets a list for scores for A, C, G, T
    node_labels = {}  # Stores the final computed label for each node

    # Initialize scores for the leaves (for a specific position in the DNA sequence)
    for leaf in dna_strings:
        leaf_char = dna_strings[leaf][pos]  # Get the character at the given position
        node_labels[leaf] = leaf_char  # Save the character for leaves at this position
        for i, symbol in enumerate(alphabet):
            if leaf_char == symbol:
                s[leaf][i] = 0  # Score is 0 if the symbol matches the leaf character
            else:
                s[leaf][i] = float('inf')  # Otherwise, it's infinity

    # Checkpoint 2: Verify the initialized leaf scores
    print(f"\nCheckpoint 2 (position {pos}): Leaf initialization")
    print("Leaf scores:", dict(s))
    print("Leaf labels:", node_labels)

    return s, node_labels

# SECTION 3: Process internal nodes and compute parsimony for a specific character position
def small_parsimony(tree, dna_strings, n, pos):
    # Set of symbols (nucleotides in this case)
    alphabet = ['A', 'C', 'G', 'T']
    
    # Initialize scores for leaf nodes for a given position in the DNA string
    s, node_labels = initialize_leaves(dna_strings, alphabet, pos)
    tags = defaultdict(int)  # Tag each node (0 for unprocessed, 1 for processed)

    # Function to calculate the mutation cost δ(i, k) (1 if i != k, 0 otherwise)
    def delta(i, k):
        return 1 if i != k else 0

    # Find a ripe node (a node with unprocessed tag, whose children are processed)
    def find_ripe_node():
        for node in tree:
            if tags[node] == 0:  # Node is unprocessed
                left, right = tree[node]
                if tags[left] == 1 and tags[right] == 1:
                    return node
        return None

    # Process all internal nodes
    while True:
        v = find_ripe_node()
        if v is None:
            break  # No ripe nodes left, we're done
        
        tags[v] = 1  # Mark node as processed
        left, right = tree[v]

        # Calculate the score for each character k at node v
        for k in range(len(alphabet)):
            min_left = float('inf')
            min_right = float('inf')

            # Calculate the minimum score for the left and right children
            for i in range(len(alphabet)):
                min_left = min(min_left, s[left][i] + delta(alphabet[i], alphabet[k]))
            for j in range(len(alphabet)):
                min_right = min(min_right, s[right][j] + delta(alphabet[j], alphabet[k]))
            
            # Update the score for character k at node v
            s[v][k] = min_left + min_right

        # Assign a label to the internal node based on the minimum score
        best_char_idx = s[v].index(min(s[v]))  # Get the index of the minimum score
        node_labels[v] = alphabet[best_char_idx]  # Assign the corresponding character to the node

    # Root is the last internal node processed
    root = max(tree)
    parsimony_score = min(s[root])  # Minimum score over all characters at the root

    # Checkpoint 3: Verify the processed internal node scores and labels
    print(f"\nCheckpoint 3 (position {pos}): Internal nodes processed")
    print("Internal node scores:", dict(s))
    print("Node labels:", node_labels)
    print("Parsimony score:", parsimony_score)

    return node_labels, parsimony_score

# SECTION 4: Format the output adjacency list with mutation costs
def format_output(tree, node_labels, dna_strings, parsimony_score, n, pos):
    output = [f"Position {pos}: Parsimony score = {parsimony_score}"]
    
    # Collect the adjacency list with mutation costs
    for parent in tree:
        for child in tree[parent]:
            parent_label = node_labels.get(parent, dna_strings.get(parent, ""))  # Ensure valid label
            child_label = node_labels.get(child, dna_strings.get(child, ""))  # Ensure valid label
            mutations = sum(1 for a, b in zip(parent_label, child_label) if a != b)
            output.append(f"{parent_label}->{child_label}:{mutations}")
    
    # Checkpoint 4: Verify the final adjacency list and mutation costs
    print(f"\nCheckpoint 4 (position {pos}): Final output formatting")
    for line in output:
        print(line)

    return output

# SECTION 5: Main function to integrate all steps across DNA positions
def small_parsimony_problem(data):
    tree, dna_strings, n = parse_input(data)
    
    # Process each position in the DNA strings independently
    full_output = []
    for pos in range(len(next(iter(dna_strings.values())))):  # DNA length
        internal_labels, parsimony_score = small_parsimony(tree, dna_strings, n, pos)
        result = format_output(tree, internal_labels, dna_strings, parsimony_score, n, pos)
        full_output.extend(result)
    
    return full_output

# Example input data
data = [
    "4",
    "4->CAAATCCC",
    "4->ATTGCGAC",
    "5->CTGCGCTG",
    "5->ATGGACGA",
    "6->4",
    "6->5"
]

# Solve the small parsimony problem
output = small_parsimony_problem(data)

# Output the result
for line in output:
    print(line)