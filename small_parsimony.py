from collections import defaultdict

# Parse the input into a tree structure
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

    return tree, dna_strings, n

# Compute the parsimony score for a given tree
def small_parsimony(tree, dna_strings, n):
    # Set of symbols (nucleotides in this case)
    alphabet = ['A', 'C', 'G', 'T']
    num_symbols = len(alphabet)
    
    # Initialize the scoring dictionary for all nodes
    s = defaultdict(lambda: [float('inf')] * num_symbols)  # Each node gets a list for scores for A, C, G, T
    tags = defaultdict(int)  # Tag each node (0 for unprocessed, 1 for processed)
    
    # Initialize scores for the leaves
    for leaf in dna_strings:
        tags[leaf] = 1  # Mark leaf as processed
        leaf_string = dna_strings[leaf]
        for i, symbol in enumerate(alphabet):
            if leaf_string == symbol:
                s[leaf][i] = 0  # Score is 0 if the symbol matches the leaf character
            else:
                s[leaf][i] = float('inf')  # Otherwise, it's infinity

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
        for k in range(num_symbols):
            min_left = float('inf')
            min_right = float('inf')

            for i in range(num_symbols):
                min_left = min(min_left, s[left][i] + delta(alphabet[i], alphabet[k]))
            for j in range(num_symbols):
                min_right = min(min_right, s[right][j] + delta(alphabet[j], alphabet[k]))
            
            # Update the score for character k at node v
            s[v][k] = min_left + min_right

    # Root is the last internal node processed
    root = max(tree)
    parsimony_score = min(s[root])  # Minimum score over all characters at the root

    return s, parsimony_score

# Format the output adjacency list with mutation costs
def format_output(tree, internal_labels, dna_strings, parsimony_score, n):
    output = [str(parsimony_score)]
    
    # Helper function to convert list to string
    def label_to_string(label):
        if isinstance(label, list):
            return ''.join(label)  # Join list of characters into a string
        return label
    
    # Collect the adjacency list with mutation costs
    for parent in tree:
        for child in tree[parent]:
            parent_label = label_to_string(internal_labels.get(parent, dna_strings.get(parent, "")))
            child_label = label_to_string(internal_labels.get(child, dna_strings.get(child, "")))
            mutations = sum(1 for a, b in zip(parent_label, child_label) if a != b)
            output.append(f"{parent_label}->{child_label}:{mutations}")
    
    return output

# Main function to integrate all steps
def small_parsimony_problem(data):
    tree, dna_strings, n = parse_input(data)
    internal_labels, parsimony_score = small_parsimony(tree, dna_strings, n)
    result = format_output(tree, internal_labels, dna_strings, parsimony_score, n)
    return result

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