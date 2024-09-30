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

# Fitch's algorithm implementation for Small Parsimony
def fitch(tree, dna_strings, n):
    # DP table to store the character states for each node
    internal_labels = {}
    parsimony_score = 0
    sequence_length = len(next(iter(dna_strings.values())))

    # Initialize internal labels for all internal nodes
    for node in range(n, max(tree) + 1):
        internal_labels[node] = [''] * sequence_length

    # Function to compute the most parsimonious label for a single site
    def fitch_pass(node):
        if node in dna_strings:
            return list(dna_strings[node])

        children = tree[node]
        left_label = fitch_pass(children[0])
        right_label = fitch_pass(children[1])
        
        result = []
        for l, r in zip(left_label, right_label):
            if l == r:
                result.append(l)
            else:
                result.append(min(l, r))  # Arbitrarily choose one of them
        return result

    # Post-order traversal to fill in internal node labels
    for node in range(n, max(tree) + 1):
        internal_labels[node] = fitch_pass(node)
    
    # Calculate total parsimony score by counting differences between nodes
    def calc_parsimony(node):
        nonlocal parsimony_score
        if node not in tree:
            return dna_strings[node]
        
        left, right = tree[node]
        left_seq = calc_parsimony(left)
        right_seq = calc_parsimony(right)
        
        node_seq = internal_labels[node]
        
        for i in range(sequence_length):
            if left_seq[i] != node_seq[i]:
                parsimony_score += 1
            if right_seq[i] != node_seq[i]:
                parsimony_score += 1
        
        return node_seq
    
    calc_parsimony(max(tree))  # Traverse starting from the root
    return internal_labels, parsimony_score

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
    internal_labels, parsimony_score = fitch(tree, dna_strings, n)
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
