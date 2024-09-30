from collections import defaultdict

# Parse the input into a tree structure
def parse_input(data):
    tree = defaultdict(list)
    dna_strings = {}
    n = int(data[0])
    
    # Fill the tree and DNA string dictionary
    for line in data[1:]:
        parent, child = line.split("->")
        parent, child = int(parent), child.strip()
        
        if child.isdigit():
            tree[parent].append(int(child))
        else:
            tree[parent].append(len(dna_strings))
            dna_strings[len(dna_strings)] = child

    return tree, dna_strings, n

# Fitch's algorithm implementation for Small Parsimony
def fitch(tree, dna_strings, n):
    # DP table to store sets of possible characters for each node
    internal_labels = {}
    parsimony_score = 0
    sequence_length = len(next(iter(dna_strings.values())))

    # Initialize internal labels for all nodes
    for node in range(n, max(tree) + 1):
        internal_labels[node] = [''] * sequence_length

    # Function to compute the parsimony for a single site
    def fitch_pass(node):
        if node in dna_strings:
            return list(dna_strings[node])

        children = tree[node]
        left = fitch_pass(children[0])
        right = fitch_pass(children[1])
        
        result = []
        for l, r in zip(left, right):
            if l == r:
                result.append(l)
            else:
                result.append(min(l, r))  # Arbitrarily choose one in case of mismatch
        return result

    # Post-order traversal to fill in the internal node labels
    for node in range(n, max(tree) + 1):
        internal_labels[node] = fitch_pass(node)
    
    # Calculate total parsimony score
    def calc_parsimony(node):
        nonlocal parsimony_score
        if node not in tree:
            return dna_strings[node]
        
        left, right = tree[node]
        left_seq = calc_parsimony(left)
        right_seq = calc_parsimony(right)

        for i in range(sequence_length):
            if left_seq[i] != internal_labels[node][i]:
                parsimony_score += 1
            if right_seq[i] != internal_labels[node][i]:
                parsimony_score += 1
        
        return internal_labels[node]
    
    calc_parsimony(max(tree))  # Traverse from the root
    return internal_labels, parsimony_score

# Format the output adjacency list with mutation costs
def format_output(tree, internal_labels, dna_strings, parsimony_score, n):
    output = [str(parsimony_score)]
    
    # Collect the adjacency list with mutation costs
    for parent in tree:
        for child in tree[parent]:
            parent_label = internal_labels.get(parent, dna_strings[parent])
            child_label = internal_labels.get(child, dna_strings.get(child))
            mutations = sum(1 for a, b in zip(parent_label, child_label) if a != b)
            output.append(f"{parent_label}->{child_label}:{mutations}")
            output.append(f"{child_label}->{parent_label}:{mutations}")
    
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
