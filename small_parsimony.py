from collections import defaultdict

# Function to parse input dynamically
def parse_input():
    # Sample input
    data = """
    4
    1->GATTACA
    2->CTTGATT
    3->ATCTGGA
    4->GTCGATG
    5->1
    5->2
    6->3
    6->4
    7->5
    7->6
    """.strip().splitlines()

    n = int(data[0].strip())  # Number of leaf nodes
    tree = defaultdict(list)
    sequences = {}

    for line in data[1:]:
        line = line.strip()
        if not line:
            continue  # Skip any empty lines
        parts = line.split("->")
        node = int(parts[0].strip())
        if parts[1].isnumeric():
            child = int(parts[1].strip())
            tree[node].append(child)
        else:
            sequences[node] = parts[1].strip()  # Leaf node sequence
    
    return n, tree, sequences

# Find the root node by identifying the node that is not a child of any other node
def find_root(tree):
    all_nodes = set(tree.keys())
    child_nodes = {child for children in tree.values() for child in children}
    root_nodes = all_nodes - child_nodes
    return root_nodes.pop()  # Return the only node that isn't a child of another

# Fitch's Algorithm Step 1: Bottom-up calculation of sets and score
def fitch_bottom_up(node, tree, sequences, sets, score, seq_len):
    if node in sequences:  # It's a leaf node
        sets[node] = [{nuc} for nuc in sequences[node]]
        return

    for child in tree[node]:
        fitch_bottom_up(child, tree, sequences, sets, score, seq_len)

    node_sets = []
    for i in range(seq_len):
        left_set = sets[tree[node][0]][i]
        right_set = sets[tree[node][1]][i]
        intersection = left_set & right_set
        if intersection:
            node_sets.append(intersection)
        else:
            node_sets.append(left_set | right_set)
            score[0] += 1  # Add 1 to the score when there's no intersection
    sets[node] = node_sets

# Fitch's Algorithm Step 2: Top-down backtracking to assign sequences
def fitch_top_down(node, parent_seq, sets, sequences, tree):
    if node not in sequences:
        seq = []
        for i, s in enumerate(sets[node]):
            if parent_seq and parent_seq[i] in s:
                seq.append(parent_seq[i])
            else:
                seq.append(next(iter(s)))
        sequences[node] = ''.join(seq)

    for child in tree[node]:
        fitch_top_down(child, sequences[node], sets, sequences, tree)

# Calculate mutations between sequences
def calculate_mutations(seq1, seq2):
    return sum(1 for a, b in zip(seq1, seq2) if a != b)

# Main function to solve the Small Parsimony problem
def small_parsimony():
    # Parse hardcoded input
    n, tree, sequences = parse_input()
    
    # Sequence length (assuming all sequences are of equal length)
    seq_len = len(next(iter(sequences.values())))
    
    # Find the root dynamically
    root = find_root(tree)
    
    # Step 1: Bottom-up phase using Fitch's algorithm
    sets = {}
    score = [0]  # Parsimony score as a list to allow modification in recursion
    fitch_bottom_up(root, tree, sequences, sets, score, seq_len)
    
    # Step 2: Top-down phase to reconstruct internal sequences
    fitch_top_down(root, None, sets, sequences, tree)
    
    # Output the total parsimony score
    print(score[0])
    
    # Step 3: Output the adjacency list with mutation counts
    adjacency_list = set()  # Use a set to avoid duplicate entries
    for node in tree:
        for child in tree[node]:
            seq1 = sequences[node]
            seq2 = sequences[child]
            mutations = calculate_mutations(seq1, seq2)
            # Avoid self-loops
            if seq1 != seq2:
                adjacency_list.add((seq1, seq2, mutations))
    
    # Print the adjacency list without duplicates
    for seq1, seq2, mutations in adjacency_list:
        print(f"{seq1}->{seq2}:{mutations}")

# Run the small parsimony solver
if __name__ == "__main__":
    small_parsimony()
