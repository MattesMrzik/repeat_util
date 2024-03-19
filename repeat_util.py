import matplotlib.pyplot as plt
import networkx as nx
from Bio.Seq import Seq

aa_full_names = {
    'A': 'Alanine',
    'R': 'Arginine',
    'N': 'Asparagine',
    'D': 'Aspartic Acid',
    'C': 'Cysteine',
    'E': 'Glutamic Acid',
    'Q': 'Glutamine',
    'G': 'Glycine',
    'H': 'Histidine',
    'I': 'Isoleucine',
    'L': 'Leucine',
    'K': 'Lysine',
    'M': 'Methionine',
    'F': 'Phenylalanine',
    'P': 'Proline',
    'S': 'Serine',
    'T': 'Threonine',
    'W': 'Tryptophan',
    'Y': 'Tyrosine',
    'V': 'Valine',
    'U': 'Selenocysteine',  # Special case for selenocysteine
    'O': 'Pyrrolysine',  # Special case for pyrrolysine
}


def codon_to_aa(codon):
    seq = Seq(codon)
    aa = seq.translate()
    Seq.translate(codon)

    return aa


def shift_sequence(seq, shift):
    seq = Seq(seq)
    shifted_seq = seq[shift:] + seq[:shift]

    return shifted_seq


def all_shifts_and_complements(codon):
    for i in range(3):
        shifted_codon = shift_sequence(codon, i)
        amino_acid = codon_to_aa(shifted_codon)
        # print also the reverse complement
        shifted_codon_rc = shift_sequence(codon, i).reverse_complement()
        amino_acid_rc = codon_to_aa(shifted_codon_rc)
        # print also the full name of the amino acid

        print(
            f"Shifted codon: {shifted_codon}, amino acid: {amino_acid} (={aa_full_names[amino_acid]}), reverse complement: "
            f"{shifted_codon_rc}, amino acid: {amino_acid_rc} (={aa_full_names[amino_acid_rc]})")


def find_atomic_pattern(pattern):
    """The notion of atomic pattern is used in the case where the
    DNA strand of the ETR is not known. For a given pattern p,
    the atomic pattern of p is the first pattern (in lexicographic
    order) obtained from p by circular permutations on p and
    on its reverse complement. Thus, several patterns share the
    same atomic pattern. See doi:10.1109/BIBM.2014.6999134
    """

    atomic_patterns = []
    rc = pattern.reverse_complement()
    for i in range(len(pattern)):
        shifted_pattern = shift_sequence(pattern, i)
        shifted_pattern_rc = shift_sequence(rc, i)
        atomic_patterns.append(min(shifted_pattern, shifted_pattern_rc))

    return min(atomic_patterns)


def de_bruijn_graph(k, sequence, with_reverse_complement = False):
    g = nx.DiGraph()
    rc = Seq(sequence).reverse_complement()
    for i in range(len(sequence) - k):
        node_from = sequence[i:i + k]
        node_to = sequence[i + 1:i + k + 1]
        rc_node_from = rc[i:i + k]
        rc_node_to = rc[i + 1:i + k + 1]
        g.add_edge(rc_node_from, rc_node_to)
        g.add_edge(node_from, node_to)
    return g


def plot_de_bruijn_graph(G):
    pos = nx.planar_layout(G)
    nx.draw(G, pos, with_labels = True, font_weight = 'bold', node_size = 700, node_color = 'skyblue', font_size = 8)
    plt.show()


def show_de_bruijn_graph(k, sequence, with_reverse_complement = False):
    g = de_bruijn_graph(k, sequence, with_reverse_complement)
    plot_de_bruijn_graph(g)


if __name__ == "__main__":
    # all_shifts_and_complements("CAG")
    # print()
    # all_shifts_and_complements("CTG")

    # print(find_atomic_pattern(Seq("CAG")))
    show_de_bruijn_graph(4, "CAG"*4+"CTG"*3+"CAG"*3)