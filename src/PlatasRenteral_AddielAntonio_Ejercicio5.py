from Bio.Seq import Seq

def find_complete_orfs(dna_sequence):
    """Encuentra ORFs completos que empiezan con ATG y terminan en un codon de paro, incluyendo '*' al final"""
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    proteins = []

    def extract_orfs(seq):
        seq_len = len(seq)
        orfs = []
        for frame in range(3):
            i = frame
            while i < seq_len - 2:
                codon = seq[i:i+3]
                if codon == start_codon:
                    protein_seq = ""
                    j = i
                    found_stop = False
                    while j < seq_len - 2:
                        current_codon = seq[j:j+3]
                        if current_codon in stop_codons:
                            found_stop = True
                            break
                        aa = Seq(current_codon).translate()
                        protein_seq += str(aa)
                        j += 3
                    # Solo agregamos ORF si encontramos un codon de paro, y agregamos '*'
                    if protein_seq and protein_seq[0] == "M" and found_stop:
                        protein_seq += "*"  # Marcar codon de stop
                        orfs.append(protein_seq)
                    i = j + 3  # Continuar despues del codon stop
                else:
                    i += 3
        return orfs

    # Hebra principal
    proteins.extend(extract_orfs(dna_sequence))
    # Hebra Reverso complementaria
    rev_comp_sequence = str(Seq(dna_sequence).reverse_complement())
    proteins.extend(extract_orfs(rev_comp_sequence))

    return proteins

if __name__ == "__main__":
    dna_sequence = "AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG"

    proteins = find_complete_orfs(dna_sequence)

    print("===== ORFs completos (empiezan con M y terminan con un codon de paro *) =====")
    if proteins:
        for i, p in enumerate(proteins, 1):
            print(f"Proteina {i}: {p}")
        longest = max(proteins, key=len)
        print("\nORF mas largo encontrado:", longest)
    else:
        print("No se encontraron ORFs...")