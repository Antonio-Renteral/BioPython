class Gene:
    # Atributos de la clase
    nucleotidos_validos = {'A', 'T', 'C', 'G'}

    # Constructor
    def __init__(self, nombre, inicio, final, secuencia=""):
        self.nombre = nombre # Nombre del gene
        self.inicio = inicio # Donde inicia el gene
        self.final = final # Donde termina el gene
        self.secuencia = secuencia.upper()

    def longitud(self): # Devuelve la longitud de la secuencia dada
        return len(self.secuencia)
    
    def contenido_gc(self): # Calculamos el contenido de gc
        conteo_gc = self.secuencia.count('G') + self.secuencia.count('C')
        return (conteo_gc / len(self.secuencia)) * 100 if self.secuencia else 0 # Se calcula el porcentaje, si la secuencia esta vacia regresa 0

    def transcribir(self): # Convierte ADN a ARN remplazando T por U
        return self.secuencia.replace('T', 'U')


class tRNA(Gene):
    def __init__(self, nombre, inicio, final, secuencia="", anticodon="", aminoacido=""):
        super().__init__(nombre, inicio, final, secuencia) # Pasamos los mismos datos del constructor de gene
        self.aminoacido = aminoacido # Guardamos el aminoacido que transporta
        self.anticodon = anticodon.upper() if anticodon else "" # Guardamos la secuencia del anticodon, si el usuario no lo da se queda vacia la cadena

    def calcular_anticodon(self): # En caso de que el usuario no nos de el anticodon podemos calcularlo, usando las reglas de apareamiento de bases
        complementos = {"A": "U", "U": "A", "G": "C", "C": "G"}
        return "".join(complementos.get(base, "?") for base in self.transcribir())

    def emparejamiento(self, codon): # Revisa que el anticodon se empareje correctamente con el codon dado
        codon = codon.upper()
        anticodon_a_usar = self.anticodon if self.anticodon else self.calcular_anticodon()
        return anticodon_a_usar == codon


class RNA_no_codificante(Gene):
    def __init__(self, nombre, inicio, final, secuencia="", funcion=""):
        super().__init__(nombre, inicio, final, secuencia) # Pasamos los mismos datos del constructor de gene
        self.funcion = funcion # Pedimos al usuario que nos diga que funcion tiene


class Proteina(tRNA):
    # Diccionario del codigo genetico estandar
    codigo_genetico = {
        "UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
        "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
        "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
        "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
        "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
        "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
        "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
        "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
        "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP", "UGA":"STOP",
        "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
        "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
        "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
        "UGU":"C", "UGC":"C", "UGG":"W",
        "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
        "AGA":"R", "AGG":"R",
        "AGU":"S", "AGC":"S",
        "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",
    }

    def __init__(self, nombre, inicio, final, secuencia="", anticodon=""):
        super().__init__(nombre, inicio, final, secuencia, anticodon) # Pasamos los mismos datos del constructor de gene
        self.secuencia_aa = self.traducir()  # Genera los aminoacidos al crear la proteina

    def traducir(self): # Toma de 3 en 3 (codon) busca el aa correspondiente, si encuentra stop se detiene, devuelve la secuencia de aa
        rna = self.transcribir()
        aminoacidos = []
        for i in range(0, len(rna), 3):
            codon = rna[i:i+3]
            if len(codon) == 3:
                aa = self.codigo_genetico.get(codon, "?")
                aminoacidos.append(aa)
                if aa == "STOP":  # se detiene en el stop
                    break
        return "".join(aminoacidos)

    def longitud(self): # Devolvemos el numero de nucleotidos y el numero de aminoacidos
        return {
            "nucleotidos": len(self.secuencia),
            "aminoacidos": len(self.secuencia_aa.replace("STOP", ""))  # no cuenta STOP como aa real
        }
    
if __name__ == "__main__":
    print("=== Resultados de todas las clases ===\n")

    # Resultados Gene
    gen = Gene("Gen1", 1, 12, "ATGCATGCATGC")
    print("Gene:", gen.nombre)
    print("Longitud:", gen.longitud())
    print("Contenido GC:", gen.contenido_gc())
    print("RNA transcrito:", gen.transcribir())
    print()

    # Resultados tRNA
    trna = tRNA("tRNA1", 1, 6, "ATGCGT", aminoacido="Met")
    print("tRNA:", trna.nombre)
    print("Anticodon calculado:", trna.calcular_anticodon())
    print("Emparejamiento con 'UAC':", trna.emparejamiento("UAC"))
    print()

    # Resultados RNA_no_codificante
    rna_nc = RNA_no_codificante("RNA1", 1, 8, "ATGCTAGC", funcion="Regulación")
    print("RNA_no_codificante:", rna_nc.nombre)
    print("Función declarada:", rna_nc.funcion)
    print()

    # Resultados Proteina
    prot = Proteina("Prot1", 1, 9, "ATGGTTTAA")
    print("Proteina:", prot.nombre)
    print("Secuencia nucleotidos:", prot.secuencia)
    print("Traducción a aminoácidos:", prot.secuencia_aa)
    print("Longitudes:", prot.longitud())