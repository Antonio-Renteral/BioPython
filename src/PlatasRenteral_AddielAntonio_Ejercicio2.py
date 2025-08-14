class Gene():
    # Atributos de clase
    nucleotidos_validos = {'A', 'T', 'C', 'G'}

    def __init__(self, nombre, inicio, final, secuencia=""):
        self.nombre = nombre
        self.inicio = inicio
        self.final = final
        self.secuencia = secuencia.upper()

    def longitud(self):
        return len(self.secuencia)
    
    def contenido_gc(self):
        conteo_gc = self.secuencia.count('G') + self.secuencia.count('C')
        return (conteo_gc / len(self.secuencia)) * 100 if self.secuencia else 0 # Verifica que haya algo en la secuencia

    def transcribir(self):
        return self.secuencia.replace('T', 'U')
    
# Probamos nuestros metodos
gene = Gene("LacZ", 34, 675, "ATGCGTACGTTAGC")

print("Longitud:", gene.longitud())
print(f"Contenido GC: {gene.contenido_gc():.2f}%")
print("Cadena de RNA:", gene.transcribir())