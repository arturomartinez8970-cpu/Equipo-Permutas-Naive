from Bio import Entrez, SeqIO
import time
import matplotlib.pyplot as plt

# CONFIGURACIÓN DE NCBI
Entrez.email = "arturo.morales7406@alumnos.udg.mx"


# FUNCIÓN PARA DESCARGAR SECUENCIA DE ADN DESDE NCBI

def obtener_secuencia_ncbi(id_ncbi):
	"""Descarga una secuencia de ADN desde NCBI y devuelve la secuencia como string."""
	print(f"Descargando secuencia {id_ncbi} desde NCBI...")
	handle = Entrez.efetch(db="nucleotide", id=id_ncbi, rettype="fasta", retmode="text")
	registro = SeqIO.read(handle, "fasta")
	handle.close()
	sec = str(registro.seq)
	print(f"Secuencia descargada: {len(sec)} bases.")
	return sec


# ALGORITMO NAIVE (Fuerza Bruta)


def naive_search(text, pattern):
	"""Fuerza Bruta: devuelve lista de posiciones donde aparece `pattern` en `text`."""
	n, m = len(text), len(pattern)
	matches = []
	for i in range(n - m + 1):
		if text[i:i + m] == pattern:
			matches.append(i)
	return matches

# MEDICIÓN Y GRÁFICA PARA NAIVE

def medir_tiempo_naive(texto, patron):
	inicio = time.perf_counter()
	naive_search(texto, patron)
	return time.perf_counter() - inicio


def comparar_naive(texto, patron):
	tamaños = [max(100, len(texto) // 20), max(100, len(texto) // 10), max(100, len(texto) // 5), max(100, len(texto) // 2), len(texto)]
	tiempos_naive = []

	for t in tamaños:
		subtexto = texto[:t]
		tiempos_naive.append(medir_tiempo_naive(subtexto, patron))

	plt.plot(tamaños, tiempos_naive, label="Naive", marker='o')
	plt.xlabel("Tamaño del texto (bases)")
	plt.ylabel("Tiempo (segundos)")
	plt.title("Rendimiento: Algoritmo Naive (Fuerza Bruta)")
	plt.legend()
	plt.grid(True)
	plt.show()


# MAIN

if __name__ == "__main__":
	id_secuencia = "NC_045512.2"  # SARS-CoV-2
	secuencia = obtener_secuencia_ncbi(id_secuencia)

	# Ejemplo de patrón a buscar
	patron = "ATGTTTGTTT"

	print("\n=== Búsqueda Naive ===")
	resultados_naive = naive_search(secuencia, patron)
	print(f"Coincidencias encontradas: {len(resultados_naive)}\nPosiciones: {resultados_naive}")

	# Comparativa (solo Naive)
	comparar_naive(secuencia, patron)


