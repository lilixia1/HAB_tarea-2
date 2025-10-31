'''
Escribe un script en Python que implemente un ejemplo de propagación en redes utilizando
algún algoritmo de GUILD y/o DIAMOnD. Usa los genes ENO1, PGK1 y HK2 como semillas.
Se proporcionan archivos de red formateada para ambos algoritmos, así como una red
filtrada obtenida de STRING y un script de ejemplo para procesarla.
'''

genes_semilla = ['ENO1', 'PGK1', 'HK2']
nodos_añadidos = 5
especie = 'human'

## 0. Librerias
import pandas as pd
import os
import requests
from io import StringIO
from mygene import MyGeneInfo
from tqdm import tqdm
import networkx as nx

## 1. Carga y preparacion de datos
def cargar_datos(fichero):
    try:
        interactions = pd.read_csv(fichero, sep=",")
        interactions.columns = [
            "protein1_entrez", "protein2_entrez"
        ]
        return interactions
    except Exception as e:
        print("Ha habido un error al cargar los datos: ", e)
    
from mygene import MyGeneInfo
from tqdm import tqdm

def hugo_to_entrez(genes_list, species, batch_size=1000, show_progress=True):
    """
    Convierte una lista de nombres de genes HUGO a sus IDs Entrez usando MyGene.info.

    Args:
        genes_list (list[str]): Lista de nombres de genes (por ejemplo ['ENO1', 'PGK1', 'HK2'])
        species (str): Especie, por defecto 'human'
        batch_size (int): Tamaño de lote para las consultas a la API
        show_progress (bool): Si True, muestra barra de progreso (requiere tqdm)

    Returns:
        list[str]: Lista de IDs Entrez (None si algún gen no se encuentra)
    """
    mg = MyGeneInfo()
    entrez_mapping = {}

    iterator = range(0, len(genes_list), batch_size)
    if show_progress:
        iterator = tqdm(iterator, desc="Mapeando HUGO → Entrez")

    for i in iterator:
        batch = genes_list[i:i + batch_size]
        results = mg.querymany(
            batch,
            scopes='symbol',
            fields='entrezgene',
            species=species
        )

        for res in results:
            gene = res.get('query')
            entrez = res.get('entrezgene')
            if entrez:
                entrez_mapping[gene.upper()] = int(entrez)

    entrez_list = [entrez_mapping.get(g, None) for g in genes_list]
    return entrez_list


def lista_id(diccionario):
    lista=[]
    for x,y in diccionario.items():
        lista.append(y)
    return lista


## 2. Construir la red
def construir_red(interactions, SEED_GENES):
    G = nx.Graph()
    for _, row in interactions.iterrows():
        G.add_edge(row["protein1_entrez"], row["protein2_entrez"])

    print("Verificando interacciones para los genes semilla...")
    seed_interactions = interactions[interactions["protein1_entrez"].isin(SEED_GENES) | interactions["protein2_entrez"].isin(SEED_GENES)]
    print(f"Interacciones encontradas para genes semilla: {seed_interactions.shape[0]}")

    if seed_interactions.empty:
        print("No se encontraron interacciones para los genes semilla.")
    return G

## 3a. Propagacion red con DIAMOND

def diamond_algorithm(graph, seed_genes, max_added_nodes):
    nodes_added = []
    candidate_scores = {}

    # Verificación de la presencia de genes semilla en el grafo
    for seed in seed_genes:
        if seed not in graph:
            print(f"El gene semilla {seed} no está en la red.")
        else:
            print(f"El gene semilla {seed} está en la red.")

    for seed in seed_genes:
        if seed not in graph:
            continue
        neighbors = graph.neighbors(seed)
        for neighbor in neighbors:
            if neighbor in seed_genes:
                continue
            candidate_scores[neighbor] = candidate_scores.get(neighbor, 0) + 1

    for _ in range(max_added_nodes):
        if not candidate_scores:
            break
        best_candidate = max(candidate_scores, key=candidate_scores.get)
        nodes_added.append(best_candidate)
        
        for neighbor in graph.neighbors(best_candidate):
            if neighbor not in candidate_scores and neighbor not in seed_genes:
                candidate_scores[neighbor] = 0
            candidate_scores[neighbor] = candidate_scores.get(neighbor, 0) + 1
        
        del candidate_scores[best_candidate]
        
    return nodes_added
## 3b. Análisis funcional
## 4a. Propagación red con GUILD
## 4b. Análisis funcional
## 5. Discusión resultados

def main():
    ##RUTAS
    # Ruta absoluta del script actual
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    # Carpeta de datos
    DATA_DIR = os.path.join(BASE_DIR, "data", "network_diamond.txt")
    
    ##GENES SEMILLA A HUGO
    genes_semilla_entrez = hugo_to_entrez(genes_semilla, especie)
    print(genes_semilla_entrez)
    
    interacciones = cargar_datos(DATA_DIR)
    print(interacciones.head())

    ##Construir grafo
    red = construir_red(interacciones, genes_semilla_entrez)
    return
    ##DIAMOND
    grafo_diamond = diamond_algorithm(red, genes_semilla_entrez, nodos_añadidos)
    

main()