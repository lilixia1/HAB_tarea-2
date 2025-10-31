'''
Escribe un script en Python que implemente un ejemplo de propagación en redes utilizando
algún algoritmo de GUILD y/o DIAMOnD. Usa los genes ENO1, PGK1 y HK2 como semillas.
Se proporcionan archivos de red formateada para ambos algoritmos, así como una red
filtrada obtenida de STRING y un script de ejemplo para procesarla.
'''

fichero = "data/network_diamond.txt"
genes_semilla = ['ENO1', 'PGK1', 'HK2']
nodos_añadidos = 5

## 0. Librerias
import pandas as pd
import os
import requests
from io import StringIO

## 1. Carga y preparacion de datos
def cargar_datos(fichero):
    try:
        interactions = pd.read_csv(fichero, sep=",")
        interactions.columns = [
            "protein1_hugo", "protein2_hugo"
        ]
        return interactions
    except Exception as e:
        print("Ha habido un error al cargar los datos: ", e)
    
def hugo_to_string_ids(genes, species=9606):
    """
    Convierte una lista de nombres de genes HUGO a IDs de STRING (e.g. 9606.ENSP...).
    Usa la API de STRING.
    """
    url = "https://string-db.org/api/tsv/get_string_ids"
    params = {
        "identifiers": "%0d".join(genes),  # separador URL
        "species": species,
        "limit": 1,
    }

    response = requests.post(url, data=params)
    if response.status_code != 200:
        raise RuntimeError(f"Error en la solicitud a STRING: {response.status_code}")

    df = pd.read_csv(StringIO(response.text), sep="\t")
    mapping = dict(zip(df["queryItem"], df["stringId"]))
    return mapping

## 2. Visualizar red
def representar_red(grafo):
    return imagen

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
    mapping = hugo_to_string_ids(genes_semilla)
    print(mapping)
    
    grafo = cargar_datos(DATA_DIR)

    ##DIAMOND
    #diamond_algorithm(grafo, genes_semilla, nodos_añadidos)
    

main()