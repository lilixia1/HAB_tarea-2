import os
import argparse
import pandas as pd
from mygene import MyGeneInfo
from tqdm import tqdm
import networkx as nx
import scipy.stats
import matplotlib.pyplot as plt
import numpy as np
import scipy.special
import requests
import operator
import importlib
import subprocess
import sys

# --- Parámetros de la Tarea ---
nodos_añadidos = 50
especie = 'human'

# Variable global para compartir los Entrez ID de las semillas
genes_semilla_entrez = []

# ----------------------------------------------------------------------
#                         FUNCIONES DE LECTURA Y CONVERSIÓN
# ----------------------------------------------------------------------
def importar_genes(file_path: str) -> list[str]:
    """
    Carga una lista de genes desde un archivo donde los genes están separados por coma (e.g., COX4I2, ND1, ATP6).

    Args:
        file_path (str): Ruta al archivo de texto.

    Returns:
        list[str]: Lista de símbolos de genes limpios y únicos.
    """
    print(f"\nCargando genes desde: {file_path}")
    
    if not os.path.exists(file_path):
        print(f"[ERROR] Archivo no encontrado: {file_path}")
        return []
        
    try:
        with open(file_path, 'r') as f:
            # 1. Leer todo el contenido (asume que los genes están en una o pocas líneas)
            content = f.read() 
            
            # 2. Reemplazar saltos de línea, dividir por comas y limpiar espacios
            raw_genes = content.replace('\n', ',').split(',')
            
            # 3. Limpiar espacios y filtrar entradas vacías
            lista_genes = [gen.strip() for gen in raw_genes if gen.strip()]
            
            # 4. Eliminar duplicados manteniendo el orden
            lista_genes = list(dict.fromkeys(lista_genes))

    except Exception as e:
        print(f"[ERROR] Error al leer el archivo {file_path}: {e}")
        return []
    
    print(f"Genes cargados ({len(lista_genes)}): {lista_genes}")

    return lista_genes

def cargar_datos(fichero):
    """Carga las interacciones de un archivo TSV y nombra las columnas como HUGO."""
    try:
        interactions = pd.read_csv(fichero, sep="\t")
        interactions.columns = [
            "protein1_hugo", "protein2_hugo", "combined_score"
        ]
        return interactions
    except Exception as e:
        print(f"Error al cargar los datos de {fichero}: {e}")
        return None

def convertir_hugo_a_entrez(genes_hugo, especie):
    """Mapea una lista de genes HUGO a sus respectivos Entrez ID."""
    mg = MyGeneInfo()
    res = mg.querymany(genes_hugo, scopes='symbol', fields='entrezgene', species=especie, as_dataframe=True, verbose=False)
    
    hugo_to_entrez = res[res['entrezgene'].notna()]['entrezgene'].astype(str).to_dict()
    final_mapping = {k: v[0] if isinstance(v, list) else str(v) for k, v in hugo_to_entrez.items()}
    
    return final_mapping

def preparar_red(interactions_df, entrez_map):
    """Añade columnas de Entrez ID al DataFrame y lo limpia."""
    print("Mapeando todos los genes de interacciones a Entrez ID...")
    
    interactions_df["protein1_entrez"] = interactions_df["protein1_hugo"].map(entrez_map)
    interactions_df["protein2_entrez"] = interactions_df["protein2_hugo"].map(entrez_map)
    
    interactions_df.dropna(subset=["protein1_entrez", "protein2_entrez"], inplace=True)
    
    return interactions_df

def construir_red(interactions_df, SEED_GENES_ENTREZ):
    """Construye un grafo NetworkX usando los Entrez ID y el score combinado como peso."""
    G = nx.Graph()
    for _, row in interactions_df.iterrows():
        # Asumiendo que 'combined_score' ya está entre 0 y 1000, lo dividimos por 1000.
        G.add_edge(row["protein1_entrez"], row["protein2_entrez"], weight=row["combined_score"] / 1000.0)

    seed_interactions = interactions_df[
        interactions_df["protein1_entrez"].isin(SEED_GENES_ENTREZ) | 
        interactions_df["protein2_entrez"].isin(SEED_GENES_ENTREZ)
    ]
    print(f"Interacciones encontradas para genes semilla: {seed_interactions.shape[0]}")
    if seed_interactions.empty:
        print("Advertencia: No se encontraron interacciones para los genes semilla.")
    return G

# ----------------------------------------------------------------------
#                         FUNCIONES DIAMOnD
# ----------------------------------------------------------------------

def compute_all_gamma_ln(N):
    gamma_ln = {i: scipy.special.gammaln(i) for i in range(1, N + 1)}
    return gamma_ln

def gauss_hypergeom(x, r, b, n, gamma_ln):
    max_index = len(gamma_ln) - 1
    if r + b > max_index or x < 0 or r < x or b < (n - x) or n < x:
        return 0 

    log_p = (gamma_ln[r] - (gamma_ln[x] + gamma_ln[r - x]) +
             gamma_ln[b] - (gamma_ln[n - x] + gamma_ln[b - n]) -
             gamma_ln[r + b])
             
    return np.exp(log_p)

def pvalue(kb, k, N, s, gamma_ln):
    """Calcula el p-valor hipergeométrico."""
    return sum(gauss_hypergeom(n, s, N - s, k, gamma_ln) for n in range(kb, k + 1))

def diamond_iteration_of_first_X_nodes(G, S, X, alpha=1):
    added_nodes = []
    neighbors = {node: set(G.neighbors(node)) for node in G.nodes}
    degrees = dict(G.degree())
    cluster_nodes = set(S)
    
    if len(G.nodes) == 0:
        return []

    gamma_ln = compute_all_gamma_ln(len(G.nodes))

    while len(added_nodes) < X:
        min_p = float('inf')
        next_node = None
        
        candidate_nodes = set()
        for seed in cluster_nodes:
            candidate_nodes.update(neighbors.get(seed, set()))
        candidate_nodes -= cluster_nodes

        if not candidate_nodes:
            break

        for node in tqdm(candidate_nodes, desc=f"DIAMOnD Iteración {len(added_nodes) + 1}/{X}"):
            k = degrees.get(node, 0)
            kb = sum((1 for neighbor in neighbors[node] if neighbor in cluster_nodes))

            if k == 0 or k < kb: 
                continue

            try:
                p = pvalue(kb, k, len(G.nodes), len(cluster_nodes), gamma_ln)
            except Exception:
                continue

            if p < min_p:
                min_p = p
                next_node = node

        if next_node:
            added_nodes.append(next_node)
            cluster_nodes.add(next_node)
        else:
            break
            
    return added_nodes

# ----------------------------------------------------------------------
#                         FUNCIÓN GUILD
# ----------------------------------------------------------------------

def guild_propagation(G, S, X, alpha=0.8):
    """
    Implementa un algoritmo de difusión de información (Random Walk with Restart 
    simplificado) similar a GUILD.
    """
    
    print(f"Iniciando propagación GUILD (alpha={alpha})...")
    
    # Inicialización de scores: 1 para semillas, 0 para el resto
    current_scores = {node: 0.0 for node in G.nodes}
    for seed in S:
        if seed in G.nodes:
            current_scores[seed] = 1.0
        else:
            print(f"⚠️ Semilla {seed} no encontrada en el grafo. Omitiendo.")

    if not any(current_scores.values()):
        return []
        
    f0 = current_scores.copy()
    
    # Iteración de difusión 
    for _ in tqdm(range(20), desc="GUILD Iteraciones de Difusión"): 
        new_scores = {}
        for node in G.nodes:
            
            neighbor_sum = 0.0
            for neighbor, data in G[node].items():
                weight = data.get('weight', 1.0)
                
                # Grado ponderado del vecino
                neighbor_weighted_degree = sum(G[neighbor][n].get('weight', 1.0) for n in G[neighbor])
                
                if neighbor_weighted_degree > 0:
                    neighbor_sum += current_scores[neighbor] * (weight / neighbor_weighted_degree)
            
            new_scores[node] = alpha * neighbor_sum + (1 - alpha) * f0[node]

        current_scores = new_scores

    # Ordenar y seleccionar los X mejores nodos (excluyendo las semillas)
    ranked_nodes = sorted(current_scores.items(), key=operator.itemgetter(1), reverse=True)
    
    added_nodes = []
    seed_set = set(S)
    
    for node, score in ranked_nodes:
        if node not in seed_set:
            added_nodes.append(node)
            if len(added_nodes) >= X:
                break
                
    return added_nodes

# ----------------------------------------------------------------------
#                         FUNCIONES DE SALIDA Y MAIN
# ----------------------------------------------------------------------

def guardar_resultados(diamond_genes, guild_genes, archivo_salida):
    """Guarda los genes semilla y los genes añadidos por ambos métodos en un archivo TSV."""
    
    entrez_ids = list(set(genes_semilla_entrez + diamond_genes + guild_genes))
    
    # Mapeo inverso de Entrez ID a HUGO
    mg = MyGeneInfo()
    res = mg.querymany(entrez_ids, scopes='entrezgene', fields='symbol', species=especie, as_dataframe=True, verbose=False)
    entrez_to_hugo = res[res['symbol'].notna()]['symbol'].to_dict()
    
    # Crear el DataFrame de resultados
    data = []
    seed_set = set(genes_semilla_entrez)
    diamond_set = set(diamond_genes)
    guild_set = set(guild_genes)
    
    for entrez_id in entrez_ids:
        hugo_symbol = entrez_to_hugo.get(entrez_id, 'N/A')
        
        tipo = []
        if entrez_id in seed_set:
            tipo.append('Seed_Gene')
        if entrez_id in diamond_set:
            tipo.append('DIAMOnD_Added')
        if entrez_id in guild_set:
            tipo.append('GUILD_Added')
            
        data.append({
            'Entrez_ID': entrez_id, 
            'HUGO_Symbol': hugo_symbol, 
            'Tipo': '|'.join(tipo)
        })
        
    results_df = pd.DataFrame(data)
    results_df.to_csv(archivo_salida, sep="\t", index=False)
    print(f"\nResultados de ambos algoritmos guardados en: {archivo_salida}")


def graficar_red_enriquecida(G, seed_genes, diamond_genes, guild_genes, output_image_file):
    """Genera la visualización de la red enriquecida comparando ambos resultados."""
    
    all_nodes = set(seed_genes) | set(diamond_genes) | set(guild_genes)
    subgraph = G.subgraph(all_nodes)
    
    if len(subgraph.nodes) < 2:
        print("Advertencia: Grafo enriquecido demasiado pequeño para dibujar.")
        return

    pos = nx.spring_layout(subgraph, k=0.15, iterations=20) 
    plt.figure(figsize=(15, 15))

    # Nodos de semillas (Azul)
    nx.draw_networkx_nodes(subgraph, pos, nodelist=seed_genes, node_color='blue', node_size=800, label="Seed Genes", alpha=0.8)

    # Nodos solo DIAMOnD (Naranja)
    diamond_only = list(set(diamond_genes) - set(guild_genes))
    nx.draw_networkx_nodes(subgraph, pos, nodelist=diamond_only, node_color='orange', node_size=500, label="DIAMOnD Only", alpha=0.7)
    
    # Nodos solo GUILD (Verde)
    guild_only = list(set(guild_genes) - set(diamond_genes))
    nx.draw_networkx_nodes(subgraph, pos, nodelist=guild_only, node_color='green', node_size=500, label="GUILD Only", alpha=0.7)
    
    # Nodos Comunes (Rojo)
    common_nodes = list(set(diamond_genes) & set(guild_genes))
    nx.draw_networkx_nodes(subgraph, pos, nodelist=common_nodes, node_color='red', node_size=600, label="Common (DIAMOnD & GUILD)", alpha=0.9)


    # Enlaces (Edges)
    nx.draw_networkx_edges(subgraph, pos, alpha=0.4, edge_color='gray')
    
    # Etiquetas (Intentar usar HUGO)
    node_labels = {n: n for n in subgraph.nodes()}
    try:
        mg = MyGeneInfo()
        res = mg.querymany(list(subgraph.nodes()), scopes='entrezgene', fields='symbol', species=especie, as_dataframe=True, verbose=False)
        entrez_to_hugo = res[res['symbol'].notna()]['symbol'].to_dict()
        node_labels = {k: entrez_to_hugo.get(k, k) for k in subgraph.nodes()}
    except Exception:
        pass 
        
    nx.draw_networkx_labels(subgraph, pos, labels=node_labels, font_size=8, font_weight='bold')

    plt.legend(loc="upper left", markerscale=0.7)
    plt.title(f"Comparación de Redes Enriquecidas: DIAMOnD vs GUILD ({len(subgraph.nodes)} Nodos)")
    plt.axis('off')
    
    try:
        plt.savefig(output_image_file, format='png', bbox_inches='tight')
        print(f"Grafo comparativo guardado como imagen en: {output_image_file}")
    except Exception as e:
        print(f"Error al guardar la imagen: {e}")
        
    plt.close()

# ----------------------------------------------------------------------
#                             MAIN CLI
# ----------------------------------------------------------------------

def main():
    global genes_semilla_entrez
    
    # 1. Configuración de argparse para CLI
    parser = argparse.ArgumentParser(
        description="Propagación de genes en redes de interacción usando DIAMOnD y GUILD."
    )
    parser.add_argument(
        '--seed-file', 
        required=True, 
        help="Ruta al archivo de texto que contiene los genes semilla (separados por coma)."
    )
    parser.add_argument(
        '--input', 
        required=True, 
        help="Ruta al archivo de entrada de la red (e.g., string_network_filtered_hugo-400.tsv)"
    )
    parser.add_argument(
        '--output', 
        default='propagation_results.tsv', 
        help="Nombre del archivo de resultados TSV (ambos métodos)."
    )
    parser.add_argument(
        '--plot',
        default='propagation_network_comparison.png', 
        help="Nombre del archivo para la imagen comparativa de la red."
    )
    args = parser.parse_args()
    
    print(f"--- 🧬 Iniciando Propagación (DIAMOnD & GUILD) ---")
    
    # 2. Creación de la Carpeta de Resultados
    RESULTS_DIR = "results"
    try:
        if not os.path.exists(RESULTS_DIR):
            os.makedirs(RESULTS_DIR)
            print(f"Carpeta de resultados creada: {RESULTS_DIR}")
    except OSError as e:
        print(f"Error al crear el directorio {RESULTS_DIR}: {e}")
        return

    output_tsv_path = os.path.join(RESULTS_DIR, args.output)
    output_plot_path = os.path.join(RESULTS_DIR, args.plot)
    print(f"Archivo de entrada: {args.input}")
    print(f"Archivos de salida se guardarán en: {os.path.abspath(RESULTS_DIR)}")
    
    ## 3. CONVERSIÓN DE IDs Y PREPARACIÓN
    seed_file_path = args.seed_file
    genes_semilla = importar_genes(seed_file_path) # Llamar a la nueva función
    print("\nConvirtiendo genes semilla HUGO a Entrez ID...")
    seed_map = convertir_hugo_a_entrez(genes_semilla, especie)
    genes_semilla_entrez = list(seed_map.values())
    
    if not genes_semilla_entrez:
        print("ERROR: No se pudieron obtener los Entrez ID para los genes semilla. Abortando.")
        return
        
    interacciones = cargar_datos(args.input)
    if interacciones is None:
        return

    all_genes_hugo = pd.concat([interacciones['protein1_hugo'], interacciones['protein2_hugo']]).unique()
    hugo_to_entrez_map = convertir_hugo_a_entrez(all_genes_hugo, especie)
    
    interacciones_entrez = preparar_red(interacciones, hugo_to_entrez_map)
    
    ## 4. CONSTRUIR GRAFO
    red = construir_red(interacciones_entrez, genes_semilla_entrez)
    
    if len(red.nodes) == 0:
        print("El grafo está vacío. No se puede ejecutar la propagación.")
        return

    # Determinamos el número real de nodos a añadir
    n = min(nodos_añadidos, len(red.nodes) - len(genes_semilla_entrez))
    if n <= 0:
        print("No hay nodos para añadir o la red es muy pequeña.")
        return

    ## 5. EJECUTAR DIAMOnD
    print(f"\n--- Ejecutando DIAMOnD ---")
    diamond_genes = diamond_iteration_of_first_X_nodes(red, genes_semilla_entrez, n)

    ## 6. EJECUTAR GUILD
    print(f"\n--- Ejecutando GUILD ---")
    guild_genes = guild_propagation(red, genes_semilla_entrez, n)


    ## 7. GUARDAR Y GRAFICAR RESULTADOS
    guardar_resultados(diamond_genes, guild_genes, output_tsv_path)
    graficar_red_enriquecida(red, genes_semilla_entrez, diamond_genes, guild_genes, output_plot_path)

if __name__ == '__main__':
    main()