
# Tarea 2: Network Propagation

*Escribe un script en Python que implemente un ejemplo de propagación en redes utilizando algún algoritmo de **GUILD** y/o **DIAMOnD**. Usa los genes **ENO1**, **PGK1** y **HK2** como semillas. Se proporcionan archivos de red formateada para ambos algoritmos, así como una red filtrada obtenida de STRING y un script de ejemplo para procesarla.*

## 1. Requisitos

Este script requiere Python 3.x y las siguientes librerías, que deben instalarse antes de la ejecución.
```
pandas
numpy
mygene
requests
networkx
scipy
matplotlib
tqdm
```

## 2. Estructura del Proyecto

```
/network\_propagation/
├── data/
│   ├── network\_guild.txt                        # Red formateada para GUILD
│   ├── network\_diamond.txt                      # Red formateada para DIAMOnD
│   ├── string\_network\_filtered\_hugo-400.tsv     # Red filtrada de STRING
│   └── genes\_seed.txt                           # Genes semilla: ENO1, PGK1, HK2
├── scripts/
│   ├── process\_STRING.py                        # Script de ejemplo para procesar la red
│   └── tu\_script.py                             # Script que debe entregar el estudiante
├── results/                                     # Carpeta para resultados generados
├── README.md                                    # Este archivo
└── requirements.txt                             # Dependencias: networkx, pandas

```

## 3. Archivo de Entrada (seed.txt, string_network_filtered_hugo-400.tsv)

El archivo debe contener los símbolos de los genes, ya sea separados por línea o separados por coma. Aparte, se nos proporcionó una red de StringDB reducida, ya que sino sería computacionalmente costoso ejecutar todo el código con la totalidad de la base de datos.
```
ENO1, PGK1, HK2
```

## 4. Ejecución desde la Línea de Comandos (CLI)

**Comando Principal**
```{Bash}
python scripts/scriptDiamond.py --seed-file [RUTA_GENES] --input [RUTA_RED] [OPCIONALES]
```
**Ejemplo Práctico**
```{Bash}
python scripts/scriptDiamond.py \
    --seed-file data/genes_seed.txt \
    --input data/string_network_filtered_hugo-400.tsv \
    --output propagacion_comparativa.tsv \
    --plot red_final_comparacion.png
```

## 5. Archivos de Resultados Generados
El script generará dos archivos dentro de la carpeta results/:

1. *Archivo TSV de Resultados* (propagation_comparativa.tsv):
    - Contiene la lista completa de genes (semillas y añadidos).
    - La columna Tipo compara la pertenencia a los clusters de DIAMOnD y GUILD (ej., Seed_Gene, DIAMOnD_Added|GUILD_Added).

2. *Imagen PNG del Grafo* (red_final_comparacion.png):

    - Visualización del subgrafo enriquecido para la inspección visual de la propagación.

    - Codificación por Colores para mostrar la comparación de resultados:

        - Azul: Genes Semilla (ENO1, PGK1, HK2).

        - Rojo: Genes añadidos por DIAMOnD y GUILD (comunes).

        - Naranja: Genes añadidos solo por DIAMOnD.

        - Verde: Genes añadidos solo por GUILD.

## 6. Metodología de Análisis
El script implementa un pipeline de propagación de genes en la red de interacción proteína-proteína (PPI):

1. Carga y Mapeo de IDs: La red de STRING se carga. Todos los genes HUGO (Human Genome Organisation) son mapeados a su respectivo Entrez ID utilizando la librería mygene, que es el formato interno utilizado por los algoritmos.

2. Construcción del Grafo: Se utiliza networkx para construir un grafo no dirigido, usando el Entrez ID como nodo y el combined_score de STRING (normalizado) como peso de la arista.

3. Algoritmos de Propagación:

    - DIAMOnD: Es un método de enriquecimiento local basado en la distribución hipergeométrica. Mide si la conexión de un gen vecino a un cluster es estadísticamente más fuerte de lo esperado por azar (calculando un p-valor). Selecciona el nodo con el p-valor más bajo en cada iteración.

    - GUILD: Utiliza un algoritmo de difusión global (similar a Random Walk with Restart) para asignar un score de prioridad a cada nodo en función de su distancia ponderada y su conectividad al conjunto semilla. Selecciona los nodos con la puntuación más alta.

4. Selección y Comparación: Ambos algoritmos añaden los 50 genes más significativos/priorizados (definido por nodos_añadidos). Los resultados se comparan por intersección y diferencia de conjuntos (sets) para identificar las coincidencias.

## 7. Funciones implementadas en el código
| Categoría | Función | Descripción |
| :--- | :--- | :--- |
| **Carga/Limpieza** | `importar_genes(file_path)` | Lee los genes semilla desde el archivo `.txt` y los limpia. |
| | `cargar_datos(fichero)` | Carga el TSV de STRING. |
| | `convertir_hugo_a_entrez(...)` | Mapea una lista de genes HUGO a Entrez ID (usando `mygene`). |
| | `preparar_red(...)` | Mapea la red completa de HUGO a Entrez ID y limpia filas nulas. |
| **Grafo** | `construir_red(...)` | Crea el objeto `networkx.Graph` con Entrez IDs y pesos normalizados. |
| **DIAMOnD** | `diamond_iteration_of_first_X_nodes(...)` | Implementa el algoritmo iterativo de DIAMOnD basado en p-valores hipergeométricos. |
| | `pvalue(...)`, `gauss_hypergeom(...)` | Funciones auxiliares para el cálculo del p-valor. |
| **GUILD** | `guild_propagation(...)` | Implementa la difusión de información (Random Walk) para priorizar genes por *score*. |
| **Salida** | `guardar_resultados(...)` | Genera el archivo TSV de resultados con la clasificación comparativa (DIAMOnD vs GUILD). |
| | `graficar_red_enriquecida(...)` | Genera el archivo PNG con el grafo enriquecido y codificado por colores. |
| **Control** | `main()` | Gestiona los argumentos CLI, crea la carpeta `results/` y coordina la ejecución de todos los pasos. |
