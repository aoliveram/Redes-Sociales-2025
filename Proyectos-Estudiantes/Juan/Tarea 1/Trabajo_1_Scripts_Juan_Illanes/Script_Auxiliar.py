#Cargado de librerias

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx

from matplotlib.patches import Patch

#Funciones auxiliares

def ordenar_matriz_por_sumas(matriz):
    """
    Ordena una matriz:
    1. Filas: de mayor a menor según la suma de cada fila.
    2. Columnas: de mayor a menor según la suma de cada columna.

    Parámetro:
    matriz (array-like 2D): matriz de entrada.

    Retorna:
    matriz_ordenada (ndarray): matriz ordenada.
    """
    M = np.array(matriz)

    # Ordenar filas por suma (descendente)
    fila_indices = np.argsort(M.sum(axis=1))[::-1]
    M = M[fila_indices, :]

    # Ordenar columnas por suma (descendente)
    col_indices = np.argsort(M.sum(axis=0))[::-1]
    M = M[:, col_indices]

    return M


def reconectar_maximum_spanning(G, k_star):
    """
    G: grafo networkx no dirigido y pesado
    k_star: grado medio objetivo

    Devuelve:
    - H: red resultante
    - edge_star: umbral mínimo de peso que logra k_star
    """
    N = len(G.nodes)

    def average_degree(G):
        return sum(dict(G.degree()).values()) / N

    # 1. Árbol generador máximo
    H = nx.maximum_spanning_tree(G, weight='weight')

    # 2. Grado medio actual
    k_actual = average_degree(H)
    if k_actual >= k_star:
        return H, None  # no se necesitó edge_star

    # 3. Aristas que no están en el MST, ordenadas por peso descendente
    edges_sorted = sorted(
        (e for e in G.edges(data=True) if not H.has_edge(e[0], e[1])),
        key=lambda x: x[2]['weight'],
        reverse=True
    )

    edge_star = None

    # 4. Añadir hasta alcanzar k_star
    for u, v, data in edges_sorted:
        H.add_edge(u, v, weight=data['weight'])
        k_actual = average_degree(H)
        edge_star = data['weight']  # el último que se agregó
        if k_actual >= k_star:
            break

    return H, edge_star


def Balassa(X):
    X_c = np.sum(X, axis=1)
    X_p = np.sum(X, axis=0)
    X_net = np.sum(X)
    M_cp = np.zeros(X.shape)
    for i in range(len(X_c)):
        for j in range(len(X_p)):
            if (X_c[i] != 0) and (X_p[j] != 0):
                M_cp[i, j] = (X[i, j] / X_c[i]) / (X_p[j] / X_net)
    return M_cp


def similarity(M):
    c, p = M.shape
    phi = np.zeros((p, p))
    M_p = np.sum(M, axis=0)[:, np.newaxis]

    SUM = M.T @ M
    DEN = M_p * (M_p > M_p.T) + M_p.T * (M_p.T >= M_p)
    DEN = (DEN != 0) * DEN + (DEN == 0) * 1
    return SUM / DEN - np.diag(np.diag(SUM / DEN))


edge_normalization = lambda x: x - 0.5 * x ** 2


def normalization(df):
    a = df.min()
    b = df.max()
    return (df - a) / (b - a)


def array_trip(comunas):
    return [comuna[10:] for comuna in comunas]


def geographic_layout(G, coordenadas_dict):
    pos = {}

    lat_central = -33.0

    for nodo in G.nodes():
        if nodo in coordenadas_dict:
            lat, lng = coordenadas_dict[nodo]

            # Proyección cilíndrica equidistante centrada en Chile
            # Mantiene distancias verdaderas en el meridiano central
            x = lng * np.cos(np.radians(lat_central))  # Escala longitud por latitud central
            y = lat  # Latitud directa para mantener Chile "recto"

            pos[nodo] = (2 * x, y)

    return pos

#Parametros utiles

eq_mon = {
    'USD' : 950,
    'CLF' : 38_500,
    'UTM': 66_000,
    'CLP': 1
}

stringing = lambda x: float(x.replace(',', '.')) if isinstance(x, str) else float(x)

Columnas = [
    'Codigo', 'Nombre', 'Tipo de Adquisicion', 'Estado', 'CodigoOrganismo', 'NombreOrganismo', 'sector', 'RutUnidad', 'ComunaUnidad', 'RegionUnidad', 'CodigoMoneda', 'FechaAdjudicacion',
    'CodigoProductoONU', 'RutProveedor', 'NombreProveedor', 'MontoUnitarioOferta', 'CantidadAdjudicada', 'MontoLineaAdjudica', 'Oferta seleccionada'
]

#Cargado y limpieza de dataframe

data_root = rf'./data/raw/lic_2024-1.csv'
df = pd.read_csv(data_root, encoding = 'latin', sep = ';').loc[:, Columnas]

df['NombreProveedor'] = '_' + df['NombreProveedor']

df['FechaAdjudicacion'] = pd.to_datetime(df['FechaAdjudicacion'])
df = df.set_index('FechaAdjudicacion')

df['MontoUnitarioOferta'] = (df['MontoUnitarioOferta'].apply(stringing))
df['EquivalenciaMoneda'] = (df['CodigoMoneda'].map(eq_mon))
df['CantidadAdjudicada'] = df['CantidadAdjudicada'].apply(stringing)

df['MontoLineaAdjudica'] = df['MontoUnitarioOferta'] * df['EquivalenciaMoneda'] * df['CantidadAdjudicada']

Condiciones = (
    (df['MontoUnitarioOferta'] > 1e1) & (df['MontoLineaAdjudica'] < 1e14) & ( df['MontoUnitarioOferta'] != df['CantidadAdjudicada']) & (df['Estado'] == 'Adjudicada')
)

df = df.loc[
    Condiciones, :
].dropna()

df_selec = df.loc[
    df['Oferta seleccionada'] == 'Seleccionada'
].copy()

df_unselec = df.loc[
    df['Oferta seleccionada'] == 'No Seleccionada'
].copy()

del df

#Grafo bipartito

B_seleccionados = nx.Graph()

proveedores_seleccionados = df_selec['NombreProveedor'].unique()
comunas_seleccionadas = df_selec['ComunaUnidad'].unique()

B_seleccionados.add_nodes_from(proveedores_seleccionados, bipartite="proveedor")
B_seleccionados.add_nodes_from(comunas_seleccionadas, bipartite="comuna")

df_counts = (
    df_selec.groupby(["ComunaUnidad", "NombreProveedor"])
    .size()
    .reset_index(name="weight")
)

for _, row in df_counts.iterrows():
    B_seleccionados.add_edge(row["ComunaUnidad"], row["NombreProveedor"], weight=row["weight"])

# Metodo HH

X = nx.bipartite.biadjacency_matrix(B_seleccionados, row_order=comunas_seleccionadas, column_order=proveedores_seleccionados).toarray()

map_comunas_selec = {i: nodo for i, nodo in enumerate(comunas_seleccionadas)}
map_proveedores_selec = {j: nodo for j, nodo in enumerate(proveedores_seleccionados)}

RCA = Balassa(X)
M_cp = RCA > 1
phi_pq = similarity(M_cp)

G = nx.from_numpy_array(phi_pq)

#Gráfico proyectado - Empresas

G_spanning= nx.maximum_spanning_tree(G)

N = G_spanning.number_of_nodes()
L = np.sum([k for (u, k) in G_spanning.degree(weight = 'weight')])

G_spanning_, a = reconectar_maximum_spanning(G, 3)

Enlaces = np.array([w['weight'] for (u, v, w) in G_spanning_.edges(data = True)])
Grado = np.array([k for (u, k) in G_spanning_.degree(weight = 'weight')])
pos = nx.spring_layout(G_spanning_)

fig, ax = plt.subplots(figsize = (10, 8), dpi = 100)

df_sector_unico = (
    df_selec.groupby("NombreProveedor")["sector"]
      .first()   # toma el primer valor encontrado
      .reset_index()
)

color = dict(
    zip(
        list(df_selec['sector'].unique()),
        [
    '#1f77b4',
    '#ff7f0e',
    '#2ca02c',
    '#d62728',
    '#9467bd',
    '#8c564b',
    '#e377c2'
]
    )
)

coloreado = [
    color[sector] for sector in list(df_sector_unico['sector'])
]

nodes = nx.draw_networkx_nodes(
    G_spanning_,
    pos = pos,
    node_size = Grado ** 0.9,
    node_color = coloreado,
    ax = ax
)

nx.draw_networkx_edges(
    G_spanning_,
    pos = pos,
    width = edge_normalization(Enlaces),
    alpha = 0.2,
    ax = ax
)

legend_elements = [Patch(facecolor=col, label=sec)
                   for sec, col in color.items()]
ax.legend(handles=legend_elements, title="Sector")

plt.axis("off")

plt.tight_layout()
plt.savefig('./figs/Public_Licitation_space.pdf', dpi = 100)

#Grafico Bipartito

pos_bip = nx.bipartite_layout(B_seleccionados, proveedores_seleccionados)

prov_degree = [B_seleccionados.degree()[nodo] for nodo in proveedores_seleccionados]
com_degree = [B_seleccionados.degree()[nodo] for nodo in comunas_seleccionadas]

fig, ax = plt.subplots(figsize = (10, 6), dpi = 100)

ax.text(-0.12, -0.47, 'Empresas', fontsize = 12)
ax.text(0.94, -0.47, 'Comunas', fontsize = 12)

nx.draw_networkx_nodes(B_seleccionados, pos = pos_bip, nodelist = proveedores_seleccionados, node_color = 'tab:blue', node_size = prov_degree, ax = ax)
nx.draw_networkx_nodes(B_seleccionados, pos = pos_bip, nodelist = comunas_seleccionadas, node_color = 'tab:orange', node_size = com_degree, ax = ax)
nx.draw_networkx_edges(B_seleccionados, pos = pos_bip, width = 0.001, ax = ax, alpha = 0.03)

ax.axis('off')

plt.savefig('./figs/biparite_licitation.pdf', dpi = 100)

#Grafico proyectado - Comunas

geografia = pd.read_csv('./data/geo_comunas.csv')

geografia['latitud'] = geografia['latitud']
geografia['longitud'] = geografia['longitud']

M = X @ X.T
M = M - np.diag(np.diag(M))

comunas_proyec =  nx.from_numpy_array(M)
comunas_proyec = nx.relabel_nodes(comunas_proyec, map_comunas_selec)

dicc_comunas = {
    row["nombre"]: np.array([row["latitud"], row["longitud"]])
    for _, row in geografia.iterrows()
}

G = comunas_proyec
coordenadas_dict = dicc_comunas

nodos_con_coords = [n for n in G.nodes() if n in coordenadas_dict]
subgrafo = G.subgraph(nodos_con_coords)

pos = geographic_layout(subgrafo, coordenadas_dict);

fig, ax = plt.subplots(figsize=(4, 10), dpi = 100)

nodos_norte = []  # lat > -30° (más al norte, menos negativo)
nodos_centro = []  # -37.5° <= lat <= -30°
nodos_sur = []    # lat < -37.5° (más al sur, más negativo)

for nodo in nodos_con_coords:
    lat = coordenadas_dict[nodo][0]
    if lat > -30:
        nodos_norte.append(nodo)
    elif lat >= -37.5:
        nodos_centro.append(nodo)
    else:
        nodos_sur.append(nodo)


# Dibujar aristas
nx.draw_networkx_edges(subgrafo, pos, edge_color='gray', alpha=0.1, width = 0.5)

# Dibujar nodos por categoría de latitud
if nodos_norte:
    nx.draw_networkx_nodes(subgrafo, pos, nodelist=nodos_norte,
                          node_color='red', node_size=20)

if nodos_centro:
    nx.draw_networkx_nodes(subgrafo, pos, nodelist=nodos_centro,
                          node_color='blue', node_size=10)

if nodos_sur:
    nx.draw_networkx_nodes(subgrafo, pos, nodelist=nodos_sur,
                          node_color='green', node_size=10)

zona = {'red':'Norte', 'blue':'Centro', 'green':'Sur'}

legend_elements = [Patch(facecolor=col, label=sec)
                   for col, sec in zona.items()]
ax.legend(handles=legend_elements, title="Sector")

ax.axis('off')
plt.savefig('./figs/comuna_proyection.pdf')