import os
import geopandas as gpd
import pandas as pd
import networkx as nx
from shapely import LineString, Polygon, MultiPolygon, distance, intersects, minimum_bounding_radius as min_radius
from shapely.geometry import box
from shapely.wkt import loads, dumps


# Função para defi
UTMCODES = {
    '17S':"EPSG:31977",
    '18S':"EPSG:31978",
    '19S':"EPSG:31979",
    '20S':"EPSG:31980",
    '21S':"EPSG:31981",
    '22S':"EPSG:31982",
    '23S':"EPSG:31983",
    '24S':"EPSG:31984",
    '25S':"EPSG:31985",
    '17N':"EPSG:31971",
    '18N':"EPSG:31972",
    '19N':"EPSG:31973",
    '20N':"EPSG:31974",
    '21N':"EPSG:31975",
    '22N':"EPSG:31976",
    '23N':"EPSG:6210",
    '24N':"EPSG:6211"
    }

# Função para identificar projeção local UTM em SIRGAS 2000
def find_utm_proj(X, Y):
    # Hemisfério
    h = 'N' if Y > 0 else 'S'
    # Fuso
    if -84 <= X < -30: # Ponto dentro do intervalo do Brasil
        f = int((X//6) + 31)
    else:
        f = ''
    fuse = f'{f}{h}'
    return UTMCODES[fuse]

# Função/regra para avaliar se geometria deve ser considerada espuria
def geomIsResidual(geom):
    ap_ratio = geom.area/geom.length
    return (ap_ratio > 0.5 or geom.area > 500)

# Função para importar malha de setores e configurá-la
def makeCensusGrid(gpkg_file, id_column, layer=None, mun_filter=[]):
    gdf = gpd.read_file(gpkg_file, layer=layer)[[id_column, 'geometry']]
    gdf = gdf.rename(columns={id_column:'CD_GEOCODI'})
    gdf['CD_GEOCODI'] = gdf['CD_GEOCODI'].astype(str)
    if mun_filter:
        gdf = gdf[gdf['CD_GEOCODI'].apply(lambda x: x[:7] in mun_filter)]

    # Remoção de sufixos literais
    gdf['CD_GEOCODI'] = gdf['CD_GEOCODI'].apply(lambda x: ''.join([i for i in str(x) if i.isdigit()]))

    # Correção das geometrias
    gdf['geometry'] = gdf['geometry'].make_valid()

    # Campo de grupo a partir da estrutura do id
    gdf['GROUP'] = gdf['CD_GEOCODI'].apply(lambda x: str(x)[:11])

    # Reprojetar camada
    Y = (gdf.total_bounds[1] + gdf.total_bounds[3])/2
    X = (gdf.total_bounds[0] + gdf.total_bounds[2])/2
    UTMCRS = find_utm_proj(X, Y)
    gdf = gdf.to_crs(UTMCRS)

    # Calcular vizinhos
    gdf['NEIGHBOR'] = gdf.apply(lambda x: getNeighbors(x, gdf), axis=1)

    return gdf

# Função para identificar vizinhos de cada geometria de uma malha
def getNeighbors(row, geogrid):
    df_neighbors = geogrid.iloc[geogrid['geometry'].sindex.query(row['geometry'], predicate='intersects')]
    df_neighbors = df_neighbors[df_neighbors['CD_GEOCODI'] != row['CD_GEOCODI']]
    if 'GROUP' in geogrid.columns:
        group = row['GROUP']
        df_neighbors = df_neighbors.query('GROUP == @group')

    return df_neighbors['CD_GEOCODI'].to_list()

# Criar grafo de vizinhança para conferência
def makeNeighborhoodGraph(gdf):
    G = nx.Graph()
    # Adicionar nós
    for i, row in gdf[['CD_GEOCODI', 'geometry']].iterrows():
        G.add_node(row['CD_GEOCODI'], center=row['geometry'].representative_point())
    # Dados das arestas
    dic = gdf[['CD_GEOCODI', 'NEIGHBOR']]\
          .explode('NEIGHBOR').dropna()\
          .to_dict(orient='records')
    # Adicionar arestas
    for i in dic:
        if (i['CD_GEOCODI'], i['NEIGHBOR']) not in G.edges:
            G.add_edge(i['CD_GEOCODI'], i['NEIGHBOR'],
                       geom=LineString([G.nodes[i['CD_GEOCODI']]['center'], G.nodes[i['NEIGHBOR']]['center']]))
    return G

# Exportar grafo de vizinhança para conferência
def exportNeighborhoodGraph(gdf, gpkg_file, layer=''):
    G = makeNeighborhoodGraph(gdf)
    df = pd.DataFrame(G.edges.data()).rename(columns={0:'node1', 1:'node2', 2:'geometry'})
    df['geometry'] = df['geometry'].apply(lambda x: x['geom'])
    geograph = gpd.GeoDataFrame(df, geometry='geometry', crs=gdf.crs)
    geograph.to_file(gpkg_file, layer=layer, driver='GPKG')

# Listar vizinhos como dicionário
def getNeighborhoods(geogrid):
    dic = geogrid[['CD_GEOCODI', 'NEIGHBOR']].to_dict(orient='records')
    dic = {k['CD_GEOCODI']:k['NEIGHBOR'] for k in dic}
    return dic

# Obter lista de geometrias indexadas pelo id dos setores
def getGeometries(geodataframe):
    geoms = geodataframe[['CD_GEOCODI', 'geometry']]
    geoms = geoms[['CD_GEOCODI', 'geometry']].to_dict(orient='records')
    geoms = {k['CD_GEOCODI']:k['geometry'] for k in geoms}
    return geoms

# Teste de interseção de bboxes
def bboxIntersects(geom1, geom2):
    bbox1 = box(*geom1.bounds)
    bbox2 = box(*geom2.bounds)
    return intersects(bbox1, bbox2)

# Classificar mudanças entre duas malhas
def classifyIdChanges(m1, m2):
    # Dicionários de relações
    dic_1 = getNeighborhoods(m1)
    dic_2 = getNeighborhoods(m2)

    # Dicionários de geometrias
    geoms_1 = getGeometries(m1)
    geoms_2 = getGeometries(m2)

    # Classificação das mudanças
    class_data = {}
    for k, neighbors in dic_1.items():
        if k not in dic_2: 
            class_data[k] = 'Exclusão'
    for k, neighbors in dic_2.items():
        neighbors_out = k in dic_1 and k in dic_2 and not any([i in dic_1 for i in neighbors])
        distance_out = k in dic_1 and k in dic_2 \
                        and distance(geoms_1[k].centroid, geoms_2[k].centroid) > (min_radius(geoms_1[k])+min_radius(geoms_2[k])) \
                        and distance(geoms_1[k].centroid, geoms_2[k].centroid) >= 200 \
                        and not bboxIntersects(geoms_1[k], geoms_2[k])
        if distance_out:
            class_data[k] = 'Desassociação'
        elif neighbors_out:
            class_data[k] = 'Desassociação'
        elif k in dic_1 and k in dic_2:
            class_data[k] = 'Manutenção'
        elif k not in dic_1 and k in dic_2:
            class_data[k] = 'Criação'
        else:
            class_data[k] = ''
    return class_data

def makeCompatGraph(m1, m2):
    G_compat = nx.Graph()

    # Adicionar todos os nós
    for _, row in m1[['CD_GEOCODI', 'GROUP', 'CLASS', 'geometry']].iterrows():
        G_compat.add_node(f"A.{row['CD_GEOCODI']}",
                        malha='A',
                        nome=row['CD_GEOCODI'],
                        group=row['GROUP'],
                        geom=row['geometry'],
                        center=row['geometry'].representative_point(),
                        classe=row['CLASS'])
    for _, row in m2[['CD_GEOCODI', 'GROUP', 'CLASS', 'geometry']].iterrows():
        G_compat.add_node(f"B.{row['CD_GEOCODI']}",
                        malha='B',
                        nome=row['CD_GEOCODI'],
                        group=row['GROUP'],
                        geom=row['geometry'],
                        center=row['geometry'].representative_point(),
                        classe=row['CLASS'])
    
    return G_compat