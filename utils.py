import os
import geopandas as gpd
import pandas as pd
import networkx as nx
from shapely import LineString, Polygon, MultiPolygon, distance, intersects, minimum_bounding_radius as min_radius
from shapely.geometry import box
from shapely.wkt import loads, dumps


# Códigos URM em SIRGAS 2000
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

def find_utm_proj(X, Y):
    '''
    Identifica projeção local UTM em SIRGAS 2000
    '''
    # Hemisfério
    h = 'N' if Y > 0 else 'S'
    # Fuso
    if -84 <= X < -30: # Ponto dentro do intervalo do Brasil
        f = int((X//6) + 31)
    else:
        f = ''
    fuse = f'{f}{h}'
    return UTMCODES[fuse]

def geomIsResidual(geom):
    '''
    Avalia se geometria é espúria
    '''
    ap_ratio = geom.area/geom.length
    return (ap_ratio < 0.2)

def removeHoles(geom, area_min=1):
    '''
    Remove buracos de uma geometria
    '''
    if isinstance(geom, Polygon):
        geom = MultiPolygon([geom])
    out_polys = []
    for part in geom.geoms:
        interiors = []
        for i in part.interiors:
            p = Polygon(i)
            if p.area > area_min:
                interiors.append(i)
        out_polys.append(Polygon(part.exterior.coords, holes=interiors))
    return MultiPolygon(out_polys) if len(out_polys)>1 else out_polys[0]

def makeCensusGrid(gpkg_file, id_column, layer=None, mun_filter=[]):
    '''
    Importa malha de setores e a configura
    '''
    gdf = gpd.read_file(gpkg_file, layer=layer)[[id_column, 'geometry']]
    gdf = gdf.rename(columns={id_column:'ID'})
    gdf['ID'] = gdf['ID'].astype(str)
    if mun_filter:
        gdf = gdf[gdf['ID'].apply(lambda x: x[:7] in mun_filter)]

    # Remoção de sufixos literais
    gdf['ID'] = gdf['ID'].apply(lambda x: ''.join([i for i in str(x) if i.isdigit()]))

    # Correção das geometrias
    gdf['geometry'] = gdf['geometry'].make_valid()

    # Campo de grupo a partir da estrutura do id
    gdf['GROUP'] = gdf['ID'].apply(lambda x: str(x)[:11])

    # Reprojetar camada
    Y = (gdf.total_bounds[1] + gdf.total_bounds[3])/2
    X = (gdf.total_bounds[0] + gdf.total_bounds[2])/2
    UTMCRS = find_utm_proj(X, Y)
    gdf = gdf.to_crs(UTMCRS)

    # Calcular vizinhos
    gdf['NEIGHBOR'] = gdf.apply(lambda x: getNeighbors(x, gdf), axis=1)

    # Calcular área
    gdf['AREA'] = gdf.area

    return gdf


def getNeighbors(row, geogrid):
    '''
    Identifica vizinhos de cada geometria de uma malha
    '''
    df_neighbors = geogrid.iloc[geogrid['geometry'].sindex.query(row['geometry'], predicate='intersects')]
    df_neighbors = df_neighbors[df_neighbors['ID'] != row['ID']]
    if 'GROUP' in geogrid.columns:
        group = row['GROUP']
        df_neighbors = df_neighbors.query('GROUP == @group')

    return df_neighbors['ID'].to_list()


def makeNeighborhoodGraph(gdf):
    '''
    Cria grafo de vizinhança para conferência
    '''
    G = nx.Graph()
    # Adicionar nós
    for i, row in gdf[['ID', 'geometry']].iterrows():
        G.add_node(row['ID'], center=row['geometry'].representative_point())
    # Dados das arestas
    dic = gdf[['ID', 'NEIGHBOR']]\
          .explode('NEIGHBOR').dropna()\
          .to_dict(orient='records')
    # Adicionar arestas
    for i in dic:
        if (i['ID'], i['NEIGHBOR']) not in G.edges:
            G.add_edge(i['ID'], i['NEIGHBOR'],
                       geom=LineString([G.nodes[i['ID']]['center'], G.nodes[i['NEIGHBOR']]['center']]))
    return G


def exportNeighborhoodGraph(gdf, gpkg_file, layer=''):
    '''
    Exporta grafo de vizinhança em gpkg
    '''
    G = makeNeighborhoodGraph(gdf)
    df = pd.DataFrame(G.edges.data()).rename(columns={0:'node1', 1:'node2', 2:'geometry'})
    df['geometry'] = df['geometry'].apply(lambda x: x['geom'])
    geograph = gpd.GeoDataFrame(df, geometry='geometry', crs=gdf.crs)
    geograph.to_file(gpkg_file, layer=layer, driver='GPKG')


def getNeighborhoods(geogrid):
    '''
    Retorna coluna de vizinhos dos setores como dicionário
    '''
    dic = geogrid[['ID', 'NEIGHBOR']].to_dict(orient='records')
    dic = {k['ID']:k['NEIGHBOR'] for k in dic}
    return dic


def getGeometries(geogrid):
    '''
    Obtém lista de geometrias indexadas pelo id dos setores
    '''
    geoms = geogrid[['ID', 'geometry']]
    geoms = geoms[['ID', 'geometry']].to_dict(orient='records')
    geoms = {k['ID']:k['geometry'] for k in geoms}
    return geoms

def evaluateUnionArea(gridUnion, ):
    gdf = gridUnion.query('`A.ID` == `B.ID` & `A.PCT_AREA')


def bboxIntersects(geom1, geom2):
    '''
    Testa interseção de bboxes
    '''
    bbox1 = box(*geom1.bounds)
    bbox2 = box(*geom2.bounds)
    return intersects(bbox1, bbox2)

def format_corresp(a, b):
    if a>1:
        if b>1:
            return 'm:n'
        else:
            return 'n:1'
    else:
        if b>1:
            return '1:n'
        else:
            return '1:1'

class compatibility_graph(nx.Graph):
    def __init__(self, m1, m2):
        '''
        Classe de grafo de compatibilização criado das malhas m1 e m2
        '''
        super().__init__()
        self.m1 = m1
        self.m2 = m2
        self.makeGridUnion()
        self.classifyGridChanges()
        self.setNodes()


    def _classifyIdChanges(self):
        '''
        Classifica mudanças entre duas malhas
        '''
        # Dicionários de relações
        dic_1 = getNeighborhoods(self.m1)
        dic_2 = getNeighborhoods(self.m2)

        # Dicionários de geometrias
        geoms_1 = getGeometries(self.m1)
        geoms_2 = getGeometries(self.m2)

        # Calcular correspondencias externas
        outermatchA = self.gridUnion[self.gridUnion['OUTER_MATCH']]['A.ID'].to_list()
        outermatchB = self.gridUnion[self.gridUnion['OUTER_MATCH']]['B.ID'].to_list()
        outermatch = outermatchA + outermatchB

        # Classificação das mudanças
        class_data = {}
        for k, neighbors in dic_1.items():
            if k not in dic_2: 
                class_data[k] = 'Exclusão'
        for k, neighbors in dic_2.items():
            id_kept = k in dic_1
            if id_kept:
                geometry_out = k in outermatch
            neighbors_out = id_kept and not any([i in dic_1 for i in neighbors])
            distance_out = id_kept and distance(geoms_1[k].centroid, geoms_2[k].centroid) > (min_radius(geoms_1[k])+min_radius(geoms_2[k]))
            
            if id_kept and (distance_out or neighbors_out or geometry_out):
                class_data[k] = 'Desassociação'
            elif id_kept:
                class_data[k] = 'Manutenção'
            elif k not in dic_1 and k in dic_2:
                class_data[k] = 'Criação'
            else:
                class_data[k] = ''
        return class_data
    
    def classifyGridChanges(self):
        '''
        Aplica o classificador nas malhas
        '''
        alteracoes = self._classifyIdChanges()
        self.m1['CLASS'] = self.m1['ID'].apply(lambda x: alteracoes[x])
        self.m2['CLASS'] = self.m2['ID'].apply(lambda x: alteracoes[x])
        self.gridUnion['A.CLASS'] = self.gridUnion['A.ID'].apply(lambda x: alteracoes[x])
        self.gridUnion['B.CLASS'] = self.gridUnion['B.ID'].apply(lambda x: alteracoes[x])
        
    def setNodes(self):
        '''
        Cria os nós do grafo
        '''
        for _, row in self.m1[['ID', 'GROUP', 'CLASS', 'geometry']].iterrows():
            self.add_node(f"A.{row['ID']}",
                        malha='A',
                        nome=row['ID'],
                        group=row['GROUP'],
                        geom=row['geometry'],
                        center=row['geometry'].representative_point(),
                        classe=row['CLASS'])
        for _, row in self.m2[['ID', 'GROUP', 'CLASS', 'geometry']].iterrows():
            self.add_node(f"B.{row['ID']}",
                        malha='B',
                        nome=row['ID'],
                        group=row['GROUP'],
                        geom=row['geometry'],
                        center=row['geometry'].representative_point(),
                        classe=row['CLASS'])


    def _prepareGridForUnion(self, geogrid, prefix=''):
        '''
        Faz tratamento das geometrias para interseção
        '''
        utmcrs = geogrid.crs
        geogrid = geogrid.rename(columns={k:f'{prefix}.{k}' for k in geogrid.columns})
        geogrid = geogrid.set_geometry(f'{prefix}.geometry', crs=utmcrs)
        # Correção contra bug GEOSException: TopologyException: found non-noded intersection
        geogrid[f'{prefix}.geometry'] = [loads(dumps(geom, rounding_precision=3)) for geom in geogrid[f'{prefix}.geometry']]
        return geogrid
    
    def makeGridUnion(self):
        '''
        Cria grid intersecionado de m1 e m2
        '''
        # Preparação para interseção
        gA = self._prepareGridForUnion(self.m1, 'A')
        gB = self._prepareGridForUnion(self.m2, 'B')

        # Interseção e tratamento para exclusão de geometrias residuais
        self.gridUnion = gpd.overlay(gA, gB, how='union')\
                        .dropna(subset=['A.ID', 'B.ID'])
        self.gridUnion['RESIDUAL'] = self.gridUnion['geometry'].apply(geomIsResidual)
        self.gridUnion['UNION_AREA'] = self.gridUnion.area

        # Indicadores de área de interseção
        area_data_A = (self.gridUnion[~self.gridUnion['RESIDUAL']]
                       .pivot_table(index='A.ID', values='UNION_AREA', aggfunc='sum')
                       .reset_index()
                       .rename(columns={'UNION_AREA':'A.ORIGINAL_AREA'}))
        self.gridUnion = self.gridUnion.merge(area_data_A, on='A.ID', how='left')
        self.gridUnion['A.PCT_AREA'] = self.gridUnion.area/self.gridUnion['A.ORIGINAL_AREA']

        area_data_B = (self.gridUnion[~self.gridUnion['RESIDUAL']]
                       .pivot_table(index='B.ID', values='UNION_AREA', aggfunc='sum')
                       .reset_index()
                       .rename(columns={'UNION_AREA':'B.ORIGINAL_AREA'}))
        self.gridUnion = self.gridUnion.merge(area_data_B, on='B.ID', how='left')
        self.gridUnion['B.PCT_AREA'] = self.gridUnion.area/self.gridUnion['B.ORIGINAL_AREA']

        # Recalcular geometria residual
        self.gridUnion['RESIDUAL'] = self.gridUnion.apply(lambda x: x['RESIDUAL'] or x['A.PCT_AREA']<0.05 or x['B.PCT_AREA']<0.05, axis=1)

        # Avalia correspondência de área entre ids distintos acima de 80%
        self.gridUnion['OUTER_MATCH'] = self.gridUnion.apply(lambda x: x['A.ID'] != x['B.ID'] and (x['A.PCT_AREA'] > 0.8 or x['B.PCT_AREA'] > 0.8), axis=1)



    def compatManutencao(self):
        '''
        Adiciona arestas de manutenção
        '''
        for s in self.m2.query('CLASS == "Manutenção"')['ID']:
            self.add_edge(f"A.{s}",
                          f"B.{s}",
                          metodo='Manutenção')

    
    def compatDivisao(self, threshold):
        '''
        Adiciona arestas de divisão
        '''
        # Selecionar setores para análise
        intersecao = self.gridUnion.query('`A.CLASS` != "Manutenção" | `B.CLASS` != "Manutenção"')

        # Caso setores desassociados com interseção tenham o mesmo id
        for s in intersecao.query('`A.ID` == `B.ID`')['B.ID']:
            self.add_edge(f"A.{s}", f"B.{s}", metodo='Manutenção desassociada')

        # Registrar divisão (>=threshold da área original de B em um único setor de A)
        for _, row in intersecao.query('RESIDUAL == False & `B.PCT_AREA` >= @threshold').iterrows():
            self.add_edge(f"A.{row['A.ID']}",
                          f"B.{row['B.ID']}",
                          metodo='Divisão')
        self.clearGroups()

    
    def compatSobreposicao(self, buffer, use_all=False):
        '''
        Adiciona arestas de sobreposição forçada
        '''
        list_isolados = list(nx.isolates(self))

        # Selecionar setores ainda não vinculados
        list_isolados_A = [i.split('.')[-1] for i in list_isolados if i.startswith('A.')]
        isolados_A = self.m1.query('ID in @list_isolados_A')
        isolados_A = isolados_A.rename(columns={k:f'A.{k}' for k in isolados_A.columns})
        isolados_A['A.geometry'] = isolados_A['A.geometry'].buffer(buffer)
        isolados_A = isolados_A.set_geometry('A.geometry', crs=self.m1.crs)
        if not use_all:
            buffered_B = self.m2.query('CLASS != "Manutenção"')
        else:
            buffered_B = self.m2.copy()
        buffered_B['B.geometry'] = buffered_B['geometry'].buffer(buffer)
        buffered_B = buffered_B.set_geometry('B.geometry', crs=self.m2.crs)
        # Interseção
        intersecao = gpd.overlay(isolados_A, buffered_B, how='union')\
                        .dropna(subset=['A.ID', 'ID'])
        intersecao['RESIDUAL'] = intersecao['geometry'].apply(geomIsResidual)
        for _, row in intersecao.query('RESIDUAL == False').iterrows():
            self.add_edge(f"A.{row['A.ID']}",
                            f"B.{row['ID']}",
                            metodo=f'Sobreposição ({buffer}m)')

        # Selecionar setores ainda não vinculados
        list_isolados_B = [i.split('.')[-1] for i in list_isolados if i.startswith('B.')]
        isolados_B = self.m2.query('ID in @list_isolados_B')
        isolados_B = isolados_B.rename(columns={k:f'B.{k}' for k in isolados_B.columns})
        isolados_B['B.geometry'] = isolados_B['B.geometry'].buffer(buffer)
        isolados_B = isolados_B.set_geometry('B.geometry', crs=self.m2.crs)
        if not use_all:
            buffered_A = self.m1.query('CLASS != "Manutenção"')
        else:
            buffered_A = self.m1.copy()
        buffered_A['A.geometry'] = buffered_A['geometry'].buffer(buffer)
        buffered_A = buffered_A.set_geometry('A.geometry', crs=self.m1.crs)
        # Interseção
        intersecao = gpd.overlay(isolados_B, buffered_A, how='union')\
                        .dropna(subset=['ID', 'B.ID'])
        intersecao['RESIDUAL'] = intersecao['geometry'].apply(geomIsResidual)
        for _, row in intersecao.query('RESIDUAL == False').iterrows():
            self.add_edge(f"A.{row['ID']}",
                            f"B.{row['B.ID']}",
                            metodo=f'Sobreposição ({buffer}m)')
        self.clearGroups()



    def clearGroups(self):
        '''
        Remove arestas que não pertencem ao mesmo grupo
        '''
        to_remove = []
        for u, v, d in self.edges(data=True):
            group_A = self.nodes[u]['group']
            group_B = self.nodes[v]['group']
            if group_A != group_B:
                to_remove.append([u, v])

        for u, v in to_remove:
            # Remove as arestas entre grupos apenas de nós que tem vínculos válidos
            if not all([e in to_remove for e in self.edges(u)]) and not all([e in to_remove for e in self.edges(v)]):
                self.remove_edge(u,v)


    def _codifyCompat(self):
        '''
        Codifica as componentes do grafo resultante
        '''
        componentes = [{'group':self.nodes[list(G)[0]]['group'], 'nodes':list(G)}\
               for G in nx.connected_components(self)]
        increment_id = {k['group']:0 for k in componentes}

        matriz_A = []
        matriz_B = []
        for c in componentes:
            increment_id[c['group']] += 1
            cod_c = f"{c['group']}{increment_id[c['group']]:05d}"

            c_A = [self.nodes[i]['nome'] for i in c['nodes'] if self.nodes[i]['malha']=='A']
            matriz_A.append({'CD_PERIMETRO':cod_c, 'ID':c_A})
            c_B = [self.nodes[i]['nome'] for i in c['nodes'] if self.nodes[i]['malha']=='B']
            matriz_B.append({'CD_PERIMETRO':cod_c, 'ID':c_B})

        # Criação dos DataFrames finais
        df_matriz_A = pd.DataFrame(matriz_A)
        self.compatTable_A = df_matriz_A.explode('ID')
        df_matriz_B = pd.DataFrame(matriz_B)
        self.compatTable_B = df_matriz_B.explode('ID')


    def _prepareExport(self):
        '''
        Cria tabelas de exportação
        '''
        self._codifyCompat()
        # Contagem de membros A dos perímetros
        data_matriz_A = self.compatTable_A.pivot_table(index='CD_PERIMETRO',
                                                values='ID',
                                                aggfunc='count').reset_index()
        data_matriz_A = data_matriz_A.rename(columns={'ID':'membros_A'})
        # Contagem de membros B dos perímetros
        data_matriz_B = self.compatTable_B.pivot_table(index='CD_PERIMETRO',
                                                values='ID',
                                                aggfunc='count').reset_index()
        data_matriz_B = data_matriz_B.rename(columns={'ID':'membros_B'})
        # Agregação dos dados
        data_matrizes = data_matriz_A.merge(data_matriz_B, on='CD_PERIMETRO')
        data_matrizes['membros'] = data_matrizes['membros_A'] + data_matrizes['membros_B']
        gdf_pc = self.compatTable_B.merge(self.m2, on='ID')
        gdf_pc = gdf_pc.merge(data_matrizes, on='CD_PERIMETRO')
        gdf_pc = gpd.GeoDataFrame(gdf_pc, geometry='geometry', crs=self.m2.crs)
        gdf_pc = gdf_pc[['CD_PERIMETRO', 'GROUP', 'membros', 'membros_A', 'membros_B', 'geometry']].dissolve(by='CD_PERIMETRO')
        gdf_pc = gdf_pc.rename(columns={'GROUP':'CD_DIST'})
        gdf_pc['CD_MUN'] = gdf_pc['CD_DIST'].apply(lambda x: x[:7])
        gdf_pc['TIPO_CORRESP'] = gdf_pc.apply(lambda x: format_corresp(x['membros_A'], x['membros_B']), axis=1)
        gdf_pc['geometry'] = gdf_pc['geometry'].apply(removeHoles)
        self.AMC = gdf_pc


    def _prepareGraphExport(self):
        '''
        Cria representação geográfica do grafo de compatibilidade
        '''
        # Arestas
        edge_data = []
        for u, v in list(self.edges):
            data_u = self.nodes[u]
            data_u = {f"{data_u['malha']}.{k}":value for k, value in data_u.items()}
            data_v = self.nodes[v]
            data_v = {f"{data_v['malha']}.{k}":value for k, value in data_v.items()}
            data_u.update(data_v)
            data_u.update(self.edges[(u, v)])
            data_u['geometry'] = LineString([data_u['A.center'], data_u['B.center']])
            edge_data.append(data_u)

        edge_gdf = gpd.GeoDataFrame(edge_data, geometry='geometry', crs=self.m2.crs)
        self.edge_gdf = edge_gdf[['A.nome', 'A.classe', 'B.nome', 'B.classe', 'metodo', 'geometry']]
       
        # Nós
        node_data = [i for _, i in list(self.nodes.data())]
        for k in node_data:
            k['grau'] = len(self[f"{k['malha']}.{k['nome']}"])
        node_gdf = gpd.GeoDataFrame(node_data, geometry='center', crs=self.m2.crs)
        self.node_gdf = node_gdf[['nome', 'malha', 'classe', 'group', 'grau', 'center']]
        

    def exportCompatFiles(self, compatName, name_C1, name_C2):
        '''
        Exporta os arquivos de compatibilização
        '''
        self._prepareExport()
        self._prepareGraphExport()
        self.compatTable_A[['ID', 'CD_PERIMETRO']].to_csv(f'malhas/{compatName}_{name_C1}.csv', sep='\t', index=False)
        self.compatTable_B[['ID', 'CD_PERIMETRO']].to_csv(f'malhas/{compatName}_{name_C2}.csv', sep='\t', index=False)
        self.AMC.to_file(f'malhas/{compatName}_AMC.gpkg',
                         layer=f'{name_C1}-{name_C2}',
                         driver='GPKG')
        # Exportar representação geográfica do grafo
        self.edge_gdf.to_file(f'malhas/{compatName}_AMC.gpkg',
                            layer=f'{name_C1}-{name_C2}_edges',
                            driver='GPKG')
        self.node_gdf.to_file(f'malhas/{compatName}_AMC.gpkg',
                            layer=f'{name_C1}-{name_C2}_nodes',
                            driver='GPKG')