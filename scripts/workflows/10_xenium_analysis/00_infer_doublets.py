import tifffile
import matplotlib.pyplot as plt
from pathlib import Path
import os
import io
from shapely.geometry import Polygon

from rtree import index
from tqdm import tqdm
import numpy as np
import pandas as pd
import argparse
import shapely
import gzip

from shapely.geometry import Polygon
from shapely import minimum_rotated_rectangle

import networkx as nx


sample_path = Path("data_raw/output-XETG00078__0018188__Region_1__20240301__143248/")


x_limit = 6000
y_limit = 3500


x_lines = [0, 5000, 11000]
y_lines = [0, 3500, 6500, 11000, 16000, 21000]



def get_long_short_sides_lengths(geometry):
    # Get the minimum rotated rectangle
    min_rect = minimum_rotated_rectangle(geometry)
    # Extract coordinates of the rectangle
    rect_coords = min_rect.exterior.coords
    # Calculate lengths of sides
    side_lengths = []
    for i in range(len(rect_coords) - 1):
        x1, y1 = rect_coords[i]
        x2, y2 = rect_coords[i + 1]
        side_length = ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5
        side_lengths.append(side_length)
    # Find the longest and shortest sides
    longest_side = max(side_lengths)
    shortest_side = min(side_lengths)
    return longest_side, shortest_side


def create_index(boundaries_df, column='cell_id', extension=1):
    idx = index.Index()
    polygons = []
    i = 0
    for _, (cell_id, polygon) in tqdm(enumerate(boundaries_df.groupby(column))):
        # print(cell_id)
        polygon = Polygon(polygon[['vertex_x', 'vertex_y']])
        long, short = get_long_short_sides_lengths(polygon)
        if long / short > 10:
            print(_)
            continue
        polygons.append({
            column: cell_id,
            "id": i,
            "polygon": polygon,
            "bounds": polygon.bounds,
            "extended_polygon": polygon.buffer(extension),
            "extended_bounds": polygon.buffer(extension).bounds
        })
        idx.insert(i, polygons[i]['extended_bounds'])
        i += 1
    return polygons, idx


with gzip.open(sample_path / "nucleus_boundaries.csv.gz", 'rb') as file:
    nucleus_boundaries_df = pd.read_csv(io.BytesIO(file.read()))


# f1_nuc = nucleus_boundaries_df[
#     (nucleus_boundaries_df["vertex_x"] > x_limit) & 
#     (nucleus_boundaries_df["vertex_y"] < y_limit )
# ]

gg = create_index(nucleus_boundaries_df, extension=0.5)

idx = index.Index()
polygons = []
i = 0
for _, (cell_id, polygon) in tqdm(enumerate(nucleus_boundaries_df.groupby('cell_id'))):
    # print(cell_id)
    polygon = Polygon(polygon[['vertex_x', 'vertex_y']])
    long, short = get_long_short_sides_lengths(polygon)
    if long / short > 10:
        print(_)
        continue
    polygons.append({
        'cell_id': cell_id,
        "id": i,
        "polygon": polygon,
        "bounds": polygon.bounds,
        "extended_polygon": polygon.buffer(1),
        "extended_bounds": polygon.buffer(1).bounds
    })
    idx.insert(i, polygons[i]['extended_bounds'])
    i += 1

f1_polygons, f1_idx = (polygons, idx)

edges = []
for polygon_a in tqdm(f1_polygons):
    for polygon_b_index in f1_idx.intersection(polygon_a['extended_bounds']):
        polygon_b = f1_polygons[polygon_b_index]
        if polygon_b['cell_id'] == polygon_a['cell_id']:
            continue
        if not polygon_a['extended_polygon'].intersects(polygon_b['extended_polygon']):
            continue
        if not polygon_a['polygon'].intersects(polygon_b['extended_polygon']):
            continue
        if not polygon_a['extended_polygon'].intersects(polygon_b['polygon']):
            continue
        intersection = polygon_a['extended_polygon'].intersection(polygon_b['extended_polygon'])
        edges.append((polygon_a['cell_id'], polygon_b['cell_id'], intersection.area))




G = nx.Graph()
G.add_nodes_from([p['cell_id'] for p in f1_polygons])
for edge in edges:
    G.add_edge(edge[0], edge[1], weight=edge[2])

positions = [[p['polygon'].centroid.x, p['polygon'].centroid.y] for p in f1_polygons]

components = list(nx.connected_components(G))

processed_components = [{cell_id} if len(component) != 2 else component for component in components for cell_id in component]

id2component_size = {}
for c in components:
    for cell_id in c:
        id2component_size[cell_id] = len(c)


df = pd.DataFrame(f1_polygons)
df['component_size'] = df['id'].apply(lambda x: id2component_size[x])

df.to_csv('data_tidy/inferred_components.csv')


import json

# Save the processed fusion list to a JSON file
with open('inferred_components.json', 'w') as f:
    json.dump([list(component) for component in processed_components], f)

# triplets_pairs = []
# for c in [c for c in components if len(c) == 3]:
#     sub_G = G.subgraph(c)
#     edges, weights = zip(*nx.get_edge_attributes(sub_G,'weight').items())
#     triplets_pairs.append(sorted(weights, reverse=True)[:2])