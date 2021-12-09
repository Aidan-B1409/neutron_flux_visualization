import argparse
import pandas as pd
import numpy as np
from scipy.spatial import Voronoi, Delaunay
from plyfile import PlyData, PlyElement

def main():
    args = parseargs()
    df = pd.read_csv(args.file, delim_whitespace=True)
    df.columns = df.columns.str.replace(" ", "")
    df.drop(['Rel', 'Error'], axis=1, inplace=True)
  
    data = list(df.itertuples(index=False, name=None))
    vertices = np.array(data, dtype = [('x', 'f4'), ('y', 'f4'), ('z', 'f4'), ('s', 'f4')])

    points = np.array(list(map(lambda a: [a[1]['X'], a[1]['Y']], df.iterrows())))
    quads = Delaunay(points)

    faces = np.array([tuple([a]) for a in filter(lambda a: len(a) == 3, quads.simplices)], dtype = [('vertex_indices', 'i4', (3,))])
    
    verts_element = PlyElement.describe(vertices, "vertex")
    face_element = PlyElement.describe(faces, "face")
    PlyData([verts_element, face_element], text=True).write('geometry.ply')


def parseargs() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="A program for converting MCNP6.2 output file into .ply for visualization. CS553 @ OSU, Fall 2021")
    parser.add_argument("-f", "--file", required=True, dest='file', help='Path to input MCNP data.')
    return parser.parse_args()

if __name__ == "__main__":
    main()