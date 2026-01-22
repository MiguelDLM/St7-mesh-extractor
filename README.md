# St7 Mesh Extractor

A Python script to extract 3D meshes from Strand7 `.st7` files and export them to OBJ format.

## Description

This tool parses binary Strand7 (.st7) structural analysis files and extracts geometric mesh data including:
- **Nodes**: 3D vertex coordinates (X, Y, Z)
- **Plate elements**: Surface meshes (triangles and quads)
- **Brick elements**: Volume meshes (tetrahedra, hexahedra, pentahedra) with automatic surface skinning

The extracted mesh is saved as a standard Wavefront OBJ file that can be imported into 3D modeling software like Blender, Maya, MeshLab, etc.

## Features

- ✅ Parses Strand7 binary format (.st7)
- ✅ Extracts all node coordinates
- ✅ Handles multiple element types:
  - Triangles and quads (plate elements)
  - Tetrahedra (4/10 nodes)
  - Hexahedra (8/20 nodes)
  - Pentahedra/Prisms (6/15 nodes)
- ✅ Automatic surface extraction from volume elements (skinning)
- ✅ Exports to standard OBJ format

## Requirements

- Python 3.6+
- No external dependencies (uses only standard library)

## Installation

Clone this repository:

```bash
git clone https://github.com/MiguelDLM/St7-mesh-extractor.git
cd St7-mesh-extractor
