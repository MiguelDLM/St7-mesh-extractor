#!/usr/bin/env python3
import struct
import sys
import math
import argparse
from pathlib import Path
from collections import defaultdict

def verify_element_block(data, offset, plate_count, brick_count, node_count):
    curr = offset
    total_to_check = min(plate_count + brick_count, 10) 
    if total_to_check == 0: return True 
    
    valid_elements = 0
    
    for _ in range(total_to_check):
        if curr + 10 > len(data): return False
        
        elem_id = struct.unpack_from('<I', data, curr)[0]
        prop_id = struct.unpack_from('<I', data, curr+4)[0]
        num_nodes = data[curr+8]
        
        if elem_id == 0 or elem_id > 500000000: return False
        if num_nodes == 0 or num_nodes > 32: return False
        
        # Check indices
        if curr + 9 + num_nodes*4 > len(data): return False
        n1 = struct.unpack_from('<I', data, curr+9)[0]
        if n1 < 1 or n1 > node_count: return False
        
        curr += 9 + num_nodes * 4
        valid_elements += 1
        
    return valid_elements > 0

def extract(filepath, outdir):
    path = Path(filepath)
    with open(path, 'rb') as f:
        data = f.read()
    
    print(f"Processing {path}...")
    
    # Header Reading
    # 0x10: Magic
    # 0x3B: NodeCount
    # 0x43: PlateCount
    # 0x47: BrickCount
    
    magic = b'STRAUS/STRAND'
    if data[0x10:0x10+len(magic)] != magic:
        print("Invalid Magic")
        return

    try:
        node_count = struct.unpack_from('<I', data, 0x3B)[0]
        plate_count = struct.unpack_from('<I', data, 0x43)[0]
        brick_count = struct.unpack_from('<I', data, 0x47)[0]
    except:
        print("Failed to read header counts")
        return

    print(f"Header: Nodes={node_count}, Plates={plate_count}, Bricks={brick_count}")
    
    # Find Node Block
    block_size = node_count * 24
    search_limit = len(data) - block_size
    node_offset = None
    
    for i in range(0, search_limit, 2): 
        try:
            v1 = struct.unpack_from('<d', data, i)[0]
            v2 = struct.unpack_from('<d', data, i+block_size-24)[0]
            if not (math.isfinite(v1) and abs(v1)<1e30 and math.isfinite(v2) and abs(v2)<1e30):
                continue
        except:
            continue
            
        elem_block_offset = i + block_size
        if verify_element_block(data, elem_block_offset, plate_count, brick_count, node_count):
            node_offset = i
            break
             
    if node_offset is None:
        print("Node block not found")
        return
        
    print(f"Node block confirmed at {node_offset}")
    
    # Extract Nodes
    nodes_xyz = []
    for i in range(node_count):
        off = node_offset + i * 24
        x, y, z = struct.unpack_from('<ddd', data, off)
        nodes_xyz.append((x,y,z))
        
    # Parse Elements
    curr = node_offset + block_size
    faces = [] 
    
    # 1. Plates (Surface)
    print(f"Parsing {plate_count} Plates...")
    for _ in range(plate_count):
        if curr + 10 > len(data): break
        num = data[curr+8]
        curr += 9
        elem_nodes = []
        for _ in range(num):
            elem_nodes.append(struct.unpack_from('<I', data, curr)[0])
            curr += 4
        
        if num >= 3:
            if num == 3: faces.append(tuple(elem_nodes))
            elif num == 4: faces.append(tuple(elem_nodes))
            elif num == 6: faces.append(tuple(elem_nodes[:3]))
            elif num >= 8: faces.append(tuple(elem_nodes[:4]))
    
    # 2. Bricks (Volume) - Skinning
    print(f"Parsing {brick_count} Bricks (Skinning)...")
    face_counts = defaultdict(int)
    step = max(1, brick_count // 10)
    
    for idx in range(brick_count):
        if idx % step == 0: print(f"  ... {int(idx/brick_count*100)}%")
        if curr + 10 > len(data): break
        num = data[curr+8]
        curr += 9
        en = []
        for _ in range(num):
            en.append(struct.unpack_from('<I', data, curr)[0])
            curr += 4
            
        sub_faces = []
        if num == 4 or num == 10: # Tet
            n = en if num==4 else en[:4]
            if len(n) >= 4:
                sub_faces = [
                    (n[0], n[2], n[1]), 
                    (n[0], n[1], n[3]),
                    (n[1], n[2], n[3]),
                    (n[2], n[0], n[3])
                ]
        elif num == 8 or num == 20: # Hex
            n = en if num==8 else en[:8]
            if len(n) >= 8:
                sub_faces = [
                    (n[0], n[3], n[2], n[1]),
                    (n[4], n[5], n[6], n[7]),
                    (n[0], n[1], n[5], n[4]),
                    (n[1], n[2], n[6], n[5]),
                    (n[2], n[3], n[7], n[6]),
                    (n[3], n[0], n[4], n[7])
                ]
        elif num == 6 or num == 15: # Penta
            n = en if num==6 else en[:6]
            if len(n) >= 6:
                sub_faces = [
                    (n[0], n[2], n[1]),
                    (n[3], n[4], n[5]),
                    (n[0], n[1], n[4], n[3]),
                    (n[1], n[2], n[5], n[4]),
                    (n[2], n[0], n[3], n[5])
                ]

        for sf in sub_faces:
            key = tuple(sorted(sf))
            face_counts[key] += 1
            
    # Write OBJ
    outname = outdir / (path.stem + ".obj")
    
    with open(outname, 'w') as f:
        f.write(f"# Extracted from {path.name}\n")
        f.write(f"# Nodes: {node_count}\n")
        
        for x,y,z in nodes_xyz:
            f.write(f"v {x} {y} {z}\n")
            
        # Write Plates
        for face in faces:
             f.write("f " + " ".join(str(n) for n in face) + "\n")
             
        # Write Skin Faces
        skin_count = 0
        for key, count in face_counts.items():
            if count == 1:
                f.write("f " + " ".join(str(n) for n in key) + "\n")
                skin_count += 1
                
    print(f"Written {outname}")
    print(f"Stats: {len(faces)} plates, {skin_count} brick-skin faces.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract OBJ mesh from Strand7 (.st7) file.')
    parser.add_argument('file', help='Path to .st7 file')
    parser.add_argument('outdir', help='Output directory')
    args = parser.parse_args()
    
    outd = Path(args.outdir)
    outd.mkdir(exist_ok=True, parents=True)
    
    extract(args.file, outd)
