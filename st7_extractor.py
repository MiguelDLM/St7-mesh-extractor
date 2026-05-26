#!/usr/bin/env python3
"""
Strand7 (.st7) extractor.

Pulls the geometric mesh and the engineering metadata we have reverse-engineered
out of a Strand7 model file and writes any combination of:

    OBJ  surface skin (brick boundary + plates as separate "attachment" groups)
    VTU  full solid volume mesh (every tetrahedron) for ParaView
    TXT  human-readable summary of materials, load cases, freedom cases, paths

See ST7_FORMAT.md for the file-format notes this is based on.
"""
import argparse
import math
import struct
import sys
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------

MAGIC = b"STRAUS/STRAND"
# Marker that begins a property record (Brick/Plate/Beam Property block).
# Bytes 0..2 identify the property family (04=brick, 03=plate, 02=beam), then
# `sub_kind` (always 0x01) and a uint16 length-of-name. The literal "Property"
# string never appears in a valid element record, so finding it is a reliable
# way to detect the end of the brick block when `brick_count` over-reads.
PROP_PROBE = b"\x01\x00\x10\x00Brick Prope"  # sub-kind 01, len=16, "Brick Prope"


def _u8(d, o):   return d[o]
def _u16(d, o):  return struct.unpack_from('<H', d, o)[0]
def _u32(d, o):  return struct.unpack_from('<I', d, o)[0]
def _f64(d, o):  return struct.unpack_from('<d', d, o)[0]


def _read_pstring(data, off):
    """Read a uint16-length-prefixed ASCII string. Returns (text, bytes_consumed)."""
    L = _u16(data, off)
    if L == 0 or off + 2 + L > len(data):
        return "", 2
    chunk = data[off + 2:off + 2 + L]
    try:
        return chunk.decode('ascii', errors='replace'), 2 + L
    except Exception:
        return "", 2 + L


# ---------------------------------------------------------------------------
# Header / node-block discovery
# ---------------------------------------------------------------------------

def _walk_plates(data, start, plate_count, node_count):
    """Walk forward through `plate_count` plate records. Returns end-offset
    on success, or None on misalignment. Plate records are small (3-8 node
    triangles / quads) so this is fast even for 100 K plates."""
    curr = start
    end = len(data)
    for _ in range(plate_count):
        if curr + 10 > end:
            return None
        type_id = _u32(data, curr)
        prop_id = _u32(data, curr + 4)
        n = _u8(data, curr + 8)
        if type_id == 0 or type_id > 256:    return None
        if prop_id == 0 or prop_id > 256:    return None
        if n == 0 or n > 32:                 return None
        rec_end = curr + 9 + n * 4
        if rec_end > end:                    return None
        # node-id range check (only the first one — fast)
        n1 = _u32(data, curr + 9)
        if n1 < 1 or n1 > node_count:        return None
        curr = rec_end
    return curr


def _check_brick_stride(data, start, brick_count, stride=25):
    """Verify that `brick_count` records of length `stride` (4-node tet)
    all have `num_nodes == 4` at byte +8. Uses byte slicing for speed."""
    end_needed = start + stride * brick_count
    if end_needed > len(data):
        return False
    # bytes at offsets start+8, start+8+stride, start+8+2*stride, ...
    slice_ = data[start + 8 : start + 8 + stride * brick_count : stride]
    return len(slice_) == brick_count and slice_.count(4) == brick_count


def _find_node_offset(data, node_count, plate_count, brick_count):
    """Locate the start of the node block.

    Strategy:
      1. Find the `Brick Property 1` marker (end of brick block).
      2. The brick block consists of 4-node tets (25 bytes each), so its
         start = marker − 25 × brick_count (or × (brick_count − 1) if the
         file is one of the off-by-one cases).
      3. Walk forward through plate records: from any candidate plate_start,
         the right one parses exactly `plate_count` records and lands
         precisely at the brick-block start.
      4. node_offset = plate_start − node_count × 24.

    This avoids scanning the 88 KB pre-mesh region byte by byte.
    """
    block_size = node_count * 24
    end = len(data)

    prop_marker = b"\x04\x01\x00\x10\x00Brick Property 1"
    prop_off = data.find(prop_marker)

    if prop_off != -1:
        for bc_guess in (brick_count, brick_count - 1):
            brick_start = prop_off - 25 * bc_guess
            if brick_start < 0:
                continue
            if not _check_brick_stride(data, brick_start, bc_guess):
                continue
            # Find plate_start: walk backwards from brick_start through
            # likely plate-region offsets. Plates are mostly 3-node tris
            # (21 bytes); allow a per-plate range of 13..100 bytes.
            # We expect: walking plate_count records from plate_start ends
            # at brick_start.
            avg = 21
            lo = brick_start - 100 * plate_count
            hi = brick_start - 13 * plate_count
            # Start at the most-likely offset and spiral out, so common
            # cases (3-node tris) hit immediately.
            est = brick_start - avg * plate_count
            tried = set()
            for off in [est] + [est + d for d in range(-2048, 2049, 1)
                                if est + d != est]:
                if off < lo or off > hi or off in tried:
                    continue
                tried.add(off)
                got = _walk_plates(data, off, plate_count, node_count)
                if got == brick_start:
                    plate_start = off
                    node_offset = plate_start - block_size
                    if node_offset < 0:
                        continue
                    # Quick sanity check on the very first node's coords.
                    try:
                        x = _f64(data, node_offset)
                        if math.isfinite(x) and abs(x) < 1e30:
                            return node_offset
                    except struct.error:
                        pass

    # ----- FALLBACK: linear scan (only reached for unrecognised files) ----
    limit = end - block_size
    for i in range(0, limit, 2):
        try:
            v1 = _f64(data, i)
            v2 = _f64(data, i + block_size - 24)
            if not (math.isfinite(v1) and abs(v1) < 1e30
                    and math.isfinite(v2) and abs(v2) < 1e30):
                continue
        except struct.error:
            continue
        plate_end = _walk_plates(data, i + block_size, plate_count, node_count)
        if plate_end is None:
            continue
        # Found plate_start; now check brick alignment.
        if not _check_brick_stride(data, plate_end, brick_count - 1):
            if not _check_brick_stride(data, plate_end, brick_count):
                continue
        return i
    return None


# ---------------------------------------------------------------------------
# Element parsing (with safe brick truncation)
# ---------------------------------------------------------------------------

def _parse_element_block(data, start, count, max_end):
    """Parse `count` element records starting at `start`, but stop early if we
    hit a property-block header (which means `count` over-states the real
    number of elements by some amount — see ST7_FORMAT.md §5).

    Each record:
        uint32 type, uint32 prop_id, uint8 num_nodes, uint32 node_id[num_nodes]
    Returns (elements, new_offset, truncated_at_idx_or_None).
    """
    elements = []
    curr = start
    for idx in range(count):
        if curr + 9 > max_end:
            return elements, curr, idx
        # Heuristic: if the next 12 bytes match a property-block header,
        # stop. The `01 00` at offset 1-2 is the property "sub-kind" and the
        # ASCII "Prope" never appears in element node IDs.
        if (curr + 16 < max_end
                and data[curr + 5:curr + 13] == b"\x00\x10\x00Brick"
                and data[curr + 1:curr + 3] == b"\x01\x00"):
            return elements, curr, idx
        type_id = _u32(data, curr)
        prop_id = _u32(data, curr + 4)
        n = _u8(data, curr + 8)
        if n == 0 or n > 32:
            return elements, curr, idx
        if curr + 9 + n * 4 > max_end:
            return elements, curr, idx
        nodes = struct.unpack_from(f'<{n}I', data, curr + 9)
        # Reject obviously bad node IDs (cheap guard against misalignment).
        if any(nid == 0 for nid in nodes):
            return elements, curr, idx
        elements.append((type_id, prop_id, n, nodes))
        curr += 9 + n * 4
    return elements, curr, None


# ---------------------------------------------------------------------------
# Materials / properties / loads / groups
# ---------------------------------------------------------------------------

def _scan_post_mesh_strings(data, start):
    """Walk the post-mesh region and return a list of every plausible
    length-prefixed ASCII string with its file offset and length."""
    out = []
    i = start
    end = len(data)
    while i < end - 2:
        L = _u16(data, i)
        if 3 <= L <= 200 and i + 2 + L <= end:
            chunk = data[i + 2:i + 2 + L]
            if all(32 <= b < 127 for b in chunk) and any(
                    c.isalpha() for c in chunk.decode('ascii', 'ignore')):
                out.append((i, L, chunk.decode('ascii')))
                i += 2 + L
                continue
        i += 1
    return out


def _parse_properties(data, start):
    """Extract the property records that follow the brick block.

    The decoded layout is approximate (see ST7_FORMAT.md §6); we just pull
    name, property-id, RGBA colour, material name, and the trailing
    (scale, Young's modulus, Poisson's ratio) doubles that follow the
    `0B 00 40 20 00 …` marker.
    """
    props = []
    i = start
    end = len(data)
    while i < end - 16:
        # Each property record begins with `KK 01 00` where KK is the family
        # byte (04/03/02). The next 2 bytes are a string length.
        kind = data[i]
        if kind not in (0x02, 0x03, 0x04):
            i += 1
            continue
        if data[i + 1] != 0x01 or data[i + 2] != 0x00:
            i += 1
            continue
        name_len = _u16(data, i + 3)
        if not (4 <= name_len <= 80) or i + 5 + name_len > end:
            i += 1
            continue
        name = data[i + 5:i + 5 + name_len].decode('ascii', 'ignore')
        if not (name.startswith(("Brick Prop", "Brick prop",
                                 "Plate Prop", "Plate prop",
                                 "Beam Prop",  "Beam prop"))):
            i += 1
            continue

        # Skip the family/version bytes after the name (heuristic — variable).
        cursor = i + 5 + name_len
        prop_id = None
        color = None
        material_name = ""
        young = None
        poisson = None
        thickness_or_scale = None

        # Look for the property_id in the next 16 bytes (uint32 1-based).
        for o in range(cursor, min(cursor + 16, end - 4)):
            v = _u32(data, o)
            if 1 <= v < 256:
                prop_id = v
                color = tuple(data[o + 4:o + 8])
                cursor = o + 8
                break

        # Find the material name (another length-prefixed string within
        # the next 64 bytes).
        for o in range(cursor, min(cursor + 80, end - 2)):
            L = _u16(data, o)
            if 4 <= L <= 60 and o + 2 + L <= end:
                chunk = data[o + 2:o + 2 + L]
                if all(32 <= b < 127 for b in chunk):
                    text = chunk.decode('ascii')
                    if any(c.isalpha() for c in text) and not text.startswith(
                            ("Brick", "Plate", "Beam")):
                        material_name = text
                        cursor = o + 2 + L
                        break

        # The doubles block always starts with a double = 8.0
        # (`40 20 00 00 00 00 00 00` little-endian), preceded by a 2-byte flag
        # field (`0B 00` for the common case, `0B C0` when temperature data is
        # also present). After the 8.0 marker the layout is:
        #   double  8.0                                 ← thickness / scale
        #   double  0.0                                 ← filler
        #   (optional) double  ref_temperature ×2       ← only in some files
        #   double  Young's modulus                     ← typically 2.0e4 (MPa)
        #   double  Poisson ratio                       ← typically 0.3
        # We find the 8.0 anchor, then scan forward for the first double in
        # the plausible Young's-modulus range and assume Poisson immediately
        # follows it.
        anchor = b"\x40\x20\x00\x00\x00\x00\x00\x00"  # double 8.0 LE
        idx = data.find(anchor, cursor, cursor + 200)
        if idx != -1:
            thickness_or_scale = 8.0
            # Doubles after the anchor are NOT necessarily 8-byte aligned to it
            # (the property record header is variable-length). Scan byte-by-byte
            # for the first double that looks like a Young's modulus, then take
            # the immediately-following double as Poisson's ratio.
            scan = idx + 8
            for o in range(scan, min(scan + 64, end - 16)):
                d = _f64(data, o)
                if math.isfinite(d) and 1.0e2 <= d <= 1.0e15:
                    p = _f64(data, o + 8)
                    if math.isfinite(p) and -1.0 <= p <= 0.6:
                        young = d
                        poisson = p
                        break

        kind_name = {0x04: "Brick", 0x03: "Plate", 0x02: "Beam"}[kind]
        props.append({
            "kind": kind_name,
            "name": name,
            "prop_id": prop_id,
            "color_rgba": color,
            "material": material_name,
            "thickness_or_scale": thickness_or_scale,
            "young_modulus": young,
            "poisson_ratio": poisson,
            "file_offset": i,
        })
        i = cursor + 1  # advance past the record we just parsed
    return props


def _parse_loads(data, start):
    """Find all nodal force records (0x42 0x01 0x07 0x03 marker)."""
    marker = b"\x42\x01\x07\x03\x00\x00\x00"
    loads = []
    i = start
    while True:
        pos = data.find(marker, i)
        if pos == -1 or pos + 42 > len(data):
            break
        nid = _u32(data, pos + 7)
        fx = _f64(data, pos + 18)
        fy = _f64(data, pos + 26)
        fz = _f64(data, pos + 34)
        # cheap sanity check on the doubles
        if all(math.isfinite(v) and abs(v) < 1e20 for v in (fx, fy, fz)):
            loads.append((nid, fx, fy, fz))
        i = pos + 1
    return loads


def _parse_paths_and_titles(data, start):
    """Pull the trailing Windows-style file references and the model title."""
    strings = _scan_post_mesh_strings(data, start)
    paths   = [s for _, _, s in strings if 'C:\\' in s or '/' in s and s.endswith(
        ('.lsa', '.lsl', '.txt', '.srf', '.drf'))]
    titles  = [s for _, _, s in strings if 'volume adjusted' in s
                                       or 'morphology' in s
                                       or 'point sampling' in s]
    return paths, titles


# ---------------------------------------------------------------------------
# Surface extraction with consistent outward normals
# ---------------------------------------------------------------------------

# Canonical tet faces — these point OUTWARD when the tet's nodes are listed
# in the standard positive-volume ordering (n[3] is the "apex" opposite face 0).
_TET_FACES = (
    (0, 2, 1),  # face opposite vertex 3
    (0, 1, 3),  # face opposite vertex 2
    (1, 2, 3),  # face opposite vertex 0
    (0, 3, 2),  # face opposite vertex 1
)


def _tet_face_outward(nodes, face_idx, coords):
    """Return the face triangle (3 node IDs) oriented so its normal points
    AWAY from the tet's opposite vertex. `nodes` is the 4 tet vertex IDs,
    `coords` is the global coordinate list (1-based, coords[i] = (x,y,z))."""
    f = _TET_FACES[face_idx]
    a, b, c = nodes[f[0]], nodes[f[1]], nodes[f[2]]
    d = nodes[6 - f[0] - f[1] - f[2]]   # the 4th vertex
    pa, pb, pc, pd = coords[a], coords[b], coords[c], coords[d]
    # face normal = (pb-pa) x (pc-pa)
    ux, uy, uz = pb[0]-pa[0], pb[1]-pa[1], pb[2]-pa[2]
    vx, vy, vz = pc[0]-pa[0], pc[1]-pa[1], pc[2]-pa[2]
    nx = uy*vz - uz*vy
    ny = uz*vx - ux*vz
    nz = ux*vy - uy*vx
    # vector from face centroid to opposite vertex
    wx, wy, wz = pd[0]-pa[0], pd[1]-pa[1], pd[2]-pa[2]
    if nx*wx + ny*wy + nz*wz > 0:
        # normal points TOWARDS the opposite vertex → triangle is wound
        # the wrong way; flip it.
        return (a, c, b)
    return (a, b, c)


def _build_brick_skin(bricks, coords):
    """Given a list of (type, prop, n, nodes) brick records, return a list of
    outward-oriented boundary triangles."""
    face_owner = {}   # sorted-tuple → (tet_index, face_index)
    duplicates = set()
    for ti, (_type, _prop, n, nodes) in enumerate(bricks):
        if n != 4 and n != 10:
            # not a linear/quadratic tet — skip volume-skinning here.
            # (no hex/penta observed in sample data; can be added later)
            continue
        tet = nodes if n == 4 else nodes[:4]
        for fi in range(4):
            face = (tet[_TET_FACES[fi][0]],
                    tet[_TET_FACES[fi][1]],
                    tet[_TET_FACES[fi][2]])
            key = tuple(sorted(face))
            if key in face_owner:
                duplicates.add(key)
            else:
                face_owner[key] = (ti, fi)
    out = []
    for key, (ti, fi) in face_owner.items():
        if key in duplicates:
            continue
        tet = bricks[ti][3][:4]
        out.append(_tet_face_outward(tet, fi, coords))
    return out


# ---------------------------------------------------------------------------
# File-level extractor
# ---------------------------------------------------------------------------

def extract(path):
    """Parse a .st7 file and return a dict with everything we know how to pull."""
    with open(path, 'rb') as f:
        data = f.read()

    if data[0x10:0x10 + len(MAGIC)] != MAGIC:
        raise ValueError(f"{path}: invalid Strand7 magic")

    header = {
        "header_len":  _u32(data, 0x00),
        "build":       _u32(data, 0x04),
        "major":       _u32(data, 0x08),
        "minor":       _u32(data, 0x0C),
        "magic":       data[0x10:0x24].decode('ascii', 'ignore'),
        "node_count":  _u32(data, 0x3B),
        "plate_count": _u32(data, 0x43),
        "brick_count": _u32(data, 0x47),
    }
    nc, pc, bc = header["node_count"], header["plate_count"], header["brick_count"]
    block_size = nc * 24

    node_offset = _find_node_offset(data, nc, pc, bc)
    if node_offset is None:
        raise RuntimeError(f"{path}: could not locate node block")

    # ----- nodes
    coords = [None]  # 1-based indexing
    for i in range(nc):
        off = node_offset + i * 24
        coords.append((_f64(data, off), _f64(data, off + 8), _f64(data, off + 16)))

    # ----- plates
    plate_start = node_offset + block_size
    plates, after_plates, _ = _parse_element_block(
        data, plate_start, pc, len(data))

    # ----- bricks (with safe truncation if header over-counts)
    bricks, after_bricks, truncated_at = _parse_element_block(
        data, after_plates, bc, len(data))

    # ----- post-mesh metadata
    properties = _parse_properties(data, after_bricks)
    loads      = _parse_loads(data, after_bricks)
    strings    = _scan_post_mesh_strings(data, after_bricks)
    paths      = [s for _, _, s in strings if 'C:\\' in s
                                            or s.endswith(('.lsa', '.lsl',
                                                           '.txt', '.srf',
                                                           '.drf'))]
    titles     = [s for _, _, s in strings if 'volume adjusted' in s
                                            or 'morphology' in s]
    skip = {"Global XYZ", "XYZ", "Model"}
    group_names = []
    for _, _, s in strings:
        if s in skip:                          continue
        if s in paths:                         continue
        if s in titles:                        continue
        if any(p["name"] == s for p in properties): continue
        if any(p["material"] == s for p in properties): continue
        if s.startswith(("Brick Prop", "Plate Prop", "Plate prop",
                         "Beam Prop",  "Brick prop", "Beam prop")):
            continue
        if not any(c.isalpha() for c in s):    continue
        # filter random binary-looking blobs (mostly punctuation)
        letters = sum(1 for c in s if c.isalpha())
        if letters < len(s) * 0.4:             continue
        group_names.append(s)

    return {
        "path":         Path(path),
        "header":       header,
        "node_offset":  node_offset,
        "coords":       coords,
        "plates":       plates,
        "bricks":       bricks,
        "brick_truncated_at": truncated_at,
        "properties":   properties,
        "loads":        loads,
        "paths":        paths,
        "titles":       titles,
        "group_names":  group_names,
        "file_size":    len(data),
        "after_bricks": after_bricks,
    }


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

def write_obj(model, out_path):
    """Write surface skin + plate attachment groups as a Wavefront OBJ."""
    coords = model["coords"]
    plates = model["plates"]
    bricks = model["bricks"]

    # Group plates by their `type` field (1..7) — these are the muscle
    # attachment patches in this dataset.
    plates_by_type = defaultdict(list)
    for type_id, _prop, n, nodes in plates:
        if n == 3:
            plates_by_type[type_id].append(tuple(nodes))
        elif n == 4:
            plates_by_type[type_id].append(tuple(nodes))
        elif n >= 6:
            plates_by_type[type_id].append(tuple(nodes[:3]))

    skin = _build_brick_skin(bricks, coords)

    with open(out_path, 'w') as f:
        f.write(f"# Extracted from {model['path'].name}\n")
        f.write(f"# Nodes: {model['header']['node_count']}\n")
        f.write(f"# Bricks: {len(bricks)}  Plates: {len(plates)}  "
                f"Attachment groups: {len(plates_by_type)}\n")

        for x, y, z in coords[1:]:
            f.write(f"v {x} {y} {z}\n")

        # Bone surface (skinned bricks) first
        f.write("o Bone_Surface\n")
        for face in skin:
            f.write("f " + " ".join(str(n) for n in face) + "\n")

        # Each plate type as its own selectable object
        for type_id in sorted(plates_by_type):
            f.write(f"o Attachment_{type_id}\n")
            for face in plates_by_type[type_id]:
                f.write("f " + " ".join(str(n) for n in face) + "\n")
    return out_path


def write_vtu(model, out_path):
    """Write the full solid mesh (all tetrahedra) as a VTK XML Unstructured Grid.

    Cell types we emit:
        VTK_TETRA   = 10  (4-node)
        VTK_QUADRATIC_TETRA = 24  (10-node, mid-side nodes preserved)
    Cell data arrays:
        property_id : per-cell property index
        type_id     : per-cell `type` field from the .st7
    """
    coords = model["coords"]
    bricks = model["bricks"]
    n_nodes = len(coords) - 1
    n_cells = len(bricks)

    # build connectivity / offsets / types
    conn = []
    offsets = []
    types = []
    prop_ids = []
    type_ids = []
    running = 0
    for type_id, prop_id, n, nodes in bricks:
        if n == 4:
            ct = 10
            conn.extend(nid - 1 for nid in nodes)
            running += 4
        elif n == 10:
            ct = 24
            conn.extend(nid - 1 for nid in nodes)
            running += 10
        else:
            # skip unsupported cells but keep arrays consistent
            continue
        offsets.append(running)
        types.append(ct)
        prop_ids.append(prop_id)
        type_ids.append(type_id)

    with open(out_path, 'w') as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="UnstructuredGrid" version="0.1" '
                'byte_order="LittleEndian">\n')
        f.write('  <UnstructuredGrid>\n')
        f.write(f'    <Piece NumberOfPoints="{n_nodes}" '
                f'NumberOfCells="{len(types)}">\n')

        # Points
        f.write('      <Points>\n')
        f.write('        <DataArray type="Float64" NumberOfComponents="3" '
                'format="ascii">\n')
        for x, y, z in coords[1:]:
            f.write(f"          {x} {y} {z}\n")
        f.write('        </DataArray>\n')
        f.write('      </Points>\n')

        # Cells
        f.write('      <Cells>\n')
        f.write('        <DataArray type="Int64" Name="connectivity" '
                'format="ascii">\n          ')
        # Wrap connectivity at a reasonable line width
        for k, v in enumerate(conn):
            f.write(f"{v} ")
            if (k + 1) % 16 == 0:
                f.write("\n          ")
        f.write('\n        </DataArray>\n')

        f.write('        <DataArray type="Int64" Name="offsets" '
                'format="ascii">\n          ')
        for k, v in enumerate(offsets):
            f.write(f"{v} ")
            if (k + 1) % 16 == 0:
                f.write("\n          ")
        f.write('\n        </DataArray>\n')

        f.write('        <DataArray type="UInt8" Name="types" '
                'format="ascii">\n          ')
        for k, v in enumerate(types):
            f.write(f"{v} ")
            if (k + 1) % 32 == 0:
                f.write("\n          ")
        f.write('\n        </DataArray>\n')
        f.write('      </Cells>\n')

        # Cell data
        f.write('      <CellData Scalars="property_id">\n')
        f.write('        <DataArray type="Int32" Name="property_id" '
                'format="ascii">\n          ')
        for k, v in enumerate(prop_ids):
            f.write(f"{v} ")
            if (k + 1) % 32 == 0:
                f.write("\n          ")
        f.write('\n        </DataArray>\n')

        f.write('        <DataArray type="Int32" Name="type_id" '
                'format="ascii">\n          ')
        for k, v in enumerate(type_ids):
            f.write(f"{v} ")
            if (k + 1) % 32 == 0:
                f.write("\n          ")
        f.write('\n        </DataArray>\n')
        f.write('      </CellData>\n')

        f.write('    </Piece>\n')
        f.write('  </UnstructuredGrid>\n')
        f.write('</VTKFile>\n')
    return out_path


def write_txt(model, out_path):
    """Write a human-readable summary of model + simulation parameters."""
    h = model["header"]
    plates_by_type = defaultdict(int)
    plates_by_prop = defaultdict(int)
    for type_id, prop_id, _n, _nodes in model["plates"]:
        plates_by_type[type_id] += 1
        plates_by_prop[prop_id] += 1
    bricks_by_prop = defaultdict(int)
    for _t, prop_id, _n, _nodes in model["bricks"]:
        bricks_by_prop[prop_id] += 1

    loads = model["loads"]
    if loads:
        sx = sum(L[1] for L in loads)
        sy = sum(L[2] for L in loads)
        sz = sum(L[3] for L in loads)
        smag = math.sqrt(sx*sx + sy*sy + sz*sz)
    else:
        sx = sy = sz = smag = 0.0

    with open(out_path, 'w') as f:
        f.write(f"Strand7 model summary — {model['path'].name}\n")
        f.write("=" * 60 + "\n\n")
        f.write("FILE\n")
        f.write(f"  size              : {model['file_size']:,} bytes\n")
        f.write(f"  magic             : {h['magic']!r}\n")
        f.write(f"  Strand7 version   : {h['major']}.{h['minor']} (build {h['build']})\n")
        f.write(f"  node block offset : 0x{model['node_offset']:X}\n")
        f.write(f"  after-bricks      : 0x{model['after_bricks']:X}\n")
        f.write("\nMESH\n")
        f.write(f"  nodes             : {h['node_count']:,}\n")
        f.write(f"  plates (header)   : {h['plate_count']:,}\n")
        f.write(f"  plates (parsed)   : {len(model['plates']):,}\n")
        f.write(f"  bricks (header)   : {h['brick_count']:,}\n")
        f.write(f"  bricks (parsed)   : {len(model['bricks']):,}")
        if model["brick_truncated_at"] is not None:
            f.write(f"   (truncated at idx {model['brick_truncated_at']} — "
                    f"header over-counts by {h['brick_count'] - len(model['bricks'])})")
        f.write("\n")
        f.write(f"  plate types       : {dict(plates_by_type)}\n")
        f.write(f"  plate prop_ids    : {dict(plates_by_prop)}\n")
        f.write(f"  brick prop_ids    : {dict(bricks_by_prop)}\n")

        f.write("\nMATERIALS / PROPERTIES\n")
        if not model["properties"]:
            f.write("  (none decoded)\n")
        for p in model["properties"]:
            f.write(f"  - [{p['kind']}] {p['name']!r}\n")
            f.write(f"      prop_id     = {p['prop_id']}\n")
            f.write(f"      colour RGBA = {p['color_rgba']}\n")
            f.write(f"      material    = {p['material']!r}\n")
            if p['thickness_or_scale'] is not None:
                f.write(f"      thickness/scale = {p['thickness_or_scale']:g}\n")
            if p['young_modulus'] is not None:
                f.write(f"      Young's E   = {p['young_modulus']:g}\n")
            if p['poisson_ratio'] is not None:
                f.write(f"      Poisson ν   = {p['poisson_ratio']:g}\n")

        f.write("\nLOAD CASE (nodal forces)\n")
        f.write(f"  records           : {len(loads):,}\n")
        if loads:
            f.write(f"  sum F             : ({sx:.4e}, {sy:.4e}, {sz:.4e}) N\n")
            f.write(f"  |sum F|           : {smag:.4e} N\n")
            f.write(f"  first 5           :\n")
            for nid, fx, fy, fz in loads[:5]:
                f.write(f"     node {nid:6d}  F=({fx:+.4e}, {fy:+.4e}, {fz:+.4e})\n")

        f.write("\nNAMED GROUPS / LOAD-CASE LABELS\n")
        if not model["group_names"]:
            f.write("  (none decoded — file may use a single anonymous case)\n")
        for g in model["group_names"]:
            f.write(f"  - {g}\n")

        f.write("\nSOURCE TITLES\n")
        for t in model["titles"]:
            f.write(f"  - {t}\n")

        f.write("\nRESULT-FILE REFERENCES (paths embedded by Strand7)\n")
        for p in model["paths"]:
            f.write(f"  - {p}\n")
    return out_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Extract mesh & metadata from a Strand7 (.st7) file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument('file', help='Path to .st7 file')
    ap.add_argument('outdir', help='Output directory')
    ap.add_argument('--obj', action='store_true',
                    help='Write surface skin OBJ (bone surface + attachments)')
    ap.add_argument('--vtu', action='store_true',
                    help='Write full solid mesh as VTU for ParaView')
    ap.add_argument('--txt', action='store_true',
                    help='Write human-readable parameter summary')
    ap.add_argument('--all', action='store_true',
                    help='Write all output formats (equivalent to --obj --vtu --txt)')
    args = ap.parse_args(argv)

    if args.all:
        args.obj = args.vtu = args.txt = True
    if not (args.obj or args.vtu or args.txt):
        ap.error("nothing to do — pass at least one of --obj / --vtu / --txt / --all")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[*] Reading {args.file}")
    model = extract(args.file)
    h = model["header"]
    stem = model["path"].stem

    print(f"    nodes={h['node_count']:,}  plates={len(model['plates']):,}  "
          f"bricks={len(model['bricks']):,}  "
          f"properties={len(model['properties'])}  loads={len(model['loads']):,}")
    if model["brick_truncated_at"] is not None:
        print(f"    [note] truncated brick block at idx {model['brick_truncated_at']}"
              f"  (header over-counts by "
              f"{h['brick_count'] - len(model['bricks'])})")

    if args.obj:
        p = write_obj(model, outdir / f"{stem}.obj")
        print(f"[+] wrote {p}")
    if args.vtu:
        p = write_vtu(model, outdir / f"{stem}.vtu")
        print(f"[+] wrote {p}")
    if args.txt:
        p = write_txt(model, outdir / f"{stem}.txt")
        print(f"[+] wrote {p}")


if __name__ == '__main__':
    main()
