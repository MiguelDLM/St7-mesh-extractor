# Strand7 (.st7) Binary File Format — Reverse-Engineered Notes

> **Status: unofficial / work-in-progress.**
> Strand7 (now Straus7) does not publish the binary layout of its native
> `.st7` model files. The notes below were produced by inspecting saved
> models with a hex editor and cross-checking the parsed values against
> what the Strand7 GUI exposes. They cover the bytes we need in order to
> extract the geometric mesh, the per-element material assignments, and
> the load / freedom / metadata blocks that follow.
>
> If you only need the *results* of a solved model, use the companion
> `.lsa` / `.lsl` / `.drf` files written by the solver — those are simpler
> to parse and are referenced from inside every `.st7`.
>
> Sections marked **(unknown)** are blocks we can locate but have not
> fully decoded. Pull requests welcome.

All multi-byte values are **little-endian**. Strings are length-prefixed
with a `uint16` length and stored without a NUL terminator.

---

## 1. File-level layout

```
┌────────────────────────────────────────────────────────────────┐
│  HEADER + saved GUI state             (~88 KB, mostly fixed)   │
├────────────────────────────────────────────────────────────────┤
│  NODE BLOCK            node_count × 24 bytes  (X, Y, Z doubles)│
├────────────────────────────────────────────────────────────────┤
│  PLATE BLOCK           plate_count × variable element records  │
├────────────────────────────────────────────────────────────────┤
│  BRICK BLOCK           brick_count × variable element records  │
├────────────────────────────────────────────────────────────────┤
│  PROPERTY DEFINITIONS  per-property name, material, E, ν, …    │
├────────────────────────────────────────────────────────────────┤
│  FREEDOM CASE / RESTRAINTS    (per-node — ~25 B/record)        │
├────────────────────────────────────────────────────────────────┤
│  LOAD CASE / NODAL FORCES     (42 B/record, X/Y/Z component)   │
├────────────────────────────────────────────────────────────────┤
│  COORDINATE SYSTEMS, GROUPS, RESULT-FILE PATHS, TITLE          │
└────────────────────────────────────────────────────────────────┘
```

Each major block is contiguous; the blocks appear in the order listed
above. The reference implementation in [`st7_extractor.py`](st7_extractor.py)
locates them by:

1. Reading the element counts at fixed offsets in the header.
2. Finding the *first* property record (`04 01 00 10 00 Brick Property 1`)
   via byte-pattern search — this is the end of the brick block.
3. Walking backwards through the 25-byte 4-node tetrahedron records to
   recover the brick-block start.
4. Walking the plate records forward from successive candidate offsets
   until exactly `plate_count` of them land at the brick-block start.
5. Subtracting `node_count × 24` to get the node-block start.

This is robust against the 88 KB of variable GUI state at the top of the
file (window sizes, font names, saved colours, last-opened result paths)
which would otherwise force a linear scan.

---

## 2. Fixed header (offsets 0x00 – 0x3A)

| Offset | Size | Field                | Notes                                |
|--------|------|----------------------|--------------------------------------|
| 0x00   | 4    | header length / version marker | Observed: `0x20` (= 32 bytes? or version tag) |
| 0x04   | 4    | build number         | Identifies the Strand7 build that wrote the file |
| 0x08   | 4    | major version        | Observed: `7`                        |
| 0x0C   | 4    | minor version        | Observed: `3`                        |
| 0x10   | 20   | **magic string**     | ASCII `STRAUS/STRAND (c)G+D` (use this as the file-type check) |
| 0x24   | 4    | constant `0x13` = 19 | Unknown purpose                      |
| 0x28   | 4    | constant `0x5A` = 90 | Unknown purpose                      |
| 0x2C   | 4    | reserved             | Zero in every sample                 |
| 0x30   | 4    | `1400`               | Looks like a saved window pixel width|
| 0x34   | 4    | `958`                | Looks like a saved window pixel height|
| 0x38   | 3    | constant `01 03 00`  | Unknown                              |

The major / minor version pair and the magic together identify a
Strand7 model file. Reject anything that does not match the magic.

---

## 3. Variable header fields (offsets 0x3B onward)

| Offset | Type    | Field         | Description                                |
|--------|---------|---------------|--------------------------------------------|
| 0x3B   | uint32  | `node_count`  | Total number of nodes                      |
| 0x3F   | uint32  | reserved (0)  |                                            |
| 0x43   | uint32  | `plate_count` | Total number of plate (2-D / shell) elements |
| 0x47   | uint32  | `brick_count` | Total number of brick (3-D solid) elements |
| 0x4B   | uint32  | reserved (0)  |                                            |
| 0x4F   | uint32  | reserved (0)  |                                            |
| 0x53   | uint32  | _count?_      | Often `plate_count + 19` but the relation is not universal — possibly group or property-set count |

### Header doubles (around 0x9B – 0xBF)

Three double-precision values appear in this region; their meaning is
not confirmed but they are constant within a file:

| Offset | Type   | Observed behaviour                              | Guess                              |
|--------|--------|-------------------------------------------------|------------------------------------|
| 0x9B   | double | Per-file constant, in the 10³–10⁵ range         | Bounding-box size or unit scale    |
| 0xA3   | double | Nearly identical across files (~`4.06e+04`)     | Possibly an engine-version stamp   |
| 0xAB   | uint32 | Constant `4`                                    | Unknown                            |
| 0xB7   | double | Duplicate of the value at 0x9B                  | Same as 0x9B                       |

### Pre-mesh region (0xC0 – start of node block)

Approximately 88 KB of mostly-zero space punctuated by:

* Length-prefixed display-font names: `Small Fonts`, `Arial`,
  `Arial Unicode MS`, `Tahoma`.
* Repeated `ColorA=FFFFFF`, `ColorB=FFFFFF`, … (GUI palette).
* Absolute Windows paths to the source model and any solved result files
  (`.lsa`, `.lsl`, `.txt`, `.srf`, `.drf`) saved alongside the model.

The size of this region varies by a few hundred bytes depending on how
many result paths Strand7 has remembered. For a robust parse, do not
hard-code the node-block start — locate it from the property block at
the end of the brick block (see §1, step 2).

---

## 4. Node block

Begins at the offset discovered from the property-block back-walk.
Layout is dense:

```
nodes[i] = (x: double, y: double, z: double)   # 24 bytes
```

`node_count` consecutive records, no per-record header. Node IDs are
implicit (1-based), matching the Strand7 GUI.

---

## 5. Element blocks (plates, then bricks)

Both element blocks share the same per-record layout. Records are
variable-length:

```
element_record:
    uint32  type_or_kind          # appears constant within a file (see notes below)
    uint32  property_id           # 1-based index into the property table
    uint8   num_nodes             # 3, 4, 6, 8 for plates; 4, 8, 10, 15, 20 for bricks
    uint32  node_id[num_nodes]    # 1-based node indices
```

Each record is `9 + 4 × num_nodes` bytes long.

### Notes / unknowns

* **`type_or_kind`** is one tiny integer (1–7 in every file we've
  examined) and is constant within a given file. It is **not** an
  element ID — element IDs are implicit positional (the *k*-th plate
  record is element `k+1`). Within a single plate block we have seen
  `type_or_kind` take 1–7 different values, and exporting one OBJ
  object per value cleanly separates the underlying "regions" (e.g.
  per-attachment patches, per-surface markings).
* In some files, the brick block consumes one more record than the
  `brick_count` header implies, *or* the last record is silently
  truncated. If you trust `brick_count`, watch for a record whose
  `num_nodes` is implausible (e.g. > 32) — that is the start of the
  property block bleeding into the loop. The reference implementation
  detects this by walking until it sees the `04 01 00 ... Brick`
  property header, then stops one record early.
* No quadratic plate / brick variants have been observed in the sample
  set, but the per-record layout naturally supports any `num_nodes`
  value.

---

## 6. Property / material block

Begins immediately after the brick block. Repeating records, one per
defined property. The decoded layout is approximate:

```
property_record:
    uint8   record_kind            # 0x04 = brick prop, 0x03 = plate prop, 0x02 = beam prop
    uint8   sub_kind               # always 0x01
    uint8   reserved               # 0x00
    uint16  name_length
    char    name[name_length]      # e.g. "Brick Property 1"
    uint8   ?                      # 0x01
    uint8   ?                      # 0x07
    uint32  property_id            # matches element_record.property_id
    uint8   colorR, colorG, colorB, colorA   # display colour
    … ~5 unknown bytes …
    uint16  material_name_length
    char    material_name[…]       # e.g. "Steel", "Bone", project-specific names
    … ~13 unknown bytes (control / flags) …
    [optional uint16 length-marker, depends on `sub_kind` flags]
    double  scale_or_thickness     # in our samples always 8.0
    double  filler / 0.0
    (optional) double × 2          # reference temperature(s) for thermal materials
    double  young_modulus          # MPa
    double  poisson_ratio
    … trailing flags / counts …
```

The first record in the block is always `Brick Property 1`; that makes
it a convenient marker for locating the end of the brick block (see §1
step 2 and §5).

For `record_kind = 0x02` (beam) the trailing doubles encode the
section profile and not material constants, but the layout up to the
material name is the same.

---

## 7. Freedom case / nodal restraints

Begins shortly after the property block. Repeating records, roughly
25 bytes each. Approximate layout:

```
restraint_record:
    uint32  node_id
    uint32  flag                       # 0x01
    uint8   mask[8]                    # 0xFF×8 sentinel for "free"
    uint16  ?                          # 0x0002
    uint8   ?, ?                       # 0x01 0x01
    uint8   dof_mask                   # bit-field of restrained DOFs
                                       #   bit 0 = DX, bit 1 = DY, bit 2 = DZ
                                       #   bit 3 = RX, bit 4 = RY, bit 5 = RZ
                                       #   (mapping not officially confirmed)
    uint8   ?                          # 0x01
    uint16  ?                          # 0x0037
    uint16  ?                          # 0x0019
    uint16  coordinate_system_id?      # 0 = global XYZ
```

Records cover only the nodes that have at least one restrained DOF;
free nodes are omitted.

---

## 8. Load case / nodal forces

Easy to identify by the recurring 7-byte prefix `42 01 07 03 00 00 00`.
Each record is 42 bytes:

```
force_record:
    uint8   marker[7] = 42 01 07 03 00 00 00
    uint32  node_id
    uint32  flag1                      # observed = 2
    uint16  flag2                      # observed = 1794 (0x0702)
    double  Fx
    double  Fy
    double  Fz
                                       # 7 + 4 + 4 + 2 + 24 = 41
                                       # in practice 42 (1 byte of padding)
```

Records are one per loaded node. The vector sum across all records
should equal the total applied load expected for the model.

---

## 9. Coordinate systems, groups, metadata, result-file paths

The tail of the file contains, in approximate order:

* **Coordinate system** definitions — the default `Global XYZ`, plus
  any user-defined systems. Each has a name and nine doubles encoding
  three unit vectors.
* **Group / named-selection tree** — the root group is always called
  `Model`. Children include user-defined groups (which may correspond
  to anatomical regions, named load configurations, etc.).
* **Load case** and **Freedom case** labels (`Load Case 1`,
  `Freedom Case 1`, …).
* A free-text **model title** (one length-prefixed string).
* Up to ~5 length-prefixed **absolute Windows paths** to companion
  result files written by the solver:
  * `.lsa` — Linear Static Analysis result archive
  * `.lsl` — solver log
  * `.txt` — text dump of selected results
  * `.srf` — surface-pressure file
  * `.drf` — displacement results

The internal structure of the group tree and of each coordinate-system
record is only partially decoded; the reference extractor just pulls
the labels and paths.

---

## 10. What the reference extractor reads vs. ignores

| Section                                | Reference extractor | Notes |
|----------------------------------------|---------------------|-------|
| Fixed header (magic, version)          | ✅ (magic + version) | Reads version for the TXT summary |
| Pre-mesh GUI state / paths             | ❌                  | Could be parsed for provenance |
| Node coordinates                       | ✅                  | All `node_count × (x,y,z)` doubles |
| Plate elements (connectivity)          | ✅                  | Quadratic mid-side nodes dropped if present |
| Brick elements (connectivity)          | ✅                  | Output as VTU cells (type 10/24) for the volume mesh; bricks are also surface-skinned for the OBJ |
| Property → material assignments        | ✅                  | Name, RGBA, material name |
| Material constants (E, ν, thickness)   | ✅                  | Heuristic — see §6 |
| Freedom case (nodal restraints)        | ⚠️ partial          | Counted but DOF mask not interpreted |
| Load case (nodal forces Fx, Fy, Fz)    | ✅                  | Sum and per-record listing emitted in TXT |
| Named groups / load-case labels        | ✅                  | Names only — group membership not yet parsed |
| Coordinate systems                     | ❌                  | Only the default is meaningful in the sample data |
| Result-file paths                      | ✅                  | Listed in the TXT summary |

---

## 11. Conventions for surface extraction

The bricks in the sample data are all linear tetrahedra. To recover
the outer surface (skin):

1. For each tet `(n0, n1, n2, n3)` emit the four candidate faces:
   `(n0, n2, n1)`, `(n0, n1, n3)`, `(n1, n2, n3)`, `(n0, n3, n2)`.
   The orderings above are the *outward-pointing* triangulation
   when the tet has positive signed volume.
2. Sort each face's vertex IDs and use the result as a key.
3. Faces that appear in exactly **one** tet are boundary faces.
   Faces that appear in two tets are interior and should be discarded.
4. For each surviving boundary face, recompute the normal and ensure
   it points away from the opposing vertex of its owning tet. Flip
   the triangle if not — this corrects for tets stored with negative
   signed volume.

The reference extractor implements this in `_build_brick_skin` /
`_tet_face_outward`.

---

## 12. Open questions

1. The 4-byte `type_or_kind` at the start of each element record is
   constant within a file but takes different values across files.
   Does Strand7 use this as an "element family" tag, a mesh-generation
   revision, or something else? It correlates with how the GUI groups
   plates into selectable regions but we have no documentation.
2. Why does the brick block in some files consume one record more
   than `brick_count`? Candidates: a per-block terminator we are not
   consuming, a padding byte we are not skipping, or a header field
   that holds the *highest* element ID rather than the count.
3. The block of ~13 control bytes between the material name and the
   doubles in §6 is not fully decoded. The first double of each
   property record is invariably `8.0` in our samples — confirming
   the field exists, but not what it actually represents (plate
   thickness for plate records seems likely, but it appears unchanged
   in brick records too).
4. The freedom-case DOF mask bits (§7) should be cross-referenced
   against the Strand7 *Restraint* dialog to confirm which bit
   corresponds to which translational / rotational DOF.
5. The two header doubles at 0x9B / 0xB7 are identical within a file
   but vary between files. They might be a model-wide length / unit
   scale, or simply a saved camera distance from the last GUI session.

---

## 13. Reference and acknowledgements

This document was produced by reverse-engineering a corpus of `.st7`
files. None of the information here came from official Strand7 / G+D
sources. If you have access to the Strand7 API documentation
(`St7API.h`) and can confirm or refute any of the field meanings
above, please open an issue or PR.

The companion `.st7` model files used to derive these notes are
available from their original publishing authors and are not included
in this repository.
