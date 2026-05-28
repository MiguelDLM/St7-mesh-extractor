# Strand7 (.st7) Binary File Format — Reverse-Engineered Notes

> **Status: unofficial.**
> Strand7 (now Straus7) does not publish the binary layout of its native
> `.st7` model files, but it does ship a text-format equivalent that the
> GUI can export from any model. We have used that text export as ground
> truth to cross-check the field meanings below: every field marked
> ✅ has been confirmed by comparing the binary record against the
> corresponding row in Strand7's own `.txt` export of the same model.
>
> If you only need post-solve results (stress, displacement, …) use the
> companion `.lsa` / `.lsl` / `.drf` files instead — those are simpler.
>
> Sections marked **(unknown)** are blocks we can locate but have not
> fully decoded. Pull requests welcome.

All multi-byte values are **little-endian**. Strings are length-prefixed
with a `uint16` length and stored without a NUL terminator.

The companion text-export reference is also documented here in the same
terms Strand7 uses (`PlateShellProp`, `BrickProp`, `NdFreedom`,
`PlGlobalLoad`, `Modulus`, `Poisson`, `MemThick`, …) — see §11.

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
│  PROPERTY DEFINITIONS  PlateShellProp / BrickProp / BeamProp   │
├────────────────────────────────────────────────────────────────┤
│  FREEDOM CASE / RESTRAINTS  (NdFreedom — per node, per case)   │
├────────────────────────────────────────────────────────────────┤
│  LOAD CASE / PLATE FACE LOADS  (PlGlobalLoad — per plate)      │
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
| 0x00   | 4    | header length / version marker | Observed: `0x20` = 32. Possibly the size of the fixed header block. |
| 0x04   | 4    | build number         | Identifies the Strand7 build that wrote the file |
| 0x08   | 4    | major version        | Observed: `7`                        |
| 0x0C   | 4    | minor version        | Observed: `3`                        |
| 0x10   | 20   | **magic string**     | ASCII `STRAUS/STRAND (c)G+D` (use this as the file-type check) |
| 0x24   | 4    | constant `0x13` = 19 | Unknown                              |
| 0x28   | 4    | constant `0x5A` = 90 | Unknown                              |
| 0x2C   | 4    | reserved             | Zero in every sample                 |
| 0x30   | 4    | `1400`               | Looks like a saved window pixel width|
| 0x34   | 4    | `958`                | Looks like a saved window pixel height|
| 0x38   | 3    | constant `01 03 00`  | Unknown                              |

The major / minor version pair and the magic together identify a
Strand7 model file. Reject anything that does not match the magic.

The corresponding text export has `FileFormat Straus7.<major>.<minor>.<patch>`
at the top, so the version pair here can be cross-referenced.

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
| 0x53   | uint32  | unknown count | Often `plate_count + 19` in sample data but the relation is not universal |

### Header doubles (around 0x9B – 0xBF)

Three double-precision values appear in this region; their meaning is
not confirmed but they are constant within a file:

| Offset | Type   | Observed behaviour                              | Guess                              |
|--------|--------|-------------------------------------------------|------------------------------------|
| 0x9B   | double | Per-file constant, in the 10³–10⁵ range         | Bounding-box size or saved camera distance |
| 0xA3   | double | Nearly identical across files (~`4.06e+04`)     | Possibly an engine-version stamp   |
| 0xAB   | uint32 | Constant `4`                                    | Unknown                            |
| 0xB7   | double | Duplicate of the value at 0x9B                  | Same as 0x9B                       |

### Units

**Units are NOT stored in the binary header** (or at least not in a way
we have recovered). Strand7's GUI lets the user pick units when opening
the file. The corresponding text export *does* carry a `UNITS` block
with explicit names:

```
LengthUnit           mm
MassUnit             kg
EnergyUnit           J
PressureUnit         MPa
ForceUnit            N
TemperatureUnit      C
```

In practice the sample data is consistently mm / kg / N / MPa and the
reference extractor reports those as the assumed units in its TXT
summary; if you process a model with different conventions you'll need
to override.

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

In the text export this maps directly to:
```
Node  <id>   <x>   <y>   <z>
```

---

## 5. Element blocks (plates, then bricks)

Both element blocks share the same per-record layout. Records are
variable-length:

```
element_record:
    uint32  group_id            ✅ matches the 3rd column of `Tri3 / Quad4 /
                                   Tetra4 / Hexa8` rows in the text export
    uint32  property_id         ✅ matches the 2nd column (PropertyId)
                                   → references PlateShellProp / BrickProp
    uint8   num_nodes           ✅ 3, 4, 6, 8 for plates;
                                   4, 8, 10, 15, 20 for bricks
    uint32  node_id[num_nodes]  ✅ 1-based node indices
```

Each record is `9 + 4 × num_nodes` bytes long. The element ID itself
is **implicit positional** (the *k*-th record is element `k+1`).

### group_id

The `group_id` field selects which group (named in the GROUP DEFINITIONS
block) the element belongs to. Strand7 organises groups as a tree
rooted at `Model` (always group 1). In the text export this looks like:

```
Group  1   16711680   "\\Model"
Group  2    3355647   "cranium"
Group  3   16757299   "frontal dome"
```

In files where the modeller has used sub-groups to flag muscle
attachment patches or other regions of interest, the OBJ writer emits
one selectable `Plate_Group_N` object per `group_id`.

A note on persistence: groups that have been deleted from a model
sometimes still appear referenced from element records but no longer
have a name in GROUP DEFINITIONS. The reference extractor preserves
the numeric `group_id` regardless.

### Observations across sample data

* All bricks observed are 4-node linear tetrahedra (25 bytes per record,
  num_nodes = 4). The fixed-stride layout is what makes the back-walk
  from the property marker in §1 step 3 work.
* No quadratic plate / brick variants have been observed in the sample
  set, but the per-record layout naturally supports any `num_nodes`
  value.
* Earlier observation: "in some files the brick block consumes one
  more record than `brick_count` implies." With the corrected node-block
  detection this is no longer seen — those failures were caused by a
  wrong `node_offset` overshooting into the plate / brick blocks, not
  by any genuine over-count in the header.

---

## 6. Property / material block

Begins immediately after the brick block. One record per defined
property. The text-export equivalents are:

```
PlateShellProp   <id>   "<name>"
  MaterialName   "<material>"
  Modulus        <E in MPa>
  Poisson        <ν>
  MemThick       <membrane thickness>
  BendThick      <bending thickness>
  NonLinType     <Elasticplastic | ...>
  YieldCriterion <VonMises | ...>
  NumLayers      <integer>

BrickProp        <id>   "<name>"
  MaterialName   "<material>"
  Modulus        <E>
  Poisson        <ν>
  NonLinType     <...>
  YieldCriterion <...>
```

Approximate binary layout:

```
property_record:
    uint8   record_kind             ✅  0x04 = BrickProp
                                        0x03 = PlateShellProp
                                        0x02 = BeamProp
    uint8   sub_kind                    always 0x01
    uint8   reserved                    0x00
    uint16  name_length
    char    name[name_length]       ✅  e.g. "Brick Property 1"
    uint8   ?                           0x01
    uint8   ?                           0x07
    uint32  property_id             ✅  matches element_record.property_id
    uint8   color_R, color_G, color_B, color_A   ✅  display colour
    … ~5 unknown bytes …
    uint16  material_name_length
    char    material_name[…]        ✅  e.g. "Steel", "Bone",
                                        project-specific material names
    … ~13 unknown bytes (control / flags) …
    [optional length-marker, 0B 00 or 0B C0 — see note]
    double  constant 8.0                purpose unknown; appears in every record
    double  filler / 0.0
    (optional) double × 2               reference temperatures, only when
                                        the 0B C0 flag is set
    double  Modulus                 ✅  Young's modulus (MPa in sample data)
    double  Poisson                 ✅  Poisson's ratio
    … 10 unknown bytes …
    uint32  NumLayers               ✅  only present in PlateShellProp records
                                        (Strand7 default = 10)
    … 2 unknown bytes …
    double  MemThick                ✅  PlateShellProp only — membrane thickness
    double  BendThick               ✅  PlateShellProp only — bending thickness
    … trailing flags …
```

For `record_kind = 0x02` (Beam) the trailing doubles encode the section
profile and not material constants, but the layout up to the material
name is the same.

The leading constant double = 8.0 has been observed in every property
record in every sample file. We do not know what it represents — it
does not match any of the values exposed in the corresponding TXT
export.

---

## 7. Freedom case / nodal restraints

Begins shortly after the property block. The text-export equivalent is:

```
NdFreedom  <freedom_case_id>  <node_id>  <kind>  DX DY DZ RX RY RZ
```

(only the restrained DOFs are listed; omitted ones are free.)

Approximate binary layout (one record per (freedom_case, node) pair):

```
restraint_record:
    uint32  node_id
    uint32  freedom_case_id / flag        observed = 0x01 in single-case files
    uint8   mask[8]                       0xFF×8 sentinel for "free"
    uint16  ?                             0x0002
    uint8   ?, ?                          0x01 0x01
    uint8   dof_mask                      bit-field of restrained DOFs
                                          bits map to DX, DY, DZ, RX, RY, RZ
                                          (mapping not yet confirmed against
                                          the GUI's bit order)
    uint8   ?                             0x01
    uint16  ?                             0x0037
    uint16  ?                             0x0019
    uint16  coordinate_system_id?         0 = global XYZ
```

Records cover only the nodes that have at least one restrained DOF;
free nodes are omitted.

---

## 8. Load case / plate face loads

Easy to identify by the recurring 7-byte prefix `42 01 07 03 00 00 00`.
Each record is 42 bytes:

```
plate_face_load_record:
    uint8   marker[7] = 42 01 07 03 00 00 00
    uint32  plate_element_id            ✅  matches the 2nd column of
                                            PlGlobalLoad rows in the text
                                            export
    uint32  flag1                           observed = 2
    uint16  flag2                           observed = 1794 (0x0702)
    double  Fx                          ✅
    double  Fy                          ✅
    double  Fz                          ✅
                                            # 7 + 4 + 4 + 2 + 24 = 41
                                            # in practice 42 (1 byte of padding)
```

In Strand7 terminology this is `PlGlobalLoad <load_case_id>
<plate_element_id> <Fx> <Fy> <Fz>` — a global-coordinate force vector
applied to the face of plate element *N*. (Plates in this dataset are
all `Tri3`, so face-1 is unambiguous; for quadrilateral or higher-order
plates Strand7 also writes a face number, but we have not yet seen one
in the binary.)

The binary records do not appear to carry the parent load-case ID
explicitly — they are grouped by position in the file (each `Load Case`
section in the text export corresponds to a contiguous run of these
records).

The reference extractor reports the sum of all force vectors so you can
cross-check it against the total load you expect from the model
(e.g. a ~40 kN bite force).

---

## 9. Coordinate systems, groups, metadata, result-file paths

The tail of the file contains, in approximate order:

* **Coordinate system** definitions — the default `Global XYZ`, plus
  any user-defined systems. Each has a name and nine doubles encoding
  three unit vectors (`CoordSys <id> "<name>" GlobalXYZ` or
  `CoordSys <id> "<name>" Cartesian / Cylindrical / Spherical`).
* **Group / named-selection tree** — the root group is always called
  `Model`. Children correspond to user-defined groups (anatomical
  regions, material zones, etc.).
* **Freedom case** and **Load case** labels (`Freedom Case 1`,
  `Load Case 1`, …). These also appear as the headers of the
  NdFreedom / PlGlobalLoad sections in the text export.
* A free-text **model title / project / author** trio (each a length-prefixed
  string). These map to `ModelName`, `Title`, `Project`, `Author`,
  `Reference`, `Comments` in the text export.
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
| Plate elements (connectivity, group, prop) | ✅              | Quadratic mid-side nodes dropped if present |
| Brick elements (connectivity, group, prop) | ✅              | Output as VTU cells (type 10 / 24) for the volume mesh; bricks are also surface-skinned for the OBJ |
| Property → material assignments        | ✅                  | Name, RGBA, material name |
| Material constants (Modulus, Poisson, MemThick, BendThick, NumLayers) | ✅ | NonLinType / YieldCriterion not yet parsed |
| Freedom case (nodal restraints)        | ⚠️ partial          | Counted but DOF mask not interpreted |
| Load case (PlGlobalLoad: Fx, Fy, Fz)   | ✅                  | Sum and per-record listing emitted in TXT |
| Named groups / load-case labels        | ✅                  | Names only — group membership not yet parsed |
| Coordinate systems                     | ❌                  | Only the default `Global XYZ` is meaningful in the sample data |
| Result-file paths                      | ✅                  | Listed in the TXT summary |

---

## 11. Strand7 text-export reference

If you have access to a Strand7 installation you can write a model out
in text format (`File → Save As → *.txt`). This is the closest thing to
official documentation for the field names used above. Section order
in the text file mirrors the binary block order, which is what we used
to cross-check our parser:

| Text section                       | Binary block        |
|------------------------------------|---------------------|
| `MODEL INFORMATION`                | Header §2 + tail metadata §9 |
| `UNITS`                            | (not in binary — assumed) |
| `GROUP DEFINITIONS`                | Tail §9 (only names) |
| `FREEDOM CASE DEFINITIONS`         | Tail §9             |
| `LOAD CASE DEFINITIONS`            | Tail §9             |
| `COORDINATE SYSTEM DEFINITIONS`    | Tail §9             |
| `NODE COORDINATES`                 | Node block §4       |
| `PLATE ELEMENTS`                   | Plate block §5      |
| `BRICK ELEMENTS`                   | Brick block §5      |
| `NODE RESTRAINTS (ROTATION AS RADIAN)` | Freedom-case block §7 |
| `PLATE FACE GLOBAL LOADS`          | Load-case block §8  |
| `PLATE PROPERTIES`                 | Property block §6   |
| `BRICK PROPERTIES`                 | Property block §6   |

---

## 12. Conventions for surface extraction

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

## 13. Open questions

1. The leading constant double = 8.0 in every property record (§6).
   It does not match any field exposed in the text export. Could be
   an internal default precision / layer count / version stamp.
2. The freedom-case DOF mask bits (§7) should be cross-referenced
   against the Strand7 *Restraint* dialog (or the `St7API.h` header
   constants) to confirm which bit corresponds to which translational /
   rotational DOF.
3. The two header doubles at 0x9B / 0xB7 (§3) are identical within a
   file but vary between files. They might be a model-wide length /
   unit scale, or simply a saved camera distance from the last GUI
   session.
4. The exact width / meaning of the 13 control bytes between the
   material name and the doubles in §6.
5. How load-case grouping is encoded for the `PlGlobalLoad` records.
   The text export groups them under headers (e.g. `/ BoneLoad`), but
   we have not located the boundary marker in the binary.
6. Group hierarchy — the GROUP DEFINITIONS lists groups but does not
   show the parent/child relationship that exists in the GUI. Likely
   encoded elsewhere in the tail.

---

## 14. Reference and acknowledgements

This document was produced by reverse-engineering a corpus of `.st7`
files alongside their corresponding Strand7-produced `.txt` exports.
None of the information here came from the official Strand7 / G+D API
documentation. If you have access to the `St7API.h` header and can
confirm or refute any of the field meanings above, please open an
issue or PR.

The `.st7` model files used to derive these notes are available from
their original publishing authors and are not included in this
repository.
