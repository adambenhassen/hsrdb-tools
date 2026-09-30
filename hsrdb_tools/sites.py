"""Symmetry expansion of the asymmetric unit into the unit cell, with anisotropic displacement tensors
rotated along with each image (gemmi's get_all_unit_cell_sites copies them unrotated)."""

from dataclasses import dataclass

import gemmi
import numpy as np

DEFAULT_B = 0.5  # Å²; HighScore's displacement parameter for atoms without one, as in the official COD databases
SAME_POSITION = 0.01  # Å; closer images of a site are one atom (coordinates of a special position, rounded)
SPECIAL_POSITION_TOLERANCE = 0.4  # Å; images this close are also one atom if the CIF says so or they cannot all be occupied


@dataclass
class Shell:
    """Atoms spread evenly over a sphere, refined as one scatterer (FullProf's SPHS/SASH: a rotating C60 molecule,
    a disordered water layer lining a cage): at each centre, form factor sum(n f_element(q)) * sin(qR)/(qR) with
    q = 2 pi / d and isotropic displacement u_iso. per_cell: its atoms per unit cell (formula, density); each
    of the len(xyz) centres carries per_cell / len(xyz) of them."""
    xyz: np.ndarray  # (n, 3) fractional coordinates of the centres in [0, 1)
    radius: float  # Å
    u_iso: float  # Å²
    per_cell: dict  # element -> atoms per unit cell


@dataclass
class Images:
    site: np.ndarray  # (n,) index into st.sites
    xyz: np.ndarray  # (n, 3) fractional coordinates in [0, 1)
    ucart: np.ndarray  # (n, 3, 3) Cartesian displacement tensor in Å² (isotropic: U·I)
    reasons: list | None = None  # per site of the structure: why its images were merged or not (merge_tolerances)
    shells: tuple = ()  # Shell scatterers besides the atoms

    def multiplicities(self, nsites):
        return np.bincount(self.site, minlength=nsites)


def ucif_to_cart(u, cell):
    """Cartesian U from a CIF (Uij against reciprocal-axis lengths) anisotropic tensor."""
    rc = cell.reciprocal()
    n = np.diag([rc.a, rc.b, rc.c])
    a = np.array(cell.orth.mat.tolist())
    ucif = np.array([[u.u11, u.u12, u.u13], [u.u12, u.u22, u.u23], [u.u13, u.u23, u.u33]])
    return a @ n @ ucif @ n @ a.T


def _occupancy_at(st, i, metric):
    """Summed occupancy of the sites at the position of site i (mixed sites are listed once per element)."""
    xyz = np.array([[s.fract.x, s.fract.y, s.fract.z] for s in st.sites])
    d = (xyz - xyz[i] + 0.5) % 1.0 - 0.5
    return sum(st.sites[k].occ for k in np.flatnonzero(np.einsum("ai,ij,aj->a", d, metric, d) < SAME_POSITION**2))


def _operations(st, ops):
    ops = list(st.spacegroup.operations()) if ops is None else ops
    return (np.array([op.rot for op in ops], dtype=float) / gemmi.Op.DEN,
            np.array([op.tran for op in ops], dtype=float) / gemmi.Op.DEN)


def _groups(x, diff, dist2, metric, tol):
    """Images that are one atom, as the index of the first image of each image's group.

    Images connected by steps shorter than tol form a group (unlike the images within tol of one image, the same
    whichever image comes first, so merged atoms keep the site's symmetry); groups whose mean positions are then
    still closer than tol merge too, until none are (e.g. images 0.39 A apart in pairs whose centres are 0.35 A
    apart: one atom on the special position, not two)."""
    n = len(dist2)
    owner = np.full(n, -1)
    for a in range(n):
        if owner[a] >= 0:
            continue
        owner[a], todo = a, [a]
        while todo:
            b = todo.pop()
            new = np.flatnonzero((owner < 0) & (dist2[b] < tol**2))
            owner[new] = a
            todo.extend(new)
    while True:
        heads = np.flatnonzero(owner == np.arange(n))
        centre = np.array([x[a] + diff[owner == a, a].mean(axis=0) for a in heads])
        d = (centre[:, None, :] - centre[None, :, :] + 0.5) % 1.0 - 0.5
        close = np.einsum("abi,ij,abj->ab", d, metric, d) < tol**2
        np.fill_diagonal(close, False)
        if not close.any():
            return owner
        i, j = np.argwhere(close)[0]  # i < j: heads are sorted
        owner[owner == heads[j]] = heads[i]


def _image_distances(st, i, rot, tran, metric):
    """Fractional images of site i under the operations, their pairwise nearest-image differences and squared
    distances (A^2)."""
    s = st.sites[i]
    x = (rot @ np.array([s.fract.x, s.fract.y, s.fract.z]) + tran) % 1.0  # (ops, 3)
    diff = (x[:, None, :] - x[None, :, :] + 0.5) % 1.0 - 0.5  # diff[a, b] = x[a] - x[b], nearest lattice image
    return x, diff, np.einsum("abi,ij,abj->ab", diff, metric, diff)


def merge_tolerances(st, ops=None, site_orders=None, separate_stated=False):
    """Per site: the distance (Å) below which its images are one atom, and why.

    Images closer than SAME_POSITION are one atom on a special position whose coordinates were rounded. Images
    up to SPECIAL_POSITION_TOLERANCE apart are ambiguous: one atom on a special position, or the separate
    positions of a split (disordered) atom whose occupancy already accounts for the split. The CIF decides when
    it states the site-symmetry order (site_orders[i]: the set of orders the CIF allows for site i, or None);
    otherwise the images are one atom only if they could not all be occupied (occupancy at the position times
    the number of distinct positions merged above 1.01). Reasons: "exact", "stated", "occupancy", "split", and
    "stated_general" for a site within SAME_POSITION of a special position that the CIF states to have a lower
    site symmetry (e.g. a general position in a negative SHELXL part, whose images count separately): its images
    are merged as the other reasons say, or, with separate_stated, not at all.
    """
    orth = np.array(st.cell.orth.mat.tolist())
    metric = orth.T @ orth
    rot, tran = _operations(st, ops)
    out = []
    for i in range(len(st.sites)):
        x, diff, dist2 = _image_distances(st, i, rot, tran, metric)
        # site-symmetry order: operations mapping the site into the group of its first image
        order_same = np.count_nonzero(_groups(x, diff, dist2, metric, SAME_POSITION) == 0)
        order_near = np.count_nonzero(_groups(x, diff, dist2, metric, SPECIAL_POSITION_TOLERANCE) == 0)
        stated = site_orders[i] if site_orders else None
        if order_near == order_same:
            tol, why = SAME_POSITION, "exact"
        elif stated and order_near in stated and order_same not in stated:
            tol, why = SPECIAL_POSITION_TOLERANCE, "stated"
        elif stated and order_same in stated and order_near not in stated:
            tol, why = SAME_POSITION, "stated"
        elif _occupancy_at(st, i, metric) * order_near / order_same > 1.01:
            tol, why = SPECIAL_POSITION_TOLERANCE, "occupancy"
        else:
            tol, why = SAME_POSITION, "split"
        if stated and order_same > 1 and order_same not in stated and min(stated) < order_same:
            tol, why = (0.0 if separate_stated else tol), "stated_general"
        out.append((tol, why))
    return out


def unit_cell_images(st, ops=None, site_orders=None, separate_stated=False):
    """Images of every site in the unit cell under ops (gemmi.Op list), by default st's space-group operations.

    Images closer than the site's merge_tolerances distance are one atom, at their mean position with their
    mean displacement tensor, which puts it exactly on its special position.
    """
    cell = st.cell
    orth = np.array(cell.orth.mat.tolist())
    metric = orth.T @ orth
    rot, tran = _operations(st, ops)
    # rotation of each operation in Cartesian space: C = A R A^-1
    cart_rot = np.einsum("ij,ojk,kl->oil", orth, rot, np.linalg.inv(orth))
    site_idx, xyz, ucart = [], [], []
    tolerances = merge_tolerances(st, ops, site_orders, separate_stated)
    for i, (s, (tol, _)) in enumerate(zip(st.sites, tolerances)):
        x, diff, dist2 = _image_distances(st, i, rot, tran, metric)
        if s.aniso.nonzero():
            us = np.einsum("oij,jk,olk->oil", cart_rot, ucif_to_cart(s.aniso, cell), cart_rot)
        else:
            us = np.broadcast_to(s.u_iso * np.eye(3), (len(rot), 3, 3))
        owner = _groups(x, diff, dist2, metric, tol)  # first image of the atom each image belongs to
        keep = np.flatnonzero(owner == np.arange(len(rot)))
        site_idx.append(np.full(len(keep), i))
        xyz.append(np.array([x[a] + diff[owner == a, a].mean(axis=0) for a in keep]) % 1.0)
        ucart.append(np.array([us[owner == a].mean(axis=0) for a in keep]))
    return Images(np.concatenate(site_idx), np.concatenate(xyz), np.concatenate(ucart), [w for _, w in tolerances])
