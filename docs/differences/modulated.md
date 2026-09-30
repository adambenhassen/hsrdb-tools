# Modulated structures

**COD24:** a modulated structure is stored as its basic structure: the pattern has main lines only, calculated as if the atoms sat at their average positions with their average occupancies.
**hsrdb-tools:** the pattern is calculated from the CIF's full superspace description: main lines as the modulation changes them, and satellite lines. The structure record holds the basic structure, the only part of the model it can represent.

## What a modulated structure is

In a modulated crystal the atoms are displaced from their average positions, or their sites are partly occupied, by waves whose period does not fit the unit cell (wave vector **q**, e.g. 0.4513 **a**\* + ½ **b**\* in COD [4001801](https://www.crystallography.net/cod/4001801.cif)). No finite cell ever repeats, so such a structure has no ordinary space group. A CIF describes it in (3+d)-dimensional superspace: the basic structure, the modulation functions of each atom along its internal coordinate x₄ = t + **q**·x̄, and superspace symmetry operations. JANA, the program that refined almost all of them, writes these to the CIF.

Its diffraction pattern holds main reflections at **h** and satellites at **h** + m**q**. The modulation changes the main reflections too: an atom that moves ±0.5 Å scatters less coherently into them, as a large displacement parameter would.

## The calculation

`hsrdb_tools/superspace.py` calculates the structure factor of each main and satellite reflection by averaging every atom's scattering over its internal coordinate, with its displacement, occupancy and displacement parameters at each point, and applies the superspace symmetry operations. Supported, as JANA defines them: Fourier waves for displacement, occupancy and displacement parameters; crenel functions (an atom present over part of x₄) with sawtooth and Legendre-polynomial displacement and displacement-parameter functions inside the interval; commensurate structures (a section t₀ of a rational modulation, whose pattern is that of its supercell).

- Satellite orders are added until two successive orders give no reflection of 0.5 on the 0–1000 scale; satellite reflections weaker than 0.5, which HighScore would show as 0, are left out. Without that limit, the weak high-order satellites of an incommensurate structure fill the 203-line list and push out the main lines above 30–45° 2θ.
- A commensurate structure is a supercell: every reflection counts, as for any other structure.
- The line list stores three Miller indices; a line whose strongest reflection is a satellite carries hkl 0 0 0.
- CIFs that list only superspace operations (and no 3D symmetry) get the space group of their basic structure from the 3D parts of those operations.

## Evidence

JANA lists the structure factors it calculated for the observed main and satellite reflections in most of these CIFs. `dev/verify_modulated.py` compares them with ours (with the CIF's own f′ and f″, neutron scattering lengths or electron form factors, and twin domains). The 30 CIFs that list them, X-ray, neutron and electron data, incommensurate and commensurate, (3+1)- and (3+2)-dimensional, agree within R = 0.86% (median 0.15%), the rounding of the listed values. A CIF whose listed values are not reproduced within 1% is stored as its basic structure (below). For the CIFs that list none, the unit-cell content is checked against the CIF's density: 171 of the 198 agree within 3%; the others disagree under any reading of their occupancies (a density or Z given for another cell, or occupancies scaled by site multiplicity), and COD24 reads them the same way.

JANA writes the occupancy of an atom confined to a crenel interval as its average over the modulation; some CIFs give the value within the interval instead. The reading that reproduces the CIF's formula × Z, else its density, is used (14 CIFs within the interval); an atom confined only by a sawtooth function is given with the value within its interval (COD 4002590, where Na and Mn alternate on one site).

Against the basic structure, COD24's model, the modulation changes a main line by 10 or more on the 0–1000 scale in 147 of the 198 entries (50 or more in 73), and the strongest satellite line reaches 10 or more in 138 (50 or more in 70). Examples:

| COD | Compound | Largest main-line change | Strongest satellite | Agreement with the CIF's structure factors |
|---|---|---|---|---|
| [1549598](https://www.crystallography.net/cod/1549598.cif) | mullite | 1.8 | 18 | none listed |
| [4001801](https://www.crystallography.net/cod/4001801.cif) | Li₀.₁₅Nd₀.₆₁₇TiO₃, (3+2)-dimensional | 94 | 6.7 | none listed |
| [1573359](https://www.crystallography.net/cod/1573359.cif) | Zn₄Si₂O₁₀ | 53 | 41 | R 0.44% |
| [1545111](https://www.crystallography.net/cod/1545111.cif) | barium Ba-IVb at 16.5 GPa | 108 | 27 | none listed |
| [1559960](https://www.crystallography.net/cod/1559960.cif) | sodium saccharinate 1.875-hydrate, commensurate | 109 | 49 | R 0.17% |
| [2310781](https://www.crystallography.net/cod/2310781.cif) | Pb₀.₇₆Bi₀.₂₀Fe₁.₀₅O₂.₆₃, crenel ordering | 101 | 299 | R 0.22% (neutron) |
| [4002590](https://www.crystallography.net/cod/4002590.cif) | NaMnO₂, commensurate, sawtooth | 393 | 936 | R 0.65% (neutron) |

## Stored as the basic structure, as in COD24

Where the superspace model cannot be evaluated, the entry is the basic structure, as COD24 stores it, and its comment and the build report say why ("modulation not included (…): pattern of the basic structure, as in COD24"):

- **Composite crystals** (30 CIFs, e.g. misfit-layer compounds): two interpenetrating subsystems with different basic lattices. As in COD24, the second subsystem's atoms are read in the first subsystem's cell.
- **Functions specific to JANA** (35 CIFs): orthogonalised harmonic functions in crenel intervals, x-harmonic functions and modulated anharmonic displacement parameters, whose exact definitions are internal to JANA.
- **CIFs whose own structure factors are not reproduced** (3 CIFs): where a CIF lists the structure factors its refinement calculated and the calculation here differs from them by more than 1% (R), a convention of the refinement is not followed; the entry is the basic structure.
- **Incomplete or ambiguous modulation data** (18 CIFs): wave vectors that do not identify their waves, a wave written without its wave number, a (3+3)-dimensional modulation.

CIFs that state a modulation but give only the average structure (25) are stored as that structure, noted "satellites missing".

Tests: `tests/test_modulated.py` (Bessel-function satellites of a displacement wave, sinc satellites of a crenel, occupancy waves, sawtooth = Legendre order 1, symmetry images and special positions against explicit atoms, commensurate structure against its explicit supercell, centred basic cell stored primitive).
