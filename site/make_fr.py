"""fr.html from index.html: every visible English string replaced by its French translation; fails if a string is
missing, so the two pages cannot drift apart silently. Also adds the language links to both pages."""
from pathlib import Path

HERE = Path(__file__).parent
NB = " "  # narrow no-break space: thousands separator and before % : ; ! ? in French
EN_LINK = ('<div class="langswitch" role="group" aria-label="Language / Langue"><a aria-current="true" lang="en">EN</a>'
           '<a href="fr.html" hreflang="fr" lang="fr" title="Français">FR</a></div>')
FR_LINK = ('<div class="langswitch" role="group" aria-label="Language / Langue"><a href="./" hreflang="en" lang="en" '
           'title="English">EN</a><a aria-current="true" lang="fr">FR</a></div>')
# the reader's own choice wins; without one, a browser whose first language is French gets the French page
LANG_JS = """<script>
(function () {
  var k = "cod-lang", chosen = null;
  try { chosen = localStorage.getItem(k); } catch (e) {}
  if (!chosen && document.documentElement.lang === "en") {
    var l = (navigator.languages && navigator.languages[0]) || navigator.language || "";
    if (/^fr\\b/i.test(l)) location.replace("fr.html" + location.hash);
  }
  document.addEventListener("click", function (e) {
    var a = e.target.closest && e.target.closest(".langswitch a[hreflang]");
    if (a) try { localStorage.setItem(k, a.getAttribute("hreflang")); } catch (err) {}
  });
})();
</script>
"""
ALT = ('<link rel="alternate" hreflang="en" href="https://cod.adambh.net/">\n'
       '<link rel="alternate" hreflang="fr" href="https://cod.adambh.net/fr.html">\n')
LANG_CSS = (".lang { margin-left: 18px; font-size: 0.9rem; color: var(--text-2); text-decoration: none; border: 1px solid "
                "var(--border); border-radius: 6px; padding: 2px 8px; white-space: nowrap; }\n.lang:hover { color: var(--text); }\n")
LANG_CSS = ("header.top .wrap > div:nth-child(2) { display: flex; align-items: center; }\n"
            ".langswitch { display: inline-flex; margin-left: 18px; border: 1px solid var(--border); border-radius: 8px; "
            "overflow: hidden; vertical-align: middle; }\n"
            ".langswitch a { padding: 3px 10px; font-size: 0.85rem; font-weight: 600; color: var(--text-2); text-decoration: none; }\n"
            ".langswitch a[href]:hover { color: var(--text); background: var(--surface-2); }\n"
            ".langswitch a[aria-current] { background: var(--accent); color: #fff; }\n")


def n(x):  # French thousands
    return x.replace(",", NB)


PAIRS = [
    ('<html lang="en">', '<html lang="fr">'),
    ('<meta name="description" content="Reference databases for Malvern Panalytical HighScore 3.x and 4.x, calculated from the Crystallography Open Database: 533,950 entries, independently verified with cctbx.">',
     f'<meta name="description" content="Bases de données de référence pour Malvern Panalytical HighScore 3.x et 4.x, calculées à partir de la Crystallography Open Database{NB}: {n("533,950")} fiches, vérifiées indépendamment avec cctbx.">'),
    ('<div class="brand">COD for HighScore <span>· reference databases</span></div>',
     '<div class="brand">Bases de référence HighScore <span>· issues de la COD</span></div>'),
    ('<a href="#download">Download</a><a href="#model">Calculation</a><a href="#verification">Verification</a><a href="#differences">Differences</a><a href="#install">Install</a><a href="#changelog">Changelog</a>',
     '<a href="#download">Téléchargement</a><a href="#model">Calcul</a><a href="#verification">Vérification</a><a href="#differences">Différences</a><a href="#install">Installation</a><a href="#changelog">Versions</a>'),
    ('<h1>The Crystallography Open Database, as HighScore reference patterns</h1>',
     '<h1>La Crystallography Open Database en diagrammes de référence pour HighScore</h1>'),
    ('<p class="lede">Calculated powder patterns and structure records for every usable entry of the <a href="https://www.crystallography.net/">Crystallography Open Database</a>, for Malvern Panalytical HighScore 3.x and 4.x and later. Built from the October 2026 COD, calculated independently of the official databases, and checked against <a href="https://cctbx.github.io/">cctbx</a>.</p>',
     '<p class="lede">Diagrammes de poudre calculés et données structurales pour chaque fiche exploitable de la <a href="https://www.crystallography.net/">Crystallography Open Database</a>, pour Malvern Panalytical HighScore 3.x et 4.x et versions ultérieures. Construits à partir de la COD de octobre 2026, calculés indépendamment des bases officielles et vérifiés avec <a href="https://cctbx.github.io/">cctbx</a>.</p>'),
    ('<ul class="badges" aria-label="Summary">', '<ul class="badges" aria-label="Résumé">'),
    ('<li>COD, October 2026</li><li>Cu Kα₁</li><li>HighScore 3.x · 4.x+</li><li>Independent recalculation</li><li>Free to use</li>',
     '<li>COD, octobre 2026</li><li>Cu Kα₁</li><li>HighScore 3.x · 4.x+</li><li>Recalcul indépendant</li><li>Libre d’utilisation</li>'),
    ('<div class="stat"><b>533,950</b><span>entries in each database</span></div>',
     f'<div class="stat"><b>{n("533,950")}</b><span>fiches dans chaque base</span></div>'),
    ('<div class="stat"><b>22,367</b><span>entries newer than the official COD 2024 database</span></div>',
     f'<div class="stat"><b>{n("22,367")}</b><span>fiches plus récentes que la base officielle COD 2024</span></div>'),
    ('<div class="stat"><b>202</b><span>modulated structures with satellite reflections</span></div>',
     '<div class="stat"><b>202</b><span>structures modulées avec réflexions satellites</span></div>'),
    ('<div class="stat"><b>0</b><span>calculation differences from cctbx</span></div>',
     '<div class="stat"><b>0</b><span>écart de calcul avec cctbx</span></div>'),
    ('<h2>Download</h2>', '<h2>Téléchargement</h2>'),
    ('<p class="muted">Each archive holds one <code>.hsrdb</code> file (about 7 GB unpacked). Checksums are read from <a href="https://cod-dl.adambh.net/SHA256SUMS.txt">SHA256SUMS.txt</a>.</p>',
     '<p class="muted">Chaque archive contient un fichier <code>.hsrdb</code> (environ 7 Go une fois décompressé). Les sommes de contrôle proviennent de <a href="https://cod-dl.adambh.net/SHA256SUMS.txt">SHA256SUMS.txt</a>.</p>'),
    ('<p class="muted">The tool that builds these databases, with the list of corrections made from publications, is open source: <a href="https://github.com/adambenhassen/hsrdb-tools">source code on GitHub</a>.</p>',
     '<p class="muted">L’outil qui construit ces bases, avec la liste des corrections faites à partir des publications, est libre : <a href="https://github.com/adambenhassen/hsrdb-tools">code source sur GitHub</a>.</p>'),
    ('<h3>HighScore 4.x and later</h3>', '<h3>HighScore 4.x et ultérieures</h3>'),
    ('<p class="meta">Structure records in the format of HighScore 4 (phase version 6).</p>',
     '<p class="meta">Données structurales au format de HighScore 4 (phase, version 6).</p>'),
    ('download>Download for HighScore 4.x</a>', 'download>Télécharger pour HighScore 4.x</a>'),
    ('<p class="meta">The same entries with structure records in the HighScore 3 format (phase version 3).</p>',
     '<p class="meta">Les mêmes fiches, avec les données structurales au format de HighScore 3 (phase, version 3).</p>'),
    ('download>Download for HighScore 3.x</a>', 'download>Télécharger pour HighScore 3.x</a>'),
    ('<summary>Show as a table</summary>', '<summary>Afficher sous forme de tableau</summary>'),
    ('<summary>How to check and unpack a download</summary>', '<summary>Vérifier et décompresser un téléchargement</summary>'),
    ('<p class="small muted" style="margin-top:12px">On Windows, check a download in PowerShell with <code>Get-FileHash</code> on the file (SHA-256 is the default) and compare the hash with <a href="https://cod-dl.adambh.net/SHA256SUMS.txt">SHA256SUMS.txt</a>, or right-click the file › 7-Zip › CRC SHA › SHA-256; then unpack it with <a href="https://www.7-zip.org/">7-Zip</a>. On Linux or macOS: <code>sha256sum -c --ignore-missing SHA256SUMS.txt</code>.</p>',
     f'<p class="small muted" style="margin-top:12px">Sous Windows, vérifiez un téléchargement dans PowerShell avec <code>Get-FileHash</code> sur le fichier (SHA-256 par défaut) et comparez l’empreinte avec <a href="https://cod-dl.adambh.net/SHA256SUMS.txt">SHA256SUMS.txt</a>, ou clic droit sur le fichier › 7-Zip › CRC SHA › SHA-256{NB}; puis décompressez-le avec <a href="https://www.7-zip.org/">7-Zip</a>. Sous Linux ou macOS{NB}: <code>sha256sum -c --ignore-missing SHA256SUMS.txt</code>.</p>'),
    ('<h2>How the patterns are calculated</h2>', '<h2>Calcul des diagrammes</h2>'),
    ('<h3>Diffraction model</h3>', '<h3>Modèle de diffraction</h3>'),
    ('<li>Kinematic structure factors with the exact International Tables (1992) form-factor coefficients and Cromer–Liberman f′ and f″ at Cu Kα₁ (1.540598 Å).</li>',
     '<li>Facteurs de structure en théorie cinématique, avec les coefficients exacts des facteurs de diffusion des International Tables (1992) et les termes f′ et f″ de Cromer–Liberman à Cu Kα₁ (1,540598 Å).</li>'),
    ('<li>|F|² averaged over Friedel pairs, so anomalous scattering is treated correctly in non-centrosymmetric structures.</li>',
     '<li>|F|² moyenné sur les paires de Friedel, pour traiter correctement la diffusion anomale dans les structures non centrosymétriques.</li>'),
    ('<li>The published displacement parameters: anisotropic U, B or β tensors, each rotated with its symmetry image; B = 0.5 Å² only where a CIF gives none, as HighScore and the official databases do.</li>',
     '<li>Les paramètres de déplacement publiés : tenseurs anisotropes U, B ou β, chacun tourné avec son image par symétrie ; B = 0,5 Å² seulement quand le CIF n’en donne pas, comme HighScore et les bases officielles.'.replace(" :", NB + ":").replace(" ;", NB + ";") + '</li>'),
    ('<li>Powder Lorentz–polarisation (1 + cos²2θ)/(sin²θ cos θ); lines to 90° 2θ, merged within 0.005°, at most 203 lines, scaled to 1000.</li>',
     f'<li>Facteur de Lorentz–polarisation de poudre (1 + cos²2θ)/(sin²θ cos θ){NB}; raies jusqu’à 90° en 2θ, fusionnées à 0,005° près, 203 raies au plus, normalisées à 1000.</li>'),
    ('<li>Reflections that fall at the same angle without being symmetry-equivalent are added, as a measured line contains all of them.</li>',
     '<li>Les réflexions de même angle qui ne sont pas équivalentes par symétrie sont additionnées, car une raie mesurée les contient toutes.</li>'),
    ('<h3>Structures as published</h3>', '<h3>Structures telles que publiées</h3>'),
    ('<li>Symmetry from the CIF\'s listed operations, else its Hall symbol, else its Hermann–Mauguin symbol; non-standard settings, shifted origins and supercells standardised with spglib only when it finds exactly the listed group.</li>',
     f'<li>Symétrie tirée des opérations listées dans le CIF, sinon de son symbole de Hall, sinon de son symbole de Hermann–Mauguin{NB}; réglages non standard, origines décalées et surstructures standardisés avec spglib seulement s’il retrouve exactement le groupe listé.</li>'),
    ('<li>Stored in the reference setting HighScore expects, in the same cell as the official databases (reduced monoclinic and triclinic cells, conventional orthorhombic axes).</li>',
     '<li>Enregistrées dans le réglage de référence attendu par HighScore, avec la même maille que les bases officielles (mailles monocliniques et tricliniques réduites, axes orthorhombiques conventionnels).</li>'),
    ('<li>Split atoms next to special positions kept or merged by the CIF\'s own site symmetry, occupancies and density; dummy atoms (ring centroids) do not scatter.</li>',
     f'<li>Atomes dédoublés près des positions spéciales conservés ou fusionnés selon la symétrie de site, les taux d’occupation et la densité donnés par le CIF{NB}; les atomes fictifs (centroïdes de cycles) ne diffusent pas.</li>'),
    ('<li>Modulated structures calculated from their superspace model, with satellites; the result is checked against the structure factors the refinement program lists in the CIF.</li>',
     f'<li>Structures modulées calculées à partir de leur modèle en superespace, avec les satellites{NB}; le résultat est comparé aux facteurs de structure que le programme d’affinement a inscrits dans le CIF.</li>'),
    ('<li>500 CIFs completed or corrected from their original publications (or from the CIF\'s own statements), each cited and marked in the entry; the CIF files themselves are never changed.</li>',
     f'<li>500 CIF complétés ou corrigés d’après leur publication d’origine (ou les indications du CIF lui-même), chacun cité et signalé dans sa fiche{NB}; les fichiers CIF eux-mêmes ne sont jamais modifiés.</li>'),
    ('<p class="small muted">Every entry carries the cell, space group, unit-cell content, calculated density, I/Ic, the structure record (atoms, occupancies, displacement parameters, Wyckoff positions) and the literature reference, and belongs to the inorganic, organic, mineral and metallic subfiles as in the official databases.</p>',
     '<p class="small muted">Chaque fiche contient la maille, le groupe d’espace, le contenu de la maille, la densité calculée, I/Ic, les données structurales (atomes, taux d’occupation, paramètres de déplacement, positions de Wyckoff) et la référence bibliographique, et appartient aux sous-fichiers inorganique, organique, minéral et métallique comme dans les bases officielles.</p>'),
    ('<h2>Verification</h2>', '<h2>Vérification</h2>'),
    ('<h3>Against cctbx</h3>', '<h3>Avec cctbx</h3>'),
    ('<p>Every pattern of COD was recalculated with <a href="https://cctbx.github.io/">cctbx</a>, an independent crystallographic library, in two ways: from the structure as read here (structure factors, form and displacement factors, multiplicities, Lorentz–polarisation), and from cctbx\'s own reading of the raw CIF (parsing, settings, site handling). The full run gave no calculation differences; each reading difference has a documented cause, such as cctbx treating dummy atoms as scatterers. The 507 CIFs corrected or newly built for this release were checked again, and so was the finished database.</p>',
     f'<p>Tous les diagrammes de la COD ont été recalculés avec <a href="https://cctbx.github.io/">cctbx</a>, une bibliothèque cristallographique indépendante, de deux manières{NB}: à partir de la structure telle que lue ici (facteurs de structure, facteurs de diffusion et de déplacement, multiplicités, Lorentz–polarisation), et à partir de la lecture du CIF brut par cctbx (analyse du fichier, réglages, traitement des sites). Le calcul complet n’a donné aucun écart de calcul{NB}; chaque écart de lecture a une cause documentée, par exemple cctbx qui traite les atomes fictifs comme des diffuseurs. Les 507 CIF corrigés ou nouvellement construits pour cette version ont été vérifiés de nouveau, de même que la base terminée.</p>'),
    ('<h3>Against the official COD databases</h3>', '<h3>Avec les bases COD officielles</h3>'),
    ('<thead><tr><th>511,583 entries shared with COD 2024</th><th class="num">Agreement</th></tr></thead>',
     f'<thead><tr><th>{n("511,583")} fiches communes avec COD 2024</th><th class="num">Accord</th></tr></thead>'),
    ('<tr><td>Space group</td><td class="num">99.47%</td></tr>', f'<tr><td>Groupe d’espace</td><td class="num">99,47{NB}%</td></tr>'),
    ('<tr><td>Cell volume within 0.5%</td><td class="num">99.99%</td></tr>', f'<tr><td>Volume de maille à 0,5{NB}% près</td><td class="num">99,99{NB}%</td></tr>'),
    ('<tr><td>Element set</td><td class="num">99.89%</td></tr>', f'<tr><td>Éléments chimiques</td><td class="num">99,89{NB}%</td></tr>'),
    ('<tr><td>I/Ic within 10%</td><td class="num">97.17%</td></tr>', f'<tr><td>I/Ic à 10{NB}% près</td><td class="num">97,17{NB}%</td></tr>'),
    ('<tr><td>All strong official lines present</td><td class="num">99.98%</td></tr>', f'<tr><td>Toutes les raies fortes officielles présentes</td><td class="num">99,98{NB}%</td></tr>'),
    ('<tr><td>Largest strong-line difference, median / 90th percentile (0–1000)</td><td class="num">12 / 111</td></tr>',
     '<tr><td>Plus grand écart sur une raie forte, médiane / 90<sup>e</sup> centile (0–1000)</td><td class="num">12 / 111</td></tr>'),
    ('<p class="small muted" style="margin-top:8px">Every entry of both databases also passes a structural audit: valid peak lists, consistent metadata and element filters, and structure records that decode and re-encode identically.</p>',
     f'<p class="small muted" style="margin-top:8px">Chaque fiche des deux bases passe aussi un audit structurel{NB}: listes de raies valides, métadonnées et filtres d’éléments cohérents, données structurales qui se décodent et se réencodent à l’identique.</p>'),
    ('<h2>Where it differs from the official databases</h2>', '<h2>Différences avec les bases officielles</h2>'),
    ('<p>Most differences come from calculation choices; each was checked on official entries.</p>',
     f'<p>La plupart des différences viennent de choix de calcul{NB}; chacune a été vérifiée sur des fiches officielles.</p>'),
    ('<thead><tr><th>Topic</th><th>Official COD 2024</th><th>Here</th></tr></thead>',
     '<thead><tr><th>Sujet</th><th>COD 2024 officielle</th><th>Ici</th></tr></thead>'),
    ('<tr><td>Coincident reflections (Laue classes −3, 6/m, 4/m, m−3)</td><td>appear to count one of them</td><td>all added (below)</td></tr>',
     '<tr><td>Réflexions coïncidentes (classes de Laue −3, 6/m, 4/m, m−3)</td><td>une seule semble comptée</td><td>toutes additionnées (ci-dessous)</td></tr>'),
    ('<tr><td>Modulated structures</td><td>basic structure, main lines only</td><td>superspace calculation with satellites</td></tr>',
     '<tr><td>Structures modulées</td><td>structure de base, raies principales seulement</td><td>calcul en superespace avec satellites</td></tr>'),
    ('<tr><td>Atom types that are not elements (<code>WO</code>, <code>ST1</code>, <code>OS1</code>)</td><td>first letters read as an element (tungsten, sulfur, osmium)</td><td>the element where certain, or the article\'s model</td></tr>',
     '<tr><td>Types d’atomes qui ne sont pas des éléments (<code>WO</code>, <code>ST1</code>, <code>OS1</code>)</td><td>premières lettres lues comme un élément (tungstène, soufre, osmium)</td><td>l’élément quand il est certain, sinon le modèle de l’article</td></tr>'),
    ('<tr><td>Dummy atoms (centroids)</td><td>counted as atoms</td><td>not scattering, not in the formula</td></tr>',
     '<tr><td>Atomes fictifs (centroïdes)</td><td>comptés comme des atomes</td><td>ne diffusent pas, absents de la formule</td></tr>'),
    ('<tr><td>Physically invalid anisotropic tensors</td><td>used as given</td><td>replaced by the isotropic value</td></tr>',
     '<tr><td>Tenseurs anisotropes physiquement impossibles</td><td>utilisés tels quels</td><td>remplacés par la valeur isotrope</td></tr>'),
    ('<tr><td>Cyanide pseudo-atom <code>CN</code></td><td>carbon only</td><td>C and N on the site</td></tr>',
     '<tr><td>Pseudo-atome cyanure <code>CN</code></td><td>carbone seul</td><td>C et N sur le site</td></tr>'),
    ('<tr><td>COD content</td><td>May 2024</td><td>October 2026, including revisions</td></tr>',
     '<tr><td>Contenu de la COD</td><td>mai 2024</td><td>octobre 2026, révisions comprises</td></tr>'),
    ('<i style="background:var(--series-1)"></i>This database</span><span><i style="background:var(--series-2)"></i>Official COD 2024</span>',
     '<i style="background:var(--series-1)"></i>Cette base</span><span><i style="background:var(--series-2)"></i>COD 2024 officielle</span>'),
    ('<title id="chartTitle">Calculated stick patterns of finnemanite, COD 9009599</title><desc id="chartDesc">This database above the axis, the official COD 2024 database mirrored below, 2θ from 10 to 45 degrees, Cu Kα1. At 26.36° the (1 2 0) and (2 1 0) reflections coincide; this database gives 353, the official one 216.</desc>',
     '<title id="chartTitle">Diagrammes en bâtons calculés de la finnemanite, COD 9009599</title><desc id="chartDesc">Cette base au-dessus de l’axe, la base officielle COD 2024 en miroir au-dessous, 2θ de 10 à 45 degrés, Cu Kα1. À 26,36° les réflexions (1 2 0) et (2 1 0) coïncident ; cette base donne 353, la base officielle 216.</desc>'),
    ('<figcaption>Finnemanite, Pb₅(AsO₃)₃Cl, P6₃/m (COD <a href="https://www.crystallography.net/cod/9009599.html">9009599</a>), Cu Kα₁. In Laue class 6/m, (1 2 0) and (2 1 0) are not symmetry-equivalent but have the same d. The official line at 3.379 Å equals (1 2 0) alone; with both it is 353, and the whole pattern rescales. Hover or tap a line for d, hkl and intensity.</figcaption>',
     f'<figcaption>Finnemanite, Pb₅(AsO₃)₃Cl, P6₃/m (COD <a href="https://www.crystallography.net/cod/9009599.html">9009599</a>), Cu Kα₁. Dans la classe de Laue 6/m, (1 2 0) et (2 1 0) ne sont pas équivalentes par symétrie mais ont le même d. La raie officielle à 3,379 Å correspond à (1 2 0) seule{NB}; avec les deux, elle vaut 353 et tout le diagramme est renormalisé. Survolez ou touchez une raie pour voir d, hkl et l’intensité.</figcaption>'),
    ('<summary>Show the lines as a table</summary>', '<summary>Afficher les raies sous forme de tableau</summary>'),
    ('<th class="num">This database</th><th class="num">Official</th>', '<th class="num">Cette base</th><th class="num">Officielle</th>'),
    ('<h2>Using it in HighScore</h2>', '<h2>Utilisation dans HighScore</h2>'),
    ('<li>Download the archive for your HighScore version and check its SHA-256.</li>',
     '<li>Téléchargez l’archive correspondant à votre version de HighScore et vérifiez son SHA-256.</li>'),
    ('<li>Unpack it with 7-Zip; keep the <code>.hsrdb</code> on a local disk.</li>',
     f'<li>Décompressez-la avec 7-Zip{NB}; gardez le fichier <code>.hsrdb</code> sur un disque local.</li>'),
    ('<li>Add the file as a reference database in HighScore.</li>',
     '<li>Ajoutez le fichier comme base de données de référence dans HighScore.</li>'),
    ('<li>Untick or remove other COD databases first: they use the same <code>96-</code> reference codes, and HighScore refuses a second database with conflicting prefixes.</li>',
     f'<li>Désactivez ou retirez d’abord les autres bases COD{NB}: elles utilisent les mêmes codes de référence <code>96-</code>, et HighScore refuse une seconde base dont les préfixes entrent en conflit.</li>'),
    ('<p class="callout">Reference codes follow the official scheme, <span class="nowrap"><code>96-</code> + (COD id + 1)</span>: COD 9009599 is <code class="nowrap">96-900-9600</code>. Each entry links back to its COD page through its id.</p>',
     f'<p class="callout">Les codes de référence suivent le schéma officiel, <span class="nowrap"><code>96-</code> + (numéro COD + 1)</span>{NB}: la COD 9009599 devient <code class="nowrap">96-900-9600</code>. Chaque fiche renvoie à sa page COD par son numéro.</p>'),
    ('<h2>What is not included</h2>', '<h2>Ce qui n’est pas inclus</h2>'),
    ('<p>1,228 CIFs cannot give a pattern: 1,017 have no atom coordinates, 152 no symmetry, 46 no unit cell, and a few have symmetry or atom types that cannot be read without guessing. None of them is in the official COD 2024 database either. Each was checked; where the original article gives what the CIF lacks, the entry is built from it.</p>',
     f'<p>{n("1,228")} CIF ne permettent pas de calculer un diagramme{NB}: {n("1,017")} n’ont pas de coordonnées atomiques, 152 pas de symétrie, 46 pas de maille, et quelques-uns ont une symétrie ou des types d’atomes illisibles sans deviner. Aucun ne figure non plus dans la base officielle COD 2024. Chacun a été examiné{NB}; lorsque l’article d’origine donne ce qui manque au CIF, la fiche est construite à partir de lui.</p>'),
    ('<p class="small muted" style="margin-bottom:0">Skipped CIFs by reason</p>', '<p class="small muted" style="margin-bottom:0">CIF non retenus, par raison</p>'),
    ('<li><span>No atom coordinates</span><span class="track"><i style="width:100%"></i></span><b>1,017</b></li>',
     f'<li><span>Pas de coordonnées atomiques</span><span class="track"><i style="width:100%"></i></span><b>{n("1,017")}</b></li>'),
    ('<li><span>No usable symmetry</span><span class="track"><i style="width:16%"></i></span><b>162</b></li>',
     f'<li><span>Pas de symétrie exploitable</span><span class="track"><i style="width:16%"></i></span><b>162</b></li>'),
    ('<li><span>No valid unit cell</span><span class="track"><i style="width:4.5%"></i></span><b>46</b></li>',
     f'<li><span>Pas de maille valide</span><span class="track"><i style="width:4.5%"></i></span><b>46</b></li>'),
    ('<li><span>Atom type that is not an element</span><span class="track"><i style="width:0.2%"></i></span><b>2</b></li>',
     f'<li><span>Type d’atome qui n’est pas un élément</span><span class="track"><i style="width:0.2%"></i></span><b>2</b></li>'),
    ('<li><span>Operations do not fit the structure</span><span class="track"><i style="width:0.1%"></i></span><b>1</b></li>',
     f'<li><span>Opérations incompatibles avec la structure</span><span class="track"><i style="width:0.1%"></i></span><b>1</b></li>'),
    ('<p class="small muted">Five entries of the official database have no usable coordinates in their current CIFs (COD 2300247, 2300248, 2300253, 2300257, 5900030). They are included as the official database has them, and their entry comments say that they are not correct structures.</p>',
     '<p class="small muted">Cinq fiches de la base officielle n’ont pas de coordonnées utilisables dans leur CIF actuel (COD 2300247, 2300248, 2300253, 2300257, 5900030). Elles sont reprises telles que la base officielle les contient, et leur commentaire indique qu’il ne s’agit pas de structures correctes.</p>'),
    ('<p><a href="https://cod-dl.adambh.net/SHA256SUMS.txt">SHA256SUMS.txt</a> · <a href="https://github.com/adambenhassen/hsrdb-tools">Source code (GitHub)</a></p>',
     '<p><a href="https://cod-dl.adambh.net/SHA256SUMS.txt">SHA256SUMS.txt</a> · <a href="https://github.com/adambenhassen/hsrdb-tools">Code source (GitHub)</a></p>'),
    ('<p>An independent build, not affiliated with or endorsed by the Crystallography Open Database or Malvern Panalytical. HighScore is a trademark of Malvern Panalytical. The patterns are calculated; compare them with your own measurements.</p>',
     f'<p>Une construction indépendante, sans lien avec la Crystallography Open Database ni Malvern Panalytical et non approuvée par elles. HighScore est une marque de Malvern Panalytical. Les diagrammes sont calculés{NB}; comparez-les à vos propres mesures.</p>'),
    ('if (n) parts.push((n / 1e9).toFixed(2) + " GB");', 'if (n) parts.push((n / 1e9).toFixed(2).replace(".", ",") + "\\u202fGo");'),
    ('if (lm) parts.push("published " + new Date(lm).toISOString().slice(0, 10));', 'if (lm) parts.push("publié le " + new Date(lm).toLocaleDateString("fr-FR", { day: "numeric", month: "long", year: "numeric", timeZone: "UTC" }));'),
    ('xl.textContent = "2θ (degrees), Cu Kα₁";', 'xl.textContent = "2θ (degrés), Cu Kα₁";'),
    ('"<b>" + (name === "ours" ? "This database" : "Official COD 2024")', '"<b>" + (name === "ours" ? "Cette base" : "COD 2024 officielle")'),
    ('"aria-label": (name === "ours" ? "This database" : "Official COD 2024") + ", 2θ " + r[0].toFixed(2) + "°, hkl "',
     '"aria-label": (name === "ours" ? "Cette base" : "COD 2024 officielle") + ", 2θ " + r[0].toFixed(2).replace(".", ",") + "°, hkl "'),
    ('lab1.textContent = "This database";', 'lab1.textContent = "Cette base";'),
    ('lab2.textContent = "Official COD 2024";', 'lab2.textContent = "COD 2024 officielle";'),
    ('"</b><br>2θ " + r[0].toFixed(2) + "°, d " + r[1].toFixed(4) + " Å<br>hkl "',
     '"</b><br>2θ " + r[0].toFixed(2).replace(".", ",") + "°, d " + r[1].toFixed(4).replace(".", ",") + " Å<br>hkl "'),
    ('<p class="release">Release <b data-release>2026-10-30</b>, built from the COD of October 2026 · <a href="#changelog">changelog</a></p>',
     '<p class="release">Version du <b data-release>2026-10-30</b>, construite à partir de la COD d’octobre 2026 · <a href="#changelog">historique des versions</a></p>'),
    ('<h2>Changelog</h2>', '<h2>Historique des versions</h2>'),
    ('<p class="muted">Each release is named by its build date, in the archive and in the database file HighScore shows: <code>COD_2026_10_30_HS4x.hsrdb</code>.</p>',
     f'<p class="muted">Chaque version porte la date de sa construction, dans le nom de l’archive et dans celui du fichier de base affiché par HighScore{NB}: <code>COD_2026_10_30_HS4x.hsrdb</code>.</p>'),
    ('        <p class="charttitle">Skipped CIFs by decade of publication (432 more have no year)</p>\n        <ul class="cols"><li><b>1</b><i style="height:2px"></i><em>1910</em></li><li><b>13</b><i style="height:5px"></i><em>1920</em></li><li><b>201</b><i style="height:82px"></i><em>1930</em></li><li><b>43</b><i style="height:17px"></i><em>1940</em></li><li><b>38</b><i style="height:15px"></i><em>1950</em></li><li><b>13</b><i style="height:5px"></i><em>1960</em></li><li><b>12</b><i style="height:5px"></i><em>1970</em></li><li><b>12</b><i style="height:5px"></i><em>1980</em></li><li><b>208</b><i style="height:85px"></i><em>1990</em></li><li><b>246</b><i style="height:100px"></i><em>2000</em></li><li><b>8</b><i style="height:3px"></i><em>2010</em></li><li><b>1</b><i style="height:2px"></i><em>2020</em></li></ul>\n',
     f'        <p class="charttitle">CIF non retenus, par décennie de publication (432 autres n’ont pas d’année)</p>\n        <ul class="cols"><li><b>1</b><i style="height:2px"></i><em>1910</em></li><li><b>13</b><i style="height:5px"></i><em>1920</em></li><li><b>201</b><i style="height:82px"></i><em>1930</em></li><li><b>43</b><i style="height:17px"></i><em>1940</em></li><li><b>38</b><i style="height:15px"></i><em>1950</em></li><li><b>13</b><i style="height:5px"></i><em>1960</em></li><li><b>12</b><i style="height:5px"></i><em>1970</em></li><li><b>12</b><i style="height:5px"></i><em>1980</em></li><li><b>208</b><i style="height:85px"></i><em>1990</em></li><li><b>246</b><i style="height:100px"></i><em>2000</em></li><li><b>8</b><i style="height:3px"></i><em>2010</em></li><li><b>1</b><i style="height:2px"></i><em>2020</em></li></ul>\n'),
    ('        <p class="charttitle">Agreement with COD 2024 (axis from 96%)</p>\n        <ul class="dots">\n          <li><span>Space group</span><span class="track"><i style="left:86.75%"></i></span><b>99.47%</b></li>\n          <li><span>Cell volume within 0.5%</span><span class="track"><i style="left:99.75%"></i></span><b>99.99%</b></li>\n          <li><span>Element set</span><span class="track"><i style="left:97.25%"></i></span><b>99.89%</b></li>\n          <li><span>I/Ic within 10%</span><span class="track"><i style="left:29.25%"></i></span><b>97.17%</b></li>\n          <li><span>All strong official lines present</span><span class="track"><i style="left:99.5%"></i></span><b>99.98%</b></li>\n          <li><span></span><span class="ticks"><em>96%</em><em>97%</em><em>98%</em><em>99%</em><em>100%</em></span><b></b></li>\n        </ul>\n',
     f'        <p class="charttitle">Accord avec la COD 2024 (axe à partir de 96 %)</p>\n        <ul class="dots">\n          <li><span>Groupe d’espace</span><span class="track"><i style="left:86.75%"></i></span><b>99,47{NB}%</b></li>\n          <li><span>Volume de maille à 0,5 % près</span><span class="track"><i style="left:99.75%"></i></span><b>99,99{NB}%</b></li>\n          <li><span>Jeu d’éléments</span><span class="track"><i style="left:97.25%"></i></span><b>99,89{NB}%</b></li>\n          <li><span>I/Ic à 10 % près</span><span class="track"><i style="left:29.25%"></i></span><b>97,17{NB}%</b></li>\n          <li><span>Toutes les raies fortes officielles présentes</span><span class="track"><i style="left:99.5%"></i></span><b>99,98{NB}%</b></li>\n          <li><span></span><span class="ticks"><em>96{NB}%</em><em>97{NB}%</em><em>98{NB}%</em><em>99{NB}%</em><em>100{NB}%</em></span><b></b></li>\n        </ul>\n'),
    ('        <p class="charttitle">Largest strong-line difference per entry, official lines (0–1000 scale)</p>\n        <ul class="bars short">\n          <li><span>Median</span><span class="track"><i style="width:1.2%"></i></span><b>12</b></li>\n          <li><span>90th percentile</span><span class="track"><i style="width:11.1%"></i></span><b>111</b></li>\n          <li><span>99th percentile</span><span class="track"><i style="width:38.1%"></i></span><b>381</b></li>\n          <li><span></span><span class="ticks"><em>0</em><em>500</em><em>1,000</em></span><b></b></li>\n        </ul>\n',
     f'        <p class="charttitle">Plus grand écart sur une raie forte, par fiche, raies officielles (échelle 0–1000)</p>\n        <ul class="bars short">\n          <li><span>Médiane</span><span class="track"><i style="width:1.2%"></i></span><b>12</b></li>\n          <li><span>90ᵉ centile</span><span class="track"><i style="width:11.1%"></i></span><b>111</b></li>\n          <li><span>99ᵉ centile</span><span class="track"><i style="width:38.1%"></i></span><b>381</b></li>\n          <li><span></span><span class="ticks"><em>0</em><em>500</em><em>1{NB}000</em></span><b></b></li>\n        </ul>\n'),
    ('    <div class="grid2">\n      <div>\n        <p class="charttitle">Corrected CIFs per release</p>\n        <ul class="bars short"><li><span>2026-10-30</span><span class="track"><i style="width:100%"></i></span><b>500</b></li><li><span>2026-10-29</span><span class="track"><i style="width:100%"></i></span><b>500</b></li><li><span>2026-10-15</span><span class="track"><i style="width:100%"></i></span><b>500</b></li><li><span>2026-10-04</span><span class="track"><i style="width:99.8%"></i></span><b>499</b></li><li><span>2026-10-03</span><span class="track"><i style="width:100%"></i></span><b>499</b></li><li><span>2026-09-30</span><span class="track"><i style="width:11.4%"></i></span><b>57</b></li></ul>\n      </div>\n      <div>\n        <p class="charttitle">Skipped CIFs per release</p>\n        <ul class="bars short"><li><span>2026-10-30</span><span class="track"><i style="width:73.6%"></i></span><b>1,228</b></li><li><span>2026-10-29</span><span class="track"><i style="width:73.6%"></i></span><b>1,228</b></li><li><span>2026-10-15</span><span class="track"><i style="width:73.5%"></i></span><b>1,227</b></li><li><span>2026-10-04</span><span class="track"><i style="width:73.4%"></i></span><b>1,225</b></li><li><span>2026-10-03</span><span class="track"><i style="width:73.3%"></i></span><b>1,224</b></li><li><span>2026-09-30</span><span class="track"><i style="width:100%"></i></span><b>1,669</b></li></ul>\n      </div>\n    </div>\n',
     f'    <div class="grid2">\n      <div>\n        <p class="charttitle">CIF corrigés, par version</p>\n        <ul class="bars short"><li><span>2026-10-30</span><span class="track"><i style="width:100%"></i></span><b>500</b></li><li><span>2026-10-29</span><span class="track"><i style="width:100%"></i></span><b>500</b></li><li><span>2026-10-15</span><span class="track"><i style="width:100%"></i></span><b>500</b></li><li><span>2026-10-04</span><span class="track"><i style="width:99.8%"></i></span><b>499</b></li><li><span>2026-10-03</span><span class="track"><i style="width:100%"></i></span><b>499</b></li><li><span>2026-09-30</span><span class="track"><i style="width:11.4%"></i></span><b>57</b></li></ul>\n      </div>\n      <div>\n        <p class="charttitle">CIF non retenus, par version</p>\n        <ul class="bars short"><li><span>2026-10-30</span><span class="track"><i style="width:73.6%"></i></span><b>1{NB}228</b></li><li><span>2026-10-29</span><span class="track"><i style="width:73.6%"></i></span><b>1{NB}228</b></li><li><span>2026-10-15</span><span class="track"><i style="width:73.5%"></i></span><b>1{NB}227</b></li><li><span>2026-10-04</span><span class="track"><i style="width:73.4%"></i></span><b>1{NB}225</b></li><li><span>2026-10-03</span><span class="track"><i style="width:73.3%"></i></span><b>1{NB}224</b></li><li><span>2026-09-30</span><span class="track"><i style="width:100%"></i></span><b>1{NB}669</b></li></ul>\n      </div>\n    </div>\n'),
    ('<h3>2026-10-30 <span>533,950 entries</span></h3>', f'<h3>2026-10-30 <span>{n("533,950")} fiches</span></h3>'),
    ("<li>Deuterium has its own element flag, as in the official databases: a HighScore element search for D finds the 1,760 deuterated entries (before, D shared hydrogen's flag). Everything else is unchanged from 2026-10-29.</li>",
     f'<li>Le deutérium a désormais son propre indicateur d’élément, comme dans les bases officielles{NB}: une recherche de D dans HighScore trouve les {n("1,760")} fiches deutérées (auparavant, D partageait l’indicateur de l’hydrogène). Le reste est inchangé par rapport à 2026-10-29.</li>'),
    ('<h3>2026-10-29 <span>533,950 entries</span></h3>', f'<h3>2026-10-29 <span>{n("533,950")} fiches</span></h3>'),
    ('<li>11 CIFs added to the COD since the previous release give entries (22,367 newer than the official COD 2024 database); one more unpublished deposit, without coordinates, is skipped (1,228 in all).</li>',
     f'<li>11 CIF ajoutés à la COD depuis la version précédente donnent des fiches ({n("22,367")} plus récentes que la base officielle COD 2024){NB}; un autre dépôt non publié, sans coordonnées, est écarté ({n("1,228")} au total).</li>'),
    ('<li>The COD revised seven CIFs in the meantime (1010542, 4100872–4100877); the entries 4100872–4100877 now carry their correct authors, title, journal and pages.</li>',
     f'<li>La COD a révisé sept CIF entre-temps (1010542, 4100872–4100877){NB}; les fiches 4100872–4100877 portent désormais leurs auteurs, titre, revue et pages corrects.</li>'),
    ('<li>A calculation error in four commensurately modulated structures is fixed (COD 2310038, 2310043, 2310046, 2310055): line intensities change by up to 5 on the 0–1000 scale, I/Ic by up to 1.4%.</li>',
     f'<li>Une erreur de calcul dans quatre structures modulées commensurables est corrigée (COD 2310038, 2310043, 2310046, 2310055){NB}: les intensités des raies changent d’au plus 5 sur l’échelle 0–1000, I/Ic d’au plus 1,4{NB}%.</li>'),
    ('<li>A line holding several reflections of equal intensity at one d-spacing now carries the same hkl label on every platform (2,849 entries); patterns are unchanged.</li>',
     f'<li>Une raie qui regroupe plusieurs réflexions de même intensité à une même distance d porte désormais le même indice hkl sur toutes les plateformes ({n("2,849")} fiches){NB}; les diagrammes sont inchangés.</li>'),
    ('<li>Rebuilt from the COD mirror of 29 October 2026.</li>', '<li>Reconstruit à partir du miroir de la COD du 29 octobre 2026.</li>'),
    ('<h3>2026-10-15 <span>533,939 entries</span></h3>', f'<h3>2026-10-15 <span>{n("533,939")} fiches</span></h3>'),
    ('<li>7 CIFs added to the COD since the previous release give entries (22,356 newer than the official COD 2024 database); three more unpublished deposits, without coordinates, are skipped.</li>',
     f'<li>7 CIF ajoutés à la COD depuis la version précédente donnent des fiches ({n("22,356")} plus récentes que la base officielle COD 2024){NB}; trois autres dépôts non publiés, sans coordonnées, sont écartés.</li>'),
    ('<li>LaONO₃ (COD 6000784), skipped before, is now built from its article (500 corrected CIFs); 1,227 CIFs are skipped in all.</li>',
     f'<li>LaONO₃ (COD 6000784), écarté jusqu’ici, est désormais construit d’après son article (500 CIF corrigés){NB}; {n("1,227")} CIF sont écartés au total.</li>'),
    ('<li>Rebuilt from the COD mirror of 15 October 2026; the calculation is unchanged.</li>',
     '<li>Reconstruit à partir du miroir de la COD du 15 octobre 2026{}; le calcul est inchangé.</li>'.format(NB)),
    ('<h3>2026-10-04 <span>533,931 entries</span></h3>', f'<h3>2026-10-04 <span>{n("533,931")} fiches</span></h3>'),
    ('<li>31 CIFs added to the COD since the previous release give entries (22,348 newer than the official COD 2024 database); one more unpublished deposit, without coordinates, is skipped (1,225 in all).</li>',
     f'<li>31 CIF ajoutés à la COD depuis la version précédente donnent des fiches ({n("22,348")} plus récentes que la base officielle COD 2024){NB}; un autre dépôt non publié, sans coordonnées, est écarté ({n("1,225")} au total).</li>'),
    ('<li>The COD revised one CIF in the meantime (9011997); its entry follows the revision.</li>',
     '<li>La COD a révisé un CIF entre-temps (9011997){}; sa fiche suit la révision.</li>'.format(NB)),
    ('<li>Rebuilt from the COD mirror of 3 October 2026; the calculation is unchanged.</li>',
     '<li>Reconstruit à partir du miroir de la COD du 3 octobre 2026{}; le calcul est inchangé.</li>'.format(NB)),
    ('<h3>2026-10-03 <span>533,900 entries</span></h3>', f'<h3>2026-10-03 <span>{n("533,900")} fiches</span></h3>'),
    ('<li>445 more CIFs give an entry (22,317 newer than the official COD 2024 database); 1,224 CIFs remain skipped, down from 1,669.</li>',
     f'<li>445 CIF de plus donnent une fiche ({n("22,317")} plus récentes que la base officielle COD 2024){NB}; {n("1,224")} CIF restent non retenus, contre {n("1,669")}.</li>'),
    ('<li>499 CIFs are now completed or corrected from their original publications (57 before), each cited and marked in the entry; each correction was checked against the article, and the CIF files themselves are never changed.</li>',
     f'<li>499 CIF sont désormais complétés ou corrigés d’après leur publication d’origine (57 auparavant), chacun cité et signalé dans sa fiche{NB}; chaque correction a été vérifiée d’après l’article, et les fichiers CIF eux-mêmes ne sont jamais modifiés.</li>'),
    ('<li>The corrected CIFs were audited again afterwards: sanity checks on all of them, and coordinate tables compared row by row with the article wherever its PDF was available; two errors were found and fixed.</li>',
     f'<li>Les CIF corrigés ont ensuite été audités à nouveau{NB}: contrôles de cohérence sur tous, et tableaux de coordonnées comparés ligne à ligne à l’article lorsque son PDF était disponible{NB}; deux erreurs ont été trouvées et corrigées.</li>'),
    ('<li>Four more modulated structures are calculated from the superspace model their article gives (202 entries in all); structures whose CIF gives only a superspace symbol are built as average structures.</li>',
     f'<li>Quatre structures modulées de plus sont calculées à partir du modèle en superespace donné par leur article (202 fiches au total){NB}; les structures dont le CIF ne donne qu’un symbole de superespace sont construites comme structures moyennes.</li>'),
    ('<li>The source code and the list of corrections are public on GitHub.</li>',
     '<li>Le code source et la liste des corrections sont publics sur GitHub.</li>'),
    ('<h3>2026-09-30 <span>533,455 entries</span></h3>', f'<h3>2026-09-30 <span>{n("533,455")} fiches</span></h3>'),
    ('<li>Every entry of the official COD 2024 database is included, with 21,872 newer ones.</li>',
     f'<li>Toutes les fiches de la base officielle COD 2024 sont incluses, ainsi que {n("21,872")} fiches plus récentes.</li>'),
    ('<li>Modulated structures are calculated from their superspace model, with satellite reflections (198 entries); where a model cannot be evaluated, the basic structure is stored, as in COD 2024, and the entry says so.</li>',
     f'<li>Les structures modulées sont calculées à partir de leur modèle en superespace, avec les réflexions satellites (198 fiches){NB}; lorsqu’un modèle ne peut pas être évalué, la structure de base est enregistrée, comme dans COD 2024, et la fiche l’indique.</li>'),
    ('<li>57 CIFs are completed or corrected from their original publications or their own statements.</li>',
     '<li>57 CIF sont complétés ou corrigés d’après leur publication d’origine ou leurs propres indications.</li>'),
    ('<li>Cells are stored as the official databases store them: reduced monoclinic and triclinic cells (P2₁/n and P2₁/a now become P2₁/c with the short cell, not β ≈ 140°), conventional orthorhombic axes; about 94,000 entries. Patterns are unchanged; hkl labels follow the cell.</li>',
     f'<li>Les mailles sont enregistrées comme dans les bases officielles{NB}: mailles monocliniques et tricliniques réduites (P2₁/n et P2₁/a deviennent P2₁/c avec la maille courte, et non β ≈ 140°), axes orthorhombiques conventionnels{NB}; environ {n("94,000")} fiches. Les diagrammes sont inchangés{NB}; les indices hkl suivent la maille.</li>'),
    ('<li>The listed symmetry operations or the Hall symbol decide over a disagreeing Hermann–Mauguin symbol, so entries whose atoms were doubled are now right.</li>',
     '<li>Les opérations de symétrie listées ou le symbole de Hall l’emportent sur un symbole de Hermann–Mauguin discordant{}: les fiches dont les atomes étaient doublés sont désormais correctes.</li>'.format(NB)),
    ('<li>Split atoms next to special positions follow the CIF\'s own site symmetry, occupancies and density; displacement parameters where U<sub>iso</sub> is missing come from B<sub>iso</sub>; tensors with an empty diagonal are replaced.</li>',
     f'<li>Les atomes dédoublés près des positions spéciales suivent la symétrie de site, les taux d’occupation et la densité donnés par le CIF{NB}; lorsque U<sub>iso</sub> manque, les paramètres de déplacement viennent de B<sub>iso</sub>{NB}; les tenseurs à diagonale vide sont remplacés.</li>'),
    ('<li>Archive and database names carry the build date.</li>', '<li>Les noms des archives et des bases portent la date de construction.</li>'),
    ('<h3>2026-09-27 <span>533,229 entries</span></h3>', f'<h3>2026-09-27 <span>{n("533,229")} fiches</span></h3>'),
    ('<li>First public release.</li>', '<li>Première version publique.</li>'),
    ('<p><a href="skipped.html">All 1,228 skipped CIFs, with their references, DOIs and reasons →</a></p>',
     f'<p><a href="skipped.fr.html">Les {n("1,228")} CIF non retenus, avec leurs références, DOI et raisons →</a></p>'),
    ('<p><a href="corrections.html">The 500 corrected CIFs, each with its problem, correction and source →</a></p>',
     '<p><a href="corrections.fr.html">Les 500 CIF corrigés, chacun avec son problème, sa correction et sa source →</a></p>'),
]


def main():
    en = (HERE / "index.html").read_text()
    # language links and alternates on the English page (idempotent)
    if EN_LINK not in en:
        en = en.replace("      </nav>\n    </div>", "      </nav>\n      " + EN_LINK + "\n    </div>", 1)
        en = en.replace("<title>", ALT + LANG_JS + "<title>", 1)
        en = en.replace("</style>", LANG_CSS + "</style>", 1)
        (HERE / "index.html").write_text(en)
    fr = en
    for a, b in [(EN_LINK, FR_LINK), ('<title>COD for HighScore</title>', '<title>Bases de référence HighScore issues de la COD</title>')]:
        if fr.count(a) != 1:
            raise SystemExit(f"not found exactly once: {a[:80]!r}")
        fr = fr.replace(a, b)
    for a, b in PAIRS:
        if fr.count(a) != 1:
            raise SystemExit(f"not found exactly once: {a[:80]!r}")
        fr = fr.replace(a, b)
    (HERE / "fr.html").write_text(fr)


main()
