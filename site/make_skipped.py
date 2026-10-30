"""skipped.html and skipped.fr.html: every CIF the build skips, with its COD link, reference, DOI and reason.

Input: skipmeta.json (one record per skipped CIF, from the build's skipped report and the CIFs themselves; see the
handoff) and docs/corrections.md (notes for the CIFs checked against their articles). The style is taken from
index.html so the pages stay alike.
"""
import html
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).parent
REPO = HERE.parent
NB = " "

# key: (English label, English explanation, French label, French explanation)
KINDS = {
    "retracted": ("Retracted or marked erroneous",
                  "COD keeps the entry as the record of a retraction or flags it as erroneous; it has no usable structure.",
                  "Rétractée ou signalée erronée",
                  "La COD conserve la fiche comme trace d’une rétractation ou la signale comme erronée ; elle ne contient pas de structure exploitable."),
    "tbp": ("Unpublished deposit",
            "Deposited before publication without atom coordinates; the structure itself was not made available.",
            "Dépôt non publié",
            "Déposée avant publication sans coordonnées atomiques ; la structure elle-même n’a pas été rendue disponible."),
    "powder": ("Cell only (powder study)",
               "The CIF records the cell and symmetry of a powder study, but no atom positions. Where the article gives the structure, the entry is built from it; the articles checked so far are noted.",
               "Maille seule (étude sur poudre)",
               "Le CIF donne la maille et la symétrie d’une étude sur poudre, mais aucune position atomique. Lorsque l’article donne la structure, la fiche est construite à partir de lui ; les articles déjà examinés sont indiqués."),
    "coords": ("No atom coordinates",
               "The CIF has a cell and symmetry but no atom positions; they may exist only in the original paper.",
               "Pas de coordonnées atomiques",
               "Le CIF donne une maille et une symétrie mais aucune position atomique ; elles n’existent peut-être que dans l’article d’origine."),
    "nosym": ("No symmetry",
              "Neither symmetry operations nor a space-group symbol: nearly all are empty entries from papers of 1926–1962, the others modulated structures whose superspace model the CIF lacks.",
              "Pas de symétrie",
              "Ni opérations de symétrie ni symbole de groupe d’espace : presque toutes sont des fiches vides d’articles de 1926 à 1962, les autres des structures modulées dont le CIF n’a pas le modèle en superespace."),
    "nocell": ("No unit cell",
               "An empty entry that gives only a formula and a reference.",
               "Pas de maille",
               "Une fiche vide qui ne donne qu’une formule et une référence."),
    "symbol": ("Unreadable space-group symbol",
               "The symbol is not a space group (for example Pbc2, ?P?, unknown) and no operations are listed; reading it would mean guessing.",
               "Symbole de groupe d’espace illisible",
               "Le symbole n’est pas un groupe d’espace (par exemple Pbc2, ?P?, unknown) et aucune opération n’est donnée ; le lire reviendrait à deviner."),
    "superspace": ("Superspace symbol without operations",
                   "A modulated structure whose superspace symmetry operations the CIF does not list.",
                   "Symbole de superespace sans opérations",
                   "Une structure modulée dont le CIF ne donne pas les opérations de symétrie en superespace."),
    "atomtype": ("Atom type that is not an element",
                 "Atoms whose element neither the CIF nor its article decides.",
                 "Type d’atome qui n’est pas un élément",
                 "Des atomes dont ni le CIF ni l’article ne permettent de fixer l’élément."),
    "ops": ("Operations do not fit the structure",
            "The listed symmetry operations do not form the space group of the atoms given.",
            "Opérations incompatibles avec la structure",
            "Les opérations de symétrie données ne forment pas le groupe d’espace des atomes fournis."),
}
ORDER = list(KINDS)


GREEK = dict(zip("abgdezhqiklmnxoprstfcyw", "αβγδεζηθικλμνξοπρστφχψω"))


def cif_text(t):
    """CIF text markup as HTML: \\m -> μ (Greek letters), ^2^ superscript, ~3~ subscript; the rest escaped."""
    t = html.escape(t)
    t = re.sub(r"\\([A-Za-z])", lambda m: (GREEK.get(m.group(1).lower(), m.group(1)).upper() if m.group(1).isupper()
                                               else GREEK.get(m.group(1), m.group(1))), t)
    t = re.sub(r"\^([^^\s]{1,12})\^", r"<sup>\1</sup>", t)
    return re.sub(r"~([^~\s]{1,12})~", r"<sub>\1</sub>", t)


def kind(r):
    reason = r["reason"]
    if r["retracted"]:
        return "retracted"
    if r["tbp"]:
        return "tbp"
    if reason.startswith("no atom coordinates"):
        return "powder" if r["id"].startswith("6") else "coords"
    if reason.startswith("no valid unit cell"):
        return "nocell"
    if "listed operations do not match" in reason:
        return "ops"
    if reason.startswith("atom type"):
        return "atomtype"
    if reason.startswith("unknown space group"):
        hm = r["hm"]
        if not hm:
            return "nosym"
        return "superspace" if re.search(r"\(|^X", hm) else "symbol"
    raise SystemExit(f"unclassified skip reason for {r['id']}: {reason}")


def notes():
    """COD id -> the reason given in docs/corrections.md, 'Skipped CIFs checked against their articles'."""
    text = (REPO / "docs" / "corrections.md").read_text()
    part = text.split("## Skipped CIFs checked against their articles", 1)[1]
    out = {}
    for line in part.splitlines():
        if not line.startswith("| ["):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split(" | ")]
        why = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", cells[-1]).replace("**", "")
        for cod in re.findall(r"\[(\d{7})\]", cells[0]):
            out[cod] = why
    return out


def reference(r):
    parts = [r["journal"]]
    if r["volume"]:
        parts.append(r["volume"] + (", " + r["page"] if r["page"] else ""))
    ref = " ".join(p for p in parts if p)
    return ref + (f" ({r['year']})" if r["year"] else "")


def style():
    src = (HERE / "index.html").read_text()
    return src[src.index("<style>"):src.index("</style>") + len("</style>")]


T = {
    "en": dict(lang="en", title="Skipped CIFs", home="index.html", other="skipped.fr.html", other_lang="fr", other_label="FR",
               brand="COD for HighScore <span>· skipped CIFs</span>", back="← Back to the databases",
               h1="Skipped CIFs",
               lede=f"The {{n}} COD CIFs from which no correct pattern can be calculated, each with its reason. None of them is in the official COD 2024 database either. Where an article gives what a CIF lacks, the entry is built from it; the CIFs checked against their articles carry a note.",
               summary="Reasons", count="CIFs", search="Search COD id, compound, formula, journal or DOI", all="All reasons",
               th=("COD", "Compound", "Reference", "Reason"), shown="{a} of {b} shown", note="Checked against the article:",
               dup="COD marks it as a duplicate of", nodoi="no DOI", foot="Generated from the build report of release {rel}."),
    "fr": dict(lang="fr", title="CIF non retenus", home="fr.html", other="skipped.html", other_lang="en", other_label="EN",
               brand="COD pour HighScore <span>· CIF non retenus</span>", back="← Retour aux bases",
               h1="CIF non retenus",
               lede=f"Les {{n}} CIF de la COD à partir desquels aucun diagramme correct ne peut être calculé, chacun avec sa raison. Aucun ne figure non plus dans la base officielle COD 2024. Lorsqu’un article donne ce qui manque au CIF, la fiche est construite à partir de lui ; les CIF examinés d’après leur article portent une note (en anglais).",
               summary="Raisons", count="CIF", search="Rechercher un numéro COD, un composé, une formule, une revue ou un DOI", all="Toutes les raisons",
               th=("COD", "Composé", "Référence", "Raison"), shown="{a} sur {b} affichés", note="Examiné d’après l’article :",
               dup="La COD la signale comme doublon de", nodoi="pas de DOI", foot="Produit à partir du rapport de construction de la version {rel}."),
}


def page(rows, lang, release):
    t = T[lang]
    li = 0 if lang == "en" else 2
    fmt = (lambda x: f"{x:,}") if lang == "en" else (lambda x: f"{x:,}".replace(",", NB))
    counts = {k: sum(1 for r in rows if r["kind"] == k) for k in ORDER}
    esc = html.escape
    out = [f"""<!doctype html>
<html lang="{t['lang']}">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{t['title']}</title>
<link rel="alternate" hreflang="en" href="https://cod.adambh.net/skipped.html">
<link rel="alternate" hreflang="fr" href="https://cod.adambh.net/skipped.fr.html">
{style()}
<style>
.wrap {{ max-width: 1440px; }}
.lede {{ max-width: 820px; }}
#rows td:nth-child(2) {{ min-width: 180px; }}
#rows td:nth-child(3) {{ min-width: 260px; }}
.filters {{ display: flex; flex-wrap: wrap; gap: 10px; margin: 14px 0 10px; }}
.filters input, .filters select {{ font: inherit; font-size: 0.95rem; padding: 7px 10px; border: 1px solid var(--border); border-radius: 8px; background: var(--surface); color: var(--text); }}
.filters input {{ flex: 1 1 320px; min-width: 0; }}
#rows td {{ font-size: 0.9rem; }}
#rows td.id {{ white-space: nowrap; }}
#rows .why {{ display: block; color: var(--text-2); font-size: 0.85rem; margin-top: 3px; }}
.kinds td:first-child {{ width: 32%; }}
.langlink {{ font-size: 0.85rem; font-weight: 600; border: 1px solid var(--border); border-radius: 8px; padding: 3px 10px; text-decoration: none; color: var(--text-2); }}
</style>
</head>
<body>
<header class="top"><div class="wrap"><div class="brand">{t['brand']}</div><div><a class="langlink" href="{t['other']}" hreflang="{t['other_lang']}" lang="{t['other_lang']}">{t['other_label']}</a></div></div></header>
<main><div class="wrap" style="padding-top:28px;padding-bottom:40px">
<p><a href="{t['home']}">{t['back']}</a></p>
<h1>{t['h1']}</h1>
<p class="lede">{t['lede'].format(n=fmt(len(rows)))}</p>
<h2 style="margin-top:24px">{t['summary']}</h2>
<div class="tablewrap"><table class="kinds"><tbody>"""]
    for k in ORDER:
        if counts[k]:
            out.append(f"<tr><td><b>{esc(KINDS[k][li])}</b></td><td class=\"num\">{fmt(counts[k])}</td>"
                       f"<td>{esc(KINDS[k][li + 1])}</td></tr>")
    out.append(f"""</tbody></table></div>
<div class="filters"><input id="q" type="search" placeholder="{esc(t['search'])}" aria-label="{esc(t['search'])}">
<select id="k" aria-label="{esc(t['summary'])}"><option value="">{esc(t['all'])}</option>""")
    for k in ORDER:
        if counts[k]:
            out.append(f'<option value="{k}">{esc(KINDS[k][li])} ({fmt(counts[k])})</option>')
    out.append(f"""</select></div>
<p class="small muted" id="shown" aria-live="polite"></p>
<div class="tablewrap"><table><thead><tr>{''.join(f'<th>{h}</th>' for h in t['th'])}</tr></thead><tbody id="rows">""")
    for r in rows:
        name = cif_text(r["name"] or "")
        formula = esc(r["formula"] or "")
        comp = name + ("<br>" if name and formula else "") + (f'<span class="mono small">{formula}</span>' if formula else "")
        doi = (f'<br><a href="https://doi.org/{esc(r["doi"])}">doi:{esc(r["doi"])}</a>' if r["doi"]
               else f'<br><span class="small muted">{t["nodoi"]}</span>')
        why = KINDS[r["kind"]][li]
        extra = ""
        if r["note"]:
            extra = f'<span class="why">{t["note"]} {esc(r["note"])}</span>'
        elif r["dup"]:
            extra = f'<span class="why">{t["dup"]} <a href="https://www.crystallography.net/cod/{esc(r["dup"])}.html">{esc(r["dup"])}</a>.</span>'
        search = " ".join([r["id"], r["name"], r["formula"], r["journal"], r["doi"], r["year"]]).lower()
        out.append(f'<tr data-k="{r["kind"]}" data-s="{esc(search)}"><td class="id"><a href="https://www.crystallography.net/cod/{r["id"]}.html">{r["id"]}</a></td>'
                   f'<td>{comp}</td><td>{esc(reference(r))}{doi}</td><td>{esc(why)}{extra}</td></tr>')
    out.append(f"""</tbody></table></div>
<p class="small muted" style="margin-top:14px">{t['foot'].format(rel=release)}</p>
</div></main>
<script>
(function () {{
  var q = document.getElementById("q"), k = document.getElementById("k"), shown = document.getElementById("shown");
  var rows = Array.prototype.slice.call(document.querySelectorAll("#rows tr")), total = rows.length;
  var fmt = function (x) {{ return x.toLocaleString("{t['lang']}"); }};
  function apply() {{
    var s = q.value.trim().toLowerCase(), kk = k.value, n = 0;
    rows.forEach(function (tr) {{
      var ok = (!kk || tr.getAttribute("data-k") === kk) && (!s || tr.getAttribute("data-s").indexOf(s) >= 0);
      tr.hidden = !ok; if (ok) n++;
    }});
    shown.textContent = {json.dumps(t['shown'])}.replace("{{a}}", fmt(n)).replace("{{b}}", fmt(total));
  }}
  q.addEventListener("input", apply); k.addEventListener("change", apply);
  var h = decodeURIComponent(location.hash.slice(1)); if (h) q.value = h;
  apply();
}})();
</script>
</body>
</html>
""")
    return "\n".join(out)


def main():
    meta, release = Path(sys.argv[1]), sys.argv[2]
    rows = json.loads(meta.read_text())
    n = notes()
    for r in rows:
        for f in ("name", "formula"):  # CIF placeholders for unknown values
            if r[f].strip("'\" ") in ("", "?", "."):
                r[f] = ""
        r["kind"] = kind(r)
        r["note"] = n.get(r["id"], "")
    rows.sort(key=lambda r: (ORDER.index(r["kind"]), r["id"]))
    (HERE / "skipped.html").write_text(page(rows, "en", release))
    (HERE / "skipped.fr.html").write_text(page(rows, "fr", release))
    print({k: sum(1 for r in rows if r["kind"] == k) for k in ORDER}, sum(1 for r in rows if r["note"]), "with notes")


main()
