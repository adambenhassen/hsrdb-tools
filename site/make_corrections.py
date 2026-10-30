"""corrections.html and corrections.fr.html from docs/corrections.md: one card per corrected CIF (COD links,
compound, problem in the CIF, correction, source, check), plus the entries kept as the official database has them.
The texts stay in English on the French page (noted there); the headings are translated. Style from index.html."""
import html
import re
from pathlib import Path

HERE = Path(__file__).parent
REPO = HERE.parent
NB = " "


def inline(md):
    """The Markdown of a table cell as HTML: links, bold, code."""
    parts, pos = [], 0
    for m in re.finditer(r"`([^`]+)`|\*\*([^*]+)\*\*|\[([^\]]+)\]\(((?:[^()\s]|\([^()\s]*\))+)\)", md):  # DOIs hold balanced ()
        parts.append(html.escape(md[pos:m.start()]))
        if m.group(1) is not None:
            parts.append(f"<code>{html.escape(m.group(1))}</code>")
        elif m.group(2) is not None:
            parts.append(f"<b>{inline(m.group(2))}</b>")
        else:
            url = m.group(4)
            if url.startswith("../") or url.endswith(".md") or url.startswith("#"):
                parts.append(inline(m.group(3)))  # repository-internal link: text only
            else:
                parts.append(f'<a href="{html.escape(url)}">{inline(m.group(3))}</a>')
        pos = m.end()
    parts.append(html.escape(md[pos:]))
    return "".join(parts)


def cells(line):
    return [c.strip() for c in line.strip().strip("|").split(" | ")]


def tables():
    text = (REPO / "docs" / "corrections.md").read_text()
    main, rest = text.split("## Skipped CIFs checked against their articles", 1)
    fixes = [cells(l) for l in main.splitlines() if l.startswith("| [")]
    kept = [cells(l) for l in rest.splitlines() if l.startswith("| [") and "Built as in COD24" in l]
    return fixes, kept


def style():
    src = (HERE / "index.html").read_text()
    return src[src.index("<style>"):src.index("</style>") + len("</style>")]


T = {
    "en": dict(lang="en", title="Corrected CIFs", home="index.html", other="corrections.fr.html", other_lang="fr",
               other_label="FR", brand="COD for HighScore <span>· corrections</span>", back="← Back to the databases",
               h1="CIFs corrected from their publications",
               lede="{n} COD CIFs cannot give a correct pattern as deposited, but the original article, or the CIF's own statements, say what the file leaves out or gets wrong. For these the build applies a correction after reading the CIF; the CIF files are never changed, and each entry says in its comment in HighScore that it was corrected and from which source. A correction is added only when the source itself gives the missing information, and each one is checked against the CIF's formula, Z and density. The exception are the corrections marked “not checked against the article”: CIFs whose transcribed structure contradicts itself (atoms on top of each other, or a content against the CIF's own formula), read with a corrected structure of the same publication.",
               search="Search COD id, compound or source", labels=("Problem in the CIF", "Correction", "Source", "Check"),
               kept_h="Kept as in the official database",
               kept_p="These entries of the official COD 2024 database have no usable coordinates in their current CIFs. They are included as the official database has them, so that nothing it contains is missing, and each entry's comment says that it is not a correct structure.",
               article="Article", why="Why", note="", shown="{a} of {b} shown"),
    "fr": dict(lang="fr", title="CIF corrigés", home="fr.html", other="corrections.html", other_lang="en",
               other_label="EN", brand="COD pour HighScore <span>· corrections</span>", back="← Retour aux bases",
               h1="CIF corrigés d’après leur publication",
               lede="{n} CIF de la COD ne donnent pas de diagramme correct tels quels, mais l’article d’origine, ou les indications du CIF lui-même, précisent ce que le fichier omet ou donne de travers. Pour ceux-ci, la construction applique une correction après la lecture du CIF ; les fichiers CIF ne sont jamais modifiés, et le commentaire de chaque fiche dans HighScore indique qu’elle a été corrigée et d’après quelle source. Une correction n’est ajoutée que si la source donne elle-même l’information manquante, et chacune est vérifiée par rapport à la formule, au Z et à la densité du CIF. Font exception les corrections marquées « not checked against the article » : des CIF dont la structure transcrite se contredit (atomes superposés, ou contenu contraire à la formule du CIF lui-même), lus avec une structure corrigée de la même publication.",
               search="Rechercher un numéro COD, un composé ou une source", labels=("Problème dans le CIF", "Correction", "Source", "Vérification"),
               kept_h="Reprises telles que dans la base officielle",
               kept_p="Ces fiches de la base officielle COD 2024 n’ont pas de coordonnées utilisables dans leur CIF actuel. Elles sont reprises telles que la base officielle les contient, afin que rien de ce qu’elle contient ne manque, et le commentaire de chaque fiche indique qu’il ne s’agit pas d’une structure correcte.",
               article="Article", why="Raison", note="Le détail de chaque correction est en anglais.", shown="{a} sur {b} affichés"),
}


def page(fixes, kept, lang, release):
    t = T[lang]
    n = sum(len(re.findall(r"\[(\d{7})\]", f[0])) for f in fixes)
    count = f"{n:,}".replace(",", NB) if lang == "fr" else f"{n:,}"
    out = [f"""<!doctype html>
<html lang="{t['lang']}">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{t['title']}</title>
<link rel="alternate" hreflang="en" href="https://cod.adambh.net/corrections.html">
<link rel="alternate" hreflang="fr" href="https://cod.adambh.net/corrections.fr.html">
{style()}
<style>
.filters {{ margin: 16px 0 12px; }}
.filters input {{ font: inherit; font-size: 0.95rem; padding: 7px 10px; border: 1px solid var(--border); border-radius: 8px; width: 100%; max-width: 520px; background: var(--surface); color: var(--text); }}
.fix {{ background: var(--surface); border: 1px solid var(--border); border-radius: 12px; padding: 14px 16px; margin: 12px 0; }}
.fix h3 {{ margin: 0 0 8px; font-size: 1.05rem; overflow-wrap: anywhere; }}
.fix h3 .ids {{ font-family: system-ui, sans-serif; font-size: 0.9rem; font-weight: 600; }}
.fix dl {{ display: grid; grid-template-columns: 170px 1fr; gap: 6px 14px; margin: 0; font-size: 0.93rem; }}
.fix dt {{ color: var(--text-2); font-weight: 600; font-size: 0.85rem; }}
.fix dd {{ margin: 0; min-width: 0; overflow-wrap: anywhere; }}
@media (max-width: 640px) {{ .fix dl {{ grid-template-columns: 1fr; }} .fix dd {{ margin-bottom: 6px; }} }}
.langlink {{ font-size: 0.85rem; font-weight: 600; border: 1px solid var(--border); border-radius: 8px; padding: 3px 10px; text-decoration: none; color: var(--text-2); }}
</style>
</head>
<body>
<header class="top"><div class="wrap"><div class="brand">{t['brand']}</div><div><a class="langlink" href="{t['other']}" hreflang="{t['other_lang']}" lang="{t['other_lang']}">{t['other_label']}</a></div></div></header>
<main><div class="wrap" style="padding-top:28px;padding-bottom:40px">
<p><a href="{t['home']}">{t['back']}</a></p>
<h1>{t['h1']}</h1>
<p class="lede">{html.escape(t['lede'].format(n=count))}</p>
{f'<p class="small muted">{t["note"]}</p>' if t['note'] else ''}
<div class="filters"><input id="q" type="search" placeholder="{html.escape(t['search'])}" aria-label="{html.escape(t['search'])}"></div>
<p class="small muted" id="shown" aria-live="polite"></p>
<div id="fixes">"""]
    for f in fixes:
        ids, compound, problem, correction, source, check = f
        search = html.escape(" ".join(re.findall(r"\d{7}", ids) + [re.sub(r"[`*\[\]()]", " ", compound + " " + source)]).lower())
        out.append(f'<article class="fix" data-s="{search}"><h3><span class="ids">{inline(ids)}</span> · {inline(compound)}</h3><dl>'
                   + "".join(f"<dt>{lab}</dt><dd>{inline(v)}</dd>" for lab, v in zip(t["labels"], (problem, correction, source, check)))
                   + "</dl></article>")
    out.append(f"""</div>
<h2 style="margin-top:32px">{t['kept_h']}</h2>
<p>{html.escape(t['kept_p'])}</p>""")
    for k in kept:
        ids, compound, article, why = k
        out.append(f'<article class="fix"><h3><span class="ids">{inline(ids)}</span> · {inline(compound)}</h3><dl>'
                   f"<dt>{t['article']}</dt><dd>{inline(article)}</dd><dt>{t['why']}</dt><dd>{inline(why)}</dd></dl></article>")
    out.append(f"""
</div></main>
<script>
(function () {{
  var q = document.getElementById("q"), shown = document.getElementById("shown");
  var cards = Array.prototype.slice.call(document.querySelectorAll("#fixes .fix")), total = cards.length;
  function apply() {{
    var s = q.value.trim().toLowerCase(), n = 0;
    cards.forEach(function (c) {{ var ok = !s || c.getAttribute("data-s").indexOf(s) >= 0 || c.textContent.toLowerCase().indexOf(s) >= 0; c.hidden = !ok; if (ok) n++; }});
    shown.textContent = {repr(t['shown'])}.replace("{{a}}", n).replace("{{b}}", total);
  }}
  q.addEventListener("input", apply);
  var h = decodeURIComponent(location.hash.slice(1)); if (h) q.value = h;
  apply();
}})();
</script>
</body>
</html>
""")
    return "\n".join(out)


def main():
    import sys
    release = sys.argv[1] if len(sys.argv) > 1 else ""
    fixes, kept = tables()
    for lang, name in (("en", "corrections.html"), ("fr", "corrections.fr.html")):
        (HERE / name).write_text(page(fixes, kept, lang, release))
    print(len(fixes), "rows,", sum(len(re.findall(r"\[(\d{7})\]", f[0])) for f in fixes), "CIFs;", len(kept), "kept as COD24")


main()
