"""Build submission versions of the two papers (Markdown + PDF).

The drafts in paper/ stay the editable sources. This script removes the draft
status lines, writes a formal byline, and turns relative links into full
GitHub URLs so they work in a PDF. Output: paper/submission/ and paper/pdf/.
Run: python3 tools/build_submission.py   (needs pandoc and Playwright/Chromium)
"""
import os, re, subprocess, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = "https://github.com/Gokyona/life-origins-traversal"
DATE = "30 October 2026"
TITLE1 = "COMPASS: mapping where life-like behaviour appears in stochastic models, with an ordered, non-repeating active search"
TITLE2 = "How far up the ladder of life can environments carry a bottom-up chemistry? An ordered search over 13 environmental variables"

PAPERS = [
    ("COMPASS_method.md", "Paper1_COMPASS_method",
     f"Companion paper: *{TITLE2}* ({REPO}/blob/main/paper/COMPASS_stage4.md)"),
    ("COMPASS_stage4.md", "Paper2_Stage4_environments",
     f"Method paper: *{TITLE1}*, Zenodo preprint, doi:[10.5281/zenodo.23825371](https://doi.org/10.5281/zenodo.23825371)"),
]

def absolutise(text):
    def fix(m):
        label, target = m.group(1), m.group(2)
        if re.match(r"(https?:|mailto:|#)", target):
            return m.group(0)
        path, _, anchor = target.partition("#")
        full = os.path.normpath(os.path.join("paper", path)).replace(os.sep, "/")
        return f"[{label}]({REPO}/blob/main/{full}{'#' + anchor if anchor else ''})"
    return re.sub(r"\[([^\]]+)\]\(([^)\s]+)\)", fix, text)

def transform(src, companion):
    lines = src.split("\n")
    assert lines[0].startswith("# "), "title line expected"
    out, i = [lines[0], ""], 1
    # skip the draft byline, links line and status blockquote up to the first rule
    while i < len(lines) and lines[i].strip() != "---":
        i += 1
    out += [
        f"**Gokyona Wahiro** (ORCID [0009-0008-1514-8856](https://orcid.org/0009-0008-1514-8856)) · Independent researcher, Toronto, Canada · {DATE}",
        "",
        f"Code, raw data and detailed stage reports: <{REPO}> · Text: CC BY 4.0 · Code: MIT",
        "",
        companion,
        "",
    ]
    body = "\n".join(lines[i:])
    return absolutise("\n".join(out) + "\n" + body)

def main():
    sub = os.path.join(ROOT, "paper", "submission")
    pdf = os.path.join(ROOT, "paper", "pdf")
    os.makedirs(sub, exist_ok=True); os.makedirs(pdf, exist_ok=True)
    css = os.path.join(ROOT, "tools", "paper.css")
    names = []
    for src, name, companion in PAPERS:
        with open(os.path.join(ROOT, "paper", src), encoding="utf-8") as f:
            md = transform(f.read(), companion)
        for bad in ("Draft 1", "Draft 2", "for review", "Status of this draft"):
            if bad in md:
                sys.exit(f"{name}: draft marker left: {bad!r}")
        mdp = os.path.join(sub, name + ".md")
        with open(mdp, "w", encoding="utf-8") as f:
            f.write(md)
        html = os.path.join(sub, name + ".html")
        subprocess.run(["pandoc", mdp, "-f", "gfm", "-t", "html5", "-s",
                        "--metadata", "pagetitle=" + name, "--css", css,
                        "--embed-resources", "-o", html], check=True)
        names.append((html, os.path.join(pdf, name + ".pdf")))
    js = os.path.join(ROOT, "tools", "html2pdf.js")
    for html, out in names:
        subprocess.run(["node", js, html, out], check=True)
        os.remove(html)
        print("wrote", os.path.relpath(out, ROOT))

if __name__ == "__main__":
    main()
