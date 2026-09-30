#!/usr/bin/env python3
"""
Génère les schémas et graphiques du README (docs/images/), en thème clair et sombre.

    python3 docs/generer_images.py

Les mesures proviennent de benchmarks/RESULTATS.md (make bench) : les mettre à jour ici
après une nouvelle campagne de mesures.
"""

import math
import os
from html import escape

DOSSIER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "images")

THEMES = {
    "clair": dict(fond="#f4f6f9", surface="#ffffff", encre="#16202e", discret="#566376", trait="#d8dee7",
                  accent="#2446b3", accent_fond="#e7ecfb", partage="#b86e0a", partage_fond="#fbefd9",
                  code="#eef1f6", bon="#1f7a4d", ginac="#a3aebd"),
    "sombre": dict(fond="#0e131a", surface="#151c26", encre="#e4e9f0", discret="#98a5b7", trait="#283241",
                   accent="#8fa8ff", accent_fond="#1a2340", partage="#efb24f", partage_fond="#2a2213",
                   code="#111821", bon="#5cc98f", ginac="#4b5667"),
}

POLICE = "Inter, 'Segoe UI', 'Helvetica Neue', Arial, sans-serif"
MONO = "'SFMono-Regular', Menlo, Consolas, 'DejaVu Sans Mono', 'Liberation Mono', monospace"


# ----------------------------------------------------------------------------
# Primitives
# ----------------------------------------------------------------------------

def texte(x, y, contenu, t, taille=13, couleur="encre", gras=False, mono=False, ancre="start", italique=False):
    poids = ' font-weight="700"' if gras else ""
    style = ' font-style="italic"' if italique else ""
    return (f'<text x="{x}" y="{y}" fill="{t[couleur]}" font-size="{taille}" font-family="{MONO if mono else POLICE}"'
            f'{poids}{style} text-anchor="{ancre}" xml:space="preserve">{escape(contenu)}</text>')


def lignes(x, y, contenus, t, interligne=16, **options):
    return "".join(texte(x, y + i * interligne, c, t, **options) for i, c in enumerate(contenus))


def rect(x, y, l, h, remplissage, contour=None, rayon=6, pointille=False, epaisseur=1):
    trait = f' stroke="{contour}" stroke-width="{epaisseur}"' if contour else ""
    tirets = ' stroke-dasharray="4 3"' if pointille else ""
    return f'<rect x="{x}" y="{y}" width="{l}" height="{h}" rx="{rayon}" fill="{remplissage}"{trait}{tirets}/>'


def fleche_bas(x, y1, y2, t, legende=None):
    s = f'<line x1="{x}" y1="{y1}" x2="{x}" y2="{y2 - 6}" stroke="{t["discret"]}" stroke-width="1.5"/>'
    s += f'<path d="M {x - 5} {y2 - 7} L {x + 5} {y2 - 7} L {x} {y2} Z" fill="{t["discret"]}"/>'
    if legende:
        s += texte(x + 12, (y1 + y2) / 2 + 4, legende, t, taille=11.5, couleur="discret", mono=True)
    return s


def fleche(x1, y1, x2, y2, t, couleur="discret"):
    angle = math.atan2(y2 - y1, x2 - x1)
    xa, ya = x2 - 7 * math.cos(angle), y2 - 7 * math.sin(angle)
    g = (xa - 4.5 * math.sin(angle), ya + 4.5 * math.cos(angle))
    d = (xa + 4.5 * math.sin(angle), ya - 4.5 * math.cos(angle))
    return (f'<line x1="{x1}" y1="{y1}" x2="{xa:.1f}" y2="{ya:.1f}" stroke="{t[couleur]}" stroke-width="1.5"/>'
            f'<path d="M {g[0]:.1f} {g[1]:.1f} L {d[0]:.1f} {d[1]:.1f} L {x2} {y2} Z" fill="{t[couleur]}"/>')


def document(largeur, hauteur, t, contenu, titre):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {largeur} {hauteur}" width="{largeur}" '
            f'height="{hauteur}" role="img"><title>{escape(titre)}</title>'
            f'{rect(0, 0, largeur, hauteur, t["fond"], rayon=12)}{contenu}</svg>\n')


def en_tete(t, surtitre, titre, sous_titres, x=40):
    s = texte(x, 42, surtitre.upper(), t, taille=11.5, couleur="accent", mono=True)
    s += texte(x, 72, titre, t, taille=24, gras=True)
    s += lignes(x, 98, sous_titres, t, interligne=19, taille=14, couleur="discret")
    return s


def ecrire(nom, generateur):
    for theme, t in THEMES.items():
        chemin = os.path.join(DOSSIER, f"{nom}-{theme}.svg")
        with open(chemin, "w", encoding="utf-8") as f:
            f.write(generateur(t))


# ----------------------------------------------------------------------------
# 1. Architecture en couches
# ----------------------------------------------------------------------------

COUCHES = [
    ("Entrée", "include/Lecture.hpp", "construit des expressions", [
        ("Lecture", ['lire("x^2 = 2")', "nombres exacts, erreurs situées"], False),
        ("Helpers C++", ['var("x"), ast_sin, ast_pow', "opérateurs + - * /"], False),
    ]),
    ("Équations", "include/Equation*.hpp", "résout avec", [
        ("Equation", ["abstraite", "eval · deriveeGenerique"], False),
        ("EquationClassique", ["dérivée · primitive · DL", "limite · résolution · tracé"], False),
        ("EquationDifferentielle", ["somme a_i y^(i) = 0", "littérale · Cauchy · RK4"], False),
    ]),
    ("Algèbre", "Polynome.hpp · Solveur.hpp", "manipule", [
        ("developper · factoriser", ["(x - 1)*(x^2 - 2)", "facteurs sur Q"], False),
        ("Polynome exact", ["suites de Sturm · Yun · PGCD", "racines réelles certifiées"], True),
        ("Solveur", ["isolement · produit nul", "familles en k · Brent"], False),
    ]),
    ("Expressions", "include/ASTNode.hpp", "calcule avec", [
        ("Forme canonique", ["Noeud.cpp : x + x = 2*x", "sommes et produits n-aires"], True),
        ("Hash-consing", ["une expression, un noeud", "égalité = même adresse"], True),
        ("Règles", ["Regles.cpp", "éval · dérivées · primitives"], False),
        ("Limites", ["Limites.cpp", "L'Hôpital · signe de l'infini"], False),
        ("Séries", ["Series.cpp", "DL par séries tronquées"], False),
        ("Affichage", ["Affichage.cpp", "x^2 + 5*x + 6, relisible"], False),
    ]),
    ("Calcul numérique", "Nombre · Evaluateur · Eigen", None, [
        ("Nombre", ["rationnel exact int64 → GMP", "ou réel double"], True),
        ("ProgrammeEvaluation", ["Evaluateur.cpp", "programme compilé, blocs"], True),
        ("Eigen", ["matrice compagnon", "valeurs propres · RK4"], False),
        ("Ref", ["Ref.hpp", "compteur intrusif"], False),
    ]),
]


def architecture(t):
    largeur, x0, xb = 1000, 40, 250
    colonnes, ecart, hb = 3, 12, 62
    lb = (largeur - 40 - xb - (colonnes - 1) * ecart) // colonnes
    s = en_tete(t, "Vue d'ensemble", "Cinq couches, chacune au service de la suivante",
                ["Les équations sont l'interface ; elles s'appuient sur l'algèbre et sur des expressions",
                 "symboliques uniques en mémoire, qui calculent avec des nombres exacts. En ambre : les idées clés."])
    y = 140
    for i, (nom, chemin, lien, boites) in enumerate(COUCHES):
        rangees = (len(boites) + colonnes - 1) // colonnes
        hauteur = 28 + rangees * hb + (rangees - 1) * ecart
        s += rect(x0, y, largeur - 2 * x0, hauteur, t["surface"], t["trait"], rayon=10)
        s += texte(x0 + 20, y + 34, nom, t, taille=16, gras=True)
        s += texte(x0 + 20, y + 54, chemin, t, taille=11, couleur="discret", mono=True)
        for j, (titre, details, cle) in enumerate(boites):
            bx = xb + (j % colonnes) * (lb + ecart)
            by = y + 14 + (j // colonnes) * (hb + ecart)
            s += rect(bx, by, lb, hb, t["partage_fond"] if cle else t["code"],
                      t["partage"] if cle else t["trait"], rayon=6)
            s += texte(bx + 12, by + 21, titre, t, taille=13, gras=True, couleur="partage" if cle else "encre")
            s += lignes(bx + 12, by + 38, details, t, interligne=15, taille=11, couleur="discret", mono=True)
        y += hauteur
        if lien:
            s += fleche_bas(largeur / 2, y + 4, y + 30, t, lien)
            y += 34
    return document(largeur, y + 30, t, s, "Architecture de SymAlgo++ en cinq couches")


# ----------------------------------------------------------------------------
# 2. Construction d'une expression
# ----------------------------------------------------------------------------

ETAPES = [
    ("Vous écrivez", ['lire("x*sin(x)', '  + x*sin(x)")'], ["Texte ou helpers C++,", "même résultat."], False),
    ("Regroupement", ["x·sin(x) : 1", "x·sin(x) : 1 + 1 = 2"],
     ["Termes identiques", "fusionnés, coefficients", "exacts."], False),
    ("Signature", ["Produit", "coefficient 2", "facteurs [x, sin(x)]"],
     ["Type, valeurs et", "adresses des enfants."], False),
    ("Recherche", ["existe déjà ?", "  oui → réutilisé", "  non → alloué"],
     ["Table de hachage", "intrusive : une", "allocation au plus."], False),
    ("Un ExprPtr", ["2*x*sin(x)"], ["3 noeuds en mémoire.", "Réécrire la même", "expression rend le", "même pointeur."], True),
]


def construction(t):
    largeur, x0, ecart = 1000, 40, 12
    lc = (largeur - 2 * x0 - 4 * ecart) / 5
    s = en_tete(t, "Construction", "Ce qui se passe quand vous écrivez x*sin(x) + x*sin(x)",
                ["Aucune étape de simplification séparée : chaque opération produit directement la forme",
                 "réduite, puis vérifie si elle existe déjà en mémoire."])
    y, h = 140, 222
    for i, (titre, code, legende, final) in enumerate(ETAPES):
        cx = x0 + i * (lc + ecart)
        s += rect(cx, y, lc, h, t["partage_fond"] if final else t["surface"],
                  t["partage"] if final else t["trait"], rayon=8)
        s += texte(cx + 14, y + 30, str(i + 1), t, taille=20, gras=True, couleur="partage" if final else "accent")
        s += texte(cx + 14, y + 56, titre, t, taille=14, gras=True)
        s += rect(cx + 12, y + 70, lc - 24, 18 + 16 * len(code), t["code"], rayon=4)
        s += lignes(cx + 20, y + 88, code, t, interligne=16, taille=11, mono=True)
        s += lignes(cx + 14, y + 110 + 16 * len(code), legende, t, interligne=17, taille=12, couleur="discret")
        if i < 4:
            s += fleche(cx + lc + 1, y + h / 2, cx + lc + ecart - 1, y + h / 2, t)
    return document(largeur, y + h + 30, t, s, "Construction canonique et hash-consing")


# ----------------------------------------------------------------------------
# 3. Graphe partagé
# ----------------------------------------------------------------------------

def noeud(x, y, l, etiquette, t, genre="normal"):
    remplissage, contour, pointille = t["surface"], t["discret"], False
    if genre == "copie":
        pointille = True
    elif genre == "partage":
        remplissage, contour = t["partage_fond"], t["partage"]
    elif genre == "racine":
        remplissage, contour = t["accent_fond"], t["accent"]
    s = rect(x, y, l, 24, remplissage, contour, rayon=4, pointille=pointille, epaisseur=1.5 if genre == "partage" else 1)
    return s + texte(x + l / 2, y + 17, etiquette, t, taille=12, mono=True, ancre="middle",
                     couleur="partage" if genre == "partage" else "encre")


def arete(x1, y1, x2, y2, t, couleur="discret"):
    return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{t[couleur]}" stroke-width="1.2"/>'


def graphe(t):
    s = en_tete(t, "Mémoire", "Un graphe partagé, pas un arbre",
                ["À gauche, l'expression telle qu'on l'écrit : un arbre qui recopie x quatre fois. À droite, ce que",
                 "SymAlgo++ garde en mémoire, avec la dérivée f′ : x et sin(x) servent aux deux sans être copiés."])
    # Arbre (x*sin(x) + x*sin(x))
    a = '<g transform="translate(60 140) scale(1.2)">'
    for x1, y1, x2, y2 in [(180, 42, 100, 88), (180, 42, 260, 88), (100, 112, 60, 158), (100, 112, 140, 158),
                           (260, 112, 220, 158), (260, 112, 300, 158), (140, 182, 140, 233), (300, 182, 300, 233)]:
        a += arete(x1, y1, x2, y2, t)
    for x, y, l, e, g in [(162, 18, 36, "+", "normal"), (82, 88, 36, "×", "normal"), (242, 88, 36, "×", "copie"),
                          (42, 158, 36, "x", "normal"), (118, 158, 44, "sin", "normal"), (202, 158, 36, "x", "copie"),
                          (278, 158, 44, "sin", "copie"), (122, 233, 36, "x", "copie"), (282, 233, 36, "x", "copie")]:
        a += noeud(x, y, l, e, t, g)
    a += "</g>"
    # Graphe partagé de f et f′
    b = '<g transform="translate(530 140) scale(1.1)">'
    b += texte(95, 22, "f = 2*x*sin(x)", t, taille=12.5, gras=True, ancre="middle")
    b += texte(300, 22, "f′ = 2*x*cos(x) + 2*sin(x)", t, taille=12.5, gras=True, ancre="middle")
    b += arete(95, 72, 150, 138, t)
    b += f'<path d="M 80 72 C 50 170, 110 255, 186 255" fill="none" stroke="{t["discret"]}" stroke-width="1.2"/>'
    b += arete(150, 162, 196, 243, t)
    for x1, y1, x2, y2 in [(300, 72, 160, 138), (300, 72, 330, 128), (330, 152, 210, 243), (330, 152, 330, 203),
                           (330, 227, 214, 252)]:
        b += arete(x1, y1, x2, y2, t, "accent")
    b += texte(222, 98, "×2", t, taille=11, couleur="accent", mono=True, ancre="middle")
    b += texte(336, 104, "×2", t, taille=11, couleur="accent", mono=True)
    for x, y, l, e, g in [(40, 48, 110, "Produit ×2", "racine"), (262, 48, 76, "Somme", "racine"),
                          (296, 128, 68, "Produit", "normal"), (308, 203, 44, "cos", "normal"),
                          (126, 138, 48, "sin", "partage"), (186, 243, 36, "x", "partage")]:
        b += noeud(x, y, l, e, t, g)
    b += "</g>"
    s += a + b
    s += f'<line x1="500" y1="150" x2="500" y2="480" stroke="{t["trait"]}"/>'
    s += texte(60, 500, "Tel qu'écrit : 9 noeuds", t, taille=15, gras=True)
    s += texte(60, 520, "En pointillés : les copies de x, sin(x) et du produit.", t, taille=12.5, couleur="discret")
    s += texte(530, 500, "En mémoire : 6 noeuds pour f et f′", t, taille=15, gras=True)
    s += lignes(530, 520, ["En ambre : les noeuds partagés. Arêtes bleues : celles de f′.",
                           "Les coefficients ×2 sont portés par les arêtes, pas par des noeuds."],
                t, interligne=17, taille=12.5, couleur="discret")
    return document(1000, 560, t, s, "Graphe partagé contre arbre")


# ----------------------------------------------------------------------------
# 4. Résolution d'équations
# ----------------------------------------------------------------------------

def resolution(t):
    s = en_tete(t, "Algèbre", "Résoudre une équation",
                ["L'inconnue est isolée en inversant les opérations ; les polynômes passent par des racines",
                 "certifiées. Chaque solution est exacte quand c'est possible, avec sa valeur approchée."])

    def etape(x, y, l, h, titre, details, genre="normal"):
        fond = {"normal": t["surface"], "cle": t["partage_fond"], "fin": t["accent_fond"]}[genre]
        bord = {"normal": t["trait"], "cle": t["partage"], "fin": t["accent"]}[genre]
        r = rect(x, y, l, h, fond, bord, rayon=8)
        r += texte(x + 14, y + 22, titre, t, taille=13, gras=True, couleur="partage" if genre == "cle" else "encre")
        return r + lignes(x + 14, y + 40, details, t, interligne=15, taille=11, couleur="discret", mono=True)

    x, l = 40, 280
    s += etape(x, 140, l, 58, "gauche = droite", ["différence développée ;", "sans x : identité, ou rien"])
    s += fleche(x + l / 2, 198, x + l / 2, 222, t)
    s += etape(x, 222, l, 58, "Polynôme en x ?", ["aussi après changement de", "variable : exp(x), sin(x)…"])
    # Branche polynôme
    s += fleche(x + l, 251, 360, 251, t)
    s += texte(x + l + 6, 244, "oui", t, taille=11, couleur="discret", mono=True)
    s += etape(360, 208, 300, 88, "Polynome exact", ["sans carré (Yun)", "isolement par suites de Sturm",
                                                    "rationnelles · radicaux (degré 2)", "sinon le double le plus proche"], "cle")
    # Branche isolement
    s += fleche(x + l / 2, 280, x + l / 2, 304, t)
    s += texte(x + l / 2 + 8, 297, "non", t, taille=11, couleur="discret", mono=True)
    s += etape(x, 304, l, 88, "Isolement de x", ["somme · produit nul · u^n", "exp · ln · asin · acos · atan",
                                                 "sin, cos, tan → familles en k"])
    # Nettoyage
    s += fleche(x + l / 2, 392, x + l / 2, 416, t)
    s += fleche(510, 296, 510, 416, t)
    s += etape(x, 416, 620, 58, "Nettoyage", ["domaine vérifié (x*ln(x) = 0 ne garde que x = 1), doublons, tri"])
    s += fleche(x + l / 2, 474, x + l / 2, 498, t)
    s += etape(x, 498, 620, 58, "Solutions { liste, complet, toutReel }",
               ["complet = faux : les solutions listées sont justes, il peut en manquer"], "fin")
    s += fleche(x + l / 2, 556, x + l / 2, 580, t)
    s += etape(x, 580, 620, 58, "resoudreSurIntervalle(a, b)",
               ["familles dépliées sur [a, b] ; méthode de Brent si la résolution est incomplète"])
    # Exemples
    ex, ly = 690, 140
    s += rect(ex, ly, 270, 498, t["surface"], t["trait"], rayon=10)
    s += texte(ex + 18, ly + 30, "Exemples réels", t, taille=14, gras=True)
    exemples = [
        ("x^2 = 2", ["-2^(1/2) ; 2^(1/2)"]),
        ("x^3 - x^2 - 2x + 2 = 0", ["-2^(1/2) ; 1 ; 2^(1/2)"]),
        ("exp(2x) - 3exp(x) + 2 = 0", ["0 ; ln(2)"]),
        ("2^x = 8", ["3"]),
        ("sin(x) = 1/2", ["pi/6 + 2*pi*k", "5*pi/6 + 2*pi*k"]),
        ("exp(x) + x = 0", ["≈ -0.567143290409784", "(numérique, Brent)"]),
    ]
    yy = ly + 62
    for equation, solutions in exemples:
        s += texte(ex + 18, yy, equation, t, taille=12, mono=True)
        s += lignes(ex + 30, yy + 20, ["→ " + solutions[0]] + ["  " + v for v in solutions[1:]], t,
                    interligne=16, taille=12, mono=True, couleur="bon")
        yy += 40 + 16 * len(solutions)
    return document(1000, 668, t, s, "Résolution d'équations")


# ----------------------------------------------------------------------------
# 5. Performances face à GiNaC
# ----------------------------------------------------------------------------

# (scénario, SymAlgo++, GiNaC, facteur) : médianes de benchmarks/RESULTATS.md
MESURES = [
    ("Grande expression, un point", "1,1 µs", "1 607 µs", 1450),
    ("Série de Taylor d'ordre 10", "4,1 µs", "4 784 µs", 1170),
    ("10 000 points (par blocs)", "0,61 ms", "302 ms", 494),
    ("Évaluation en un point", "123 ns", "30 820 ns", 250),
    ("Dérivée 10e de exp(sin x)·x²", "1,38 ms", "7,45 ms", 5.4),
    ("Lecture d'un texte (73 car.)", "15,2 µs", "54,3 µs", 3.6),
    ("Dérivée première", "8,7 µs", "28,1 µs", 3.2),
    ("Collecte de 100 termes", "146 µs", "255 µs", 1.7),
    ("Mémoire du résultat (dérivée 6e)", "4,2 Ko", "4,9 Ko", 4.9 / 4.2),
]


def performances(t):
    largeur, xl, xb, xf, xv = 1000, 40, 300, 700, 760
    s = en_tete(t, "Mesures", "SymAlgo++ face à GiNaC, sur les mêmes expressions",
                ["Facteur d'avance de SymAlgo++, en échelle logarithmique. Médiane de 3 répétitions,",
                 "Intel Celeron N4120, GCC 16, -O3 ; mémoire comptée de la même façon pour les deux."])
    y0, pas = 150, 38
    maximum = math.log10(2000)
    echelle = lambda f: xb + (xf - xb) * math.log10(f) / maximum
    for graduation in (1, 10, 100, 1000):
        gx = echelle(graduation)
        s += f'<line x1="{gx:.1f}" y1="{y0 - 8}" x2="{gx:.1f}" y2="{y0 + pas * len(MESURES)}" stroke="{t["trait"]}" stroke-dasharray="3 3"/>'
        s += texte(gx, y0 + pas * len(MESURES) + 18, f"×{graduation}", t, taille=11, couleur="discret", mono=True, ancre="middle")
    s += texte(xv, y0 - 16, "SymAlgo++ · GiNaC", t, taille=11, couleur="discret", mono=True)
    for i, (nom, nous, eux, facteur) in enumerate(MESURES):
        y = y0 + i * pas
        s += texte(xl, y + 17, nom, t, taille=13)
        fin = max(echelle(facteur), xb + 3)
        s += rect(xb, y + 4, fin - xb, 18, t["accent"], rayon=3)
        etiquette = f"×{facteur:,.0f}".replace(",", " ") if facteur >= 10 else (
            "−14 %" if facteur < 1.5 and "Mémoire" in nom else f"×{facteur:.1f}".replace(".", ","))
        s += texte(fin + 8, y + 18, etiquette, t, taille=12.5, gras=True, couleur="accent")
        s += texte(xv, y + 17, nous, t, taille=12, mono=True, couleur="encre")
        s += texte(xv + 80, y + 17, eux, t, taille=12, mono=True, couleur="discret")
    return document(largeur, y0 + pas * len(MESURES) + 46, t, s, "Performances de SymAlgo++ face à GiNaC")


# ----------------------------------------------------------------------------
# 6. Matrice compagnon des EDO
# ----------------------------------------------------------------------------

# Ordre de l'EDO -> temps médian (ns) de getMatriceCompagnon()
COMPAGNON = [(1, 40.5), (3, 61.9), (5, 97.4), (7, 137), (9, 178), (11, 223), (13, 342), (15, 384),
             (17, 452), (19, 487), (21, 550), (23, 614), (25, 671), (27, 753), (29, 841)]


def compagnon(t):
    largeur, hauteur = 1000, 430
    gx, dx, hy, by = 100, 950, 140, 370
    s = en_tete(t, "Équations différentielles", "Génération de la matrice compagnon",
                ["Temps médian selon l'ordre de l'EDO : moins d'une microseconde jusqu'à l'ordre 29."])
    px = lambda o: gx + (dx - gx) * (o - 1) / 28
    py = lambda v: by - (by - hy) * v / 900
    for v in range(0, 901, 300):
        s += f'<line x1="{gx}" y1="{py(v):.1f}" x2="{dx}" y2="{py(v):.1f}" stroke="{t["trait"]}" stroke-dasharray="3 3"/>'
        s += texte(gx - 10, py(v) + 4, f"{v} ns", t, taille=11, couleur="discret", mono=True, ancre="end")
    for o in range(1, 30, 4):
        s += texte(px(o), by + 22, f"ordre {o}", t, taille=11, couleur="discret", mono=True, ancre="middle")
    chemin = " ".join(f"{'M' if i == 0 else 'L'} {px(o):.1f} {py(v):.1f}" for i, (o, v) in enumerate(COMPAGNON))
    s += f'<path d="{chemin}" fill="none" stroke="{t["accent"]}" stroke-width="2.5" stroke-linejoin="round"/>'
    for o, v in COMPAGNON:
        s += f'<circle cx="{px(o):.1f}" cy="{py(v):.1f}" r="4" fill="{t["accent"]}" stroke="{t["fond"]}" stroke-width="2"/>'
    for o, v in (COMPAGNON[0], COMPAGNON[7], COMPAGNON[-1]):
        s += texte(px(o), py(v) - 12, f"{v:g} ns".replace(".", ","), t, taille=11.5, gras=True, couleur="accent", mono=True, ancre="middle")
    return document(largeur, hauteur, t, s, "Temps de génération de la matrice compagnon")


if __name__ == "__main__":
    os.makedirs(DOSSIER, exist_ok=True)
    for nom, generateur in [("architecture", architecture), ("construction", construction), ("graphe-partage", graphe),
                            ("resolution", resolution), ("performances", performances), ("matrice-compagnon", compagnon)]:
        ecrire(nom, generateur)
    print("Images générées dans", DOSSIER)
