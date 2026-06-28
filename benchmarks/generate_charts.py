import os

def create_comparison_svg(output_path):
    # Data
    # Eval: SymAlgo++ 260 ns, GiNaC 30376 ns
    # Deriv: SymAlgo++ 50486 ns, GiNaC 28534 ns
    # We will use a split or logarithmic visual style, or two distinct sub-charts because of the 116x difference.
    # Let's make a beautiful side-by-side or two-panel SVG.
    svg = """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 800 450" width="100%" height="100%" style="background-color: #1e1e2e; font-family: system-ui, -apple-system, sans-serif; border-radius: 12px; box-shadow: 0 8px 30px rgba(0,0,0,0.3);">
  <defs>
    <!-- Gradients -->
    <linearGradient id="symalgo-grad" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#89b4fa" />
      <stop offset="100%" stop-color="#b4befe" />
    </linearGradient>
    <linearGradient id="ginac-grad" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="#f38ba8" />
      <stop offset="100%" stop-color="#cba6f7" />
    </linearGradient>
    <linearGradient id="bg-grad" x1="0%" y1="0%" x2="0%" y2="100%">
      <stop offset="0%" stop-color="#1e1e2e" />
      <stop offset="100%" stop-color="#11111b" />
    </linearGradient>
  </defs>

  <!-- Background -->
  <rect width="100%" height="100%" fill="url(#bg-grad)" rx="12" />

  <!-- Title -->
  <text x="40" y="45" fill="#cdd6f4" font-size="22" font-weight="bold" letter-spacing="0.5">Benchmark 1 : SymAlgo++ vs GiNaC</text>
  <text x="40" y="65" fill="#a6adc8" font-size="13">Temps d'exécution en nanosecondes (inférieur est meilleur)</text>

  <!-- --- PARTIE 1: EVALUATION (Échelle Logarithmique ou indicateur de facteur) --- -->
  <g transform="translate(0, 0)">
    <!-- Subtitle -->
    <text x="40" y="115" fill="#cdd6f4" font-size="16" font-weight="semibold">1. Évaluation Numérique (x = 5)</text>
    
    <!-- SymAlgo++ Bar -->
    <text x="40" y="150" fill="#a6adc8" font-size="13">SymAlgo++</text>
    <rect x="150" y="136" width="6" height="20" rx="3" fill="url(#symalgo-grad)" />
    <text x="166" y="151" fill="#89b4fa" font-size="13" font-weight="bold">260 ns</text>

    <!-- GiNaC Bar -->
    <text x="40" y="195" fill="#a6adc8" font-size="13">GiNaC</text>
    <rect x="150" y="181" width="500" height="20" rx="3" fill="url(#ginac-grad)" />
    <text x="660" y="196" fill="#f38ba8" font-size="13" font-weight="bold">30 376 ns</text>
    
    <!-- Highlight Badge -->
    <rect x="40" y="220" width="220" height="28" rx="6" fill="#252538" stroke="#313244" stroke-width="1"/>
    <text x="50" y="238" fill="#a6e3a1" font-size="12" font-weight="bold">🚀 SymAlgo++ est ~116x plus rapide</text>
  </g>

  <!-- --- PARTIE 2: DERIVATION --- -->
  <g transform="translate(0, 160)">
    <!-- Subtitle -->
    <text x="40" y="145" fill="#cdd6f4" font-size="16" font-weight="semibold">2. Dérivation Symbolique f'(x)</text>
    
    <!-- SymAlgo++ Bar -->
    <text x="40" y="180" fill="#a6adc8" font-size="13">SymAlgo++</text>
    <rect x="150" y="166" width="379" height="20" rx="3" fill="url(#symalgo-grad)" />
    <text x="539" y="181" fill="#89b4fa" font-size="13" font-weight="bold">50 486 ns</text>

    <!-- GiNaC Bar -->
    <text x="40" y="225" fill="#a6adc8" font-size="13">GiNaC</text>
    <rect x="150" y="211" width="214" height="20" rx="3" fill="url(#ginac-grad)" />
    <text x="374" y="226" fill="#f38ba8" font-size="13" font-weight="bold">28 534 ns</text>
    
    <!-- Highlight Badge -->
    <rect x="40" y="250" width="200" height="28" rx="6" fill="#252538" stroke="#313244" stroke-width="1"/>
    <text x="50" y="268" fill="#f38ba8" font-size="12" font-weight="bold">⚠️ GiNaC est ~1.7x plus rapide</text>
  </g>
</svg>"""
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(svg)

def create_scaling_svg(output_path):
    # Data for ODE scaling:
    # 1: 1283
    # 3: 2632
    # 5: 4034
    # 7: 5720
    # 9: 7587
    # 11: 9158
    # 13: 11017
    # 15: 12805
    points = [
        (1, 1283), (3, 2632), (5, 4034), (7, 5720), 
        (9, 7587), (11, 9158), (13, 11017), (15, 12805)
    ]
    
    # SVG Dimensions
    w, h = 800, 450
    margin_l, margin_r = 80, 50
    margin_t, margin_b = 80, 60
    
    # Scale mappings
    # X: 1 to 15 -> margin_l to w - margin_r
    # Y: 0 to 14000 -> h - margin_b to margin_t
    def get_x(val):
        return margin_l + (val - 1) * (w - margin_l - margin_r) / 14.0
        
    def get_y(val):
        return (h - margin_b) - (val * (h - margin_t - margin_b) / 14000.0)

    # Generate grid lines
    grid_y = ""
    for y_val in range(0, 15000, 2000):
        y_pos = get_y(y_val)
        grid_y += f'<line x1="{margin_l}" y1="{y_pos}" x2="{w - margin_r}" y2="{y_pos}" stroke="#313244" stroke-dasharray="4" />'
        grid_y += f'<text x="{margin_l - 12}" y="{y_pos + 4}" fill="#a6adc8" font-size="11" text-anchor="end">{y_val} ns</text>'

    grid_x = ""
    for x_val in range(1, 16, 2):
        x_pos = get_x(x_val)
        grid_x += f'<line x1="{x_pos}" y1="{h - margin_b}" x2="{x_pos}" y2="{margin_t}" stroke="#313244" stroke-dasharray="4" />'
        grid_x += f'<text x="{x_pos}" y="{h - margin_b + 20}" fill="#a6adc8" font-size="11" text-anchor="middle">Ordre {x_val}</text>'

    # Generate line path and dots
    path_d = ""
    dots = ""
    for i, (x_val, y_val) in enumerate(points):
        px, py = get_x(x_val), get_y(y_val)
        if i == 0:
            path_d = f"M {px} {py}"
        else:
            path_d += f" L {px} {py}"
        dots += f'<circle cx="{px}" cy="{py}" r="5" fill="#a6e3a1" stroke="#1e1e2e" stroke-width="2" />'
        dots += f'<text x="{px}" y="{py - 12}" fill="#a6e3a1" font-size="10" font-weight="semibold" text-anchor="middle">{y_val} ns</text>'

    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="100%" height="100%" style="background-color: #1e1e2e; font-family: system-ui, -apple-system, sans-serif; border-radius: 12px; box-shadow: 0 8px 30px rgba(0,0,0,0.3);">
  <defs>
    <linearGradient id="bg-grad" x1="0%" y1="0%" x2="0%" y2="100%">
      <stop offset="0%" stop-color="#1e1e2e" />
      <stop offset="100%" stop-color="#11111b" />
    </linearGradient>
  </defs>

  <!-- Background -->
  <rect width="100%" height="100%" fill="url(#bg-grad)" rx="12" />

  <!-- Title -->
  <text x="40" y="45" fill="#cdd6f4" font-size="22" font-weight="bold" letter-spacing="0.5">Scalabilité de l'Équation Différentielle</text>
  <text x="40" y="65" fill="#a6adc8" font-size="13">Temps de génération de la matrice compagnon (SymAlgo++)</text>

  <!-- Grid -->
  {grid_y}
  {grid_x}

  <!-- Axis lines -->
  <line x1="{margin_l}" y1="{h - margin_b}" x2="{w - margin_r}" y2="{h - margin_b}" stroke="#45475a" stroke-width="2" />
  <line x1="{margin_l}" y1="{h - margin_b}" x2="{margin_l}" y2="{margin_t}" stroke="#45475a" stroke-width="2" />

  <!-- Line plot -->
  <path d="{path_d}" fill="none" stroke="#a6e3a1" stroke-width="3" stroke-linecap="round" stroke-linejoin="round" />

  <!-- Points and Values -->
  {dots}
</svg>"""
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(svg)

if __name__ == "__main__":
    os.makedirs("docs/images", exist_ok=True)
    create_comparison_svg("docs/images/bench_comparison.svg")
    create_scaling_svg("docs/images/bench_ode_scaling.svg")
    print("Graphiques SVG générés avec succès dans docs/images/")
