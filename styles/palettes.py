"""
Colorblind-accessible scientific color palettes and utilities.

Standardized color schemes derived from:
- Wong, B. (2011) Nature Methods 8, 441 (Okabe-Ito palette)
- Tol, P. (2021) Colour Schemes (SRON Technical Note)
- ColorBrewer 2.0 (Brewer, C.A.)
"""

from typing import Dict, List, Optional
import matplotlib.colors as mcolors


OKABE_ITO: List[str] = [
    "#0072B2",  # Blue
    "#D55E00",  # Vermilion
    "#009E73",  # Bluish green
    "#CC79A7",  # Reddish purple
    "#E69F00",  # Orange
    "#56B4E9",  # Sky blue
    "#F0E442",  # Yellow
    "#000000",  # Black
]

TOL_BRIGHT: List[str] = [
    "#4477AA",  # Blue
    "#EE6677",  # Red
    "#228833",  # Green
    "#CCBB44",  # Yellow
    "#66CCEE",  # Cyan
    "#AA3377",  # Purple
    "#BBBBBB",  # Grey
]

TOL_MUTED: List[str] = [
    "#332288",  # Indigo
    "#88CCEE",  # Cyan
    "#44AA99",  # Teal
    "#117733",  # Green
    "#999933",  # Olive
    "#DDCC77",  # Sand
    "#CC6677",  # Rose
    "#882255",  # Wine
    "#AA4499",  # Purple
    "#DDDDDD",  # Pale grey
]

NATURE_REVIEWS: List[str] = [
    "#203864",  # Deep Navy
    "#C00000",  # Crimson
    "#385723",  # Forest Green
    "#7030A0",  # Royal Violet
    "#ED7D31",  # Burnt Orange
    "#41719C",  # Steel Blue
    "#7F7F7F",  # Charcoal
]

# Orbital specific colors commonly used in DFT literature (s, p, d, f)
ORBITAL_COLORS: Dict[str, str] = {
    "s": "#0072B2",        # Blue
    "p": "#009E73",        # Green
    "d": "#D55E00",        # Vermilion / Red-Orange
    "f": "#CC79A7",        # Purple
    "t2g": "#E69F00",      # Orange
    "eg": "#56B4E9",       # Sky Blue
    "total": "#222222",    # Dark Grey / Black
    "up": "#0072B2",       # Spin Up
    "down": "#D55E00",     # Spin Down
}

PALETTES: Dict[str, List[str]] = {
    "okabe_ito": OKABE_ITO,
    "tol_bright": TOL_BRIGHT,
    "tol_muted": TOL_MUTED,
    "nature": NATURE_REVIEWS,
}


def get_palette(name: str = "okabe_ito", n: Optional[int] = None) -> List[str]:
    """
    Retrieve a color palette by name, optionally truncated or cycled to length n.

    Parameters
    ----------
    name : str
        Palette identifier ('okabe_ito', 'tol_bright', 'tol_muted', 'nature').
    n : int, optional
        Number of colors requested. If greater than palette length, colors cycle.

    Returns
    -------
    List[str]
        List of hex color codes.
    """
    key = name.lower()
    if key not in PALETTES:
        raise ValueError(f"Unknown palette '{name}'. Available: {list(PALETTES.keys())}")
    
    colors = PALETTES[key]
    if n is None:
        return list(colors)
    if n <= len(colors):
        return colors[:n]
    # Cycle if requesting more than length
    return [colors[i % len(colors)] for i in range(n)]


def get_color(index: int, palette: str = "okabe_ito") -> str:
    """Get color at specific zero-based index with cyclic wraparound."""
    colors = get_palette(palette)
    return colors[index % len(colors)]


def get_orbital_color(orbital: str) -> str:
    """Return standard publication color for atomic orbital projection."""
    orb = orbital.strip().lower()
    return ORBITAL_COLORS.get(orb, "#555555")
