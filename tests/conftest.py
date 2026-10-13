"""
Pytest configuration ensuring headless Agg backend for matplotlib across all test sessions.
"""

import matplotlib
matplotlib.use("Agg")
