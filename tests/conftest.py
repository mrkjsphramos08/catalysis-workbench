"""
Pytest configuration ensuring headless Agg backend for matplotlib across all test sessions.
"""

import matplotlib
matplotlib.use("Agg")


import matplotlib.pyplot as plt
import pytest


@pytest.fixture(autouse=True)
def close_test_figures():
    """Release figures created by each test without touching pre-existing figures."""
    before = set(plt.get_fignums())
    yield
    for number in set(plt.get_fignums()) - before:
        plt.close(number)
