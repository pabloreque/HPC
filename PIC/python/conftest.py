"""Pytest setup: imports pic_prequeijo from the CMake build directory.

The build directory is ``$PIC_BUILD_DIR`` (set by ctest), or ``PIC/build``.
"""

import os
import sys
from pathlib import Path

import pytest

PIC_DIR = Path(__file__).resolve().parents[1]
BUILD_DIR = Path(os.environ.get("PIC_BUILD_DIR", PIC_DIR / "build"))
sys.path.insert(0, str(BUILD_DIR))


@pytest.fixture
def driver():
    """Path of the C++ driver executable built next to the module."""
    return BUILD_DIR / "pic_prequeijo"
