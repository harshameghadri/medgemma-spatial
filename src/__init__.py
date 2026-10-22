"""MedGemma spatial transcriptomics pipeline."""

try:
    import numba

    # numba's OpenMP layer and torch's bundled libomp segfault together on
    # macOS/arm64 (scanpy neighbors after `import torch`). workqueue avoids it.
    numba.config.THREADING_LAYER = "workqueue"
except ImportError:
    pass
