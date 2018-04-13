from __future__ import absolute_import

from .FastaIO import *
from .base    import *

try:
    from calign import aligner, score_alignment
except ImportError:
    from . import matrix
    from .align import AlignmentResult, aligner