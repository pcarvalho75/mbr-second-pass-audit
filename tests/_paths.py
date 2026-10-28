"""Make the package importable from the tests and keep them independent of any local configuration."""
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).resolve().parent / 'fixtures'
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
os.environ.setdefault('MBRAUDIT_CONFIG', str(REPO / 'config.template.yaml'))
sys.dont_write_bytecode = True
