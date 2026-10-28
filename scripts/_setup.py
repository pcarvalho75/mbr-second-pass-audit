"""Shared start-up for the entry scripts: parse --config, point MBRAUDIT_CONFIG at it, make the package importable."""
import argparse
import json
import os
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def setup(description, add_arguments=None):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--config', type=Path, default=None,
                        help='configuration file (default: config.local.yaml, else config.template.yaml)')
    if add_arguments:
        add_arguments(parser)
    args = parser.parse_args()
    if args.config:
        os.environ['MBRAUDIT_CONFIG'] = str(args.config.resolve())
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    sys.dont_write_bytecode = True
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    return args


def note(path, payload):
    """Write a small JSON record of what a step read and when."""
    payload = {'recorded': time.strftime('%Y-%m-%d %H:%M:%S %z'), **payload}
    Path(path).write_text(json.dumps(payload, indent=2, default=str) + '\n', encoding='utf-8')
