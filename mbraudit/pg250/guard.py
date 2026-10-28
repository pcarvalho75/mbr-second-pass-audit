"""Input identity and output checks for the 250 pg steps (plumbing; replaces the old validation_v5/io.py).

guard_inputs hashes the scientific inputs and returns a sealed object; the fold-change functions refuse to run
without it (require_guard). For synthetic fixtures every sample name must start with SYN_.
"""
from dataclasses import dataclass
from pathlib import Path
import hashlib
import json

_SEAL = object()


@dataclass(frozen=True)
class Verified:
    producer: Path
    mode: str
    inputs: tuple
    sample_names: tuple
    _seal: object


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    text = value if isinstance(value, str) else json.dumps(value, indent=2, allow_nan=False) + '\n'
    path.write_bytes(text.replace('\r\n', '\n').replace('\r', '\n').encode('utf-8'))


def verify_inputs(inputs, expected=None):
    paths = {name: Path(path).resolve() for name, path in inputs.items() if path is not None}
    if len(set(paths.values())) != len(paths):
        raise ValueError('Scientific input roles must name distinct files')
    hashes = {}
    for name, path in paths.items():
        hashes[name] = sha(path)
        if expected and (name not in expected or hashes[name] != expected[name]):
            raise ValueError('Input identity differs: ' + name)
    if expected and set(expected) != set(paths):
        raise ValueError('Expected identities must cover every input role')
    return paths, hashes


def guard_inputs(inputs, mode, *, producer, sample_names=()):
    if mode not in {'real', 'synthetic'}:
        raise ValueError('Unknown input mode')
    paths, _ = verify_inputs(inputs)
    names = tuple(sample_names)
    if mode == 'synthetic' and (not names or any(not n.startswith('SYN_') for n in names)):
        raise ValueError('Synthetic inputs need explicit SYN_ sample names')
    return Verified(Path(producer).resolve(), mode, tuple(paths.items()), names, _SEAL)


def require_guard(result, *, producer, frames=()):
    if not isinstance(result, Verified) or result._seal is not _SEAL or result.producer != Path(producer).resolve():
        raise ValueError('Verified endpoint guard result required')
    if result.mode == 'synthetic':
        allowed = set(result.sample_names)
        for frame in frames:
            if 'sample' not in frame or any(not isinstance(n, str) or not n.startswith('SYN_') or n not in allowed for n in frame['sample']):
                raise ValueError('Synthetic table has an undeclared sample name')
    return result


def fresh_output(path, inputs):
    path = Path(path).resolve()
    if path.exists() and any(path.iterdir()):
        raise ValueError('Use a new empty output directory: ' + str(path))
    if any(path == p or path.is_relative_to(p.parent) for p in inputs.values()):
        raise ValueError('Output must be separate from scientific inputs')
    path.mkdir(parents=True, exist_ok=True)
    return path
