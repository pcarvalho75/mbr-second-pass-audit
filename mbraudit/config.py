"""Configuration: one YAML file of input paths, resolved relative to the folder that holds the file.

The file is chosen in this order: the MBRAUDIT_CONFIG environment variable, config.local.yaml at the repository
root (ignored by git), then config.template.yaml. Entry scripts set MBRAUDIT_CONFIG from their --config option
before importing the scientific modules, because the single-cell species registry reads its FASTA paths at import.

The single-cell code refers to the 64-cell search as 'all_cells_64' (its original alias). Which of the three search
configurations that alias means is set by MBRAUDIT_SINGLE_CELL (on, off or standard); each configuration is
analysed in its own process.
"""
import os
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
_CACHE = {}


def config_path():
    env = os.environ.get('MBRAUDIT_CONFIG')
    if env:
        return Path(env).resolve()
    local = REPO / 'config.local.yaml'
    return local if local.exists() else REPO / 'config.template.yaml'


def read():
    path = config_path()
    if path not in _CACHE:
        _CACHE[path] = yaml.safe_load(path.read_text(encoding='utf-8')) or {}
    return _CACHE[path]


def resolve(raw):
    path = Path(str(raw)).expanduser()
    return path if path.is_absolute() else (config_path().parent / path)


def value(key, default=KeyError):
    result = read()
    for part in key.split('.'):
        if not isinstance(result, dict) or part not in result:
            if default is KeyError:
                raise KeyError('Configuration key missing: ' + key)
            return default
        result = result[part]
    return result


def input_path(key):
    return resolve(value(key))


def optional_input(key):
    """String path of a configured input, or an empty string when the key is absent."""
    raw = value(key, None)
    return '' if raw is None else str(resolve(raw))


def output_root():
    path = resolve(value('output_dir', 'outputs'))
    path.mkdir(parents=True, exist_ok=True)
    return path


def output_dir(*parts):
    path = output_root().joinpath(*parts)
    path.mkdir(parents=True, exist_ok=True)
    return path


def active_configuration():
    name = os.environ.get('MBRAUDIT_SINGLE_CELL', 'on')
    if name not in ('on', 'off', 'standard'):
        raise ValueError('MBRAUDIT_SINGLE_CELL must be on, off or standard')
    return name


def single_cell_report_dir(alias):
    block = value('single_cell')
    if alias == 'all_cells_64':
        return resolve(block['configurations'][active_configuration()])
    return resolve(block['pools'][alias])
