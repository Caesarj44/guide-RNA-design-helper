

import json
from pathlib import Path

CONFIGS_FILE = Path(__file__).parent.parent / 'configs.json'

DEFAULT_CONFIGS = {
    'main_configs': {
        'dark_theme': True,
        'doench_threshold': 50.0,
        'min_gc': 40,
        'max_gc': 70,
        'selected_tail' : 'SpCas9',
    },
    'preset_gRNA_tails':{
        'SpCas9':{
            'tail' : 'GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTT',
            'pam' : 'NGG',
            'spacer_len' : 20,
            'cut_distance' : '3'
        },
    },
    'custom_gRNA_tails':{}
}

#Jinek et al. (Science, 2013)
#GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTT


def load_configs(section: str | None = None, key : str | None = None):

    if not CONFIGS_FILE.exists():
        save_configs(DEFAULT_CONFIGS)
        data = DEFAULT_CONFIGS.copy()
    else:
        try:
            with open(CONFIGS_FILE, 'r', encoding='utf-8') as f:
                data = json.load(f)
        except (json.JSONDecodeError, OSError):
            data = DEFAULT_CONFIGS.copy()
    
    for sec, values in DEFAULT_CONFIGS.items(): # sec это section, секция данных конфига
        if sec not in data:
            data[sec] = values.copy()
        elif isinstance(values, dict):
            for k,v in values.items():
                if k not in data[sec]:
                    data[sec][k] = v

    if section is None:
        return data
    if key is None:
        return data.get(section, {})
    return data.get(section, {}).get(key)

def save_configs(data: dict):
    with open(CONFIGS_FILE, 'w', encoding='utf-8') as f:
        json.dump(data,f,indent=2, ensure_ascii=False)


def update_config (section: str, key: str, value: str):
    current = load_configs()
    if section not in current:
        current[section] = {}
    current[section][key] = value
    save_configs(current)

def load_style(filename: str) ->str:
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            return f.read()
    except(FileNotFoundError, OSError):
        return ''
    
def delete_config_key(section: str, key: str):
    current = load_configs()
    if section in current and key in current[section]:
        del current[section][key]
        save_configs(current)