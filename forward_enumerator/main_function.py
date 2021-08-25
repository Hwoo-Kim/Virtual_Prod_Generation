import argparse
import os
import sys
import yaml
from scripts.forwardEnumeration.forward_enumerator import forwardEnumerator
original_dir=os.getcwd()
virtual_lib_root= os.getcwd()+'/'
print(f'forward-enumeration root: {virtual_lib_root}')
project_root = os.path.dirname(os.getcwd())+'/'
os.environ['PYTHONPATH'] += f':{project_root}'

def forward_data_generation(args):
    with open(args.config,'r') as fr:
        config=yaml.safe_load(fr)
        Rset_1_path=config['Rset_1_path']
        Rset_2_path=config['Rset_2_path']
        r1_is_sm=bool(config['r1_is_sm'])
        r2_is_sm=bool(config['r2_is_sm'])
        numb_of_rand_pairs_per_core=config['numb_of_rand_pairs_per_core']
        rxn_temp_path=config['rxn_temp_path']
        dir_name=config['dir_name']
        numb_cores=config['numb_cores']

    # for each pair, product cannont exceed 5.
    # If you want to change, you should modify it in 'forwardEnumerator' code.
    # find 'list_of_products' in 'onestep_by_reactions' function, change it to be an integer you want.
    forwardEnumerator(Rset_1_path,
        r1_is_sm,
        Rset_2_path,
        r2_is_sm,
        numb_of_rand_pairs_per_core,
        rxn_temp_path,
        dir_name,
        numb_cores)

# main operation:
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", help='config_file', type=str)
    args = parser.parse_args()
    forward_data_generation(args)
