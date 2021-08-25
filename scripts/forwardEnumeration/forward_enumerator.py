from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MolFromSmiles as Mol
from rdkit.Chem import MolToSmiles as Smiles
from rdkit.Chem.AllChem import ReactionFromSmarts as Rxn
#from rdkit.Chem.FragmentMatcher import FragmentMatcher
#from scripts.retrosynAnalysis.utils import reactant_frags_generator
#from scripts.retrosynAnalysis.utils import working_dir_setting
#from scripts.retrosynAnalysis.utils import timeout, TimeoutError
from multiprocessing import Lock, Process, Queue, current_process 
from datetime import datetime
import queue # imported for using queue.Empty exception
import pickle
import shutil
import json
import sys, time
import os
import random
from copy import copy, deepcopy
#from bisect import bisect_left
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

class SynthesisTree:
    '''
    Tree structure to represent a synthesis tree to the given target smiles.
    '''
    def __init__(self, target_smi):
        self.target = target_smi
        self.tree = [[target_smi]]
        self.notFinished = [target_smi]
        self.lastRxnIdx = None
        self.lastPrecursors = None
        self.lastRxnLoc = None

    def getCopiedTree(self):
        #new_tree = SynthesisTree(self.target)
        #new_tree.tree = [copy(pair) for pair in self.tree]
        #new_tree.notFinished = copy(self.notFinished)
        #new_tree.lastRxnIdx = copy(self.lastRxnIdx)
        #new_tree.lastPrecursors = copy(self.lastPrecursors)
        #new_tree.lastRxnLoc = copy(self.lastRxnLoc)
        #return new_tree
        return deepcopy(self)
        #return copy.copy(self)

    def getNumNotFin(self):
        return len(self.notFinished)

    def getNumRxn(self):
        return len(self.tree) -1

    def getTarget(self):
        return self.target

    def getTree(self):
        return self.tree

    def getNotFin(self):
        return self.notFinished

    def getLastRxnInform(self):
        return self.lastRxnIdx, self.lastPrecursors

    def getLastRxnLoc(self):
        return self.lastRxnLoc

    def getNumOthers(self):
        return self.numOthers

    def setLastRxnInform(self, idx:int, result:list, rxn_position:int):
        self.lastRxnIdx = idx
        self.lastPrecursors = result
        self.lastRxnLoc = rxn_position

    def removeNotFinLoc(self, loc_removed:int):
        del self.notFinished[loc_removed]

    def removeNotFinElem(self, elem_removed):
        self.notFinished.remove(elem_removed)

    def insertList(self, loc, l):
        #self.tree.append(copy.deepcopy(self.tree[-1]))
        self.tree.append(copy(self.tree[-1]))
        last = self.tree[-1]
        if len(last) > 1:
            del last[-1]
        del last[loc]
        for idx, elem in enumerate(l):
            last.insert(loc+idx, elem)

    def insertToNotFinished(self, loc, l):
        del self.notFinished[loc]
        self.numOthers = len(self.notFinished)
        for idx, elem in enumerate(l):
            self.notFinished.insert(loc+idx, elem)

    def getExpandedTrees(self, rxn_position, rxn_results):
        '''
        rxn_position is an index for notFinished list!!
        '''
        expanded_trees =[]
        elem = self.notFinished[rxn_position]
        loc = self.tree[-1].index(elem)
        for rxn_idx, result in rxn_results:
            copied_tree = self.getCopiedTree()
            copied_tree.insertList(loc, result)
            copied_tree.tree[-1].append([loc, rxn_idx])
            copied_tree.insertToNotFinished(rxn_position, result)
            copied_tree.setLastRxnInform(rxn_idx, result, rxn_position)
            expanded_trees.append(copied_tree)
        return expanded_trees

    def updateRbagInform(self, classified_reactant_bag:list, diff:int):
        '''
        This method updates the R_bag inform of the synthesis_tree, and returns boolean for the test to yield True or False.
        If the test is not determined in the current level, this methods returns None obj.
        If the result for the test is False, the tree might be not updated properly. (Bcs not necessary)
        '''
        exit = False
        precursor_Rbag_check = []
        num_others = self.getNumOthers()
        limit_numb = diff - num_others
        if limit_numb < 0:
            return False
        rxn_idx, last_precursors = self.getLastRxnInform()
        target_Rbag = classified_reactant_bag[rxn_idx]
        for prec in last_precursors:
            check = prec in target_Rbag
            precursor_Rbag_check.append(check)
            if precursor_Rbag_check.count(False) > limit_numb:
                exit = True
                break
            if check: 
                self.removeNotFinElem(prec)

        if exit:
            return False
        elif precursor_Rbag_check.count(True) == len(last_precursors) and num_others == 0:
            return True
        else:
            return None

# Forward enumeration
def duplicate_remove(list_of_list):
    if list_of_list == []:
        return []
    result = []
    for s in list_of_list:
        if not s[1] in result and not s[1][::-1] in result:
            result.append(s)
    return result

"""
def rxn_result_dup_remove(list_of_smi):
    '''
    list_of_list: list of sorted lists
    '''
    if len(list_of_smi)==0:
        return None
    ret_list = [list_of_smi[0]]
    for l in list_of_smi[1:]:
        if not l in list_of_smi and not None in l:
            ret_list.append(l)
    return ret_list
    """

    
def onestep_by_reactions(
        reactant_pair_in_mol: list,
        rxn_objs: list,
        first: bool,
        only_regio_not_problematic: bool,
        r1_is_sm: bool,
        r2_is_sm: bool
        ):
    '''
    Args:
      reactant_pair_in_mol: A pair of mol objects of molecules to which the reaction will be applied.
      rxn_objs: List of reaction objects.
      only_regio_not_problematic: Whether considers reactions only not having regio-selectivity or not
      only_uni: Whether applies only unitary reactions or not.
      only_bin: Whether applies only binary reactions or not.
    Returns:
      list_of_products: List of generated products from both uni and bi.
    '''
    list_of_products = []
    indice_list = list(range(len(rxn_objs)))
    max = 1
    random.shuffle(indice_list)
    for rxn_idx in indice_list:
        if len(list_of_products) >= max:
            break
        rxn = rxn_objs[rxn_idx]
        if  rxn.GetNumReactantTemplates() == 1:
        # 1. unitary reactions
            for idx, mol in enumerate(reactant_pair_in_mol):
                if idx == 0:
                    if not first and r1_is_sm:
                        continue
                if idx == 1:
                    if not first and r2_is_sm:
                        continue
                try:
                    rxn_results = rxn.RunReactants([mol])
                except:
                    #print(Smiles(mol))
                    continue
                rxn_results = list(set([Smiles(prod) for pair in rxn_results for prod in pair]))
                #rxn_results = rxn_result_dup_remove(rxn_results)
                if rxn_results == None:
                    continue
                if only_regio_not_problematic:
                    if len(rxn_results) > 1:        # about regio-selectivity 
                        continue
                for p_mol in rxn_results:
                    # precursor, rxn idx, product
                    #to_add = [[[Smiles(mol)], rxn_idx], p_mol]

                    # product, rxn idx, precursor_idx
                    to_add = [p_mol, rxn_idx, 'uni', idx]
                    list_of_products.append(to_add)

        elif rxn.GetNumReactantTemplates() == 2:
        # 2. binary reactions
            try:
                rxn_results = rxn.RunReactants(reactant_pair_in_mol)
            except:
                #print(Smiles(mol))
                continue
            rxn_results = list(set([Smiles(prod) for pair in rxn_results for prod in pair]))
            #rxn_results = rxn_result_dup_remove(rxn_results)
            if rxn_results == None:
                continue
            if only_regio_not_problematic:
                if len(rxn_results) > 1:        # about regio-selectivity 
                    continue
            for p_mol in rxn_results:
                # precursor, rxn idx, product
                #to_add = [[[Smiles(reactant_pair_in_mol[0]), Smiles(reactant_pair_in_mol[1])],\
                #        rxn_idx], p_mol]

                # product, rxn idx, None
                to_add = [p_mol, rxn_idx, 'bin', None]
                list_of_products.append(to_add)
        else:
            continue

    return list_of_products

def first_reaction(random_pairs, rxn_objs, only_regio_not_problematic):
    syn_trees = []
    for pair in random_pairs:
        if pair[0][0] == pair[1][0]:
            continue
        pair_in_mol = [Mol(pair[0][0]), Mol(pair[1][0])]
        reaction_result = onestep_by_reactions(pair_in_mol, rxn_objs, first=True,
            only_regio_not_problematic=only_regio_not_problematic, r1_is_sm=True, r2_is_sm=True)

        if reaction_result == []:
            continue
        pair = [pair[0][0], pair[1][0]]
        for result in reaction_result:
            if result[2] == 'uni':
                result[3] = [pair[result[3]]]
            elif result[2] == 'bin':
                result[3] = pair
            syn_trees.append(result)
    return syn_trees

def further_reaction(random_pairs, rxn_objs, only_regio_not_problematic, r1_is_sm, r2_is_sm):
    syn_trees = []
    for pair in random_pairs:
        if pair[0][0] == pair[1][0]:
            continue
        pair_in_mol = [Mol(pair[0][0]), Mol(pair[1][0])]
        reaction_result = onestep_by_reactions(pair_in_mol, rxn_objs, first=False,
            only_regio_not_problematic=only_regio_not_problematic, r1_is_sm=r1_is_sm, r2_is_sm=r2_is_sm)

        if reaction_result == []:
            continue
        if r1_is_sm == True:
            pair[0] = pair[0][0]
        if r2_is_sm == True:
            pair[1] = pair[1][0]
        for result in reaction_result:
            if result[2] == 'uni':
            # in further reaction, there is no case that uni reaction is applied to starting material.
                result[3] = pair[result[3]]
            elif result[2] == 'bin':
                result[3] = pair
            syn_trees.append(result)
    return syn_trees

# FROM HERE (21.08.12)
# first_reaction and further_reaction also.
# see memo for data structure.
def forward_enumerator(random_pairs:list, rxn_templates:list, only_regio_not_problematic,
        r1_is_sm:bool, r2_is_sm:bool):
    rxn_objects = []
    for rxn in rxn_templates:
        try:
            rxn_objects.append(Rxn(rxn))
        except:
            rxn_objects.append(None)

    if r1_is_sm and r2_is_sm:
        return first_reaction(random_pairs, rxn_objects,
                only_regio_not_problematic=only_regio_not_problematic)
    
    else:
        return further_reaction(random_pairs, rxn_objects,
                only_regio_not_problematic=only_regio_not_problematic,\
                r1_is_sm=r1_is_sm, r2_is_sm=r2_is_sm)

def do_forward_enumeration(tasks, r1_is_sm:bool, r2_is_sm:bool):
    while True:
        try:
            args = tasks.get(timeout=1)
        except queue.Empty:
            break
        else:
            random_pairs, rxn_templates, task_idx = args
            since=time.time()
            print(f'task_idx: {task_idx}')
            syn_trees = forward_enumerator(random_pairs,rxn_templates, only_regio_not_problematic=True,\
                    r1_is_sm=r1_is_sm, r2_is_sm=r2_is_sm)
            with open(f'syn_trees_{task_idx}.json', 'w') as fw:
                json.dump(syn_trees, fw)
            print(f'  {task_idx}th task time: {(time.time()-since):.2f}')
    return True

def forwardEnumerator(
        Rset_1_path:str,
        r1_is_sm:bool,
        Rset_2_path:str,
        r2_is_sm:bool,
        numb_of_rand_pairs_per_core:int,
        rxn_temp_path:list,
        dir_name,
        numb_cores):
    '''
    Main function. This conducts multiprocessing of 'retrosynthetic_analysis_single_batch'.
    '''
    with open(Rset_1_path, 'r') as fr:
        if r1_is_sm == True:
            Rset_1 = fr.read().splitlines()
        else:
            Rset_1 = json.load(fr)
            #Rset_1 = json.load(fr)[40000:80000]
    with open(Rset_2_path, 'r') as fr:
        if r2_is_sm == True:
            Rset_2 = fr.read().splitlines()
        else:
            Rset_2 = json.load(fr)
            #Rset_2 = json.load(fr)[240000:320000]

    with open(rxn_temp_path, 'r') as fr:
        rxn_templates = json.load(fr)

    working_dir = f'{os.getcwd()}/{dir_name}'
    try:
        os.mkdir(working_dir)
    except:
        i=2
        while True:
            try:
                os.mkdir(f'{working_dir}{i}')
            except:
                i+=1
            else:
                break
        working_dir = f'{working_dir}{i}'
    os.chdir(working_dir)
    print(f'Current working directory is:\n  {working_dir}')

    Rset_1_indice = [random.randint(0,len(Rset_1)-1) for i in range(numb_of_rand_pairs_per_core*numb_cores)]
    Rset_2_indice = [random.randint(0,len(Rset_2)-1) for i in range(numb_of_rand_pairs_per_core*numb_cores)]
    random_pairs = []
    for idx, num in enumerate(Rset_1_indice):
        to_append = []
        if r1_is_sm: to_append.append([Rset_1[num]])
        else: to_append.append(Rset_1[num])
        if r2_is_sm: to_append.append([Rset_2[Rset_2_indice[idx]]])
        else: to_append.append(Rset_2[Rset_2_indice[idx]])
        random_pairs.append(to_append)

    # multiprocessing of do_forward_enumeration
    numb_of_tasks = int(numb_cores)
    numb_of_procs = int(numb_cores)

    tasks = Queue()
    procs = []

    since = time.time()
    # creating tasks
    for task_idx in range(numb_of_tasks):
        args = (random_pairs[task_idx*numb_of_rand_pairs_per_core:(task_idx+1)*numb_of_rand_pairs_per_core],\
                rxn_templates, task_idx)
        tasks.put(args)

    # creating processes
    for worker in range(numb_of_procs):
        p = Process(target = do_forward_enumeration, args = (tasks, r1_is_sm, r2_is_sm))
        procs.append(p)
        p.start()
        time.sleep(0.5)

    # completing processes
    for p in procs:
        p.join()

    # join the results
    print('-----'*4)
    print('Forward enumeration step finished.\n  Joining the results...')

    syn_trees_list = []
    for task_idx in range(numb_of_tasks):
        with open(f'syn_trees_{task_idx}.json', 'r') as fr:
            syn_trees_list += json.load(fr)

        #os.remove(f'positive_set_depth_{current_depth}_{task_idx}.smi')

    # save the result
    time_passed = int(time.time()-since)
    now = datetime.now()
    finished_at = now.strftime('%Y. %m. %d (%a) %H:%M:%S')
    
    #result_report = [f'----- Config information -----\n',
    #            f'  Retro analysis started at: {since_inform}\n', \
    #            f'  Target data path: {retrosynthetic_analysis_config["retro_analysis_target"]}\n',\
    #            f'  Target data name: {target_data_name}\n', \
    #            f'  Uni template data: {uni_temp_data}\n', f'  Bi template data: {bi_temp_data}\n', \
    #            f'  Start index: {start_index}\n', f'  Reactant bag path: {conversed_classified_reactant_data}\n', f'  Depth: {depth}\n', \
    #            f'  Number of target molecules: {numb_molecules}\n', f'  Number of cores: {numb_cores}\n', \
    #            f'  Batch size: {batch_size}\n', f'  With path search: {with_path_search}\n', \
    #            f'  Max time: {max_time}\n\n', '----- Generation result -----\n']
    #result_report += [f'  Positive set depth_{i+1}: {numb_of_mols_in_each_pos[i]}\n' for i in range(depth)]
    #result_report += [f'  Negative set depth_{depth}: {numb_of_mols_in_neg}\n',\
    #        f'\n  finished_at: {finished_at}', \
    #        '\n   time passed: [%dh:%dm:%ds]' %(time_passed//3600, (time_passed%3600)//60, time_passed%60)]
    with open('generation_result.json', 'w') as fw:
        json.dump(syn_trees_list, fw)
    print('-----'*4)
    print(f'Forward enumeration finished at:\n  {finished_at}')
    print('  time passed: [%dh:%dm:%ds]' %(time_passed//3600, (time_passed%3600)//60, time_passed%60))
    print('-----'*4)
    return True


if __name__=='__main__':
    rxn = Rxn('CC>>OO')
    rxn2 = Rxn('CC>>NN')
    mol = Mol('CC')
    mol2 = Mol('CCC')
    print(rxn.RunReactants([mol,mol2]))
