import pandas as pd
import numpy as np
import argparse
import itertools
import os
import pickle
import yaml
from utils.EvidentialClassifier import EvidentialClassifier
from utils.utils import parse_material_structure_data, get_candidate_space
from lib import *

def load_data(cfg):
    df_observations = None
    data_path = get_uspex_dir(job=job, uspex_file="results1/Individuals_magmom_ene.csv")

    for micro_prompts in cfg["input_path"]["structures"]:
        df_ms_data = pd.read_csv(data_path, index_col=0)
        df_tmp = parse_material_structure_data(
            df_ms_data, micro_prompts["micro-prompt"], cfg["learn"]["list_of_elements"],
            cfg["input_path"]["fe_threshold"]
        )
        if df_observations is None:
            df_observations = df_tmp
        else:
            df_observations = pd.concat([df_tmp, df_observations])

    # for fix space reasoning
    # df_candidates = pd.read_csv(cfg["input_path"]["candidates"], index_col=0)

    df_candidates = get_candidate_space(
        material_composition=micro_prompts["micro-prompt"], list_of_elements=cfg["learn"]["list_of_elements"]
            )

    return df_observations, df_candidates 

def acquisition_function(df_candidates, acquisition_function="belief"):
    """Acquisition functions are mathematical techniques that guide how the material space should be explored during optimization"""
    if acquisition_function == "belief":
        return df_candidates["m_High"].values
    elif acquisition_function == "belief_utility":
        return df_candidates["m_High"].values + df_candidates["m_Unk"].values/2
    elif acquisition_function == "plausibility":
        return df_candidates["m_High"].values + df_candidates["m_Unk"].values
    elif acquisition_function == "epistemic":
        return df_candidates["m_Unk"].values
    elif acquisition_function == "aleatoric":
        return np.min(df_candidates[["m_High", "m_Low"]].values, axis=1)
    elif acquisition_function == "classification_uncertainty":
        return 1 - np.abs(df_candidates["m_High"].values - df_candidates["m_Low"].values)
    elif acquisition_function == "total":
        return 1 - np.max(df_candidates[["m_High", "m_Low"]].values, axis=1)
    elif acquisition_function == "plausibility_belief":
        return 2*df_candidates["m_High"].values + df_candidates["m_Unk"].values

def train_recommender(df_data_train, cfg_learn):

    ec = EvidentialClassifier(
        core_set=cfg_learn["list_of_elements"], frame_of_discernment=frozenset({"High", "Low"}), n_gram_evidence=2, 
        alpha=cfg_learn["alpha"], version=cfg_learn["ers_version"]
    )
    X = df_data_train["set_name"].values
    y = df_data_train[cfg_learn["target_variable"]].values
    ec.fit(X, y)
    return ec

def main(cfg):

    # Load dataset of materials structures
    df_observations, df_candidates = load_data(cfg)
    
    # Select top structures with lowest formation energy for each micro configuration to learn the recommender for micro configuration
    df_data_train = df_observations.sort_values(["energy_substance_pa"], ascending=True).groupby("set_name").head(1)

    # Train recommender
    recommender = train_recommender(df_data_train, cfg["learn"])
    
    # Evaluate candidates
    y_pred, final_decisions = recommender.predict(df_candidates["set_name"].values, show_decision=True)
    
    ## Parse results of candidates estimation to dataframe
    m_High, m_Low, m_Unk = [], [], []
    for final_decision in final_decisions: 
        m_High.append(final_decision[frozenset({"High"})])
        m_Low.append(final_decision[frozenset({"Low"})])
        m_Unk.append(final_decision[frozenset({"High", "Low"})])
    df_candidates["m_High"] = m_High
    df_candidates["m_Unk"] = m_Unk
    df_candidates["m_Low"] = m_Low
    df_candidates["y_pred"] = y_pred
    
    # Ranking candidates
    
    ## Measure acquisition score
    acquisition_scores = acquisition_function(df_candidates, cfg["acquisition_function"])
    df_candidates["acquisition_scores"] = acquisition_scores
    
    ## Sorting candidates
    df_candidates = df_candidates.sort_values(by=["acquisition_scores"], ascending=False)
    makedirs(cfg["output_path"] + "/tmp.txt")
    df_candidates.to_csv("{}/ranking.csv".format(cfg["output_path"]))

if __name__ == '__main__':
    
    # Example of setting parameter for micro-prompt recommender
    # parser = argparse.ArgumentParser(description='VASP-Viz: Main')
    # parser.add_argument('--config-file', type=str)
    
    # args = parser.parse_args()  
    for job in jobs:
        config_file = get_micro_config(job=job)

        with open(config_file, "r") as yaml_file:
            cfg = yaml.safe_load(yaml_file)
        
        output_path = get_uspex_dir(job=job, uspex_file="ML/micro-prompt")
        cfg["output_path"] = output_path

        main(cfg)