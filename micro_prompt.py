import pandas as pd
import numpy as np
from utils.EvidentialClassifier import EvidentialClassifier

def load_data(cfg_input):
    df_observations = pd.read_csv(cfg_input["structures"], index_col=0)
    df_candidates = pd.read_csv(cfg_input["candidates"], index_col=0)
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
        core_set=cfg_learn["predicting_variables"], frame_of_discernment=frozenset({"High", "Low"}), n_gram_evidence=2, 
        alpha=cfg_learn["alpha"], version=cfg_learn["ers_version"]
    )
    X = df_data_train[cfg_learn["predicting_variables"]].values
    y = df_data_train[cfg_learn["target_variable"]].values
    ec.fit(X, y)
    return ec

def main(cfg):
    # Load dataset of materials structures
    df_observations, df_candidates = load_data(cfg["input_path"])
    
    # Select top 100th structures with lowest formation energy for each micro configuration to learn the recommender for micro configuration
    df_data_train = df_observations.sort_values(["energy_substance_pa"], ascending=True).groupby("set_name").head(100)
    
    # Train recommender
    recommender = train_recommender(df_data_train, cfg["learn"])
    
    # Evaluate candidates
    y_pred, final_decisions = recommender.predict(df_candidates[cfg["learn"]["predicting_variables"]].values, show_decision=True)
    
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
    df_candidates.to_csv("{}/ranking.csv".format(cfg["output_path"]))

if __name__ == '__main__':
    
    # Example of setting parameter for micro-prompt recommender
    cfg = {
        "input_path": {
            "structures": "example/Sm1Fe12/micro-prompt/structures.csv",
            "candidates": "example/Sm1Fe12/micro-prompt/candidates.csv"
        },
        "output_path": "example/Sm1Fe12/micro-prompt/output",
        "acquisition_function": "belief_utility",
        "learn": {
            "alpha": 0.001,
            "tuned_parameters": [
                {
                    "alpha": [0.001]
                }
            ],
            "features_type": "categorical",
            "predicting_variables": [
                'Sm', 'Fe', 'B', 'C', 'N'
                'SYMM_1', 'SYMM_2', 'SYMM_3', 'SYMM_4', 'SYMM_5', 'SYMM_6', 'SYMM_7', 'SYMM_8', 'SYMM_9', 'SYMM_10', 
                'SYMM_11', 'SYMM_12', 'SYMM_13', 'SYMM_14', 'SYMM_15', 'SYMM_16', 'SYMM_17', 'SYMM_18', 'SYMM_19', 'SYMM_20', 'SYMM_21', 
                'SYMM_22', 'SYMM_23', 'SYMM_24', 'SYMM_25', 'SYMM_26', 'SYMM_27', 'SYMM_28', 'SYMM_29', 'SYMM_30', 'SYMM_31', 'SYMM_32', 
                'SYMM_33', 'SYMM_34', 'SYMM_35', 'SYMM_36', 'SYMM_37', 'SYMM_38', 'SYMM_39', 'SYMM_40', 'SYMM_41', 'SYMM_42', 'SYMM_43', 
                'SYMM_44', 'SYMM_45', 'SYMM_46', 'SYMM_47', 'SYMM_48', 'SYMM_49', 'SYMM_50', 'SYMM_51', 'SYMM_52', 'SYMM_53', 'SYMM_54', 
                'SYMM_55', 'SYMM_56', 'SYMM_57', 'SYMM_58', 'SYMM_59', 'SYMM_60', 'SYMM_61', 'SYMM_62', 'SYMM_63', 'SYMM_64', 'SYMM_65', 
                'SYMM_66', 'SYMM_67', 'SYMM_68', 'SYMM_69', 'SYMM_70', 'SYMM_71', 'SYMM_72', 'SYMM_73', 'SYMM_74', 'SYMM_75', 'SYMM_76', 
                'SYMM_77', 'SYMM_78', 'SYMM_79', 'SYMM_80', 'SYMM_81', 'SYMM_82', 'SYMM_83', 'SYMM_84', 'SYMM_85', 'SYMM_86', 'SYMM_87', 
                'SYMM_88', 'SYMM_89', 'SYMM_90', 'SYMM_91', 'SYMM_92', 'SYMM_93', 'SYMM_94', 'SYMM_95', 'SYMM_96', 'SYMM_97', 'SYMM_98', 
                'SYMM_99', 'SYMM_100', 'SYMM_101', 'SYMM_102', 'SYMM_103', 'SYMM_104', 'SYMM_105', 'SYMM_106', 'SYMM_107', 'SYMM_108', 'SYMM_109', 
                'SYMM_110', 'SYMM_111', 'SYMM_112', 'SYMM_113', 'SYMM_114', 'SYMM_115', 'SYMM_116', 'SYMM_117', 'SYMM_118', 'SYMM_119', 'SYMM_120', 
                'SYMM_121', 'SYMM_122', 'SYMM_123', 'SYMM_124', 'SYMM_125', 'SYMM_126', 'SYMM_127', 'SYMM_128', 'SYMM_129', 'SYMM_130', 'SYMM_131', 
                'SYMM_132', 'SYMM_133', 'SYMM_134', 'SYMM_135', 'SYMM_136', 'SYMM_137', 'SYMM_138', 'SYMM_139', 'SYMM_140', 'SYMM_141', 'SYMM_142', 
                'SYMM_143', 'SYMM_144', 'SYMM_145', 'SYMM_146', 'SYMM_147', 'SYMM_148', 'SYMM_149', 'SYMM_150', 'SYMM_151', 'SYMM_152', 'SYMM_153', 
                'SYMM_154', 'SYMM_155', 'SYMM_156', 'SYMM_157', 'SYMM_158', 'SYMM_159', 'SYMM_160', 'SYMM_161', 'SYMM_162', 'SYMM_163', 'SYMM_164', 
                'SYMM_165', 'SYMM_166', 'SYMM_167', 'SYMM_168', 'SYMM_169', 'SYMM_170', 'SYMM_171', 'SYMM_172', 'SYMM_173', 'SYMM_174', 'SYMM_175', 
                'SYMM_176', 'SYMM_177', 'SYMM_178', 'SYMM_179', 'SYMM_180', 'SYMM_181', 'SYMM_182', 'SYMM_183', 'SYMM_184', 'SYMM_185', 'SYMM_186', 
                'SYMM_187', 'SYMM_188', 'SYMM_189', 'SYMM_190', 'SYMM_191', 'SYMM_192', 'SYMM_193', 'SYMM_194', 'SYMM_195', 'SYMM_196', 'SYMM_197', 
                'SYMM_198', 'SYMM_199', 'SYMM_200', 'SYMM_201', 'SYMM_202', 'SYMM_203', 'SYMM_204', 'SYMM_205', 'SYMM_206', 'SYMM_207', 'SYMM_208', 
                'SYMM_209', 'SYMM_210', 'SYMM_211', 'SYMM_212', 'SYMM_213', 'SYMM_214', 'SYMM_215', 'SYMM_216', 'SYMM_217', 'SYMM_218', 'SYMM_219', 
                'SYMM_220', 'SYMM_221', 'SYMM_222', 'SYMM_223', 'SYMM_224', 'SYMM_225', 'SYMM_226', 'SYMM_227', 'SYMM_228', 'SYMM_229', 'SYMM_230'
            ],
            "target_variable": "Label"
        }
    }
    
    main(cfg)