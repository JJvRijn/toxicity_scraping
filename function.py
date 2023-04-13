# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 21:00:57 2023

@author: JJvRijn
"""

def scraping_hem(search_querry, selection_vector=[1,0,0,0,0,1,0,0,1,0,0,0,0], remove_unnatural=False, **kwargs):
    #import required packages
    import pandas as pd
    from mechanize import Browser
    #standardize the input for my script
    if "https://webs.iiitd.edu.in/raghava/hemolytik/submitkey_target_search.php?ran=" not in search_querry:
        search_code = "https://webs.iiitd.edu.in/raghava/hemolytik/submitkey_target_search.php?ran=" + str(
            search_querry)
    else:
        search_code = search_querry
    #start the browser and grab the search results
    browser = Browser()
    response = browser.open(search_code)
    content = response.read()
    #start processing the results
    #seperating the data into the correct labels
    split = str(content, 'utf-8').split(sep='</td><td')
    clean_split = split[7:]
    PMID = []
    year = []
    name = []
    Nter_mod = []
    Cter_mod = []
    seq =[]
    length = []
    RBCs_source = []
    toxic = []
    mod = []
    stereo = []
    func = []
    DSSP = []
    for i in range(int((len(clean_split)/14))):
        PMID_pos = ((0) + 13*i)
        year_pos = ((1) + 13*i)
        name_pos = ((2) + 13*i)
        Nter_mod_pos = ((3)+ 13*i)
        Cter_mod_pos = ((4)+ 13*i)
        seq_pos = ((5) + 13*i)
        length_pos = ((6) + 13*i)
        RBCs_source_pos = ((7)+ 13*i)
        toxic_pos = ((8) + 13*i)
        mod_pos = ((9) + 13*i)
        stereo_pos = ((10)+ 13*i)
        func_pos = ((11) + 13*i)    
        DSSP_pos = ((12) + 13*i)
        
        PMID.append((clean_split[PMID_pos]))
        year.append(clean_split[year_pos])
        name.append(clean_split[name_pos])
        Nter_mod.append((clean_split[Nter_mod_pos]))
        Cter_mod.append(clean_split[Cter_mod_pos])
        seq.append((clean_split[seq_pos]).upper())
        length.append((clean_split[length_pos]))
        RBCs_source.append(clean_split[RBCs_source_pos])
        toxic.append(clean_split[toxic_pos])
        mod.append((clean_split[mod_pos]))
        stereo.append(clean_split[stereo_pos])
        func.append(clean_split[func_pos])
        DSSP.append(clean_split[DSSP_pos])
    #cleaning the individual pieces
    for n in range(len(seq)):
        PMID[n] = PMID[n].split(sep = "\"")[1]
        year[n] = year[n].split(sep=">")[1]
        name[n] = name[n].split(sep=">")[1]
        Nter_mod[n] = Nter_mod[n].split(sep=">")[1]
        Cter_mod[n] = Cter_mod[n].split(sep=">")[1]
        seq[n] = seq[n].split(sep=">")[1]
        length[n] = length[n].split(sep=">")[1]
        RBCs_source[n] = RBCs_source[n].split(sep=">")[1]
        toxic[n] = toxic[n].split(sep = "center>")[1]
        mod[n] = mod[n].split(sep=">")[1]
        stereo[n] = stereo[n].split(sep=">")[1]
        func[n] = func[n].split(sep=">")[1]
        DSSP[n] = DSSP[n].split(sep=">")[1]
    #selecting the data we want to have in our dataframe
    selectable = [PMID, year, name, Nter_mod, Cter_mod, seq, length, RBCs_source, toxic, mod, stereo, func, DSSP]
    boolvector = [bool(elem) for elem in selection_vector]
    selected_data = [val for selected, val in zip(boolvector, selectable) if selected]
    selectable_str = ["PMID", "year", "name", "Nter_mod", "Cter_mod", "seq", "length", "RBCs_source", "toxic", 
                      "mod", "stereo", "func", "DSSP"]
    selected_data_str = [val for selected, val in zip(boolvector, selectable_str) if selected]
    #creating the dataframe
    df_hemolytik = pd.DataFrame(selected_data, index = selected_data_str)
    df_hemolytik = df_hemolytik.transpose()
    if remove_unnatural is False:
        return df_hemolytik
    if remove_unnatural is True:
        #removing the unnatural amino acids
        NAA = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "U",
               "T", "W", "Y", "V", " "]
        unvalid =[]
        for j in df_hemolytik.index:
            if (not all ([NAA_seq in NAA for NAA_seq in df_hemolytik.seq[j]])):
                unvalid.append(j)
        df_hemolytik = df_hemolytik.drop(unvalid)
        df_hemolytik.reset_index(drop=True)
        return df_hemolytik
    
def string_to_act_conc (df_hemolytik):
    #input is a dataframe with a column that is called toxic containing the string that explains the toxicity    
    #automatic reading of the activity from string to floats
        import re
        import numpy as np
        from modlamp.descriptors import  GlobalDescriptor
        df_hemolytik['toxic2'] = np.nan
        df_hemolytik['toxic3'] = "other"
        df_hemolytik['toxic4'] = "other"
        df_hemolytik['toxic5'] = "other"
        for n in range(len(df_hemolytik)):
            if "-" in df_hemolytik.toxic[n]:
                df_hemolytik.loc[n, 'toxic']= re.sub('-[^>]+[%μµnm ]', '', df_hemolytik.toxic[n])
            if "±" in df_hemolytik.toxic[n]:
                df_hemolytik.loc[n, 'toxic']= re.sub('±[^>]+[%μµnm ]', '', df_hemolytik.toxic[n])
            if "MHC" in df_hemolytik.toxic[n]:
                df_hemolytik.loc[n, 'toxic'] = df_hemolytik.toxic[n].replace('MHC', '15% hem')
            if "No hemolysis" in df_hemolytik.toxic[n]:
                df_hemolytik.loc[n, 'toxic'] = df_hemolytik.toxic[n].replace('No hemolysis', '0% hem')
            if "Lethal" in df_hemolytik.toxic[n]:
                df_hemolytik.loc[n, 'toxic'] = df_hemolytik.toxic[n].replace('Lethal', '100% hem')
            if "not detected" in df_hemolytik.toxic[n]:
                df_hemolytik.loc[n, 'toxic'] = df_hemolytik.toxic[n].replace('not detected', '0% hem')
            if "M" in df_hemolytik.toxic[n]:
                df_hemolytik.loc[n, 'toxic3'] = "M"
            if "μM" in df_hemolytik.toxic[n]:
                df_hemolytik.loc[n, 'toxic3'] = "μM"
            if "µM" in df_hemolytik.toxic[n]:
                df_hemolytik.loc[n, 'toxic3'] = "μM"
            if "µg/ml" in df_hemolytik.toxic[n]:
                df_hemolytik.loc[n, 'toxic3'] = "μg/ml"
            if "μg/ml" in df_hemolytik.toxic[n]:
                df_hemolytik.loc[n, 'toxic3'] = "μg/ml"
            if "μg/mL" in df_hemolytik.toxic[n]:
                df_hemolytik.loc[n, 'toxic3'] = "μg/ml"
            if "µg/mL" in df_hemolytik.toxic[n]:
                df_hemolytik.loc[n, 'toxic3'] = "μg/ml"        
            if "nM" in df_hemolytik.toxic[n]:
                df_hemolytik.loc[n, 'toxic3'] = "nM"
            if "mg/L" in df_hemolytik.toxic[n]:
                df_hemolytik.loc[n, 'toxic3'] = "mg/L"
            if "mg/ml" in df_hemolytik.toxic[n]:
                df_hemolytik.loc[n, 'toxic3'] = "mg/ml"        
            if "mM" in df_hemolytik.toxic[n]:
                df_hemolytik.loc[n, 'toxic3'] = "mM"
            df_hemolytik.toxic2[n] = re.findall(r"[-+]?\d*\.\d+|\d+", str(df_hemolytik.toxic[n]))
            if (len(df_hemolytik.toxic2[n])==2):
                df_hemolytik.loc[n, 'toxic4'] = (df_hemolytik.toxic2[n])[0]
                df_hemolytik.loc[n, 'toxic5'] = (df_hemolytik.toxic2[n])[1]
        #removing peptides that could not have their activity added correctly
        unvalid1= []
        unvalid2= []
        for n in range(len(df_hemolytik)):
            if ((df_hemolytik.toxic3[n]== 'other')):
                unvalid1.append(n)
            if ((df_hemolytik.toxic4[n]== 'other')):
                unvalid2.append(n)
        unvalid = unvalid1 + list(set(unvalid2) - set(unvalid1))
        df_hemolytik = df_hemolytik.drop(unvalid)
        df_hemolytik.reset_index(drop=True)        
        #standardize the activity
        df_hemolytik = df_hemolytik.reset_index(drop=True)
        #types of units that need conversion "μg/ml", "nM", "M", "mM", "mg/L", "mg/ml"
        glob = GlobalDescriptor(df_hemolytik.seq.tolist())
        glob.calculate_MW(amide=True)
        df_hemolytik['MW'] = glob.descriptor
        #see what concentrations are in µg/ml instead of µM
        mask0 = df_hemolytik["toxic3"] == 'µg/ml'
        mask1 = df_hemolytik["toxic3"] == 'mg/L'
        mask2 = df_hemolytik["toxic3"] == 'mg/ml'
        mask3 = df_hemolytik["toxic3"] == 'μM'
        mask4 = df_hemolytik["toxic3"] == 'nM'
        mask5 = df_hemolytik["toxic3"] == 'mM'
        mask6 = df_hemolytik["toxic3"] == 'M'
        #calculate the µM
        df_hemolytik["CONCENTRATION_µM"] = 0
        df_hemolytik.loc[mask0, "CONCENTRATION_µM"] = (
            (df_hemolytik["toxic5"][mask0])*1000) / df_hemolytik[mask0]['MW']
        df_hemolytik.loc[mask1, "CONCENTRATION_µM"] = (
            (df_hemolytik["toxic5"][mask1]).astype(float)*1000) / df_hemolytik[mask1]['MW']
        df_hemolytik.loc[mask2, "CONCENTRATION_µM"] = (
            (df_hemolytik["toxic5"][mask2]).astype(float)*1000000) / df_hemolytik[mask2]['MW']
        df_hemolytik.loc[mask4, "CONCENTRATION_µM"] = (
            (df_hemolytik["toxic5"][mask4]).astype(float)/1000)
        df_hemolytik.loc[mask5, "CONCENTRATION_µM"] = (
            (df_hemolytik["toxic5"][mask5]).astype(float)*1000)
        df_hemolytik.loc[mask6, "CONCENTRATION_µM"] = (
            (df_hemolytik["toxic5"][mask6]).astype(float)*1000000)
        #fill in values already in µM
        df_hemolytik.loc[mask3, "CONCENTRATION_µM"] = df_hemolytik["toxic5"][mask3].astype(float)
        df_hemolytik.drop(columns = ["toxic2", "toxic3", "toxic5", "MW"])
        df_hemolytik.rename(columns={"toxic4": "hemolytic_activity"})
        return df_hemolytik