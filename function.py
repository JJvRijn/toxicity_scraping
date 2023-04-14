# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 21:00:57 2023

@author: JJvRijn
"""

def scraping_hem(search_querry, selection_vector=[1,0,0,0,0,1,0,0,1,0,0,0,0], remove_unnatural=True, **kwargs):
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
    
def Scraping_DRAMP (input_excel_location, remove_unnatural=True,):
    import pandas as pd
    from modlamp.descriptors import GlobalDescriptor
    import re
    #import the data
    DRAMP_data = pd.read_excel(input_excel_location)
    #extracting the usefull data
    seq = DRAMP_data.Sequence
    toxic = DRAMP_data.Hemolytic_activity
    df_DRAMP = pd.DataFrame([seq, toxic], index = ['seq', 'toxic'])
    df_DRAMP = df_DRAMP.transpose()
    df_DRAMP.seq = list(map(lambda x: x.upper(),df_DRAMP.seq))
    #dropping unknown activities
    df_DRAMP = df_DRAMP[df_DRAMP["toxic"] != "Not found"]
    df_DRAMP = df_DRAMP[df_DRAMP["toxic"] != " No hemolytic activity information found."]
    df_DRAMP = df_DRAMP[df_DRAMP["toxic"] != ' Not found in the literature']
    df_DRAMP = df_DRAMP.dropna(subset=['toxic'])
    df_DRAMP.toxic=df_DRAMP.toxic.str.replace(r'\[[^]]*\]', '')
    df_DRAMP = df_DRAMP.reset_index(drop=True)
    #extracting activities
    df_DRAMP['toxic2'] = 'NaN'
    df_DRAMP['toxic3'] = 'NaN'
    df_DRAMP['toxic4'] = 'NaN'
    for n in range(len(df_DRAMP)):
        if "-" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic']= re.sub('-[0-9. ]+[%%μµnm ]', '', df_DRAMP.toxic[n])
        if "±" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic']= re.sub('±[0-9. ]+[%%μµnm ]', '', df_DRAMP.toxic[n])    
        if "no hemolytic activity" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic'] = df_DRAMP.toxic[n].replace('no hemolytic activity', '0% hem')
        if "no detectable hemo" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic'] = df_DRAMP.toxic[n].replace('no detectable '
                                                                 , '0%')
        if "MHC" in df_DRAMP.toxic[n]:
            if "MHC10" not in df_DRAMP.toxic[n]:
                df_DRAMP.loc[n, 'toxic'] = df_DRAMP.toxic[n].replace('MHC', '15% hem')
        if "Non-hemoly" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic'] = df_DRAMP.toxic[n].replace('Non-', '0%')        
        df_DRAMP.toxic2[n] = re.findall(r"[-+]?\d*\.\d+|\d+", str(df_DRAMP.toxic[n]))
        if (len(df_DRAMP.toxic2[n])==2):
            df_DRAMP.loc[n, 'toxic3'] = (df_DRAMP.toxic2[n])[0]
            df_DRAMP.loc[n, 'toxic4'] = (df_DRAMP.toxic2[n])[1]
        if (len(df_DRAMP.toxic2[n])!= 0):
            if (len(df_DRAMP.toxic2[n])%2 == 0):
                df_DRAMP.loc[n, 'toxic3'] = (df_DRAMP.toxic2[n])[0]
                df_DRAMP.loc[n, 'toxic4'] = (df_DRAMP.toxic2[n])[1]            
                for i in range((int(len(df_DRAMP.toxic2[n])/2))):
                    df_DRAMP.loc[len(df_DRAMP.index)] = [df_DRAMP.seq[n], df_DRAMP.toxic[n],
                                                   df_DRAMP.toxic2[n], (df_DRAMP.toxic2[n])[0+2*i],
                                                   (df_DRAMP.toxic2[n])[1+2*i]]
    df_DRAMP = df_DRAMP[df_DRAMP.toxic2.str.len() != 0]
    df_DRAMP = df_DRAMP[df_DRAMP.toxic3 != "NaN"]
    df_DRAMP = df_DRAMP.reset_index(drop=True)
    #extract the unit of concentration
    df_DRAMP['toxic5'] = 'NaN'
    for n in range(len(df_DRAMP)):
        if "M" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "M"
        if "g/ml" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "g/ml"
        if "g/L" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "mg/ml" 
        if "mol/L" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "M"   
        if "μM" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "μM"
        if "µM" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "μM"
        if "µmol/L" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "μM"
        if "μmol/L" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "μM"        
        if "µg/ml" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "μg/ml"
        if "μg/ml" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "μg/ml"
        if "μg/mL" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "μg/ml"
        if "µg/mL" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "μg/ml"        
        if "nM" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "nM"       
        if "mg/L" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "mg/L"
        if "mg/ml" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "mg/ml"
        if "ng/ml" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "ng/ml" 
        if "nmol/ml" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "μM" 
        if "mM" in df_DRAMP.toxic[n]:
            df_DRAMP.loc[n, 'toxic5'] = "mM"              
    df_DRAMP= df_DRAMP[df_DRAMP.toxic5 != "NaN"]
    df_DRAMP = df_DRAMP.reset_index(drop=True)
    #remove the unnatural AA or return current info
    if remove_unnatural is True:
        NAA = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "U",
               "T", "W", "Y", "V", " "]
        unvalid =[]
        for j in df_DRAMP.index:
            if (not all ([NAA_seq in NAA for NAA_seq in df_DRAMP.seq[j]])):
                unvalid.append(j)
        df_DRAMP = df_DRAMP.drop(unvalid)
        df_DRAMP = df_DRAMP.reset_index(drop=True)
    else:
         df_final = df_DRAMP["seq", "toxic3", "toxic4", "toxic5"]
         df_final = df_final.rename( columns={"toxic3":"hem_activity", "toxic4":"concentration", "toxic5":"unit_concentration"})
         return df_final
    glob = GlobalDescriptor((df_DRAMP.seq.tolist()))
    glob.calculate_MW(amide=True)
    df_DRAMP['MW'] = glob.descriptor
    #see what concentrations are not in µM
    mask0 = df_DRAMP["toxic5"] == 'μg/ml'
    mask1 = df_DRAMP["toxic5"] == 'mg/L'
    mask2 = df_DRAMP["toxic5"] == 'mg/ml'
    mask3 = df_DRAMP["toxic5"] == 'μM'
    mask4 = df_DRAMP["toxic5"] == 'nM'
    mask5 = df_DRAMP["toxic5"] == 'mM'
    mask6 = df_DRAMP["toxic5"] == 'M'
    mask7 = df_DRAMP["toxic5"] == 'g/ml'
    mask8 = df_DRAMP["toxic5"] == 'ng/ml'
    #calculate the µM
    df_DRAMP["CONCENTRATION_µM"] = "NaN"
    df_DRAMP.loc[mask0, "CONCENTRATION_µM"] = (
        (df_DRAMP["toxic4"][mask0]).astype(float)*1000) / df_DRAMP[mask0]['MW']
    df_DRAMP.loc[mask1, "CONCENTRATION_µM"] = (
        (df_DRAMP["toxic4"][mask1]).astype(float)*1000) / df_DRAMP[mask1]['MW']
    df_DRAMP.loc[mask2, "CONCENTRATION_µM"] = (
        (df_DRAMP["toxic4"][mask2]).astype(float)*1000000) / df_DRAMP[mask2]['MW']
    df_DRAMP.loc[mask4, "CONCENTRATION_µM"] = (
        (df_DRAMP["toxic4"][mask4]).astype(float)/1000)
    df_DRAMP.loc[mask5, "CONCENTRATION_µM"] = (
        (df_DRAMP["toxic4"][mask5]).astype(float)*1000)
    df_DRAMP.loc[mask6, "CONCENTRATION_µM"] = (
        (df_DRAMP["toxic4"][mask6]).astype(float)*1000000)
    df_DRAMP.loc[mask7, "CONCENTRATION_µM"] = (
        (df_DRAMP["toxic4"][mask7]).astype(float)*1000000000) / df_DRAMP[mask7]['MW']
    df_DRAMP.loc[mask8, "CONCENTRATION_µM"] = (
        (df_DRAMP["toxic4"][mask8]).astype(float)) / df_DRAMP[mask8]['MW']
    #fill in values already in µM
    df_DRAMP.loc[mask3, "CONCENTRATION_µM"] = df_DRAMP["toxic4"][mask3].astype(float)
    df_final = df_DRAMP["seq", "toxic3", "CONCENTRATION_µM"]
    df_final = df_final.rename(columns={"toxic3":"hem_activity"})
    return df_final

def Scraping_DBAASP (input_excel_location,cell_types =["human"], remove_unnatural=True):
    import pandas as pd
    from modlamp.descriptors import GlobalDescriptor
    import regex as re
    #import the excel, this can be retrieved from https://dbaasp.org/search by pressing download the dataset
    DBAASP_data = pd.read_csv(input_excel_location)
    #extract the usefull columns
    DBAASP_usefull = DBAASP_data[["ID", "SEQUENCE", "HEMOLITIC CYTOTOXIC ACTIVITY - LYSIS VALUE",
                                  "HEMOLITIC CYTOTOXIC ACTIVITY - CONCENTRATION",
                                 "HEMOLITIC CYTOTOXIC ACTIVITY - UNIT", 'HEMOLITIC CYTOTOXIC ACTIVITY - TARGET CELL' 
                                  ]]
    #make the data workable
    #by removing unknown values and empty values
    DBAASP_usefull = DBAASP_usefull.dropna(thresh=3)
    DBAASP_usefull.reset_index(drop=True)
    #removing the unselected cell lines
    human = []
    for j in DBAASP_usefull.index:
        if cell_types in DBAASP_usefull[
        'HEMOLITIC CYTOTOXIC ACTIVITY - TARGET CELL'][j].lower():
            human.append(j)
    human = list( dict.fromkeys(human) )
    DBAASP_usefull = DBAASP_usefull.loc[human]
    DBAASP_usefull.reset_index(drop=True)
    #MEC isnt universally defined so is removed as correct input is near impossible
    #without processing the indvidual papers
    DBAASP_usefull = DBAASP_usefull[DBAASP_usefull[
        "HEMOLITIC CYTOTOXIC ACTIVITY - LYSIS VALUE"] != 'MEC']
    DBAASP_usefull = DBAASP_usefull[DBAASP_usefull[
        "HEMOLITIC CYTOTOXIC ACTIVITY - LYSIS VALUE"] != '-']
    DBAASP_usefull = DBAASP_usefull[DBAASP_usefull[
        "HEMOLITIC CYTOTOXIC ACTIVITY - CONCENTRATION"] != 'nan']
    DBAASP_usefull["HEMOLITIC CYTOTOXIC ACTIVITY - LYSIS VALUE"
                  ] = DBAASP_usefull["HEMOLITIC CYTOTOXIC ACTIVITY - LYSIS VALUE"].replace([
        "EC50", "IC50", 'LD50', 'CC50', 'LC50', 'MHC', 'LC90', 'IC90', 'ED50', 'LD90']
        , ["50%",'50%', '50%', '50%', '50%', '15%', '90%', '90%', '50%', '90%'])
    #removes rows with NA in specific cols
    DBAASP_usefull = DBAASP_usefull.dropna(subset=["HEMOLITIC CYTOTOXIC ACTIVITY - LYSIS VALUE"])
    DBAASP_usefull = DBAASP_usefull.dropna(subset=["HEMOLITIC CYTOTOXIC ACTIVITY - CONCENTRATION"])
    DBAASP_usefull = DBAASP_usefull.reset_index(drop=True)
    #Extracting int for Lysis_value
    mask3 = []
    for i in DBAASP_usefull.index:
        try:
            if "-" in str(DBAASP_usefull["HEMOLITIC CYTOTOXIC ACTIVITY - LYSIS VALUE"][i]):
                DBAASP_usefull.loc[i, "Lysis_value"] = re.findall(r"[-+]?\d*\.\d+|\d+",str(DBAASP_usefull[
                    "HEMOLITIC CYTOTOXIC ACTIVITY - LYSIS VALUE"][i]))[1]       
            else:
                DBAASP_usefull.loc[i, "Lysis_value"] = re.findall(r"[-+]?\d*\.\d+|\d+",str(DBAASP_usefull[
                    "HEMOLITIC CYTOTOXIC ACTIVITY - LYSIS VALUE"][i]))[0]
        except:
            mask3.append(i)
    DBAASP_usefull = DBAASP_usefull.drop(mask3)
    DBAASP_usefull = DBAASP_usefull.reset_index(drop=True)
    if remove_unnatural is True:
        #remove unnatural Amino acids
        NAA = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "U",
               "T", "W", "Y", "V", " "]
        unvalid =[]
        for j in DBAASP_usefull.index:
            if (not all ([NAA_seq in NAA for NAA_seq in DBAASP_usefull.SEQUENCE[j]])):
                unvalid.append(j)
        DBAASP_usefull = DBAASP_usefull.drop(unvalid)
        DBAASP_usefull = DBAASP_usefull.reset_index(drop=True)
        #add molweight to covert ug/ml to uM
        glob = GlobalDescriptor(DBAASP_usefull.SEQUENCE.tolist())
        glob.calculate_MW(amide=True)
        DBAASP_usefull['MW'] = glob.descriptor
        #see what concentrations are in µg/ml instead of µM
        mask = DBAASP_usefull["HEMOLITIC CYTOTOXIC ACTIVITY - UNIT"] == 'µg/ml'
        #extract numerical to remove >, <, <= etc.
        DBAASP_usefull["CONCENTRATION_µM"] = 0
        for o in range(len(DBAASP_usefull)):
            if "-" in str(DBAASP_usefull["HEMOLITIC CYTOTOXIC ACTIVITY - CONCENTRATION"][o]):
                 DBAASP_usefull.loc[o, "conc_num"] = re.findall(r"[-+]?\d*\.\d+|\d+",str(DBAASP_usefull[
                "HEMOLITIC CYTOTOXIC ACTIVITY - CONCENTRATION"][o]))[1]
            else:
                DBAASP_usefull.loc[o, "conc_num"] = re.findall(r"[-+]?\d*\.\d+|\d+",str(DBAASP_usefull[
                    "HEMOLITIC CYTOTOXIC ACTIVITY - CONCENTRATION"][o]))[0]
        #calculate the µM
        DBAASP_usefull["CONCENTRATION_µM"] = 0
        DBAASP_usefull.loc[mask, "CONCENTRATION_µM"] = ((DBAASP_usefull[
            "conc_num"][mask]).astype(float)*1000) / DBAASP_usefull[mask]['MW']
        #fill in values already in µM
        mask2 = ~mask
        DBAASP_usefull.loc[mask2, "CONCENTRATION_µM"] = DBAASP_usefull["conc_num"][mask2].astype(float)
        DBAASP_usefull = DBAASP_usefull.reset_index(drop=True)
        DBAASP_final = DBAASP_usefull["ID", "SEQUENCE", "Lysis_value", "CONCENTRATION_µM"]
        return DBAASP_final
    else:
        DBAASP_final = DBAASP_usefull["ID", "SEQUENCE", "Lysis_value", "conc_num", "HEMOLITIC CYTOTOXIC ACTIVITY - UNIT"]
        return DBAASP_final

def Scraping_DADP (input_excel_location, remove_unnatural=True):
    import pandas as pd
    import numpy as np
    import re
    from mechanize import Browser
    #extract the availible peptides
    browser = Browser()
    browser.open("http://split4.pmfst.hr/dadp/?a=list")
    response = browser.reload()
    content = response.read()
    names = re.findall(r'\bSP_[0-9, a-z, A-Z]*', str(content))
    #Scraping the individual peptides
    seq = []
    toxic = []
    toxic_type =[]
    for r in names:
        URL = "http://split4.pmfst.hr/dadp/?a=kartica&id=" + str(r)
        browser = Browser()
        browser.open(URL)
        response = browser.reload()
        content = response.read()
        table = content.split(b"</table>")[1]
        results = table.split(b"</td>")
        toxic_type.append(results[5].split(b"<td>",1)[1])
        seq.append(str(results[-5]).split("b'<td>",1)[1])
        if len(re.findall(r"[-+]?\d*\.\d+|\d+", str(results[-4]))) >= 2:
            toxic.append(re.findall(r"[-+]?\d*\.\d+|\d+", str(results[-4]))[1])
        else:
            toxic.append("Nan")
    #cleaning and formatting the data
    DADP_data = pd.DataFrame()
    DADP_data["ID"] = names
    DADP_data["seq"] = seq
    DADP_data['CONCENTRATION_µM'] = toxic
    DADP_data = DADP_data[DADP_data["CONCENTRATION_µM"] != "Nan"]
    DADP_data = DADP_data.reset_index(drop=True)
    DADP_data["seq"] = DADP_data["seq"].str[:-1]
    #because evrything is in HC50
    DADP_data["hem_activity"] = "HC50"
    return DADP_data
