"""
This file contains all the functions that apply calculations or conduct 
estimates based on an assumptions or specific value. 
It currently contains:
    
    General Scripts: 
        units
        filter_distance
        calculateLoads (! incomplete)
        
    Specific Scripts: 
        estimate_flows_wwtp
        concentration_wwtp
        calculate_loads_wwtp
"""
import pandas as pd
import numpy as np
import geopy.distance
# Modules
import dictionaries

# VARIABLES #

# MWWTP:
sewEst = 300 # L/day from "...\MAIN\MWWTP\Sources\WWTPSources_WORKINGCOPY_rev5.xlsx" tab: "SewageEst"
avgPerson = 2.4 # Average Number of people / house according to national statistics
flow_estimate_average_file = "MWWTP/Sources/CampRVResortOther_WWTPsFlows_July122019.xlsx"
# Concentration files
wwtp_concentration_info_file = "MWWTP\\Sources\\ConcentrationData.xlsx"

# Conversion Factors 
# ex. ng_mg: ng to mg conversion factor. 
ng_mg = 1e-6 # mg/ng
pg_mg = 1e-9 # mg/pg                   
ug_mg = 1e-3 # mg/ug
ugg_mgkg = 1 # [ug/g] * [g/1e6ug] * [1e6mg/kg] = 1 [mg/kg] / [ug/g]
mgkg_mgmL = 1 # Density conversion factor for pure water 1kg/L
mgkg_mgmL_saltWater = 1.029 # Average density of surface salt water in kg/L
ppt_pg = 1 # 1 ppt = 1pg/g
# Radium 226 only, mass 1 Curie (Ci) = 1.01g, 1 Bq = 2.703e-11 Ci
# 2.703e-11[Ci/L] * 1.01 g/Ci = 2.73e-11g/L = 2.73e-8mg/L
Bq_mg = 2.73e-8 # [mg/Bq]
d_yr = 365 # days/year
s_yr = 3600 * 24 * 365 # [s/hr] * [hr/day] * [day/yr]
L_m3 = 1e-3 # 1m3/1e3L
m3_L = 1e3 # 1e3L/m3
mg_kg = 1e-6 # 1e-6kg/mg

###############################################################################
#
def units(df, unitColumn="Units", ValueColumn="Value", convertedValueName="Converted_Value", convertedUnitName="Converted_Units", replace=True, EMS=False, p=True):
    if p: print("\t\tIn Calculations: units...")
    # make sure float/numeric
    df[ValueColumn] = df[ValueColumn].astype(float)
    # Make converted columns 
    df[convertedValueName] = np.NaN
    df[convertedUnitName] = np.NaN
        
    # Basic Unit Conversions
    condition = df[unitColumn] == "ng/L"
    df.loc[condition, convertedValueName] = df.loc[condition, ValueColumn] * ng_mg 
    df.loc[condition, convertedUnitName] = "mg/L"
    
    condition = df[unitColumn] == "pg/L"
    df.loc[condition, convertedValueName] = df.loc[condition, ValueColumn] * pg_mg
    df.loc[condition, convertedUnitName] = "mg/L"
    
    condition = df[unitColumn] == "ug/L"
    df.loc[condition, convertedValueName] = df.loc[condition, ValueColumn] * ug_mg
    df.loc[condition, convertedUnitName] = "mg/L"
    
    condition = df[unitColumn] == "mg/L"
    df.loc[condition, convertedValueName] = df.loc[condition, ValueColumn]
    df.loc[condition, convertedUnitName] = "mg/L"
    
    # Other Unit Conversions
    #TEQ
    condition = df[unitColumn] == "ppt"
    df.loc[condition, convertedValueName] = df.loc[condition, ValueColumn] * ppt_pg * pg_mg
    df.loc[condition, convertedUnitName] = "mg/L"

    condition = df[unitColumn] == "ug/g"
    df.loc[condition, convertedValueName] = df.loc[condition, ValueColumn] * ugg_mgkg * mgkg_mgmL
    df.loc[condition, convertedUnitName] = "mg/L"

    condition = df[unitColumn] == "Bq/L"
    df.loc[condition, convertedValueName] = df.loc[condition, ValueColumn] * Bq_mg
    df.loc[condition, convertedUnitName] = "mg/L"

    # Flow Conversions
    condition = df[unitColumn] == "m3/d"
    df.loc[condition, convertedValueName] = df.loc[condition, ValueColumn] * d_yr    
    df.loc[condition, convertedUnitName] = "m3/yr"

    condition = df[unitColumn] == "m3/s W"
    df.loc[condition, convertedValueName] = df.loc[condition, ValueColumn] * s_yr
    df.loc[condition, convertedUnitName] = "m3/yr"

    if EMS:
        # Now that all the units are collected, correct for the different density of salt water for the samples not originally in mg/L
        units = np.logical_or(df[unitColumn] == "pg/L", np.logical_or(df[unitColumn] == "ppt", np.logical_or(df[unitColumn] == "ug/g", df[unitColumn] =="Bq/L")))
        crit = np.logical_and(df["SAMPLE_STATE"] == "Marine Water", units)
        df.loc[crit, convertedValueName] = df.loc[crit, convertedValueName] * mgkg_mgmL_saltWater
    
    # If reaplace, rewrite original value and unit columns, replacing unconverted 
    # values with nulls and delete new columns. 
    if replace:
        df[ValueColumn] = df[convertedValueName]
        df[unitColumn] = df[convertedUnitName]
        df.drop(columns=[convertedValueName, convertedUnitName], inplace=True)
        
    return df

###############################################################################
# Calculate the distance between each point (within km/10), and drop rows from 
# df2 closer than km kilometres. 
def filter_distance(df1, df2, 
                    lat_label1="Latitude", lat_label2="Latitude", 
                    lon_label1="Longitude", lon_label2="Longitude", 
                    dropping="", km=1, p=True):
    if p: print("\tIn Calculations: filter_distance...")
    df1.reset_index(inplace=True, drop=True)
    df2.reset_index(inplace=True, drop=True)

    df1[[lat_label1, lon_label1]] = df1[[lat_label1, lon_label1]].astype(float)
    df2[[lat_label2, lon_label2]] = df2[[lat_label2, lon_label2]].astype(float)
    
    # Check if anything is in Longitude [West] - consistency is key 
    # (and the geopy library prefers -West)
    if any(df2[lon_label2] > 0):
        Longitude_west = df2[lon_label2] > 0
        df2.loc[Longitude_west,lon_label2] = df2.loc[Longitude_west,lon_label2] * -1
        
    if any(df1[lon_label1] > 0):
        Longitude_west = df1[lon_label1] > 0
        df1.loc[Longitude_west,lon_label1] = df1.loc[Longitude_west,lon_label1] * -1
    
    # Prioritize keeping df1 values so drop df2 values where they overlap
    drop = []
    
    # For each df2 value, for every pari in df1 where a df2 pair is within 0.1 
    # degrees (to reduce calculation time), calculate the distance. If it is
    # closer than the specified km, add to drop list. 
    for i in range(len(df2)): # the one to drop
        #print("\t\t{}% Complete".format(int(i/len(df2)*100)))
        lat2 = df2.loc[i, lat_label2]
        lon2 = df2.loc[i, lon_label2]
        for n in range(len(df1)):
            lat1 = df1.loc[n, lat_label1]
            lon1 = df1.loc[n, lon_label1]
             # Null values throw an error, skip missing values. 
            if any([np.isnan(lat1), np.isnan(lat2), np.isnan(lon1), np.isnan(lon2)]):
                next
            else:
                # Too time consuming to check each pair, if lat longs closer 
                # than km/10, do the calculation.
                if any([abs(lat1 - lat2) < (km/10), abs(lon1 - lon2) < (km/10)]):
                    dist = geopy.distance.geodesic((lat1, lon1), (lat2, lon2)).km
                    if dist < km: 
                        drop.append(i)
                else:
                    next

    keep_df2 = df2.drop(drop)

    dropped = df2.iloc[drop]
    # print or return some stats here
    if p: 
        print("\tDropping {} complete. Number of rows to drop: {}, Number of unique"
              " latitudes dropped: {}, Number of unique longitudes dropped: {}"
              .format(dropping, len(drop), len(set(dropped[lon_label2])), 
                      len(set(dropped[lat_label2]))))

    return df1, keep_df2, dropped

###############################################################################
# !
# GENERAL - not reviewed yet
def calculateLoads(flows, concentration, source="", medium="", p=True):
    if p: print("\tIn Calculations: calculateLoads for {} in {}".format(source, medium))
    flows["Final Flow [m3/yr]"] = flows["Final Flow [m3/yr]"].astype(float)
    flows["Index"] = flows.index
    concentration.set_index("Unnamed: 0", inplace=True)

    info = ["Index", "Latitude", "Longitude", "Final Flow [m3/yr]"]
    temp = info.copy()
    temp.extend(concentration.columns)
    cols = temp
    df = pd.DataFrame(columns=cols)
    ###########
    # add flows information:
    df[info] = df[info].append(flows.loc[:,info])
    # calculate PAWP information:
    for i in concentration.columns:
        df[i] = concentration.loc["Concentration [mg/L]",i] * flows.loc[:, "Final Flow [m3/yr]"]  * m3_L * mg_kg


    df.insert(1, "Source", value=source)
    df.insert(2, "Medium", value=medium)
    df = df.fillna(0)
    
    df.to_csv("PAWP_Loads_{}.csv".format(source), index=None)
    if p: print("Finished Loads.")
    return df

"""    
###############################################################################
# Municipal Waste Water Treatment Plants (MWWTP) #
###############################################################################
""" 
###############################################################################
# Estimate wwtp flows and treatment levels

def estimate_flows_wwtp(df, path, flow_estimate_average_file=flow_estimate_average_file, p=True):
    if p: print("\tIn Calculations: estimate_flows_wwtp...")
    flow_estimate_average_file = "{}/{}".format(path, flow_estimate_average_file)
    flow_xls = pd.ExcelFile(flow_estimate_average_file)
    flow_avg = pd.read_excel(flow_xls, sheet_name="AssumptionsMWWTPList", skiprows=1).set_index("For MWWTP List, assume")
    removefromFurtherAnalysis = ["unknown or out of scope", "marina", "ferry", "washroom", "mall", "house", "houses", "marina floating home", "Chilliwack prison - to remove for some reason"] # few values or lack of information

    # Remove text and convert to float
    text = ['Data Unsubmitted', "No Data Available", 'Data Not Submitted']
    df.loc[df["Value"].isin(text)] = None
    df["Value"] = df["Value"].str.replace("[- ]", "")
    df["Value"] = df["Value"].astype(float)
    
    # Estimate missing flows based on facility type and population where
    # information does not already exist
    recreational = (df["SubType"].isin(["rv", "campground"])) & (df["Total Flow [m3/yr]"].isnull())
    df.loc[recreational, "Total Flow [m3/yr]"] = flow_avg.loc["Resorts, Campgrounds, RV Parks w STPs", "Flow per Facility"] * d_yr 
    
    mobHome = (df["SubType"] == "mobile home") & (df["Total Flow [m3/yr]"].isnull())
    df.loc[mobHome, "Total Flow [m3/yr]"] = flow_avg.loc["Mobile Home Parks", "Flow per Facility"] * d_yr 
    
    resort = (df["SubType"] == "resort") & (df["Total Flow [m3/yr]"].isnull())
    df.loc[resort, "Total Flow [m3/yr]"] = flow_avg.loc["Resorts, Apartment and Hotels", 'Flow per Facility'] * d_yr 

    minCamp = (df["SubType"] == "camp") & (df["Total Flow [m3/yr]"].isnull())
    df.loc[minCamp, "Total Flow [m3/yr]"] = flow_avg.loc["Mining camps w STPs or Septic", "Flow per Facility"] * d_yr
    
    home= (df["SubType"] == "home") & (df["Total Flow [m3/yr]"].isnull())
    df.loc[home, "Total Flow [m3/yr]"] = avgPerson * sewEst * L_m3 * d_yr# [1 house] * [avg Person/House] * [L/person*day] * [m3/L] * [d/yr] = [m3/yr]
    
    development = (df["SubType"] == "development") & (df["Total Flow [m3/yr]"].isnull())
    df.loc[development, "Total Flow [m3/yr]"] = flow_avg.loc["Development", 'Flow per Facility'] * d_yr
    
    FirstNations = (df["SubType"] == "First Nations") & (df["Total Flow [m3/yr]"].isnull())
    df.loc[FirstNations, "Total Flow [m3/yr]"] = flow_avg.loc["First Nation Reserve", 'Flow per Facility'] * d_yr
    
    prison = (df["SubType"] == "prison") & (df["Total Flow [m3/yr]"].isnull())
    df.loc[prison, "Total Flow [m3/yr]"] = flow_avg.loc["Prison", 'Flow per Facility'] * d_yr
    
    school = (df["SubType"] == "school") & (df["Total Flow [m3/yr]"].isnull())
    df.loc[school, "Total Flow [m3/yr]"] = flow_avg.loc["School - Day", 'Flow per Facility'] * d_yr
    
    marineTerminal = (df["SubType"] == "marine terminal") & (df["Total Flow [m3/yr]"].isnull())
    df.loc[marineTerminal, "Total Flow [m3/yr]"] = flow_avg.loc["Marine Terminal STP", 'Flow per Facility'] * d_yr
    
    # Remove specific types from further analysis
    drop = df[df["SubType"].isin(removefromFurtherAnalysis)].index
    df.drop(drop, inplace=True)
    
    # Do a rudimentary calculation from population to estimate municipal STPs 
    # from municipal population
    # ! Many STP/municipality? Large are contained in WSER so 
    # only smaller towns should be filled in here. 
    # "Population"
    municipal = (df["SubType"] == "municipal") & (df["Total Flow [m3/yr]"].isnull())
    df.loc[municipal, "Total Flow [m3/yr]"] = df.loc[municipal, "Population"] * sewEst * L_m3 * d_yr # [people] * [L/per*day] * [m3/L] * [d/yr] = [m3/yr]

# !
    # Estimate Treatment Type based on facility type (or other)
    
    return df
 
###############################################################################
# Roll up wwtp concentration data
def concentration_wwtp(path, wwtp_concentration_info_file=wwtp_concentration_info_file, p=True):
    if p: print("\tIn Calculations: concentration_wwtp...")
    # Read the excel file, then pull out Iona, Secondary Removal Efficiencies
    wwtp_concentration_info_file = "{}/{}".format(path, wwtp_concentration_info_file)
    xls_wwtp_concentration = pd.ExcelFile(wwtp_concentration_info_file)
    iona = pd.read_excel(xls_wwtp_concentration, "Iona_2014", skiprows=7)
    secondary_efficiencies = pd.read_excel(xls_wwtp_concentration, "Removal Efficiencies Secondary")
    
    # Make a removal efficiency dictionary
    # Assign PAWP classes
    secondary_efficiencies = dictionaries.pawp(secondary_efficiencies, path, parameterName="Contaminant")
    # Roll up %s
    secondary_efficiencies_avg = pd.DataFrame(columns=["PAWP_class", "Percent Removal"])
    secondary_efficiencies_avg["PAWP_class"] = [np.NaN] * len(set(secondary_efficiencies["PAWP_class"]))
    n = 0
    for pawp in list(set(secondary_efficiencies["PAWP_class"])):    
        avg_efficiency = np.mean(secondary_efficiencies.loc[secondary_efficiencies["PAWP_class"]==pawp, "Percent Removal"])
        secondary_efficiencies_avg.loc[n, ["PAWP_class", "Percent Removal"]] = pawp, avg_efficiency
        n += 1
    efficiency_dict = secondary_efficiencies_avg.set_index("PAWP_class").to_dict()
    
# ! not refining this code yet -_-
    # Process Iona into something useable    
    # Map PAWP catgeories
    iona = dictionaries.pawp(iona, path, parameterName="ANALYTE_UPDATED")
    # Fill in blank medians with non-detect mean limits by assigning a mean ND limit, and a separate column combining them
    iona.insert(loc=iona.columns.get_loc("Median"), column="ND Median Limit", value=np.NaN)
    iona["ND Median Limit"] = iona["Mean"].str.replace("<", "")
    iona["ND Median Limit"] = iona["ND Median Limit"]
    iona["Median else ND Mean"] = iona["Median"]
    iona.loc[iona["Median"].str.contains("-").fillna(False), "Median else ND Mean"] = iona.loc[iona["Median"].str.contains("-").fillna(False), "ND Median Limit"]

    iona.insert(loc=iona.columns.get_loc("Median"), column="ND Median.1 Limit", value=np.NaN)
    iona["ND Median.1 Limit"] = iona["Mean.1"].str.replace("<", "")
    iona["ND Median.1 Limit"] = iona["ND Median.1 Limit"]
    iona["Median.1 else ND Mean.1"] = iona["Median.1"]
    iona.loc[iona["Median.1"].str.contains("-").fillna(False), "Median.1 else ND Mean.1"] = iona.loc[iona["Median.1"].str.contains("-").fillna(False), "ND Median.1 Limit"]
    
    # Remove - so can do calculations
    iona["Influent Median [mg/L]"] = np.NaN
    iona.loc[:,"Median else ND Mean"] = iona.loc[:,"Median else ND Mean"].str.strip("-").fillna(iona.loc[:,"Median else ND Mean"])
    iona['Median else ND Mean'] = iona['Median else ND Mean'].replace('',np.nan).astype(float) 
    iona.loc[:,"Median.1 else ND Mean.1"] = iona.loc[:,"Median.1 else ND Mean.1"].str.strip("-").fillna(iona.loc[:,"Median.1 else ND Mean.1"])
    iona['Median.1 else ND Mean.1'] = iona['Median.1 else ND Mean.1'].replace('',np.nan).astype(float)
    
    # Convert Influent then Effluent units
    iona = units(iona, unitColumn="Converted_UNITS", ValueColumn="Median else ND Mean", convertedValueName="Influent Median [mg/L]", convertedUnitName="Influent_Units", replace=False)
    iona = units(iona, unitColumn="Converted_UNITS", ValueColumn="Median.1 else ND Mean.1", convertedValueName="Effluent Median [mg/L]", convertedUnitName="Effluent_Units", replace=False)

    # Make a concentration profile
    cols = list(set(iona["PAWP_class"]))
    cols.sort()
    concentration_profile = pd.DataFrame(columns=(["ConcentratioInfo"]+cols))
    concentration_profile["ConcentratioInfo"] = ["Total Influent Concentration [mg/L]", "Total Effluent Concentration [mg/L]", "Secondary Removal Efficiency [%]"]
    concentration_profile.set_index("ConcentratioInfo", inplace=True)

    for pawp in list(set(iona["PAWP_class"])):
        concentration_profile.loc["Total Influent Concentration [mg/L]", pawp] = iona.loc[iona["PAWP_class"]==pawp, "Influent Median [mg/L]"].sum()
        concentration_profile.loc["Total Effluent Concentration [mg/L]", pawp] = iona.loc[iona["PAWP_class"]==pawp, "Effluent Median [mg/L]"].sum()
    concentration_profile.loc["Secondary Removal Efficiency [%]"] = concentration_profile.columns.map(efficiency_dict["Percent Removal"]).fillna(0)

    return concentration_profile



###############################################################################
# 
def calculate_loads_wwtp(flows, concentrations, p=True):
    if p: print("\tIn Calculations: calculate_loads_wwtp...")
    sector = "MWWTP"
    #cols = ["Index", "Sector", "SubType", "Facility", "Company", "Medium", "DataSource", "DataSourceID", "Latitude", "Longitude", "Basin", "ParameterName", "Value", "Units", "PAWP_class"]
    #cols = cols + ["TreatmentLevel", "City", "Total Flow [m3/yr]", "CSO [m3/yr]", "Population"]

    # Make sure flows index is still 0, 1, 2....
    flows["Index"] = sector + flows.index.astype(str)#flows.index
    flows.set_index("Index", drop=False, inplace=True)
#    flows_dict = flows.set_index("Index").to_dict() # cannot just assign because they will be in a different order
#    flows["Index"] = flows.index

    # Shouls all be numeric by now, otherwise: .str.replace(",","").fillna() & convert to float/numeric/whatever works best & flows.loc[flows["Final Flow [m3/yr]"] == "", "Final Flow [m3/yr]"] = np.NaN
    flows = pd.concat([flows,pd.DataFrame(columns=concentrations.columns)], sort=False)
    
    # Primary: Assume the same effluent as Iona
    condition1 = (flows["TreatmentLevel"] == "Primary") | (flows["TreatmentLevel"] == "Preliminary")
    for i in concentrations.columns: 
        flows.loc[condition1, i] = concentrations.loc["Total Effluent Concentration [mg/L]",i] * flows.loc[condition1, "Total Flow [m3/yr]"] * m3_L * mg_kg # [kg/yr] = [mg/L] * [m3/yr] * [L/m3] * [kg/mg]
        
    # Secondary: Use what removal efficiencies we have....what, raw for the others? Currently using influent with 0% removal efficiency, so, raw. 
    # Treat the same as secondary due to lack of data
    condition2 = (flows["TreatmentLevel"] == "Secondary") | (flows["TreatmentLevel"] == "Tertiary")
    for i in concentrations.columns: 
        flows.loc[condition2, i] = concentrations.loc["Total Influent Concentration [mg/L]",i] * ((100.0 - concentrations.loc["Secondary Removal Efficiency [%]",i])/100.0) * flows.loc[condition2, "Total Flow [m3/yr]"]  * m3_L * mg_kg 

    # No treatment: use Iona Influent
    condition3 = (flows["TreatmentLevel"] == "Secondary") | (flows["TreatmentLevel"] == "Tertiary")
    for i in concentrations.columns: 
        flows.loc[condition3, i] = concentrations.loc["Total Influent Concentration [mg/L]",i] * flows.loc[condition3, "Total Flow [m3/yr]"] * m3_L * mg_kg 
       
    # No info/Unknown: Assume primary
    otherCondition = np.logical_or(~np.logical_or(np.logical_or(condition1, condition2), condition3), flows["TreatmentLevel"].isnull())
    for i in concentrations.columns: 
        flows.loc[otherCondition, i] = concentrations.loc["Total Effluent Concentration [mg/L]",i] * flows.loc[otherCondition, "Total Flow [m3/yr]"] * m3_L * mg_kg 

    flows.drop(columns=["ParameterName", "Value", "Units", "PAWP_class"], inplace=True)

    return flows







