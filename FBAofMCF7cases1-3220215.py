#
# This is a python code for the FBA of MCF7 cells
#
# All functions for FBA are available from pyfba module
# Please install pulp package before execution
# conda install -c conda-forge
#
import pyfba
#
# Preparation of FBA model
#
model = pyfba.MetabolicModel("RECON2_210317_2.txt", reversible = "<=>" , onedirection = "->", plus = "+", output = "normal")

#
# Cell producible metabolites (Consumption is not allowed)
#
model.add_inconstant_onlyincrease("M_nh4_b")
model.add_inconstant_onlyincrease("M_ala_L_b")
model.add_inconstant_onlyincrease("M_asn_L_b")
model.add_inconstant_onlyincrease("M_asp_L_b")
model.add_inconstant_onlyincrease("M_co2tot_b")
model.add_inconstant_onlyincrease("M_glu_L_b")
model.add_inconstant_onlyincrease("M_lac_L_b")
model.add_inconstant_onlyincrease("M_orn_b")
model.add_inconstant_onlyincrease("M_pro_L_b")
model.add_inconstant_onlyincrease("M_urea_b")
model.add_inconstant_onlyincrease("M_gudac_b")
model.add_inconstant_onlyincrease("M_biomass_b")
#
# Cell consumable metabolites in medium
#
model.add_inconstant("M_arg_L_b")
model.add_inconstant("M_glc_D_b")
model.add_inconstant("M_gln_L_b")
model.add_inconstant("M_ile_L_b")
model.add_inconstant("M_leu_L_b")
model.add_inconstant("M_lys_L_b")
model.add_inconstant("M_phe_L_b")
model.add_inconstant("M_thr_L_b")
model.add_inconstant("M_trp_L_b")
model.add_inconstant("M_tyr_L_b")
model.add_inconstant("M_val_L_b")
model.add_inconstant("M_his_L_b")
model.add_inconstant("M_met_L_b")
model.add_inconstant("M_chol_b")
model.add_inconstant("M_h_b")
model.add_inconstant("M_h2o_b")
model.add_inconstant("M_pi_b")
model.add_inconstant("M_o2_b")
model.add_inconstant("M_biomass_other_c")
#
# ent_ex_, gfe_ex_, and gfe0_ex_ represent net enthalpy change, net Gibbs free energy change,
# and net Gibbs free energy change at standard condtion levels, respectively.
#
model.add_inconstant("ent_ex_")
model.add_inconstant("gfe_ex_")
model.add_inconstant("gfe0_ex_")

#
# Uptake and production of amino acids are fixed at measured metabolic flux levels of MCF7 cells
#
r77_SubsArg  = 37 # Arginine specific uptake rate (nnol (10^6 cells h)-1)
r51_Glu_ex   = 16 # Glutamate specific production rate (nnol (10^6 cells h)-1)
r79_Orn_ex   = 43 # Ornitine specific production rate (nnol (10^6 cells h)-1)
r53_Pro_ex   = 4 # Proline specific production rate (nnol (10^6 cells h)-1)
r59_Ala_ex   = 22 # Alanine specific production rate (nnol (10^6 cells h)-1)


model.add_constraint("R_EX_arg_L_LPAREN_e_RPAREN_", r77_SubsArg * -1)
model.add_constraint("R_EX_orn_LPAREN_e_RPAREN_",  r79_Orn_ex)
model.add_constraint("R_EX_glu_L_LPAREN_e_RPAREN_",  r51_Glu_ex)
model.add_constraint("R_EX_ala_L_LPAREN_e_RPAREN_",  r59_Ala_ex)
model.add_constraint("R_EX_pro_L_LPAREN_e_RPAREN_",  r53_Pro_ex)
model.add_constraint("R_EX_asp_L_LPAREN_e_RPAREN_",  0)
model.add_constraint("R_EX_asn_L_LPAREN_e_RPAREN_",  0)
#
# Calclation of biomass production rate
#
mu = 0.039 # Specific cell proliferation rate (h-1)
biomass =  0.392 * 1000 * mu # microg (h)-1. Using the relationship of 0.392 mg DCW/ 10^6 cells
#
# Biomass production rate is fixed at measured metabolic flux levels of MCF7 cells
#
model.add_constraint("R_biomass_reaction", biomass) #biomass

#
# Objective function: Maximization of ATP consumption
#
model.add_optstep("R_DM_atp_c_", 1.0)
#
#
# Case 1, No additional constrains
#
print("Case 1")
status,objective = model.solve()
print("ATP consumption rate", model.get_value("R_DM_atp_c_"), "nmol (10^6 cells h)-1")
print("GAPDH", model.get_value("R_GAPD"), "nmol (10^6 cells h)-1")
print("ICDH", model.get_value("R_ICDHtotal"), "nmol (10^6 cells h)-1")
print("ETC", model.get_value("R_CYOOm3")* 2, "nmol (10^6 cells h)-1")
print("Enthalpy change", model.get_value("R_ent"), "mJ (10^6 cells h)-1")
model.show_result("MCF7_case1.txt")

#
# Case 2, Limilation in metabolic heat dicipation level
#
print("Case 2")
model.add_constraint("R_ent", -105) # Net enthalpy change level
status,objective = model.solve()
print("ATP consumption rate", model.get_value("R_DM_atp_c_"), "nmol (10^6 cells h)-1")
print("GAPDH", model.get_value("R_GAPD"), "nmol (10^6 cells h)-1")
print("ICDH", model.get_value("R_ICDHtotal") , "nmol (10^6 cells h)-1")
print("ETC", model.get_value("R_CYOOm3")* 2, "nmol (10^6 cells h)-1")
print("Enthalpy change", model.get_value("R_ent"), "mJ (10^6 cells h)-1")
model.show_result("MCF7_case2.txt")

#
# Case 3, Limilation in metabolic heat dicipation and electrain transport chain levels
#
print("Case 3")
model.add_constraint("R_ent", -105) # Net enthalpy change level
model.add_constraint("R_CYOOm2", 0) # Metabolic flux level of electron transport chain (Complex iV)
model.add_constraint("R_CYOOm3", 141) # Metabolic flux level of Complex IV divided by 2 (282/2)

status,objective = model.solve()
print("ATP consumption rate", model.get_value("R_DM_atp_c_"), "nmol (10^6 cells h)-1")
print("GAPDH", model.get_value("R_GAPD"), "nmol (10^6 cells h)-1")
print("ICDH", model.get_value("R_ICDHtotal"), "nmol (10^6 cells h)-1")
print("ETC", model.get_value("R_CYOOm3")* 2, "nmol (10^6 cells h)-1")
print("Enthalpy change", model.get_value("R_ent"), "mJ (10^6 cells h)-1")
model.show_result("MCF7_case3.txt")


