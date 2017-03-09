families = {}
families['a3a'] = ['4xxo_A','4xxo_B','2m65_A','5sww_AE','5sww_BF','5sww_CG','5sww_DH','5td5_AC']
families['a3b_ntd'] = ['hom_a3b_ntd_shiv']
families['a3b_ctd'] = ['5cqd_A','5cqh_A','5cqi_A','5cqk_A','5cqd_C','hom_a3b_ctd_shiv','2nbq_A']
families['a3c'] = ['3vm8_A','3vm8_B','3vow_A','3vow_B']
families['a3d_ntd'] = ['hom_a3d_ntd_shiv']
families['a3d_ctd'] = ['hom_a3d_ctd_shiv']
families['a3f_ntd'] = ['hom_a3f_ntd_shiv']
families['a3f_ctd'] = ['3wus_A','3wus_B','4j4j_A','4j4j_B', '4iou_A','4iou_B','4iou_C','4iou_D','5hx4_A','5hx4_B','5hx5_A','5hx5_B'] 
families['a3g_ntd'] = ['2mzz_A','hom_a3g_ntd_shiv','5k81_A','5k81_B','5k81_C','5k81_D','5k81_E','5k81_F','5k82_A','5k82_B','5k82_C','5k82_D', '5k83_AH','5k83_B','5k83_CI','5k83_D','5k83_EJ','5k83_F'] 
families['a3g_ctd'] = ['2jyw_A','3e1u_A','2kbo_A','2kem_A','3iqs_A','3ir2_A','3ir2_B','3v4k_A','3v4k_B','3v4j_A','3v4j_B','4rov_A','4rov_B','4row_A']
families['a3h'] = ['hom_a3h_shiv']

families['aid'] = ['5jj4_A','5jj4_B','5jj4_C']

families['z1'] = families['a3a'] + families['a3b_ctd'] + families['a3g_ctd']
families['z2'] = families['a3b_ntd'] + families['a3c'] + families['a3d_ntd'] + families['a3d_ctd'] + families['a3f_ntd'] + families['a3f_ctd'] + families['a3g_ntd']
families['z3'] = families['a3h']

families['apobec2'] = ['2nyt_A','2nyt_B','2nyt_C','2nyt_D']

families['dna_bound'] = ['5sww_AE','5sww_BF','5sww_CG','5sww_DH','5td5_AC','5k83_AH','5k83_CI','5k83_EJ']

families['all_a3'] = families['z1'] + families['z2'] + families['z3']
families['all'] = families['all_a3'] + families['aid'] + families['apobec2']

## FAMILIES TO ADD
# vif_binding
# vif_nonbinding
