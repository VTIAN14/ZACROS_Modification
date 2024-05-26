import parsing_zacros
import modifing_zacros
import printing_zacros
import linecache

linecache.clearcache()

num_of_data,all_tot_steps,adsorbate,temperature,temperature_interval,all_specnum_result = parsing_zacros.Parsing_Specnum_in_Different_Folder()

ads_aver_list = ['CO2_Cu','H2']
a_list = modifing_zacros.Average_Specnum(ads_aver_list,adsorbate,all_tot_steps,all_specnum_result)

#index_list = [1,2]
#added_string_list,added_list = modifing_zacros.Adding_Specnum(adsorbate,specnum_result,index_list)

#ideal_temperature_interval = 5
#index_list = [7,8]
#temperature_list, all_TPD_list, TPD_list = modifing_zacros.Generate_TPD_List(added_list,ideal_temperature_interval,temperature_interval,temperature,index_list)

#name_of_list = ['CO2_Cu','H2']
#printing_zacros.Print_TPD_figure_data(name_of_list,temperature_list,TPD_list)

#printing_zacros.Print_TPD_figure_data(name_of_list,temperature_list, all_TPD_list)

# ads_aver_list = ['CO2_Cu','CO2_Pd','H2']

# index_list = [1,2]

# ideal_temperature_interval = 1

# num_of_data,all_tot_steps,adsorbate,temperature,temperature_interval,all_specnum_result = parsing_zacros.Parsing_Specnum_in_Different_Folder()

# a_list = modifing_zacros.Average_Specnum(ads_aver_list,adsorbate,all_tot_steps,all_specnum_result)

# name_of_list,a_list = modifing_zacros.Adding_Specnum(ads_aver_list,a_list,index_list)

# temperature_list, all_TPD_list = modifing_zacros.Generate_TPD_List(a_list,ideal_temperature_interval,temperature_interval,temperature)

# #name_of_list = ['CO2_Cu','H2']
# printing_zacros.Print_TPD_figure_data(name_of_list,temperature_list,TPD_list)

linecache.clearcache()
        