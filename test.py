import parsing_zacros
import modifing_zacros
import printing_zacros
import linecache

linecache.clearcache()

tot_steps,events,time,temperature_interval,temperature,energy,adsorbate,specnum_result = parsing_zacros.Parsing_Specnum()

index_list = [1,2]
added_string_list,added_list = modifing_zacros.Adding_Specnum(adsorbate,specnum_result,index_list)

ideal_temperature_interval = 5
index_list = [7,8]
temperature_list, all_TPD_list, TPD_list = modifing_zacros.Generate_TPD_List(added_list,ideal_temperature_interval,temperature_interval,temperature,index_list)

name_of_list = ['CO2_Cu','H2']
printing_zacros.Print_TPD_figure_data(name_of_list,temperature_list,TPD_list)

linecache.clearcache()
        