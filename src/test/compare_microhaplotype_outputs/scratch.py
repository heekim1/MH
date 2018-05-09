# from __future__ import print_function
# from deepdiff import DeepDiff
# from txt_file import TxtFile
# from json_file import JsonFile
#
# my_out = TxtFile('./test_compare_microhaplotypes/inputs_test_compare_microhaplotypes/IonXpress_003.txt')
# my_out.parse()
# json_out = JsonFile('./test_compare_microhaplotypes/inputs_test_compare_microhaplotypes/converge_results.json')
# json_out.parse()
# diff_dict = DeepDiff(my_out.mh_dict, json_out.mh_dict)
# print (diff_dict)
#
# exclusion_list = ['avg_mapping_quality']
# if 'type_changes' in diff_dict:
#     filtered_dict = {}
#     for changed_type in diff_dict['type_changes']:
#         changed_type_name = changed_type.split('.')[-1]
#         if changed_type_name not in exclusion_list:
#             filtered_dict[changed_type] = diff_dict['type_changes'][changed_type]
#
#     if len(filtered_dict) > 0:
#         diff_dict['type_changes'] = filtered_dict
#     else:
#         del diff_dict['type_changes']
#
# print ("--------------------------------")
# print (diff_dict)
#
