import numpy as np
import math

class cluster_info:
    def __init__(self, small_adj_matrix, inter_initialise, result_interaction, site_type, angle_list, adsorbate_list, cluster_energy):
        self.small_adj_matrix = small_adj_matrix
        self.inter_initialise = False
        self.result_interaction = []
        self.site_type = site_type
        self.angle_list = angle_list
        self.adsorbate_list = adsorbate_list
        self.cluster_energy = cluster_energy
                
    def recursive_small(self, list_interaction=None, y_traj=None, y=0):
        if list_interaction is None:
            list_interaction = []
        if y_traj is None:
            y_traj = []
        for i in range(self.small_adj_matrix.shape[0]):  # 逐个检测site i
            # print(i)
            if self.small_adj_matrix[y][i] == 1: # 找到大表中一个connection, 第一次一定成功
                if i not in list_interaction: # 第一次出现site i
                    y_traj.append(y)
                    list_interaction.append(y)
                    list_interaction.append(i)
                    self.small_adj_matrix[y][i] = 0
                    self.small_adj_matrix[i][y] = 0
                    # print(list_interaction)
                    if np.sum(self.small_adj_matrix) == 0:
                        self.result_interaction = [0, 1] + [x + 1 for x in list_interaction]
                        self.initialise = True
                        return self.result_interaction, self.inter_initialise
                    else:
                        return self.recursive_small(list_interaction, y_traj, i) # 换行搜索
                else:  # 找到和原来某site 之间的关系（要求）了，只记录，继续在这行找
                # elif i != list_interaction[-2]: # 找到和原来某site 之间的关系（要求）了，只记录，继续在这行找
                    list_interaction.append(y)
                    list_interaction.append(i)
                    self.small_adj_matrix[y][i] = 0
                    self.small_adj_matrix[i][y] = 0
                    if np.sum(self.small_adj_matrix) == 0:
                        self.result_interaction = [0, 1] + [x + 1 for x in list_interaction]
                        # print('else', list_interaction)
                        self.initialise = True
                        return self.result_interaction, self.inter_initialise
        
        # 此时已查询过所有i，但没有一个符合要求，需要退回上一步
        # print('no match')
        if y_traj != []: 
            y_new = y_traj[-1]
        else: # monodentate
            self.result_interaction = [0,1]
            self.initialise = True
            return self.result_interaction, self.inter_initialise
        del y_traj[-1:] # 回到 y 的前一次
        return self.recursive_small(list_interaction, y_traj, y_new)
        
class surface_info:
    def __init__(self, big_adj_matrix, coordinate_list, site_type, adsorbate_list, lattice_energy, energy_initialise):
        self.big_adj_matrix = big_adj_matrix
        self.coordinate_list = coordinate_list
        self.site_type = site_type
        self.adsorbate_list = adsorbate_list
        self.lattice_energy = lattice_energy
        self.energy_initialise = False
    
    def angle_ABC(self, a, b, c):
        # 提取坐标
        x1, y1 = self.coordinate_list[a]
        x2, y2 = self.coordinate_list[b]
        x3, y3 = self.coordinate_list[c]

        # 计算向量 BA 和 BC
        ba = (x1 - x2, y1 - y2)
        bc = (x3 - x2, y3 - y2)

        # 计算点积
        dot_product = ba[0] * bc[0] + ba[1] * bc[1]

        # 计算向量长度
        magnitude_ba = math.sqrt(ba[0] ** 2 + ba[1] ** 2)
        if magnitude_ba == 0:
            magnitude_ba = 0.00001
        magnitude_bc = math.sqrt(bc[0] ** 2 + bc[1] ** 2)
        if magnitude_bc == 0:
            magnitude_bc = 0.00001

        # 计算 cos(θ)
        cos_theta = dot_product / (magnitude_ba * magnitude_bc)

        # 计算角度（防止浮点数误差导致 cos 超过 [-1,1]）
        cos_theta = max(-1, min(1, cos_theta))
        theta = math.acos(cos_theta)  # 反余弦（弧度）
        
        # 转换为度
        angle = math.degrees(theta)
        
        return angle
    
    def update_surface(self, cluster_info_obj, result_subgraph_after_choose): # update adsorbate_list
        if len(result_subgraph_after_choose[0]) == 0:
            return 0
        for i in range(len(result_subgraph_after_choose[0])):
            for j in range(len(self.adsorbate_list)): # check which adsorbate needs to be deleted
                if result_subgraph_after_choose[0][i] in self.adsorbate_list[j]:
                    self.adsorbate_list[j].insert(0,-1)
        for j in reversed(range(len(self.adsorbate_list))): # delete target adsorbate
            if self.adsorbate_list[j][0] == -1:
                del self.adsorbate_list[j]
        check_repeat = []
        for i in range(len(result_subgraph_after_choose[1])): # 对每一个要加进去的 site number
            for j in range(len(cluster_info_obj.adsorbate_list)): # 对每一个新cluster上的adsorbate
                if result_subgraph_after_choose[1][i] in cluster_info_obj.adsorbate_list[j]: # 查找这个site number对应新cluster的哪一个adsorbate                 
                    add_ads = [cluster_info_obj.adsorbate_list[j][0]]
                    for k in range(len(cluster_info_obj.adsorbate_list[j])-1): # 替换site在cluster和surface上的编号： [X**, 1, 2] -> [X**, 45, 67]
                        add_ads.append(result_subgraph_after_choose[0][result_subgraph_after_choose[1].index(cluster_info_obj.adsorbate_list[j][k+1])])
                    if add_ads[1] not in check_repeat:
                        check_repeat.append(add_ads[1])
                        self.adsorbate_list.append(add_ads)
                        break
        self.energy_initialise = False
        return self.adsorbate_list, self.energy_initialise
        
def find_subgraph(surface_info_obj, cluster_info_obj, result_subgraph=None, dont_search=None, dont_search2=None, pre=None, y=0, i=1):
    # 目的：遍历第 y 行的第 j (1-L) 个 site 是否符合要求 i (从1开始)
    # print('start', dont_search, dont_search2, pre, y, i)
    if result_subgraph == None:
        result_subgraph = []
    if dont_search == None:
        dont_search =  [[] for _ in range(len(cluster_info_obj.result_interaction)//2)]
    if dont_search2 == None:
        dont_search2 = []
    if pre == None:
        pre = []

    if pre != []:
        pre_index = cluster_info_obj.result_interaction.index(cluster_info_obj.result_interaction[2*i-2])
        # print(i,pre_index, result_interaction)
        if y != pre[pre_index]: # 如果当前搜索行与要求中的搜索行不一致 (搜索行是搜索对[a,b]中第一个数 a)
            # 去搜要求的搜索行
            return find_subgraph(surface_info_obj, cluster_info_obj, result_subgraph, dont_search, dont_search2, pre, pre[pre_index], i) 

    # case = 0
    for j in range(surface_info_obj.big_adj_matrix.shape[0]):  # 逐个检测site j
        index_in_subgraph = cluster_info_obj.result_interaction[2*i-1]-1
        # 找到大表中一个符合要求的connection
        if surface_info_obj.big_adj_matrix[y][j] == 1 and surface_info_obj.site_type[j] == cluster_info_obj.site_type[index_in_subgraph]:

            # 判断对应 site 上所成的角度与要求是不是一致
            angle_criteria = True
            if j not in pre:
                pre_tem = pre + [y,j]
                lattice_angle_site = [[-1,-1,-1] for _ in range(len(cluster_info_obj.angle_list))]
                for k in range(len(cluster_info_obj.angle_list)): # 对于每个角度条件k
                    for l in range(3):
                        if len(pre_tem) > cluster_info_obj.result_interaction.index(cluster_info_obj.angle_list[k][l+1]):
                            # 完成建立 lattice_angle_site 和 angle_list 之间的映射
                            lattice_angle_site[k][l] = pre_tem[cluster_info_obj.result_interaction.index(cluster_info_obj.angle_list[k][l+1])]
                for k in range(len(cluster_info_obj.angle_list)):
                    if all(x > 0 for x in lattice_angle_site[k]):
                        angle = surface_info_obj.angle_ABC(lattice_angle_site[k][0], lattice_angle_site[k][1], lattice_angle_site[k][2])
                        if angle-1 > cluster_info_obj.angle_list[k][0] or cluster_info_obj.angle_list[k][0] > angle+1: # 角度不对
                            angle_criteria = False

            # 判断对应 site 上的 adsorbate 和 lattice 上的 adsorbate 与该 adsorbate 的 dentate number 是不是一致
            if angle_criteria == True:
                pre_tem = pre + [y,j]
                lattice_adsorbate_type = [[-1] * (cluster_info_obj.adsorbate_list[k] - 1) for k in range(len(cluster_info_obj.adsorbate_list))]
                for k in range(len(cluster_info_obj.adsorbate_list)): # 对于每个adsorbate的条件
                    for l in range(len(lattice_adsorbate_type[k])): # 一个条件中的每个dendate 
                        if len(pre_tem) > cluster_info_obj.result_interaction.index(cluster_info_obj.adsorbate_list[k][l+1]):
                            # 完成建立 lattice_adsorbate_site 和 adsorbate_list 之间的映射
                            lattice_adsorbate_type[k][l] = pre_tem[cluster_info_obj.result_interaction.index(cluster_info_obj.adsorbate_list[k][l+1])]
                    lattice_adsorbate_type[k].insert(0, cluster_info_obj.adsorbate_list[k][0]) # 补上 每个adsorbate_list_small 的名字
                for k in range(len(cluster_info_obj.adsorbate_list)):
                    if all(x > 0 for x in lattice_adsorbate_type[k][1:]):
                        if lattice_adsorbate_type[k] not in surface_info_obj.adsorbate_list:
                            adsorbate_criteria = False

            if adsorbate_criteria == True:
                if ((not dont_search[i-1]) or all([y,j] != sublist for sublist in dont_search[i-1])) and j not in dont_search2:
                    if cluster_info_obj.result_interaction[2*i-1] in cluster_info_obj.result_interaction[:2*i-2]: # "要求"不涉及新site
                        pre_index = pre[cluster_info_obj.result_interaction.index(cluster_info_obj.result_interaction[2*i-1])] 
                        if j == pre[pre_index]: # 找到的这个 j 必须符合"要求"中指定的那一个site, (对应搜索对[a,b]中第二个数 b)
                            if 2 * i == len(cluster_info_obj.result_interaction): # 所有要求都检测完成
                                pre.append(y)
                                pre.append(j)
                                result_subgraph.append(pre[2:]) # 记录该映射 (从2开始)
                                del pre[-2:]
                                dont_search2.append(j) # 重新搜索本行，但不要再搜到这个 j
                                return find_subgraph(surface_info_obj, cluster_info_obj, result_subgraph, dont_search, dont_search2, pre, y, i)
                            else: # 通过当前要求，但尚未通过所有要求，继续检测剩余要求
                                pre.append(y)
                                pre.append(j)
                                # 这里“要求”被符合，开始判断下一要求（迭代）,注意因为不涉及新site,这行没有搜完
                                return find_subgraph(surface_info_obj, cluster_info_obj, result_subgraph, dont_search, dont_search2, pre, y, i+1)

                    else: # “要求” 涉及新 site， 第一次一定在这个分支中
                        if j not in pre: # 新找到的 j 必须不能和原来找到的任何 big_site index 一致
                            if 2 * i == len(cluster_info_obj.result_interaction): # 所有要求都检测完成
                                pre.append(y)
                                pre.append(j)
                                if pre[2:] != []:
                                    result_subgraph.append(pre[2:]) # 记录该映射
                                else: # monodentate
                                    result_subgraph.append([j])
                                del pre[-2:]
                                dont_search2.append(j)
                                return find_subgraph(surface_info_obj, cluster_info_obj, result_subgraph, dont_search, dont_search2, pre, y, i)
                            else: # 通过当前要求，但尚未通过所有要求，继续检测剩余要求
                                pre.append(y)
                                pre.append(j)
                                # 符合“要求” 才有继续查询是否满足其他要求的必要，这里要求被符合，进入迭代，搜下一行
                                return find_subgraph(surface_info_obj, cluster_info_obj, result_subgraph, dont_search, dont_search2, pre, j, i+1)
    
    # 此时已查询过所有j，但没有一个符合要求，需要退回上一步,继续搜索符合上一个要求的下一组yj
    # print('preif',i)
    if i == 1: # 再也搜不到一个了
        # print('afterif',i)
        for k in range(len(result_subgraph)):
            result_subgraph[k] = list(dict.fromkeys(result_subgraph[k]))
        cluster_info_obj.result_interaction = list(dict.fromkeys(cluster_info_obj.result_interaction))
        if len(cluster_info_obj.result_interaction) == 2: # monodentate
            result_subgraph.append(cluster_info_obj.result_interaction[1:])
        else:
            result_subgraph.append(cluster_info_obj.result_interaction[2:])
        return result_subgraph
    # print(dont_search[0])
    # print(dont_search2)
    # print(pre)
    # print(i)
    # print(result)
    # print(pre[-1])
    dont_search[i-2].append([pre[-2], pre[-1]].copy())
    dont_search2 = []
    del pre[-2:]
    if len(pre) == 0:
        return find_subgraph(surface_info_obj, cluster_info_obj, result_subgraph, dont_search, dont_search2, pre, 0, i-1)
    return find_subgraph(surface_info_obj, cluster_info_obj, result_subgraph, dont_search, dont_search2, pre, pre[-1], i-1)

def perform_graph_isomorphism(surface_info_obj, cluster_info_obj):

    cluster_info.recursive_small()
    result_subgraph = find_subgraph(surface_info_obj, cluster_info_obj)

    return result_subgraph
