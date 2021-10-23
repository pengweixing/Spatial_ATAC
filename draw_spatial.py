# -*- coding: utf-8 -*-
#################################################
#  File Name:draw_spatial_info.py
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Fri 27 Aug 2021 02:35:54 PM UTC
#################################################

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib import cm
import numpy as np
import pandas as pd
import matplotlib as mpl
import collections
from matplotlib.lines import Line2D


def get_path(data):
    all_path = []
    cell_index = []
    for row in range(data.shape[0]):
        for col in range(data.shape[1]):
            left_bottom = (row,col)
            left_top = (row,col+1)
            right_top = (row+1,col+1)
            right_bottom = (row+1,col)   
            close = left_bottom
            all_path.append([left_bottom,left_top,right_top,right_bottom,close])
            cell_index.append([row,col])
        
    return all_path,cell_index

def draw_spatial(data,output,vmin):
#    data = pd.read_table(data)
    data_num = data.iloc[:,0:]
    data_num = np.array(data_num)
    
    y_label = np.array(list(data.columns)[0:])
    x_label = np.array(list(data.index))
    
    colormap = cm.get_cmap('jet')
    colormap2 = cm.colors.ListedColormap(colormap(np.linspace(0, 1, 1024)))
    
    norm = plt.Normalize(data_num[data_num>-10].min(), data_num.max())
    data_norm = data_num/data_num.max()
    for_colorbar = plt.cm.ScalarMappable(cmap=colormap2, norm=norm)
    for_colorbar.set_array([])
    all_verts,cell_index = get_path(data_norm)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=1)
    
    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]
    
    fig, ax = plt.subplots(figsize=(12,10),frameon=False)
    for each_vert,each_index in zip(all_verts,cell_index):
        path = Path(each_vert, codes)
        if data_num[each_index[0],data_num.shape[0]-1-each_index[1]] > -10 :
            patch = patches.PathPatch(path, edgecolor='w',facecolor=colormap2(norm(data_norm[each_index[0],data_num.shape[0]-1-each_index[1]])),lw=2)
            ax.add_patch(patch)
        else:
            patch = patches.PathPatch(path, edgecolor='w',facecolor='w',lw=2)
            ax.add_patch(patch)
            
    fig.colorbar(for_colorbar, ax=ax, label='deviation',shrink=0.50)
    
    ax.set_xlim(0, data_num.shape[0])
    ax.set_ylim(0, data_num.shape[1])

    xtick_index = np.arange(len(x_label))[0:len(x_label):5]
    ytick_index = np.arange(len(y_label))[0:len(y_label):5]
    
    ax.set_xticks(xtick_index+0.5)
    ax.set_yticks(ytick_index+0.5)
    
    ax.set_xticklabels(x_label[xtick_index])
    ax.set_yticklabels(y_label[data_num.shape[0]-1-ytick_index])
    
    plt.setp(ax.get_xticklabels(), ha="center")
    plt.setp(ax.get_yticklabels(), va="center")
    
    #plt.show()
    plot_img = output + '.motif.pdf'
    fig.savefig(plot_img,format='pdf')
    plt.cla()
    plt.close(fig)

def build_legend(data):
    """
    Build a legend for matplotlib plt from dict
    """
    legend_elements = []
    for key in data:
        legend_elements.append(Line2D([0], [0], marker='s', color='w', label=key,
                                        markerfacecolor=data[key], markersize=20))
    return legend_elements


def draw_spatial_cluster(data,color):
    """
    project the clusters back to spatial iamgee
    """
    color_dict = dict(zip(color['group'],color['colour'])) 
    color_dict = collections.OrderedDict(sorted(color_dict.items()))
    x_label = np.array(list(data.columns)[0:])
    y_label = np.array(list(data.index))
    data_num = data.iloc[:,0:]
    data_num = np.array(data_num)
    all_verts,cell_index = get_path(data_num)
    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]

    fig, ax = plt.subplots(figsize=(10,10),frameon=False)
    for each_vert,each_index in zip(all_verts,cell_index):
        path = Path(each_vert, codes)
        temp_value = data_num[each_index[0],data_num.shape[0]-each_index[1]-1]
        if not temp_value == -100 :
            patch = patches.PathPatch(path, edgecolor='w',facecolor=color_dict[temp_value],lw=2)
            ax.add_patch(patch)
        else:
            patch = patches.PathPatch(path, edgecolor='w',facecolor='w',lw=2)
            ax.add_patch(patch)
    ax.set_xlim(0, data_num.shape[0])
    ax.set_ylim(0, data_num.shape[1])
    xtick_index = np.arange(len(x_label))[0:len(x_label):5]
    ytick_index = np.arange(len(y_label))[0:len(y_label):5]
    ax.set_xticks(xtick_index+0.5)
    ax.set_yticks(ytick_index+0.5)
    ax.set_xticklabels(y_label[ytick_index])
    ax.set_yticklabels(x_label[data_num.shape[0]-xtick_index-1])
    plt.setp(ax.get_yticklabels(), ha="right")
    plt.setp(ax.get_xticklabels(), va="top")
    legend_elements = build_legend(color_dict)
    ax.legend(handles=legend_elements, loc='upper right',
              bbox_to_anchor=(1.2, 1),fontsize=20)
    fig.savefig("project_cluster_to_spatial.pdf",format='pdf',bbox_inches='tight')
