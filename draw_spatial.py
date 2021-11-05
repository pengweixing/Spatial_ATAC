#################################################
#  File Name:draw_spatial.py
#  Author: xingpengwei
#  Mail: xingwei421@qq.com, pengwei.xing@igp.uu.se
#  Created Time: Thu 21 Oct 2021 05:55:53 PM UTC
#################################################
"""
This is the core module for projecting gene or cluster or motif spatial matrix 
to 2-D graph or staining tissue image
"""
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
import matplotlib.image as mpimg
import cv2

class SpatialPlot:
    def __init__(self,offset_r = 1, alpha=1, offset_l=0,dpi=300, started_column = 1, oformat = 'pdf' , colormap = 'jet', data = [],width=12,height=10):
        """
        started_column: the column of the number in the matrix started
        offset : control the pixel position
        """
        self.bins = 1024
        self.offset_r = offset_r
        self.offset_l = offset_l
        self.cell_index = []
        self.all_path = []
        self.legend_elements = [] 
        self.color = colormap
        self.y_label = np.array(list(data.columns)[started_column-1:])
        self.x_label = np.array(list(data.index))
        self.started_column = started_column - 1
        self.data_num = data.iloc[:,self.started_column:]
        self.data_num = np.array(self.data_num)
        self.colormap = cm.get_cmap(self.color)
        self.colormap2 = cm.colors.ListedColormap(self.colormap(np.linspace(0, 1, self.bins)))
        self.linewidth = 2
        self.edgecolor ='w'
        self.figsize_w = width
        self.figsize_h = height
        self.alpha = alpha
        self.format = oformat
        self.dpi = dpi
        self.codes = [
            Path.MOVETO,
            Path.LINETO,
            Path.LINETO,
            Path.LINETO,
            Path.CLOSEPOLY,]

    def Get_Path(self,data):
        """
        Build the path for all pixels
        """
        all_path = self.all_path
        cell_index = self.cell_index
        offset_l = self.offset_l
        offset_r = self.offset_r
        for row in range(data.shape[0]):
            for col in range(data.shape[1]):
                left_bottom = (row+offset_l,col+offset_l)
                left_top = (row+offset_l,col+offset_r)
                right_top = (row+offset_r,col+offset_r)
                right_bottom = (row+offset_r,col+offset_l)
                close = left_bottom
                all_path.append([left_bottom,left_top,right_top,right_bottom,close])
                cell_index.append([row,col])
        return all_path,cell_index


    def Build_Legend(self,data,which_cluster=''):
        """
        Build a legend for matplotlib plt from dict
        """
        legend_elements = self.legend_elements
        for key in data:
            if key in which_cluster:
                legend_elements.append(Line2D([0], [0], marker='s', color='w', label=key,
                                            markerfacecolor=data[key], markersize=20))
        return legend_elements
   
    def add_axes_attr(self,fig='',ax='',forcolorbar=''):
        if forcolorbar == '':
            pass
        else:
            fig.colorbar(forcolorbar, ax=ax, label='Score',shrink=0.50)
        ax.set_xlim(0, self.data_num.shape[0])
        ax.set_ylim(0, self.data_num.shape[1])
        xtick_index = np.arange(len(self.x_label))[0:len(self.x_label):5]
        ytick_index = np.arange(len(self.y_label))[0:len(self.y_label):5]
        ax.set_xticks(xtick_index+0.5)
        ax.set_yticks(ytick_index+0.5)
        ax.set_xticklabels(self.x_label[xtick_index])
        ax.set_yticklabels(self.y_label[self.data_num.shape[0]-1-ytick_index])
        plt.setp(ax.get_xticklabels(), ha="center")
        plt.setp(ax.get_yticklabels(), va="center")

    def project_GeneMatrix_to_spatial(self,data=[],output=''):
        """
        plot the gene score or motif deviation to spatial grpah by the coordination
        """
        norm = plt.Normalize(np.nanmin(self.data_num), np.nanmax(self.data_num))
        data_norm = (self.data_num - np.nanmin(self.data_num))/(np.nanmax(self.data_num) - np.nanmin(self.data_num))
        for_colorbar = plt.cm.ScalarMappable(cmap=self.colormap2, norm=norm)
        for_colorbar.set_array([])
        all_verts,cell_index = self.Get_Path(data_norm)     
        norm = mpl.colors.Normalize(vmin=0, vmax=1)
        width = self.figsize_w
        height = self.figsize_h
        fig, ax = plt.subplots(figsize=(self.figsize_w,self.figsize_h),frameon=False,dpi=self.dpi)    
        for each_vert,each_index in zip(all_verts,cell_index):
            path = Path(each_vert, self.codes)
            if not pd.isna(self.data_num[each_index[0],self.data_num.shape[0]-1-each_index[1]]) :
                patch = patches.PathPatch(path, edgecolor=self.edgecolor,facecolor=self.colormap2(norm(data_norm[each_index[0],self.data_num.shape[0]-1-each_index[1]])),lw=self.linewidth)
                ax.add_patch(patch)
            else:
                patch = patches.PathPatch(path, edgecolor=self.edgecolor,facecolor='w',lw=self.linewidth)
                ax.add_patch(patch)
        self.add_axes_attr(fig=fig,ax=ax,forcolorbar=for_colorbar)
        ax.set_title(output,fontsize=20)

        if self.format == 'pdf':
            plot_img = output + '.spatial.pdf'
        elif self.format == 'png':
            plot_img = output + '.spatial.png'
        elif self.format == 'svg':
            plot_img = output + '.spatial.svg'
        elif self.format == 'eps':
            plot_img = output + '.spatial.eps'
        fig.savefig(plot_img,format=self.format,bbox_inches='tight')
        plt.cla()
        plt.close(fig)       
 
    
    def project_GeneMatrix_to_tissue(self,data=[],output='',image='',color='red'):
        """
        plot the gene score or motif deviation to spatial grpah by the coordination
        """
        img = cv2.imread(image)
        norm = plt.Normalize(np.nanmin(self.data_num), np.nanmax(self.data_num))
        data_norm = (self.data_num - np.nanmin(self.data_num))/(np.nanmax(self.data_num) - np.nanmin(self.data_num))

        color_temp = list(matplotlib.colors.to_rgba(color))
        color_temp = np.repeat(np.array(color_temp,ndmin=2),self.bins,axis=0)
        color_temp[:,3] = np.linspace(0,1,self.bins)
        colormap = cm.colors.ListedColormap(color_temp)
        norm = plt.Normalize(np.nanmin(self.data_num), np.nanmax(self.data_num))
        for_colorbar = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
        for_colorbar.set_array([])
        norm2 = mpl.colors.Normalize(vmin=0, vmax=1)

        all_verts,cell_index = self.Get_Path(data_norm)     
        norm = mpl.colors.Normalize(vmin=0, vmax=1)
        width = self.figsize_w
        height = self.figsize_h
        fig, ax = plt.subplots(figsize=(self.figsize_w,self.figsize_h),frameon=False,dpi=self.dpi)    
        ax.imshow(img,extent=(0,self.data_num.shape[0],0,self.data_num.shape[0]))
        for each_vert,each_index in zip(all_verts,cell_index):
            path = Path(each_vert, self.codes)
            if not pd.isna(self.data_num[each_index[0],self.data_num.shape[0]-1-each_index[1]]) :
                patch = patches.PathPatch(path, edgecolor=(0,0,0,0),facecolor=colormap(norm2(data_norm[each_index[0],self.data_num.shape[0]-1-each_index[1]])),lw=self.linewidth)
                ax.add_patch(patch)
            else:
                pass
                #patch = patches.PathPatch(path, edgecolor=self.edgecolor,facecolor='w',lw=self.linewidth)
               # ax.add_patch(patch)
        self.add_axes_attr(fig=fig,ax=ax,forcolorbar=for_colorbar)
        ax.set_title(output,fontsize=20)
 
        if self.format == 'pdf':
            plot_img = output + '.tissue.pdf'
        elif self.format == 'png':
            plot_img = output + '.tissue.png'
        elif self.format == 'svg':
            plot_img = output + '.tissue.svg'
        elif self.format == 'eps':
            plot_img = output + '.tissue.eps'
        fig.savefig(plot_img,format=self.format,bbox_inches='tight')
        plt.cla()
        plt.close(fig)
    
    def project_cluster_to_spatial(self,data,color,output,which_cluster):
        """
        project all clusters to the spatial graph
        """    
        color_dict = dict(zip(color['group'],color['colour']))
        color_dict = collections.OrderedDict(sorted(color_dict.items()))

        all_verts,cell_index = self.Get_Path(self.data_num)
        fig, ax = plt.subplots(figsize=(self.figsize_w,self.figsize_h),frameon=False,dpi=self.dpi)    
        for each_vert,each_index in zip(all_verts,cell_index):
            path = Path(each_vert, self.codes)
            temp_value = self.data_num[each_index[0],self.data_num.shape[0]-each_index[1]-1]
            if not pd.isna(temp_value) and temp_value in which_cluster:
                patch = patches.PathPatch(path, edgecolor='w',facecolor=color_dict[temp_value],lw=2)
                ax.add_patch(patch)
            else:
                patch = patches.PathPatch(path, edgecolor='w',facecolor='w',lw=2)
                ax.add_patch(patch)

        self.add_axes_attr(fig=fig,ax=ax)
        plt.setp(ax.get_yticklabels(), ha="right")
        plt.setp(ax.get_xticklabels(), va="top")
        legend_elements = self.Build_Legend(color_dict,which_cluster=which_cluster)
        ax.legend(handles=legend_elements, loc='upper right',
                bbox_to_anchor=(1.2, 1),fontsize=20)
        if self.format == 'pdf':
            plot_img = output + '.spatial.pdf'
        elif self.format == 'png':
            plot_img = output + '.spatial.png'
        elif self.format == 'svg':
            plot_img = output + '.spatial.svg'
        elif self.format == 'eps':
            plot_img = output + '.spatial.eps'
        fig.savefig(plot_img,format=self.format,bbox_inches='tight')
        plt.cla()
        plt.close(fig)

    def project_cluster_to_tissue(self,data=[],color_list=[],output='',image='',which_cluster=''):
        
        img = cv2.imread(image)
        color_dict = dict(zip(color_list['group'],color_list['colour']))
        color_dict = collections.OrderedDict(sorted(color_dict.items()))      
        all_verts,cell_index = self.Get_Path(self.data_num)
        fig, ax = plt.subplots(figsize=(self.figsize_w,self.figsize_h),dpi=self.dpi,frameon=False)    
        ax.imshow(img,extent=(0,self.data_num.shape[0],0,self.data_num.shape[0]))
        
        for each_vert,each_index in zip(all_verts,cell_index):
            path = Path(each_vert, self.codes)
            temp_value = self.data_num[each_index[0],self.data_num.shape[0]-each_index[1]-1]
            if not pd.isna(temp_value) and temp_value in which_cluster:
                color_temp = list(matplotlib.colors.to_rgba(color_dict[temp_value]))
                color_temp[3] = self.alpha
                color_temp = tuple(color_temp)
                patch = patches.PathPatch(path, edgecolor=(0,0,0,0),facecolor=color_temp,lw=5)
                ax.add_patch(patch)
            else:
                pass
        
        self.add_axes_attr(fig=fig,ax=ax)
       
        plt.setp(ax.get_yticklabels(), ha="right")
        plt.setp(ax.get_xticklabels(), va="top")
        legend_elements = self.Build_Legend(color_dict,which_cluster=which_cluster)
        ax.legend(handles=legend_elements, loc='upper right',
                bbox_to_anchor=(1.2, 1),fontsize=20)
        if self.format == 'pdf':
            plot_img = output + '.tissue.pdf'
        elif self.format == 'png':
            plot_img = output + '.tissue.png'
        elif self.format == 'svg':
            plot_img = output + '.tissue.svg'
        elif self.format == 'eps':
            plot_img = output + '.tissue.eps'
        fig.savefig(plot_img,format=self.format,bbox_inches='tight')
        plt.cla()
        plt.close(fig)
