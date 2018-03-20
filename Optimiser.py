#!/usr/bin/env python
# This Python file uses the following encoding: utf-8
#
#
#
#    This file is used to read the dataList_New

import os         as         os
import numpy    as         np
import shutil as         sh
from     getopt import getopt
import sys        as         sys 
from scipy.optimize import minimize
from scipy.spatial import distance
from numpy.random import random
from scipy.interpolate import griddata
from scipy import integrate
import matplotlib.pyplot as plt 


# Golable index

GLOBAL_INDEX = 0
GLOBAL_DIR   = ""
GLOBAL_COUNT = 0



def load_data():

    ##############################
    # Generate data base list:

    DATA_LIST = []

    ROOT_DIR  = "/media/uqjqi/Janry_Research/Stator_Optimization/"
    if os.path.isfile("dataList"):
        print "The file \"dataList\" has been located correctly."


        f = open("dataList",'r')
        for line in f:
            if line.startswith("["):
                temp = line.replace(" ", '').replace('[','').replace(']','').replace('\r','').replace('\n','').split(",")
                DATA_LIST.append(temp)     

        # Convert the first letter in to float   
        for i in range(len(DATA_LIST)):
            for j in range(21):
                DATA_LIST[i][j] = float(DATA_LIST[i][j])
    else:
        print "WRONG! Please check the dataList is exist in this directory."


    # 将dataList中的已有的数据加载进来
    
    Index     = [] # List to store the index
    Status    = [] # List to store the Status
    Geometry  = [] # List to store the Geometry
    Post      = [] # List to store the results of post processing
    Cost      = [] # List to store the cost functions
    Directory = [] # List to store the case directory

    #[1, 1 , 0.1, 0.0, -0.1, 0.07, 0.0682, 0.0012, 0.002, 0.00008,  70., 0.2,  0.02,  nan, nan, nan, nan, nan, nan, nan, nan, folder address for simaultion ]

    for i in range(len(DATA_LIST)):
        temp_a = []
        temp_b = []
        temp_c = []
        Index.append(DATA_LIST[i][0])
        Status.append(DATA_LIST[i][1])

        #将Geometry信息读取到Geometry Lsit中
        for a in range(2,13):
            temp_a.append(DATA_LIST[i][a])
        #print temp_a
        Geometry.append(temp_a)
 
        #将后处理信息读取到Post Lsit中
        for b in range(13,17):
            temp_b.append(DATA_LIST[i][b])
        Post.append(temp_b)

        #将后COST信息读取到Post Lsit中
        for c in range(17,21):
            temp_c.append(DATA_LIST[i][c])
        Cost.append(temp_c)

        #将后Directory信息读取到Directory Lsit中
        Directory.append(DATA_LIST[i][21])  

        #print DATA_LIST

    return DATA_LIST, Index, Status, Geometry, Post, Cost, Directory 

##############################
def EXECUTE(Geometry_list,Case_directory):

    
    #############################
    #Prepare base case
    print "...Creating case folder"
    ROOT_DIR = "/media/uqjqi/Janry_Research/Stator_Optimization/"
    BASE_DIR = "/media/uqjqi/Janry_Research/Stator_Optimization/Basic_Setting/."

    CASE_DIR = Case_directory
    os.chdir(ROOT_DIR)

    print Case_directory, os.path.isfile(Case_directory)
    #判断一下是否已经完成计算
    if os.path.exists(Case_directory):
        print Case_directory, "is already exist"

    else:
        os.mkdir(Case_directory)
        os.system("cp -r "+ BASE_DIR + " " + CASE_DIR )


        #############################
        #create updated Stator_job.py
        


        print "...Creating Stator_job.py"
        job = open("Stator_job.py",'r')
        newjob = open(CASE_DIR + '/e3prep/Stator_job.py','w')
        for line in job:
            if "OP = [" in line:
                line = 'OP = ' + str(Geometry_list) + '\r\n'
            
            newjob.write("%s" % line)
        job.close()
        newjob.close()
        #print os.getcwd()
        os.chdir(CASE_DIR + '/e3prep/')
        #print os.getcwd()
        os.system('e3prep.py --job=Stator_job.py ')
        #os.system('e3post.py --job=Stator_job.py --vtk-xml')
        #os.system(' paraview plot/Stator_job.pvd')


        
        #############################
        # create OF mesh
        print os.getcwd()
        os.system('e3prepToFoam.py --job=Stator_job.py --version=fe')


        #############################
        # Set cyclic boundary conditions
        os.chdir(CASE_DIR)
        os.system("rm ./system/createPatchDict")
        os.system("mv ./system/createPatchDict.bak ./system/createPatchDict")
        os.system("createPatch -overwrite")
        os.system("rm *.obj")


        
        #############################
        # 这一步是将初始化的文档拷贝过去

        os.system("cp -r ../Initial_Value/. ./0/.")
        os.system("decomposePar")
        #os.system("paraFoam")

        # excute with 4 cores
        os.system("mpirun -np 4 transonicMRFDyMFoam -parallel")

        # reconstruct the partial files
        os.system("reconstructPar")

        # remove the processor files to save the memory
        os.system("rm -r processor*")


        # enhance the co-number and rerun the simulation
        #os.system("mv ./system/controlDict.bak ./system/controlDict")
        #os.system("decomposePar")
        #os.system("mpirun -np 4 transonicMRFDyMFoam -parallel")
        #os.system("reconstructPar")
        #os.system("rm -r processor*")



        #############################
        #将最后一步结果拷贝到Initial_value中，方便下一个算例的初始化
        os.system("cp -r ./15000/.  ../Initial_Value/." )
        os.system("rm -r ../Initial_value/uniform")
        



        # change the root directory
        os.chdir(ROOT_DIR)

    return

def POSTPROCESSOR(Case_directory):
    
    #Doing postprocess with processors of OpenFOAM
    CASE_DIR = Case_directory
    ROOT_DIR = "/media/uqjqi/Janry_Research/Stator_Optimization/"
    os.chdir(CASE_DIR)
    
    os.system("wallCellArea -latestTime > wallCellArea")
    os.system("wallCellVelocity -latestTime")

    #接下来这一部分是把T的进口边界条件改成ZeroGradient, 以用来使用Mach
    f_latest = open('wallCellArea', 'r')
    for line in f_latest:
        if 'Create mesh for time =' in line:
            timeString = line

    TimeList = timeString.replace('\n','').split(' ')
    latestTime = TimeList[5]
    print latestTime

    os.chdir(CASE_DIR+ '/' +latestTime)
    f_T = open('T','r')
    f_T_new= open('Tnew','w')
    for line in f_T:
        if 'isentropic' in line:
            line = "        type            zeroGradient;"
        f_T_new.write("%s" % line + '\r\n')
    f_T.close()
    f_T_new.close()
    os.system("mv T OLDT")
    os.system("mv Tnew T")
    os.chdir(CASE_DIR)


    #os.system("entropyIdeal -latestTime")
    os.system("Mach -latestTime")
    os.system("massFlowRate -latestTime > mdot")

    
    os.chdir(ROOT_DIR)

    return latestTime




def File_Reader(Filed_data_type,Field_name,Case_directory,time):


    scalar_list = []
    vector_list = []

    file_path = Case_directory + "/" + str(time) + "/" + str(Field_name)

    f = open(file_path,'r')

    line_count = 0; start_line = -1000

    if Field_name == 'p':
        for line in f:
            line_count += 1
            if "OF_outlet_00" in line:

                # This returns the starting line number 
                start_line = line_count

            if line_count == (start_line + 8):
                list_length = int(line.replace('\r','').replace('\n',''))
    else:
        for line in f:
            line_count += 1
            if "OF_outlet_00" in line:

                # This returns the starting line number 
                start_line = line_count

            if line_count == (start_line + 4):
                list_length = int(line.replace('\r','').replace('\n',''))        
    f.close()


    # read fields file and return a list
    if ("Scalar" == Filed_data_type) and (Field_name == 'p'):

        f_scalar = open(file_path,'r')
        line_count = 0
        for line in f_scalar:
            line_count += 1
            if (line_count >= (start_line + 10)) and (line_count <= (start_line + 10 + list_length -1)):
                scalar_list.append( float(line.replace('\r','').replace('\n','')))
        field_data_list = scalar_list
        #print field_data_list

    elif "Vector" == Filed_data_type:
        f_vector = open(file_path,'r')
        line_count = 0
        for line in f_vector:
            line_count += 1
            if (line_count >= (start_line + 6)) and (line_count <= (start_line + 6 + list_length -1)):
                temp = line.replace('(','').replace(')','').replace('\r','').replace('\n','').split(' ')
                for i in range(len(temp)):
                    temp[i] = float(temp[i])
                vector_list.append(temp)
        field_data_list = vector_list
        #print field_data_list
    elif ("Scalar" == Filed_data_type) and (Field_name != 'p'):
        f_scalar = open(file_path,'r')
        line_count = 0
        for line in f_scalar:
            line_count += 1
            if (line_count >= (start_line + 6)) and (line_count <= (start_line + 6 + list_length -1)):
                scalar_list.append( float(line.replace('\r','').replace('\n','')))
        field_data_list = scalar_list


    else:
        print "Please double check your "
    

    return  field_data_list

def Flux_Weighted_Average_Scalar(density_list,velocity_list,area_list,face_normal_list,scalar_list):
    
    # First determin the length of each list is equal or not:

    len1 = len(density_list)
    len2 = len(velocity_list)
    len3 = len(area_list)
    len4 = len(face_normal_list)
    len5 = len(scalar_list)

    numerator_flux_weighted     = 0
    denominator_flux_weighted   = 0

    if (len1 == len2) and (len2 == len3) and (len3 == len4) and (len4 == len5):
        for i in range(len1):
            # Calculate the U component on the surface normal direction for each cell face:
            cell_absolute_velocity_meridional_component = (float(velocity_list[i][0])*float(face_normal_list[i][0]) +\
                                                         float(velocity_list[i][1])*float(face_normal_list[i][1]) +\
                                                         float(velocity_list[i][2])*float(face_normal_list[i][2]))/ \
                                                         area_list[i]

            # Sum of cell_velocity_meridional_component times the cell area, that is the area weighted avarage
            numerator_flux_weighted    = numerator_flux_weighted + cell_absolute_velocity_meridional_component* float(scalar_list[i]) * area_list[i] * float(density_list[i])
            denominator_flux_weighted  = denominator_flux_weighted + abs(cell_absolute_velocity_meridional_component) * area_list[i] * float(density_list[i])

        flux_weighted_average_scalar = numerator_flux_weighted/denominator_flux_weighted

    return flux_weighted_average_scalar



def COST_EVALUATION(Mach_target, alpha_target, mdot_target, p0_target, Mach, alpha, mdot, p0):

    cost_list = []

    Mach_cost    =  ((Mach_target - Mach)/Mach_target)**2

    alpha_cost   =  ((alpha_target - alpha)/alpha_target)**2

    mdot_cost    =  ((mdot_target - mdot)/mdot_target)**2

    p0_cost      =  ((p0_target - p0)/p0_target)**2

    cost_list.append(Mach_cost)
    cost_list.append(alpha_cost)
    cost_list.append(mdot_cost)
    cost_list.append(p0_cost)

    return cost_list



def COMPARE_VECTOR(List0,List1,VT):
    # Use Mahala-Nobis method to evaluate the distance
    # To be remind that the two list shoud have the same length
    return  distance.mahalanobis(List0,List1,VT)




def Intergrand(t, mu, sigma_sqaure):
    # Cumulative distribution function for the normal distribution
    return np.exp(- (t - mu)**2. / (2.*sigma_sqaure) )/ np.sqrt(2.*np.pi*sigma_sqaure)




def Calculation(x):

 
    print "------------------------------------------------"
    print "Starting calculation, the current geometry is:"
    print x.tolist()
    print ""


    #加上一些flag，来控制之后是否运行新的case
    flag = 'None'
    redo_flag = 'None'

    #在这儿定义目标值
    Mach_target  = 0.8 
    alpha_target = 69.  # degree
    mdot_target  = 0.036
    p0_target    = 20.e6


    #定义weighting factor
    W_ma    = 20.
    W_alpha = 40.
    W_mdot  = 100.
    W_p0    = 40.





    # 将目标的 INDEX 数变成当前的路经值
    global GLOBAL_INDEX , GLOBAL_DIR , GLOBAL_COUNT 

    global ITERATION, RESIDUAL

    GLOBAL_DIR    = "/media/uqjqi/Janry_Research/Stator_Optimization/case_" + str(GLOBAL_INDEX)

    print "The currenty Index is: ", GLOBAL_INDEX
    

    # 对比当前的尺寸値的已经有列表值，如果已经存在，那么就用已经存的值替代当前值
    DATA_LIST, Index, Status, Geometry, Post, Cost, Directory  = load_data()


    # The inverse of covariance matrix
    # 斜方差矩阵
    global COV_MATRIX


    #这个地方使用不同的flag来标记最新的数据点是否根已有数据点的马式距离足够接近。如果足够接近的话，使用不同的策略进行计算

    print "Heading to distance evaluation now..."

    #这一步将现在的Geometry与之前所有的进行比较，把得到的马式距离放到一个list里
    distance_list = []
    for i in range(len(Geometry)):
        distance_list.append( COMPARE_VECTOR(x.tolist(),Geometry[i],COV_MATRIX.T) )#注意此处使用转置的斜方差矩阵





    #得到所有的马式距离之后，就要看是否有点落到相应的区间了

    print "After distance evaluation..."

    #如果马氏距离小于1e-07，直接使用临近点的Cost，或者说100%的进行插值计算
    #for j in range(len(distance_list)):
    #    if distance_list[j] < 1e-7:
    #        COST_FUNC = sum(Cost[i])   

    #找到最小的马式距离
    minimum_distance = min(distance_list)  

   #如果马氏距离在0到0.05之间，使用概率方程进行evaluation 
    if minimum_distance < 0.0005: 
        
        print "minimum distance is:", minimum_distance
        #这一步是让距离从0～0.005 整合为 0 ～0.1,以适应0～0.1的概率累计方程
        a = minimum_distance  * 200. 
        

        # cumulative distrubution function
        # we use the cumulative distrubution fuction for the normal distrubution.
        possibility_criteria =  (integrate.quad(Intergrand, -np.inf, a , args = (0.085 , 0.00005) ) )[0]

        # 开始扔硬币！
        b = random()
        print "The coin is:",b," and the possibility criteria is:", possibility_criteria

        if b >= possibility_criteria:
            print "-----Do the interpolation!"
            flag = 'interpolation'
        else:
            print "-----Run the case!"
            flag = 'evaluation'

    else:
        flag = 'evaluation'





    if 'interpolation' == flag:
            temp_Ma = []
            temp_P  = []
            temp_Alpha = []
            temp_Mdot  = []

            for i in range(len(Post)):
                temp_Ma.append((Post[i][0]))
                temp_Alpha.append((Post[i][1]))
                temp_Mdot.append((Post[i][2]))
                temp_P.append((Post[i][3]))
            #print temp_Ma

            #print Geometry

            points = np.array(Geometry)#; print points

            value_Ma = np.array(temp_Ma)#; #print values
            value_Alpha  = np.array(temp_Alpha)
            value_Mdot = np.array(temp_Mdot)
            value_P    = np.array(temp_P)

            xi  = x
            #print xi.shape


            Ma_interpolate    = griddata(Geometry, value_Ma, xi, method='linear',rescale=True)
            Alpha_interpolate = griddata(Geometry,value_Alpha,xi, method='linear',rescale=True)
            Mdot_interpolate  = griddata(Geometry,value_Mdot,xi,method='linear',rescale=True )
            P_interpolate     = griddata(Geometry,value_P,xi, method='linear',rescale=True)


            #print "!!!!!Ma=", Ma_interpolate,Alpha_interpolate,Mdot_interpolate,P_interpolate
            Cost_list = []
            Cost_list = COST_EVALUATION(Mach_target, alpha_target, mdot_target, p0_target, Ma_interpolate, Alpha_interpolate, Mdot_interpolate, P_interpolate)

            
      

            COST_FUNC = (Cost_list[0])* W_ma + (Cost_list[1]) * W_alpha + (Cost_list[2])*W_mdot + (Cost_list[3])*W_p0

            #如果插值在hall之外，那么就把flag变成evaluation

            if np.isnan( COST_FUNC ):
                print "TRY to extrapolation!!!!!! rerun the case."
                flag = 'evaluation'#'interpolation'
            else:
                print "The total cost is:", COST_FUNC, " and cost is, Ma:", Cost_list[0], " alpha:", Cost_list[1] , " mdot:", Cost_list[2], " p0:", Cost_list[3]
                GLOBAL_INDEX += 1
                ITERATION.append(GLOBAL_INDEX)
                RESIDUAL.append(COST_FUNC)
           


    print "flag is :",flag

    if ('evaluation' == flag):#  or ('evaluation' == redo_flag ):
        print "flag:", flag

        # execute the simulation
        EXECUTE(x.tolist(),GLOBAL_DIR)


        print "Heading here for evaluation"

        # do post processing
        latestTime = POSTPROCESSOR(GLOBAL_DIR)


        # Do post processing evalutation            
        density_list     = File_Reader("Scalar","rho",GLOBAL_DIR,latestTime)
        velocity_list    = File_Reader("Vector","wallCellVelocity",GLOBAL_DIR,latestTime)
        area_list        = File_Reader("Scalar","wallCellArea",GLOBAL_DIR,latestTime)
        face_normal_list = File_Reader("Vector","wallSurfaceNormal",GLOBAL_DIR,latestTime)
        mass_list        = File_Reader("Scalar","phi",GLOBAL_DIR,latestTime)
        Mach_list        = File_Reader("Scalar","Ma",GLOBAL_DIR,latestTime)
        pressure_list    = File_Reader("Scalar","p",GLOBAL_DIR,latestTime)

        alpha_list = []
        for j in range(len(velocity_list)):
            numerator   = float(velocity_list[j][0])*float(face_normal_list[j][0]) +\
                          float(velocity_list[j][1])*float(face_normal_list[j][1]) +\
                          float(velocity_list[j][2])*float(face_normal_list[j][2])
            denominator = np.sqrt(float(velocity_list[j][0])**2 + float(velocity_list[j][1])**2 + float(velocity_list[j][2])**2) *\
                          np.sqrt(float(face_normal_list[j][0])**2 + float(face_normal_list[j][1])**2 + float(face_normal_list[j][2])**2)
            alpha_list.append(np.degrees(np.arccos(numerator/denominator))) # change radians to degree



        total_pressure_list = []
        for k in range(len(velocity_list)):
            total_pressure_list.append(pressure_list[k] + 0.5*density_list[k]*((velocity_list[k][0])**2 + (velocity_list[k][1])**2 +(velocity_list[k][2])**2))



        average_Mach     = Flux_Weighted_Average_Scalar(density_list,velocity_list,area_list,face_normal_list,Mach_list)
        average_alpha    = Flux_Weighted_Average_Scalar(density_list,velocity_list,area_list,face_normal_list,alpha_list)
        mass_outlet      = sum(mass_list)
        average_total_pressure = Flux_Weighted_Average_Scalar(density_list,velocity_list,area_list,face_normal_list,total_pressure_list) # set it here



        Post_list = []
        Post_list.append(average_Mach)
        Post_list.append(average_alpha)
        Post_list.append(mass_outlet)
        Post_list.append(average_total_pressure)


        print Post_list


        Cost_list = []
        Cost_list = COST_EVALUATION(Mach_target, alpha_target, mdot_target, p0_target, average_Mach, average_alpha, mass_outlet, average_total_pressure)

        # Updating the cost function
        #Cost[i] = Cost_list

        print Cost_list

        print("index = ", GLOBAL_INDEX)
        print("dir =", GLOBAL_DIR) 
        print("count = ", GLOBAL_COUNT)


        write = [GLOBAL_INDEX] + [1] + x.tolist() + Post_list + Cost_list + [GLOBAL_DIR]
        print write



        # 写文件，算一步写一步

        f = open("dataList_New",'a+')
        f.write("%s" % write + '\r\n')
        f.close()


        COST_FUNC = float(Cost_list[0])* W_ma + float(Cost_list[1]) * W_alpha + float(Cost_list[2])*W_mdot + float(Cost_list[3])*W_p0
        
        GLOBAL_INDEX += 1

        ITERATION.append(GLOBAL_INDEX)
        RESIDUAL.append(COST_FUNC)

    else:
        print " you use wrong entry for the flag"
        pass

    #print "The cost is:", COST_FUNC
    return COST_FUNC


###
if __name__ == "__main__":


    global ITERATION, RESIDUAL
    ITERATION = []
    RESIDUAL  = []

    initial_array_type = 'latest'

    #ROOT_DIR  = "/media/uqjqi/Janry_Research/Stator_Optimization"

    DATA_LIST, Index, Status, Geometry, Post, Cost, Directory  = load_data()


    # Find the starting point, and give it an index
    global GLOBAL_INDEX 
    GLOBAL_INDEX = int(Index[-1] +1)



    first_Geom = Geometry[-1]
    x0 = first_Geom

    # 这里定义初始化的矩阵，要求是(N+1 , N)
    # 使用最新的list做为初始化，或者加一个选项，是要用最新的还是用老的

    if initial_array_type == 'latest':

        #选取最新的12个尺寸列表，用来初始化optimiser
        initial_array = np.array(Geometry[-12:])
        print initial_array

    else:
        initial_array = np.array(
            [[ 0.1, 0.0, -0.1, 0.05, 0.0682, 0.0012, 0.002, 8e-05, 70.0, 0.2, 0.05],
            [ 0.1, -0.02, -0.12, 0.07, 0.0682, 0.0012, 0.002, 8e-05, 70.0, 0.2, 0.02],
            [ 0.1, 0.0, -0.1, 0.06, 0.0682, 0.00122, 0.002, 8e-05, 70.0, 0.2, 0.02],
            [ 0.14, 0.0, -0.1, 0.07, 0.0682, 0.0012, 0.002, 8e-05, 70.0, 0.2, 0.02 ],
            [ 0.1, -0.03, -0.1, 0.07, 0.0682, 0.0012, 0.002, 8e-05, 70.0, 0.2, 0.02 ],
            [ 0.1, 0.0, -0.14, 0.07, 0.0682, 0.0012, 0.002, 8e-05, 70.0, 0.2, 0.02] ,
            [ 0.1, 0.0, -0.1, 0.08, 0.0682, 0.0012, 0.002, 8e-05, 70.0, 0.2, 0.02 ],
            [ 0.1, 0.0, -0.1, 0.07, 0.0682, 0.0012, 0.002, 8e-05, 69.8, 0.2, 0.02],
            [ 0.1, 0.0, -0.1, 0.07, 0.0682, 0.0012, 0.002, 8e-05, 70.0, 0.2, 0.022 ],
            [ 0.08, 0.0, -0.1, 0.07, 0.0682, 0.0012, 0.002, 8e-05, 70.0, 0.2, 0.02] ,
            [ 0.1, -0.01, -0.1, 0.07, 0.0682, 0.0012, 0.002, 8e-05, 70.0, 0.2, 0.02 ],
            [ 0.1, 0.0, -0.16, 0.07, 0.0682, 0.0012, 0.002, 8e-05, 70.0, 0.2, 0.02 ]] )

    global COV_MATRIX
    COV_MATRIX = np.cov(initial_array.T)

    #print COV_MATRIX

    #fun = lambda x: Himmelblau_func(x) # turns euqation(3-inputs) into fun(1-input)
    #fun = lambda x: Calculation(x)
    res =  minimize(Calculation,x0,method='Nelder-Mead',options={'initial_simplex': initial_array})
    print res
        # Goal of optimiser find x that results in min(fun(x))



    # create the final figure for the optimisation progress
    plt.figure(figsize=(8,4))
    plt.grid()
    plt.plot(ITERATION,RESIDUAL,'k-',markersize=9,markerfacecolor='w',markeredgewidth=1.5,linewidth=1.5)
    plt.legend(loc = 'best')
    plt.tight_layout()
    plt.show()
