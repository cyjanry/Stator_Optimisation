#!/usr/bin/env python
# This Python file uses the following encoding: utf-8
#
#   Script for automatic stator optimisation

import os     as     os
import numpy  as     np
import shutil as     sh
from   getopt import getopt
import sys    as     sys 



##############################
# Set the limitation for the optimisation parameters

# Preprocessing:
# Run Baseline simulation
# Run postprocessing on Baseline Simulation
# Write Baseline simulation geometry definition (parameter list) and postprocessign result to DataList

# Input:
# Folder with completed simualtion - Baseline Case
# List simualtion parameters, Index 0 refers to Baseline Case.  

# Step 1a:
# Manually pre-populate list with N new geometries
# Step 1b:
# Automatically expand list by M new geometries

# Step 2:
# From DataList select next incomplete geometry.
# Generate geomtry
# initialise from Baseline or appropriate previous converged result 

# Step 3:
# run simualtion until converged or abort.

# Step 4:
# do postprocessing
# write simualtion outcomes + postprocessign results back to DataList

# Return to Step 2 or 1b:

# Optionally perform Steps 2-4 in parallel threads.
##############################
def EXECUTE(Geometry_list,Case_directory):

    
    #############################
    #Prepare base case
    print "...Creating case folder"
    ROOT_DIR = "/home/uqjqi/Desktop/work2/2D_stage/Stator_Optimization/"
    BASE_DIR = "/home/uqjqi/Desktop/work2/2D_stage/Stator_Optimization/Basic_Setting/."
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
        # Copy the initial value to the case folder

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
        

        # change the root directory
        os.chdir(ROOT_DIR)


        '''
        ############################
        # Do postprocess and return cost function
        os.system("entropyIdeal -latestTime > entropy.log")


        # list if output parameter
        # re-openfile
        # write results to line defined by index.
        # close file. 
        '''
    return




def POSTPROCESSOR(Case_directory):
    
    #Doing postprocess with processors of OpenFOAM
    CASE_DIR = Case_directory
    ROOT_DIR = "/home/uqjqi/Desktop/work2/2D_stage/Stator_Optimization/"
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


#def COMPARE_VECTOR(List0,List1):
    # To be remind that the two list shoud have the same length

#    temp = []
#    for i in range(len(List0)):
#        temp.append((List0[i]- List1[i])**2)
#    
#    distance = np.sum(temp)
#    return distance




if __name__ == "__main__":

    # 定义要达到的最终目标
    Mach_target  = 0.8 
    alpha_target = 69.  # degree
    mdot_target  = 0.036
    p0_target    = 20.e6



    # Running case from datalist or use optimiser

    case_flag  = "file"  # or "optimiser"


   

    ##############################
    # Generate data base list:

    DATA_LIST = []

    ROOT_DIR  = "/home/uqjqi/Desktop/work2/2D_stage/Stator_Optimization/"
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


    # 从文档中执行运行

    if    "file" ==  case_flag:
        for i in range(len(DATA_LIST)):
            
            # 确保真的完成了，如果计算完成但是由于各种原因没能保存，再重新计算
            if (0 == Status[i]) :
                print "Simulation index =", int(Index[i]),"folder address:", Directory[i], "is successfully completed"


            # make sure that the 
            elif (1 == Status[i]):
                print "CASE index=", int(Index[i]), " now starting job"
                Directory[i] = "/home/uqjqi/Desktop/work2/2D_stage/Stator_Optimization/case_"+ str(i)
                
                # Create the mesh, run the simulation
                EXECUTE(Geometry[i],Directory[i])

                latestTime = POSTPROCESSOR(Directory[i])

                # Do post processing evalutation            
                density_list     = File_Reader("Scalar","rho",Directory[i],latestTime)
                velocity_list    = File_Reader("Vector","wallCellVelocity",Directory[i],latestTime)
                area_list        = File_Reader("Scalar","wallCellArea",Directory[i],latestTime)
                face_normal_list = File_Reader("Vector","wallSurfaceNormal",Directory[i],latestTime)
                mass_list        = File_Reader("Scalar","phi",Directory[i],latestTime)
                Mach_list        = File_Reader("Scalar","Ma",Directory[i],latestTime)
                pressure_list    = File_Reader("Scalar","p",Directory[i],latestTime)
           
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
                


                Post[i] = Post_list

                print("For case:", Index[i], "the average Mach number is:", average_Mach)
                print("                       the average alpha angel is:", average_alpha)
                print("                       the outlet mass flow rate is:", mass_outlet)
                print("                       the average total pressure is:", average_total_pressure)


                Cost_list = []
                Cost_list = COST_EVALUATION(Mach_target, alpha_target, mdot_target, p0_target, average_Mach, average_alpha, mass_outlet, average_total_pressure)

                # Updating the cost function
                Cost[i] = Cost_list


                # This step is to change the Status to 0, show that the case is successfully completed
                Status[i] = 0   
                # Rewrtie "nan" to average_entropy
                #Entropy[i] = average_entropy


        #print Post, Cost
        for i in range(len(DATA_LIST)):
            #print DATA_LIST[i],Cost[i]
            DATA_LIST[i][0]  = int(DATA_LIST[i][0])
            DATA_LIST[i][1]  = Status[i]
            DATA_LIST[i][13] = Post[i][0]
            DATA_LIST[i][14] = Post[i][1]
            DATA_LIST[i][15] = Post[i][2]
            DATA_LIST[i][16] = Post[i][3]
            DATA_LIST[i][17] = Cost[i][0]
            DATA_LIST[i][18] = Cost[i][1]
            DATA_LIST[i][19] = Cost[i][2]
            #print "----",DATA_LIST[i][20], Cost[i][2]
            DATA_LIST[i][20] = Cost[i][3]
            DATA_LIST[i][21] = Directory[i]

        # The last step is re-write the dataList

        f_old = open("dataList",'r')
        f_new = open("dataList_New",'w')

        i = 0
        for line in f_old:
            if line.startswith("["):
                line = DATA_LIST[i]
                i += 1
            f_new.write("%s" % line + '\r\n')

        f_old.close()
        f_new.close()


    
    #elif "optimiser" == case_flag:
    #    print "haha I'm heading here:", Geometry

    #    Init_Geom = Geometry[-1]
    #    print last_Geom


    #    res = minimize(fun, Init_Geom, method='nelder-mead',options={'xtol': 1e-8, 'disp': True})




    #    job = open("Stator_job.py",'r')
    #    newjob = open(CASE_DIR + '/e3prep/Stator_job.py','w')
    #    for line in job:
    #        if "OP = [" in line:
    #            line = 'OP = ' + str(Geometry_list) + '\r\n'
    #        
    #        newjob.write("%s" % line)
    #    job.close()
    #    newjob.close()





