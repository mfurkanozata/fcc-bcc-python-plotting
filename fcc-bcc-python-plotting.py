# -*- coding: utf-8 -*-
"""
Created on Thu May 28 08:57:11 2020
data.extend([str(x) for x in line.split()])
@author: MUHAMMETFURKANOZATA
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import itertools
data=list() #Crete a list
with open("MetalDatabase.txt", "r") as input_file: #open the input file that the data will be read
 with open("Output2.txt", "w") as output_file:    #open the output file that the data will be written
    print('This program is a useful tool to define how many atoms are present in a given volume','',
          'Please Choose an atom','',sep='\n')
    for line in input_file:    # read the text line by line
         print(line)
         for x in line.split():   # read the text line string by string
          data.extend([str(x)])     #seperate the every single string into a list
         
    N=str(input("Atom : "))    #create a input variable to get the name of element 
    N2=N.capitalize()
    Dx=float(input("Dimension(x) : "))
    Dy=float(input("Dimension(y) : "))
    Dz=float(input("Dimension(z) : "))
    
    R=data[data.index(N2)+2]    # Every radii in the source file is the second string after the sign of atom
    
    
    S=" "     
    
    for i in S:              #continue the cycle as long as i is satisfied with conditions following
        if data.index('BCC')<data.index(N2)<data.index('FCC'):
             S="BCC"
             apf=0.68
             a=4*float(R)/(3**(1/2))
             
        elif data.index('FCC')<data.index(N2)<data.index('HCP'):
             S="FCC"
             apf=0.74
             a=4*float(R)/(2**(1/2))
        elif data.index('HCP')<data.index(N2):
             S="HCP"
             apf=0.74
             a=2*float(R)      
    Cx =list()
    Cy =list()
    Cz =list()
    Cxx=list()
    Cxxx=list()
    Cxxxx=list()
    Cyy =list()
    Cyyy=list()
    Cyyyy=list()
    Czz =list()
    Czzz=list() 
    Czzzz=list()   
         
    
    
   ###### COORDINATES ######
    if S=='FCC' or S=='BCC':
        
        ax=0
        ay=0
        az=0
        while float(ax) <= Dx:     
         Cx.append(ax)  
         ax=ax+a
        while float(ay) <= Dy:     
         Cy.append(ay)  
         ay=ay+a
        while float(az) <= Dz:     
         Cz.append(az)  
         az=az+a
       
        
        
        def main(Cxx,Cyy,Czz,bx,by,bz):
          while float(bx) <= Dx:
              Cxx.append(bx)
              bx+=a
          while float(by) <= Dy:     
              Cyy.append(by)  
              by=by+(a)
          while float(bz) <= Dz:     
              Czz.append(bz)  
              bz=bz+(a) 
      #### BCC ADDITIONAL COORDINATES ###  
        if S=='BCC':
          main(Cxx,Cyy,Czz,a/2,a/2,a/2)
      #### BCC ADDITIONAL COORDINATES ###
        
        #### FCC ADDITIONAL COORDINATES ###  
        if S=='FCC':
    
          main(Cxx,Cyy,Czz,0,a/2,a/2)
        
          main(Cxxx,Cyyy,Czzz,a/2,0,a/2)
                
          main(Cxxxx,Cyyyy,Czzzz,a/2,a/2,0)
   
   ### FCC ADDITIONALCOORDINATES #### 
   
   ### HCP ADDITIONALCOORDINATES ####
    else:
        c_hcp=a*1.63
        bx=0
        bx1=a/2
        bx2=a/2
        bx3=a
        by=(3**(1/2)*a)/2
        by1=0
        by2=(1.15)*a
        by3=(3**(1/2)*a)/6
        bz=0
        bz1=0
        bz2=c_hcp/2
        bz3=c_hcp/2
        
        while float(bx) <= Dx:
              Cx.append(bx)
              bx+=a
        while float(bx1) <= Dx:
              Cxx.append(bx1)
              bx1+=a
        while float(by) <= Dy:     
              Cy.append(by)  
              by=by+(3**(1/2)*a)
        while float(by1) <= Dy:     
              Cyy.append(by1)  
              by1=by1+(3**(1/2)*a)
        while float(bz) <= Dz:     
              Cz.append(bz)  
              bz=bz+(c_hcp)
        while float(bz1) <= Dz:     
              Czz.append(bz1)  
              bz1+=(c_hcp)
              
        while float(bx2) <= Dx:
              Cxxx.append(bx2)
              bx2+=a  
        while float(by2) <= Dy:     
              Cyyy.append(by2)  
              by2=by2+(3**(1/2)*a)
        while float(bz2) <= Dz:     
              Czzz.append(bz2)  
              bz2=bz2+(c_hcp)
              
        while float(bx3) <= Dx:
              Cxxxx.append(bx3)
              bx3+=2*a  
        while float(by3) <= Dy:     
              Cyyyy.append(by3)  
              by3=by3+(3**(1/2)*a)
        while float(bz3) <= Dz:     
              Czzzz.append(bz3)  
              bz3=bz3+(c_hcp)
  
    C=[Cx,Cy,Cz]
    CC=[Cxx,Cyy,Czz]
    CCC=[Cxxx,Cyyy,Czzz]
    CCCC=[Cxxxx,Cyyyy,Czzzz]   
    comb=list(itertools.product(*C))
    comb2=list(itertools.product(*CC))
    comb3=list(itertools.product(*CCC))
    comb4=list(itertools.product(*CCCC))
    list2 = []
    CombTOT=comb+comb2+comb3+comb4
    CombTOT.sort()
    for i in CombTOT:
      if i not in list2:
        list2.append(i)
    xc=[]
    
    yc=[]
    
    zc=[]
    
    for c in list2:
          aa=c[0]          
          xc.append(aa)    
          aa=c[1]          
          yc.append(aa)
          aa=c[2]          
          zc.append(aa)
    CaX1=np.array(xc, ndmin=2)
    CaY1=np.array(yc, ndmin=2)
    CaZ1=np.array(zc, ndmin=2)
    
    T=len(list2) 
    
    print(f"Atom : "+N2, "Atomic Radius : "+str(R), "Structure : "+str(S), "Dimensions : " 
              +'x='+str(Dx) +', y='+str(Dy) +', z='+str(Dz), "Total No. of Atoms : "
              +"{:.2f}".format(T),'', sep='\n', file=output_file)
    n=0
    print('{:<19}'.format('atom'),'{:<9}'.format('x'),'{:<9}'.format('y'),'{:<10}'.format('z'),file=output_file)
    while n<=(T-1):
        aaa=list2[n]
        x=aaa[0]
        xx=round(x,3)
        xxx=str(xx)
        
        y=aaa[1]
        yy=round(y,3)
        yyy=str(yy)
        
        z=aaa[2]
        zz=round(z,3)
        zzz=str(zz)
        
        print('{:<20}'.format(N2)+'{:<10}'.format(xxx)+'{:<10}'.format(yyy)+'{:<10}'.format(zzz),
          sep='\n', file=output_file)
        n+=1
    if S=='FCC' or S=='BCC':    
        def make_ax(grid=False):
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_zlabel("z")
            ax.grid(grid)
            return ax
        
        filled = np.array([
            [[0, 0, a], [0, a, 0]],
            [[0, a, a], [a, 0, 0], [a, 0, a]],
            [[a, a, 0], [a, a, a], [0, 0, 0]]
        ])
    
        
        ax = make_ax()
        ax.voxels(np.ones((int(Dx), int(Dy), int(Dz))), facecolors='white', edgecolors='green', shade=False)
        ax.scatter(CaX1,CaY1,CaZ1, c='r', marker='o')
    
        plt.show()

    
    
    
    input_file.close
    output_file.close