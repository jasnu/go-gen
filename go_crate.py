#!/usr/bin/python 

# Script:  go.py
# Purpose: create a random graphene oxide structure
# Example: python test.py
# Author:  Javier Rojas (USACH)

######################################################################################################################################################################################
# Para esta rutina necesito establecer una serie de procesos para la generación del oxido de grafeno (GO por sus siglas en inglés).                                                  #
# Esta se comone de poder generar las posiciones de una lámina de grafeno con unas dimensiones específicas, una cantidad definida de impurezas y una cantidad definida de vacancias: #
#                                                                                                                                                                                    #
# Generar grafeno en un cuadrado con dimensiones específicas                                                                                                                         #
# Agregar impurezas de grupos epoxi (Oxigenos entre carbonos de la red) y hidroxi (OH's aderidos a un carbono)                                                                       #
# Poder generar vacancias en la red con diversos tamaños y concentraciones                                                                                                           #
######################################################################################################################################################################################

import numpy as np
import random as rnd

#--Default distances
# aCC=1.42 #Carbon-Carbon distance

#--Ouput file options
def wr_atom(out,ele="dummy",pos=(0,0,0),first=False,edge=False):
  # Escritura de archivo .xyz para la visualización
  # Escribe una linea de un archivo xyz o el inicio del mismo
  if first:
    out.write('%i\nLattice="%.4f 0.0 0.0 0.0 %.4f 0.0 0.0 0.0 %.4f" Properties=species:S:1:pos:R:3\n' % (wr_atom.natoms,wr_atom.border[0],wr_atom.border[1],wr_atom.border[2]))
    return()
  
  if edge:
    final=-wr_atom.border*np.floor_divide(pos,wr_atom.border)*np.asarray([1,1,0])+np.fmod(pos,wr_atom.border)
  else:
    final=pos
  out.write("%3s %13.6f %13.6f %13.6f\n" % (ele,final[0],final[1],final[2]))

def wr_atom(out,ele="dummy",pos=(0,0,0),first=False,edge=False):
  # Escritura de un archivo .data para la entrada a LAMMPS
  if first: #Escribe la linea de descripción del .data
    out.write('# GO generado\n%i atoms\n3 atom types\n0 %.2f xlo xhi\n0 %.2f ylo yhi\n-15 15 zlo zhi\n\nAtoms #Atomic\n\n' % (wr_atom.natoms,wr_atom.border[0],wr_atom.border[1]))
    return()

  wr_atom.counter+=1
  if edge:
    final=-wr_atom.border*np.floor_divide(pos,wr_atom.border)*np.asarray([1,1,0])+np.fmod(pos,wr_atom.border)
  else:
    final=pos
    
  if ele=="C":
    out.write("%8i 1 %13.6f %13.6f %13.6f\n" % (wr_atom.counter,final[0],final[1],final[2]))
    return()
  if ele=="O":
    out.write("%8i 2 %13.6f %13.6f %13.6f\n" % (wr_atom.counter,final[0],final[1],final[2]))
    return()
  if ele=="H":
    out.write("%8i 3 %13.6f %13.6f %13.6f\n" % (wr_atom.counter,final[0],final[1],final[2]))
  
wr_atom.natoms=0
wr_atom.counter=0
wr_atom.border=np.array((1,1,1))
  
#--Function section
def inner_form(L,typ=(0,0),watch=False):
  # Escribe una serie de atomos dependiendo de que estructura sale,
  # Partiendo de su posición en la grid (L).
  # Hay dos columnas en typ que definen dos parte de las impurezas
  # La primera define los grupos epoxy (hasta 1 por celda).
  # La sedunda define los grupos hidroxi del sistema (hasta 2 por
  # celda).
  
  wr_atom(inner_form.wo,'C',L,edge=watch) #No hay vacancias implementadas todavía
  wr_atom(inner_form.wo,'C',L+inner_form.a,edge=watch)
  
  if np.linalg.norm(typ)==0: #Aqui es donde se agregan las impurezas
    return()

  # Poniendo epoxys
  if np.abs(typ[0])==1:
    wr_atom(inner_form.wo,'O',L+inner_form.a*np.array((-0.5,0.5,0))+np.array((0,0,np.sign(typ[0])*inner_form.poxz)),edge=watch)
  
  if np.abs(typ[0])==2:
    wr_atom(inner_form.wo,'O',L+np.array((0,-0.5*inner_form.absa,np.sign(typ[0])*inner_form.poxz)),edge=watch)

  if np.abs(typ[0])==3:
    wr_atom(inner_form.wo,'O',L+inner_form.a*np.array((0.5,0.5,0))+np.array((0,0,np.sign(typ[0])*inner_form.poxz)),edge=watch)

  # Poniendo hydroxsis
  if np.abs(typ[1])==1:
    wr_atom(inner_form.wo,'O',L+np.sign(typ[1])*np.array((0,0,inner_form.ohz)),edge=watch)
    wr_atom(inner_form.wo,'H',L+np.sign(typ[1])*np.array((0,0,inner_form.hhz)),edge=watch)
    return()

  if np.abs(typ[1])==2:
    wr_atom(inner_form.wo,'O',L+inner_form.a+np.sign(typ[1])*np.array((0,0,inner_form.ohz)),edge=watch)
    wr_atom(inner_form.wo,'H',L+inner_form.a+np.sign(typ[1])*np.array((0,0,inner_form.hhz)),edge=watch)
    return()

  if np.abs(typ[1])==3:
    wr_atom(inner_form.wo,'O',L+np.sign(typ[1])*np.array((0,0,inner_form.ohz)),edge=watch)
    wr_atom(inner_form.wo,'H',L+np.sign(typ[1])*np.array((0,0,inner_form.hhz)),edge=watch)
    wr_atom(inner_form.wo,'O',L+inner_form.a+np.sign(typ[1])*np.array((0,0,inner_form.ohz)),edge=watch)
    wr_atom(inner_form.wo,'H',L+inner_form.a+np.sign(typ[1])*np.array((0,0,inner_form.hhz)),edge=watch)
    return()

  if np.abs(typ[1])==4:
    wr_atom(inner_form.wo,'O',L-np.sign(typ[1])*np.array((0,0,inner_form.ohz)),edge=watch)
    wr_atom(inner_form.wo,'H',L-np.sign(typ[1])*np.array((0,0,inner_form.hhz)),edge=watch)
    wr_atom(inner_form.wo,'O',L+inner_form.a+np.sign(typ[1])*np.array((0,0,inner_form.ohz)),edge=watch)
    wr_atom(inner_form.wo,'H',L+inner_form.a+np.sign(typ[1])*np.array((0,0,inner_form.hhz)),edge=watch)
    return()
  

def tk_grid(name,grid):
  with open(name,'w') as o:
    shp=grid.shape
    inner_form.wo=o
    wr_atom.natoms=tk_grid.atoms
    wr_atom.border=np.array((2.46*shp[0],2.13*shp[1],10))
    wr_atom(o,first=True)

    #Bottom edge first
    for i in range(shp[0]):
      L=i*tk_grid.L1
      inner_form(L,grid[i,0,:],watch=True)
    
    for j in range(1,shp[1]):
      Lj=j*tk_grid.L2-np.ceil(j/2)*tk_grid.L1

      #Carefull with that edge
      inner_form(Lj,grid[0,j,:],watch=True)
      
      #Center atoms are always fine
      for i in range(1,shp[0]):
        L=i*tk_grid.L1+Lj
        inner_form(L,grid[i,j,:])
        
#Using tipical C-C bond 1.42A
inner_form.poxz=1.4
inner_form.ohz=1.5
inner_form.hhz=inner_form.ohz+1.11
inner_form.a=np.array((1.23,0.71,0))
inner_form.absa=np.linalg.norm(inner_form.a)
tk_grid.L1=np.array((2.46,0,0))
tk_grid.L2=np.array((1.23,2.13,0))

def put_epoxy(grid):
  #Select an availible random site from a grid array
  for i in range(put_epoxy.stop):
    xsel=rnd.randrange(grid.shape[0])
    ysel=rnd.randrange(grid.shape[1])
    if grid[xsel,ysel]!=0: continue

    site=rnd.randrange(3)

    if ysel%2==1:
      xp=(xsel+1)%grid.shape[0]
      yp=(ysel+1)%grid.shape[1]
      if site==2 and np.abs(grid[xsel,yp])!=2 and np.abs(grid[xp,ysel])!=1:return((xsel,ysel),3)
      xl=(xsel-1)%grid.shape[0]
      yl=(ysel-1)%grid.shape[1]
      if site==1 and np.abs(grid[xl,yl])!=3 and np.abs(grid[xsel,yl])!=1: return((xsel,ysel),2)
      if site==0 and np.abs(grid[xl,ysel])!=3 and np.abs(grid[xl,yp])!=2: return((xsel,ysel),1)
    else:
      xp=(xsel+1)%grid.shape[0]
      yp=(ysel+1)%grid.shape[1]
      if site==2 and np.abs(grid[xp,yp])!=2 and np.abs(grid[xp,ysel])!=1:return((xsel,ysel),3)
      xl=(xsel-1)%grid.shape[0]
      yl=(ysel-1)%grid.shape[1]
      if site==1 and np.abs(grid[xsel,yl])!=3 and np.abs(grid[xp,yl])!=1: return((xsel,ysel),2)
      if site==0 and np.abs(grid[xl,ysel])!=3 and np.abs(grid[xsel,yp])!=2: return((xsel,ysel),1)
      
    
  print('ERROR: Too many attemps in put epoxy. NO OUTPUT!')
  exit(1)
put_epoxy.stop=1000

def put_hydroxi(grid):
  for i in range(put_hydroxi.stop):
    xsel=rnd.randrange(grid.shape[0])
    ysel=rnd.randrange(grid.shape[1])    
    site=rnd.randrange(2)

    if np.abs(grid[xsel,ysel])==3: continue #Not a valid place
    
    #Check site
    if np.abs(put_hydroxi.conf[xsel,ysel])==site+1 or np.abs(put_hydroxi.conf[xsel,ysel])==4 or np.abs(put_hydroxi.conf[xsel,ysel])==3: continue #Already occupied
    if site==0: #Check for position 1 (site=0)
      if grid[xsel,ysel]!=0: continue
      return((xsel,ysel),1*rnd.randrange(-1,3,2)) 

    xp=(xsel+1)%grid.shape[0]
    yp=(ysel+1)%grid.shape[1]
    if ysel%2==1: #Checking posibilites for position 2 (site 1)
      if np.abs(grid[xsel,yp])!=2 and np.abs(grid[xp,ysel])!=1: return((xsel,ysel),2*rnd.randrange(-1,3,2))
    else:  
      if np.abs(grid[xp,yp])!=2 and np.abs(grid[xp,ysel])!=1: return((xsel,ysel),2*rnd.randrange(-1,3,2))

  print('ERROR: Too many attemps in put hydroxi. NO OUTPUT!')
  exit(1)
put_hydroxi.stop=1000
    
def epoxificate(Nx,Ny,ratio):
  conf=np.zeros((Nx,Ny))
  add=int(np.around(Nx*Ny*2*ratio))
  for i in range(add):
    k,val=put_epoxy(conf)
    conf[k]=val*rnd.randrange(-1,3,2)
    # print('[%i,%i] = %i' % (k[0],k[1],val))

  return(conf,add)

def hydroxination(epox_conf,ratio):
  print(epox_conf)
  put_hydroxi.conf=np.zeros(epox_conf.shape)
  add=int(np.around(epox_conf.size*2*ratio))  
  for i in range(add):
    k,val=put_hydroxi(epox_conf)

    #Adding the new hydroxi
    if put_hydroxi.conf[k]==0:
      put_hydroxi.conf[k]=val
      continue
    if put_hydroxi.conf[k]==2:
      if val==1:
        put_hydroxi.conf[k]=3
        continue
      if val==-1:
        put_hydroxi.conf[k]=4
        continue
      print('ERROR: Coliding Hidrogens %i %i' % (put_hydroxi.conf[k],val))
    if put_hydroxi.conf[k]==-2:
      if val==1:
        put_hydroxi.conf[k]=-4
        continue
      if val==-1:
        put_hydroxi.conf[k]=-3
        continue
      print('ERROR: Coliding Hidrogens %i %i' % (put_hydroxi.conf[k],val))
      
    if put_hydroxi.conf[k]==1:
      if val==2:
        put_hydroxi.conf[k]=3
        continue
      if val==-2:
        put_hydroxi.conf[k]=-4
        continue
      print('ERROR: Coliding Hidrogens %i %i' % (put_hydroxi.conf[k],val))
    if put_hydroxi.conf[k]==-1:
      if val==2:
        put_hydroxi.conf[k]=4
        continue
      if val==-2:
        put_hydroxi.conf[k]=-3
        continue
      print('ERROR: Coliding Hidrogens %i %i' % (put_hydroxi.conf[k],val))

    print('ERROR: Configuration unexpected %i %i' % (put_hydroxi.conf[k],val))
      
  return(put_hydroxi.conf,add)
#-----------------
#-- Mini tests
#-----------------
#if __name__ == '__main__':
if True:
  rnd.seed()

  # N=np.zeros((2,2,2))
  # N[0,0,:]=np.array([0,-4])
  # tk_grid.atoms=N[:,:,0].size*2+4
  # tk_grid('example.xyz',N)
  # exit()
  N=np.zeros((82,94,2))
  N[:,:,0],added=epoxificate(N.shape[0],N.shape[1],0.34)
  tk_grid.atoms=N[:,:,0].size*2+added
  N[:,:,1],added=hydroxination(N[:,:,0],0.17)
  tk_grid.atoms+=added*2
  tk_grid('example.data',N)
  exit()
#-----------------
#-- Test script
#-----------------
if __name__ == '__main__':
  import sys
  if not globals().has_key("argv"): argv = sys.argv
  import time
  
  
