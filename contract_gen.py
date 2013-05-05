from __future__ import print_function
import os,sys,string
import numpy as np
from contract import contract
from Cheetah.Template import Template

#   Puropse:
#       Generate tensor contraction subroutine from the contract.tmpl template.
#       Cheetah Python Template Engine is used.
#   
#   Useage:
#       Set up parameters: ash, bsh, Ain, Bin, ABorder, and shape in the following.
#       Make sure Python and Cheetah are installed.
#       Make sure contract.tmpl file is present.
#       Save this file and run: pythoon contract_gen.py.
#       A new file(subroutine) named 'tensorcontractxxxx' will be generated automatically.
#   where xxxx stands for the initial rank of A and B tensor, the intermediate rank
#   of the tensor after contraction, and the final rank of the tensor.
#       One may want to copy the content of the subroutine to the generic tensor contract
#   module 'module_contract.f90' and add the subroutine name to the interface block.
#   By Yuzhi Liu on 04/28/2013. V1.
#   

# input parameters
ash = 5            # rank of tensor A
bsh = 3             # rank of tensor B
Ain = [1,4]         # contraction index in A
Bin = [1,2]         # contraction index in B
ABorder = [4,2,1,3] # the reorder sequence
#shape = [2,2]
shape = []          # the sequeeze index

# sub placeholders with ture values
rankA = ash
rankB = bsh
rankT = ash + bsh - len(Ain) - len(Bin)

if len(shape) > 0:
    rankoutT = len(shape)
else:
    rankoutT = rankT

dimA = ':'
for x in range(rankA-1):
    dimA = dimA + ',:'

dimB = ':'
for x in range(rankB-1):
    dimB = dimB + ',:'

dimoutT = ':'
for x in range(rankoutT-1):
    dimoutT = dimoutT + ',:'

dimAnew = ':'
for x in range(rankA-1):
    dimAnew = dimAnew + ',:'

dimBnew = ':'
for x in range(rankB-1):
    dimBnew = dimBnew + ',:'

dimT = ':'
for x in range(rankT-1):
    dimT = dimT + ',:'

dimSsize = rankoutT

Aother = list( set(range(1,ash+1)) - set(Ain) )
Bother = list( set(range(1,bsh+1)) - set(Bin) )


shapeAnew =''
for x in range(1,len(Aother)+1):
    shapeAnew = shapeAnew + '&\n                     SIZE(A,Aother(' + str(x) + ')),'
for x in range(1,len(Ain)+1):
    shapeAnew = shapeAnew + '&\n                     SIZE(A,ain(' + str(x) + ')),'
shapeAnew=shapeAnew[22:-1]


shapeBnew =''
for x in range(1,len(Bin)+1):
    shapeBnew = shapeBnew + '&\n                     SIZE(B,bin(' + str(x) + ')),'
for x in range(1,len(Bother)+1):
    shapeBnew = shapeBnew + '&\n                     SIZE(B,Bother(' + str(x) + ')),'
shapeBnew=shapeBnew[22:-1]

shapeT = ''
for x in range(1,ash+bsh-len(Ain)-len(Ain)+1):
    shapeT = shapeT +       '&\n                     Tsize(' + str(x) + '),'
shapeT=shapeT[22:-1]

if rankoutT == rankT:
    outTReshape = ' '
    outTNoReshape = 'outT = T'
else:
    outTReshape = 'outT = RESHAPE( T, Ssize )'
    outTNoReshape = ' '




names=[{'rankA':rankA,'rankB':rankB,'rankT':rankT,'rankoutT':rankoutT,
'dimA':dimA,'dimB':dimB,'dimoutT':dimoutT,'dimAnew':dimAnew,'dimBnew':dimBnew,
'dimT':dimT,'dimSsize':dimSsize,'shapeAnew':shapeAnew,'shapeBnew':shapeBnew,
'shapeT':shapeT,'outTReshape':outTReshape,'outTNoReshape':outTNoReshape}]

#print contract(file='contract.tmpl',searchList=names)
#print Template(file='contract.tmpl',searchList=names)

f1 = open('./tensorcontract%s%s%s%s.f90'%(rankA,rankB,rankT,rankoutT), "w")
print(Template(file='contract.tmpl',searchList=names), file = f1) 
