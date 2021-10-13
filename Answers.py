'''
This module can be run to produce all of the results quoted in the write-up
and to perform the validation to ensure the code is working properly.
It can be run question by question or all in one go.
'''
import matplotlib.pyplot as plt

params = {
   'axes.labelsize': 12,
   'font.size': 15,
   'legend.fontsize': 10,
   'xtick.labelsize': 10,
   'ytick.labelsize': 10,
   'figure.figsize': [15, 10]
   } 
plt.rcParams.update(params)

#%%
'''
Q1
'''
print('Question 1 \n \n')
exec(open("Assignment_Q1.py").read())
print('\n \n')
#%%
'''
Q2
'''
print('Question 2 \n \n ')
exec(open("Assignment_Q2.py").read())
'''
The following code is here so the answers to question 2 are not reprinted when 
question 3 calls the matrix solve method for the cubic spline.
'''
print('Result for (a):\n',Res)
print('L:\n',L)
print('U:\n',U)
print('Det:\n', Det)
print('x:\n',x)
print('Inv:\n',Inv)
print('\n \n')
#%%
'''
Q3
'''
print('Question 3 \n \n ')
exec(open("Assignment_Q3.py").read())
print('\n \n')
#%%
'''
Q4
'''
print('Question 4 \n \n')
exec(open("Assignment_Q4.py").read())
print('\n \n')
#%%
'''
Q5
'''
print('Question 5 \n \n')
exec(open("Assignment_Q5.py").read())
print('\n \n')