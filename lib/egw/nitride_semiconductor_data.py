import pandas as pd

"""
Propaties1 
from Panalytical Datacollector

GaN |a,c,ν,Psp,e31,e33,C13,C33
AlN |
InN |
"""
columns1 = ['a', 'c', 'v', 'Psp', 'e31', 'e33', 'C13', 'C33']
properties1 = [[3.1893, 5.1851, 0.203, -0.029, -0.49, 0.73, 100, 392],
               [3.1130, 4.9816, 0.225, -0.081, -0.60, 1.46, 127, 382],
               [3.5380, 5.7020, 0.291, -0.032, -0.57, 0.97, 94, 200]]

induces = ['GaN', 'AlN', 'InN']
properties = pd.DataFrame(properties1, index=induces, columns=columns1)
