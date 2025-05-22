from typing import List
import dataclasses
import ast

l = "['a', 'b', 'c']"
ls = ast.literal_eval(l)
for i in range(len(ls)):
    print(ls[i])
# @dataclasses.dataclass
# class InputSmiles:
#     reactionSMILES: str

# @dataclasses.dataclass
# class Conditions:
#     Reagents: str
#     Solvents: str
#     Temperature: str


# @dataclasses.dataclass
# class Predictions:
#     predictions: List[Conditions]

# cond_1 = {'Reagents':'1', 'Solvents':'2', 'Temperature':'3'}
# cond_2 = {'Reagents':'4', 'Solvents':'5', 'Temperature':'6'}
# conditions = Conditions(Reagents=cond_1['Reagents'], Solvents=cond_1['Solvents'], Temperature=cond_1['Temperature'])
# conditions_1 = Conditions(Reagents=cond_2['Reagents'], Solvents=cond_2['Solvents'], Temperature=cond_2['Temperature'])
# print(conditions_1.Reagents)
# pred = Predictions([conditions, conditions_1])
# # for i in range(2):
# #     pred = Predictions([conditions])
# print(pred)