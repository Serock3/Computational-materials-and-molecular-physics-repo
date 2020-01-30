#%%
import numpy as np

alpha = (0.298073,1.242567,5.782948,38.474970)

C = [1,1,1,1]

def normalize(asd):
    return [n/np.sqrt(np.dot(asd,asd)) for n in asd]

C = normalize(C)

print("hej")
# %%
