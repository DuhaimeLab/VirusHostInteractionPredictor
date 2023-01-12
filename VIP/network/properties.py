import numpy as np




'''
Nestedness calculations
'''

test_matrix = np.array([[1,0,1,1,1], 
                        [1,1,1,0,0],
                        [0,1,1,1,0],
                        [1,1,0,0,0]])


print(test_matrix.shape)


def nestedness(adj):

    pass



def rearrange(matrix):
    pass




def pairs(n):
    lst = []
    for i in range(0,n):
        for j in range(i+1,n):
            lst.append((i, j))
    return lst






'''
Modularity
'''



def modularity(adj):

    pass
