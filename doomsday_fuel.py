==================================================================
Doomsday Fuel
========================================================================================
Making fuel for the LAMBCHOP's reactor core is a tricky process because of the exotic matter involved. It starts as raw ore, then during processing, begins randomly changing between forms, eventually reaching a stable form. There may be multiple stable forms that a sample could ultimately reach, not all of which are useful as fuel.

Commander Lambda has tasked you to help the scientists increase fuel creation efficiency by predicting the end state of a given ore sample. You have carefully studied the different structures that the ore can take and which transitions it undergoes. It appears that, while random, the probability of each structure transforming is fixed. That is, each time the ore is in 1 state, it has the same probabilities of entering the next state (which might be the same state). You have recorded the observed transitions in a matrix. The others in the lab have hypothesized more exotic forms that the ore can become, but you haven't seen all of them.

Write a function answer(m) that takes an array of array of nonnegative ints representing how many times that state has gone to the next state and return an array of ints for each terminal state giving the exact probabilities of each terminal state, represented as the numerator for each state, then the denominator for all of them at the end and in simplest form. The matrix is at most 10 by 10. It is guaranteed that no matter which state the ore is in, there is a path from that state to a terminal state. That is, the processing will always eventually end in a stable state. The ore starts in state 0. The denominator will fit within a signed 32-bit integer during the calculation, as long as the fraction is simplified regularly.

For example, consider the matrix m:

[
  [0,1,0,0,0,1],  # s0, the initial state, goes to s1 and s5 with equal probability
  [4,0,0,3,2,0],  # s1 can become s0, s3, or s4, but with different probabilities
  [0,0,0,0,0,0],  # s2 is terminal, and unreachable (never observed in practice)
  [0,0,0,0,0,0],  # s3 is terminal
  [0,0,0,0,0,0],  # s4 is terminal
  [0,0,0,0,0,0],  # s5 is terminal
]
So, we can consider different paths to terminal states, such as:

s0 -> s1 -> s3
s0 -> s1 -> s0 -> s1 -> s0 -> s1 -> s4
s0 -> s1 -> s0 -> s5
Tracing the probabilities of each, we find that:

s2 has probability 0
s3 has probability 3/14
s4 has probability 1/7
s5 has probability 9/14
So, putting that together, and making a common denominator, gives an answer in the form of [s2.numerator, s3.numerator, s4.numerator, s5.numerator, denominator] which is [0, 3, 2, 9, 14].

Languages
To provide a Python solution, edit solution.py
To provide a Java solution, edit solution.java

Test cases
Inputs:

(int) m = [[0, 2, 1, 0, 0], [0, 0, 0, 3, 4], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]
Output:

(int list) [7, 6, 8, 21]
Inputs:

(int) m = [[0, 1, 0, 0, 0, 1], [4, 0, 0, 3, 2, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]
Output:

(int list) [0, 3, 2, 9, 14]
====================================================================================================

from fractions import Fraction

def matTranspose(m):
    return map(list,zip(*m))

def matMinor(m,i,j):
    return [row[:j] + row[j+1:] for row in (m[:i]+m[i+1:])]

def matDeterminant(m):
    if len(m) == 2:
        return m[0][0]*m[1][1]-m[0][1]*m[1][0]

    determinant = 0
    for c in range(len(m)):
        determinant += ((-1)**c)*m[0][c]*matDeterminant(matMinor(m,0,c))
    return determinant

def matInv(m):
    assert(len(m) == len(m[0]))
    determinant = matDeterminant(m)
    
#   2X2 matrice
    if len(m) == 2:
        return [[m[1][1]/determinant, -1*m[0][1]/determinant],
                [-1*m[1][0]/determinant, m[0][0]/determinant]]

    cofactors = []
    for r in range(len(m)):
        cofactorRow = []
        for c in range(len(m)):
            minor = matMinor(m,r,c)
            cofactorRow.append(((-1)**(r+c)) * matDeterminant(minor))
        
        cofactors.append(cofactorRow)
    cofactors = matTranspose(cofactors)
    
    for r in range(len(cofactors)):
        for c in range(len(cofactors)):
            cofactors[r][c] = cofactors[r][c]/determinant
            
    return cofactors

def findLCM(a, b):
    maxm = max(a, b) 
    minm = min(a, b) 
    i = maxm
    while(1) : 
        if (i % minm == 0): 
            return i 
        i += maxm

def matMul(a, b):
    # confirm dimensions
    assert(len(a[0]) == len(b))
    
    rows = len(a)
    cols = len(b[0])
    # create the result matrix c = a*b
    res = [[0]*cols for _ in range(rows)]

    for row in range(rows):
        for col in range(cols):
            dot_product = Fraction(0, 1)
            for i in range(len(b)):
                dot_product += a[row][i]*b[i][col]
            res[row][col] = dot_product
    return res

def solution(arr):
    
    zero_lst = list(len(arr[0])*[0])
    
    terminal = [row_index for row_index in range(len(arr)) if arr[row_index] == zero_lst]
    
    non_terminal = [row_index for row_index in range(len(arr)) if arr[row_index] != zero_lst]
    
    if len(terminal) == 1:
        return [1, 1]
    
    lst = terminal + non_terminal

# convering this to Absorbing Markov Chain matrice Format.
#  ________
# |_I_|_O_|
# |_R_|_Q_|

    res = []
    for a in non_terminal:
        temp=[]
        count = 0
        for x in lst:
            count += arr[a][x]
            temp.append(arr[a][x])
        
        for i in range(len(temp)):
            temp[i] = Fraction(temp[i], count)
            
        res.append(temp)

#   Calculating R, Q matrices
    R,Q = [],[]
    for entry in res:
        R.append(entry[:len(terminal)])
        Q.append(entry[len(terminal):])
        
#   Calculating I-Q matrice
    I_Q = []
    entry_with_one = 0
    for row in range(len(Q)):
        temp = []
        for col in range(len(Q[0])):
            if col == entry_with_one:
                temp.append(1 - Q[row][col])
            else:
                temp.append(-Q[row][col])

        I_Q.append(temp)
        entry_with_one += 1
    
#   Calculating Inverse of (I-Q)
    F = matInv(I_Q)
    
#   F*R
    FR = matMul(F, R)
    
    answer = []
    
    lcm = 1
    for row1 in FR[0]:
        lcm = findLCM(lcm, row1.denominator)
    
    answer = [int(row.numerator*(lcm/row.denominator)) for row in FR[0]]
    
    answer.append(lcm)
    
    return answer
