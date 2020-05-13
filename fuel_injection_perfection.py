# ================================================================
# Fuel Injection Perfection
# ================================================================
# Commander Lambda has asked for your help to refine the automatic quantum antimatter fuel injection system for her LAMBCHOP doomsday device. It's a great chance for you to get a closer look at the LAMBCHOP - and maybe sneak in a bit of sabotage while you're at it - so you took the job gladly.

# Quantum antimatter fuel comes in small pellets, which is convenient since the many moving parts of the LAMBCHOP each need to be fed fuel one pellet at a time. However, minions dump pellets in bulk into the fuel intake. You need to figure out the most efficient way to sort and shift the pellets down to a single pellet at a time.

# The fuel control mechanisms have three operations:

# Add one fuel pellet
# Remove one fuel pellet
# Divide the entire group of fuel pellets by 2 (due to the destructive energy released when a quantum antimatter pellet is cut in half, the safety controls will only allow this to happen if there is an even number of pellets)
# Write a function called answer(n) which takes a positive integer as a string and returns the minimum number of operations needed to transform the number of pellets to 1. The fuel intake control panel can only display a number up to 309 digits long, so there won't ever be more pellets than you can express in that many digits.

# For example:

# answer(4) returns 2: 4 -> 2 -> 1  
# answer(15) returns 5: 15 -> 16 -> 8 -> 4 -> 2 -> 1
# Languages
# To provide a Python solution, edit solution.py To provide a Java solution, edit solution.java

# Test cases
# Inputs:

# (string) n = "4"
# Output:

# (int) 2
# Inputs:

# (string) n = "15"
# Output:

# (int) 5
# ========================================================================================================

def solution(n):
    n = int(n)
    count_op=0
    while n is not 1:
        if not n&1:
            n>>=1
            
        else:
            greater_even=n+1
            smaller_even=n-1
            
            counta=countb=0
            
            while not greater_even&1:
                greater_even>>=1
                counta+=1
            
            while not smaller_even&1:
                smaller_even>>=1
                countb+=1
            
            if counta>countb and n is not 3:
                n+=1
            else:
                n-=1
        
        count_op+=1
    
    return count_op
