#!/usr/bin/env python3
"""
radix.py:
    Contains functions for converting between radix 10 and other radixes.
    Only considers integers as floats are not gaurenteed to have finite
    forms.
"""

import numpy as np


def base10_to_baseB(n, B, dtype=int):
    """
    Description:
        Takes an integer in radix 10 and converts it into it's radix B form.
        
    Arguments:
        n: number in radix 10
        B: radix to convert to
        
    Returns:
        List of digits of the number in radix B: n = ([d1, d2, d3])_B 
    """
    # Type check
    if type(B) is not int:
        print("ERROR: radix must be integer.")
        return
    try:
        n[0]
        n = np.asarray(n, dtype=dtype)
    except:
        n = np.asarray([n], dtype=dtype)
    if not issubclass(n.dtype.type, np.integer):
        print("ERROR: Can only transform radix for integers.")
        return
    number = [None,]*len(n)
    for j in range(len(n)):
        num = n[j]
        if num == 0:
            number[j] = [0]
            continue
        digits = []
        while num:
            digits.append(int(num%B))
            num //= B
        number[j] = digits[::-1]
    return number


def baseB_to_base10(n, B, dtype=int):
    """
    Description:
        Takes a number in radix B and converts it to radix 10.

    Arguments:
        n: list of digits of number in radix B
        B: radix being converted from

    Returns:
        Number in radix 10 (integer)
    """
    # Type check
    if type(B) is not int:
        print("ERROR: radix must be integer.")
        return
    try:
        n[0][0]
        n = np.asarray(n, dtype=dtype)
    except:
        n = np.asarray([n], dtype=dtype)
    if not issubclass(n.dtype.type, np.integer):
        print("ERROR: Can only transform radix for integers.")
        return
    number = np.zeros(len(n), dtype=dtype)
    for j in range(len(n)):
        for i, d in enumerate(n[j]):
            number[j] += d*B**(len(n[j])-i-1)
    return number


def baseB_to_str(n, B):
    """
    Description:
        Formats list of digits into a string of radix B

    Arguments:
        n: list of digits of number in radix B
        B: radix of number
        
    Returns:
        Number of radix B in string form: '(d1 d2 d3)_B'
    """
    # Type check
    if type(B) is not int:
        print("ERROR: radix must be integer.")
        return
    if type(n) is not list:
        print("ERROR: n must be list of digits.")
        return
    # Select seperator
    if B <= 10:
        seperator = ''
    else:
        seperator = ','
    # Build string
    string = ''
    for i in n:
        string += str(i)+seperator
    # Toss tailing seperator
    if seperator != '':
        string = string[:-1]
    # Return string, if radix 10 return standard form
    if B == 10:
        return rf"${string:s}$"
    else:
        return rf"$({string:s})_{{{B:d}}}$"