# -------------------------------------Functions Implementations -------------------------------------

'''
Function: InitializeMatrix
Description: Initialize the Matrix with zeros on all cells 
             thus, keeping Main diagonal and Bottom side of zeros
             and filling the rest.
Input: length of sequence
Returns: Matrix of zeros
'''
def InitializeMatrix(ln):

    Matrix = [[0 for x in range(ln)] for y in range(ln)]
    return Matrix
# __________________________________________________________________________________

'''
Function: Bifurcation
Description: Iterates around the Matirx in K range, getting maximum value according 
             to Bifurcation Equation
Input: Matrix, index i of current row, index j of current Column
Returns: Score which carries maximum value from K range 
         idx which carries the value of K with maximum value
'''
def Bifurcation(i, j, matrix):
    BiScore = 0
    BiIdx = -1

    for k in range(i + 1, j, 1):
        Sum = matrix[i][k] + matrix[k + 1][j]
        if (Sum > BiScore):
            BiScore = Sum
            BiIdx = k

    return BiScore, BiIdx
# __________________________________________________________________________________

'''
Function: Pairing
Description:Boolean Function that Checks if 2 nucleotide pairs together or no.
Input: 2 Nucleotides
Returns: Binary Decision ( True for Pairs / False for unpaired )
'''
def Pairing(nuc1, nuc2):

    if Pairs[nuc1] == nuc2:
        return 1
    return 0
# __________________________________________________________________________________
'''
Function: FillMatrix
Description: Iterates around the Matirx to Fill it according to "Nussinov Algorithm"
             Computes the Cells values of Diaganol, Down, left and K-Range (Bifurcation) of current cell[i,j] and takes 
             maximum as matrix[i,j] value.
             Also, if maximum come from K-range, it Stores the Values of K in dictionary for traceback.
Invokes: Pairing, Bifurcation
Input: Matrix, User input Sequence that we are Predicting it's structure
Returns: Matrix (Filled), Bifurcation Dictionary (illustrated in description section) 
'''
def FillMatrix(Seq, matrix):

    BifurcationDic = {}
    dummy = 0

    for j in range(1, len(Seq)):
        for i in range(dummy, -1, -1):
            BiScore, Bi_idx = Bifurcation(i, j, matrix)
            paired = Pairing(Seq[i], Seq[j])
            matrix[i][j] = max(matrix[i + 1][j - 1] + paired, matrix[i][j - 1], matrix[i + 1][j], BiScore)
            if matrix[i][j] == BiScore:
                BifurcationDic[(i, j)] = Bi_idx

        dummy += 1

    return matrix, BifurcationDic
# __________________________________________________________________________________

'''
Function: DotBracket_TraceBack
Description: Recursive Function that traces back the matrix then show the structure prediction by dot bracket notation
             changes dot to brackets if sequence element i pairs with seq element j as illustrated in the code below.
Invokes : Pairing Function
Input: Matrix,
       Current Indices i and j (For recursion),
       string DotBracket (thus, to be globally/statically Change),
       Bifurcation dictionary from Filling matrix Function.
Returns: list that carries the Dot bracket Notation of the structure prediction
'''
def DotBracket_TraceBack(matrix, i, j, DotBracket, BiDic):

    while matrix[i][j] != 0:

        Paired = Pairing(Seq[i], Seq[j])
        if matrix[i][j] == (matrix[i + 1][j - 1] + Paired):
            if Paired:
                DotBracket[i] = '('
                DotBracket[j] = ')'

            i += 1
            j -= 1

        elif matrix[i][j] == matrix[i][j - 1]:
            j -= 1

        elif matrix[i][j] == matrix[i + 1][j]:
            i += 1

        else:

            DotBracket = DotBracket_TraceBack(FilledMatrix, i, BifurcationDictionary[(i, j)], DotBracket, BiDic)
            DotBracket = DotBracket_TraceBack(FilledMatrix, BifurcationDictionary[(i, j)] + 1, j, DotBracket, BiDic)
            break

    return DotBracket

# __________________________________________________________________________________

# -------------------------------------------- Main ---------------------------------------------
#                                    ************* Global Declarations *************
'''
User Input, RNA seq required for Structure Prediction using Nussionov Algorithm
#  Applied Test Cases 
# Seq = 'GGGAAAUCC'
# Seq = 'CGGACCCAGACUUUC'
'''

Seq = input("Enter the Sequence : ")
'''
Dictionary carries the pairing nucleotides where  Key K  pairs with Value V { K : V }
'''

Pairs = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
# __________________________________________________________________________________
'''
String of Dots for TraceBak Dot Bracket Notation structure Prediction
'''
DBNotation = ['.' for x in range(len(Seq))]
# __________________________________________________________________________________
#                                    ************* Functions Calls *************

InitialMatrix = InitializeMatrix(len(Seq))
FilledMatrix, BifurcationDictionary = FillMatrix(Seq, InitialMatrix )
# __________________________________________________________________________________
#                           ************* Algorithm Output "Predicted structure" *************
''' 
 .join: Built in function function used to convert list to string 
'''
PredictedStructure = ''
PredictedStructure = (PredictedStructure.join(DotBracket_TraceBack(FilledMatrix, 0, len(Seq) - 1, DBNotation,BifurcationDictionary)))

print(PredictedStructure)