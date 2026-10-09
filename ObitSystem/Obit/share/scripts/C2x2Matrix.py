# Simple python 2x2 complex functions
#exec(open("C2x2Matrix.py").read())
Make2x2=None; matx2x2Mult=None; matx2x2CTrans=None;  matx2x2Inv=None
matx2x2Add=None; matx2x2Sub=None; matx2x2Scl=None; matx2x2Print=None

# Form a 2x2 complex matrix, forming complex values from pairs of reals
del Make2x2 
def Make2x2 (a11r, a11i, a12r, a12i, a21r, a21i, a22r, a22i):
    return [[complex(a11r, a11i), complex(a12r, a12i)],\
            [complex(a21r, a21i), complex(a22r, a22i)]]

# Inner (dot) Multiply 2 x2 complex matrices
del matx2x2Mult
def matx2x2Mult (in1, in2):
    out=[[complex(0.,0.),complex(0.,0.)],[complex(0.,0.),complex(0.,0.)]]
    for i in range(0,2):
        for j in range(0,2):
            out[i][j] = in1[i][0] * in2[0][j] + in1[i][1] * in2[1][j]
    return out

# end matx2x2mult 

# Complex transpose of a 2x2 complex matrix
del matx2x2CTrans
def matx2x2CTrans (in1):
    out=[ \
          [complex(in1[0][0].real, -in1[0][0].imag), complex(in1[1][0].real, -in1[1][0].imag)],
          [complex(in1[0][1].real, -in1[0][1].imag), complex(in1[1][1].real, -in1[1][1].imag)]]
    return out

# end matx2x2CTrans

# Inverse of a 2x2 complex matrix
del  matx2x2Inv 
def matx2x2Inv (in1):
    idet = 1.0/(in1[0][0]*in1[1][1] - in1[0][1]*in1[1][0])
    out=[[idet*in1[1][1],-idet*in1[0][1]], [-idet*in1[1][0], idet*in1[0][0]]]
    return out
    
# end matx2x2Inv

# Add two 2x2 complex matrices, element by element (in1+in2)
del matx2x2Add
def matx2x2Add (in1, in2):
    return [[in1[0][0]+in2[0][0], in1[0][1]+in2[0][1]], \
            [in1[1][0]+in2[1][0], in1[1][1]+in2[1][1]]]
    
# end matx2x2Add

# Subtract two 2x2 complex matrices, element by element (in1-in2)
del matx2x2Sub 
def matx2x2Sub (in1, in2):
    return [[in1[0][0]-in2[0][0], in1[0][1]-in2[0][1]], \
            [in1[1][0]-in2[1][0], in1[1][1]-in2[1][1]]]
    
# end matx2x2Sub

# Multiply a 2x2 complex matrix by a scalar
del matx2x2Scl
def matx2x2Scl (in1, scalar):
    return [[scalar*in1[0][0],scalar*in1[0][1]], [scalar*in1[1][0], scalar*in1[1][1]]]
    
# end matx2x2Scl

# Print a 2x2 complex
del matx2x2Print 
def matx2x2Print (in1):
    print(" %10.4g %10.4g   %10.4g %10.4g"%(in1[0][0].real,in1[0][0].imag,in1[0][1].real,in1[0][1].imag))
    print(" %10.4g %10.4g   %10.4g %10.4g"%(in1[1][0].real,in1[1][0].imag,in1[1][1].real,in1[1][1].imag))

# end matx2x2Print
