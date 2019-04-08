from __future__ import print_function
import random

#This script is based on 
#https://stackoverflow.com/questions/1024754/how-to-compute-a-3d-morton-number-interleave-the-bits-of-3-ints

def prettyBinString(x,d=32,steps=4,sep=".",emptyChar="0"):
    b = bin(x)[2:]
    zeros = d - len(b)


    if zeros <= 0: 
        zeros = 0
        k = steps - (len(b) % steps)
    else:
        k = steps - (d % steps)

    s = ""
    #print("zeros" , zeros)
    #print("k" , k)
    for i in list(range(zeros)): 
        #print("k:",k)
        if(k%steps==0 and i!= 0):
            s+=sep
        s += emptyChar
        k+=1

    for i in list(range(len(b))):
        if( (k%steps==0 and i!=0 and zeros == 0) or  (k%steps==0 and zeros != 0) ):
            s+=sep
        s += b[i]
        k+=1
    return s    

def binStr(x): return prettyBinString(x,32,4," ","0")


def computeBitMaskPatternAndCode(numberOfBits, numberOfEmptyBits):
    bitDistances=[ i*numberOfEmptyBits for i in list(range(numberOfBits)) ]
    print("Bit Distances: " + str(bitDistances))
    bitDistancesB = [bin(dist)[2:] for dist in  bitDistances]
    print("Bit Distances (binary): " + str(bitDistancesB))
    moveBits=[]
    
    maxLength = len(max(bitDistancesB, key=len))
    abort = False
    for i in list(range(maxLength)):
        moveBits.append([])
        for idx,bits in enumerate(bitDistancesB):
            if not len(bits) - 1 < i:
                if(bits[len(bits)-i-1] == "1"):
                    moveBits[i].append(idx)

    for i in list(range(len(moveBits))):
        print("Shifting bits by " + str(2**i) + "\t for bits idx: " + str(moveBits[i]))

    bitPositions = list(range(numberOfBits));
    print("BitPositions: " + str(bitPositions))
    maskOld = (1 << numberOfBits) -1

    codeString = "x &= " + hex(maskOld) + "\n"
    for idx in list(range(len(moveBits)-1, -1, -1)):
        if len(moveBits[idx]):
           shifted = 0
           for bitIdxToMove in moveBits[idx]:
                shifted |= 1<<bitPositions[bitIdxToMove];
                bitPositions[bitIdxToMove] += 2**idx; # keep track where the actual bit stands! might get moved several times
           # Get the non shifted part!     
           nonshifted = ~shifted & maskOld

           print("Shifted bef.:\t" + binStr(shifted) + " hex: " + hex(shifted))
           shifted = shifted << 2**idx
           print("Shifted:\t" + binStr(shifted)+ " hex: " + hex(shifted))

           print("NonShifted:\t" + binStr(nonshifted) + " hex: " + hex(nonshifted))
           maskNew =  shifted | nonshifted
           print("Bitmask is now:\t" + binStr(maskNew) + " hex: " + hex(maskNew) +"\n")
           #print("Code: " + "x = x | x << " +str(2**idx)+ " & " +hex(maskNew))

           codeString += "x = (x | x << " +str(2**idx)+ ") & " +hex(maskNew) + "\n"
           maskOld = maskNew
    return codeString


kdhash_size=4
numberOfBits = int(64/kdhash_size);
numberOfEmptyBits = kdhash_size-1;
codeString = computeBitMaskPatternAndCode(numberOfBits,numberOfEmptyBits);
print(codeString)

def partitionBy2(x):
    exec(codeString)
    return x

def checkPartition(x):
    print("Check partition for: \t" + binStr(x))
    part = partitionBy2(x);
    print("Partition is : \t\t" + binStr(part))
    #make the pattern manualy
    partC = long(0);
    for bitIdx in range(numberOfBits):
        partC  = partC | (x & (1<<bitIdx)) << numberOfEmptyBits*bitIdx
    print("Partition check is :\t" + binStr(partC))
    if(partC == part):
        return True
    else:
        return False

checkError = False        
for i in list(range(20)):
    x = random.getrandbits(numberOfBits);
    if(checkPartition(x) == False):
        checkError = True
        break
if not checkError:
    print("CHECK PARTITION SUCCESSFUL!!!!!!!!!!!!!!!!...")
else:
    print("checkPartition has ERROR!!!!")
