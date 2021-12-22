
# Python v3.4+ needed.

import functools
import json

import openssl
import subprocess
import base64
import operator
import time

import random as rd
import math
import numpy as np
from functools import reduce
import ademy


    

    ############################################################################
    #                          ARITHMETICA                                     #
    ############################################################################


"""
discrete-log given-order 'g'
"""
def gOrder(a,b,q):
    
    test = False
    while not test:
        inf = (a-1) // q
        sup = (b-1) // q

        k = rd.randint(inf, sup)
        p = 1+k*q
        if fermat(p, 2):
            test = True
            print('ok : ')
        else:
            print('nop...')

    while True:
        x = rd.randint(a, b)
        g = pow(x, k, p)
        if g != 1:
            break

    
    return p, g


def decrypt_elgamal(a, b, alea, p):

    h = pow(a, alea, p)
    invh = my_inverse(h, p)
    plain = (b*invh) % p
    tmp = str(plain)
    return int(plain.to_bytes(length=len(tmp), byteorder='big'))


def malleable_elgamal_attack(p, g, h, ciphertexts):
    

    a = ciphertexts[0]
    b = ciphertexts[1] 

    m = oracle('elgamal/malleability', a=a, b=2*b%p) # oracle to code returning dict #
    plain = m['m']
    plain //= 2 # ***** #
    
    tmp = str(plain)
    print(plain.to_bytes(length=len(tmp), byteorder='big'))

def RSA_NED_RECOVER(n, e, d):

    k = (e*d - 1) // 2

    numTest = 350
    boole = True
    while boole:
        l = []
        for i in range(numTest):
            print(i)
            x = rd.randint(1, n)
            test = pow(x,k,n)
            if test == 1:
                l.append('bad')
        if len(l) > numTest // 2:
            k = k//2
            print('div done')
        else:
            boole = False

    x = rd.randint(1, n)
    y = pow(x,k,n)
    while y == 1 or y == n-1 : 
        print('nop')
        x = rd.randint(1, n)
        y = pow(x,k,n)

    res = my_gcd(y-1, n) #use math.gcd if pb
    print(res)
    tmp = str(res)
    print(res.to_bytes(length=len(tmp), byteorder='big'))



def primes():
    """infinite generator of all prime numbers.
    >>> for p in primes():
    ...     print(p)
    ...     if p > 100:
    ...         break
    """

    yield 2 
    D = {}  
    q = 3   
    while True:
        two_p = D.pop(q, None)  
        if two_p:               
            x = q + two_p           
            while x in D:           
                x += two_p          
            D[x] = two_p            
             
        else:                   
            D[q*q] = 2*q        
            yield q
        q += 2

"""
n : pubkey
"""
def brutRSA(n):
    for p in primes():
      print(p)
      if n%p == 0:
        print('yes :')
        print(p)
      if p > 1e9:
        break
"""
n : list of modulus
a : list of a
"""
def chinese_remainder(n, a):
    sum = 0
    prod = reduce(lambda a, b: a*b, n)
    for n_i, a_i in zip(n, a):
        p = prod // n_i
        sum += a_i * mul_inv(p, n_i) * p
    return sum % prod
 
def mul_inv(a, b):
    b0 = b
    x0, x1 = 0, 1
    if b == 1: return 1
    while a > 1:
        q = a // b
        a, b = b, a%b
        x0, x1 = x1 - q * x0, x0
    if x1 < 0: x1 += b0
    return x1

def euclide_etendu(a, b):
    x,y, u,v = 0,1, 1,0 
    while a != 0: 
        q, r = b//a, b%a 
        m, n = x-u*q, y-v*q 
        b,a, x,y, u,v = a,r, u,v, m,n 
    return b,x 
    # b == gcd
    # x == coeff 'u'

def my_inverse(a, N):
    g, x = euclide_etendu(a, N)
    if g != 1: 
        return None # a is not inversible mod N
    else: 
        return x % N


"""
Solve Ax + b = 0 (n)
"""
def equaLinear(a, b, n):
    invA = my_inverse(a, n)
    x = (-b * invA) % n
    return x

"""
Binary exponentiation
"""
def my_expo_mod(N, g, n):
    h = 1
    nb = bin(n) # n in binary
    a = nb[2:len(nb)]
    l = len(a)
    
    avr = np.ones(l)
    
    for i in range(l-1):
        avr[i] = a[l-i-1] # little endian

    for i in reversed(range(l)): # i decreases
        h = (h*h) % N
        if (int(avr[i]) == 1):
            h = (h*g) % N
    return h

"""
Returns (a * b) % mod 
"""
def moduloMultiplication(a, b, mod): 
  
    res = 0 # Initialize result 
  
    # Update a if it is more than 
    # or equal to mod 
    a = a % mod
  
    while (b): 
      
        # If b is odd, add a with result 
        if (b & 1): 
            res = (res + a) % mod
              
        # Here we assume that doing 2*a 
        # doesn't cause overflow 
        a = (2 * a) % mod
  
        b >>= 1 # b = b / 2 
      
    return res

"""
extended gcd
"""
def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)
        
"""
Modular inversion
"""
def modinv(a, m):
    g, x, y = egcd(a, m)
    y = y
    if g != 1:
        raise Exception("modular inverse does not exist")
    else:
        return x % m

"""
Modular exponentiation
"""
def exp_mod(a,e,N):
    res = 1
    b = a
    i = 0
    while e>=2**i: # Invariant: b=a**(2**i)
        if e & 2**i != 0:
            e -= 2**i
            res = (res*b) % N
        b=(b*b) % N
        i += 1
    assert e==0
    return res

############## RSA Description ##############
'''
p = 166791685058524004244193349796967891057582375628222554421900760376139810612751853884085132130103413680376691955792163225140586439318149982468024594762580427199808711677295493369555968525961648250870739683078487979155223356437642580651193789994399718240684970872639527769790020308699057535288400803506636227539
q = 155306097993810867777202570386375784402875700164433188409305616187397405722848899623585922544727618171926020489589511325429592173518794776528979937779469260643672008948585098494527357280510062366000122891314529708565861044202068343351207199638739547588569317627723760685638378487164497330735729621939236963179
e = 206307031074795478380943527976463189513
N = 25903765784251968946010716815116022923475848759245099485477983704504350393487391342418948449749839345544334394331723836576048894766225616626825582832403587124559663231323356648965696965412643947351212915096416971674530138376343367331188055022849036485675154769198086350976962897916358668091536845980370850386027241628464031791649612760193724245306639717929539053247087615414155415166974892328858603159177303849571709314635860290090492471979707838594133442153817571700655120431758786600940513246274936525349586160445359470711645662198164809497580662539136651907984653461593842495599014728145112171094737246195208786481

phiN = (p-1)*(q-1)
pk=(N,e)
sk=(N, modinv(e,phiN))

# https://courses.cs.ut.ee/all/MTAT.07.002/2015_spring/uploads/homework/sheet-03.pdf
'''

"""
sk : secret key (int)
c : cipher (int)
"""
def rsa_dec(sk, c):
    (N,d) = sk
    return exp_mod(c,d,N)


"""
Returns gcd(a,b)
"""
def my_gcd(a, b): 
    while(b != 0):
        a, b = b, a % b
    return a

"""
Searching for a divisor using Fermat thm
"""
def temoinFermat(n):
    for i in range(1, 10**9):
        res = my_expo_mod(n, i, n-1)
        if res != 1:
            print(i)
            break
        else:
            print('nop')
            continue
    return res

def fermat(n, a):
    b = my_expo_mod(n, a, n-1)
    if b != 1:
        return False
    return True

"""
Miller-Rabin probabilistic test to check if
n is a prime number (stringer than Fermat test, cf : Carmichael numbers) 
T = 64 is a good choice
"""
def test_miller_rabin(n, T):
    if (T > 0):
        
        x = n-1 # x = (2^h)*m
        h = 1
        
        if (n > 2): 
            while( ((x / 2) % 2) == 0 ):
                h += 1
                x /= 2

        m = int((n-1) / (2**h))

        i = 0

        # Miller-Rabin test
        while (i < T): # Number of rounds
            i += 1
            a = rd.randint(1, n-1)
            b = my_expo_mod(n, a, m)

            if (b != 1 and b != (n-1)):
                j = 0
                while (j < h):
                    j += 1
                    if(b != n-1) and ((b*b) % n == 1):
                        return False
                    elif(b == n-1):
                        break    
                    b = (b*b) % n
                if(b != n-1):
                    return False
        return True

    return False

def testPrimeInterval(a, b):
    l = []
    for i in range(a, b, 2):
        #if(test_miller_rabin(i, 64)):
        if(fermat(i, 2)):
            print('yes')
            time.sleep(2)
            l.append(i)
        else:
            print('nop...')
    return(l)


    ############################################################################
    #                          SOME OTHER ATTACKS                              #
    ############################################################################


def bruteDico():    

    f = open('dico.txt', 'r')
    Lines = f.readlines()
    Lines = list(map(lambda s: s.strip(), Lines))

    for passphrase in Lines:
        '''
        As encryption contains randomness, two encryptions won't give same ciphertext
        '''
        myFile = open('../test/brute/cipherout/'+passphrase+'.txt', 'r') # Test all ciphertexts
        try:
            test = myFile.read()
            print(test)
            print("GOOOOOD PSWD : ", passphrase)
            myFile.close()
            time.sleep(15)
        except UnicodeDecodeError:
            print("Nop...")
    f.close()

"""
Exploits OTP using same mask twice
"""
def one_time_pad(a,b):
    a=openssl.decode_base64(a)
    b=openssl.decode_base64(b)

    xored64 = tmp(a,b) # Remove the mask

    bPart=b'0'*len(xored64) 

    xoredutf8 = tmp(bPart, xored64).decode()
    print(xoredutf8)
    
"""
Basic XOR of bytes strings
"""
def XOR(*seqs):
    return bytearray([functools.reduce(operator.xor, t, 0) for t in zip(*seqs)])

#Another XOR
def tmp(a,b):
    p=[]
    for a,b in zip(a,b):
        p.append(bytes([a ^ b]))
    return b''.join(p)



def basicRSA(n):
    out = subprocess.run(['openssl', 'pkeyutl', '-encrypt', '-pubin', '-inkey', 'pbkey.txt', '-in', 'plain.txt'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    res = out.stdout

    response = base64.b64encode(res)
    cipher = response.decode()


"""
Checks validity of card certificats
"""
def forensics():
        
    statuses = []
    CA = c.get('/bin/banks/CA') # CA certif

    dico = c.get('/bin/banks/forensics') # cards to check

    identifier = dico['identifier'] 

    card_numbers = dico['card-numbers']

    for i in range(0, len(card_numbers)):
        c1 = c.get('/bin/banks/card-data/'+card_numbers[i])

        card_certificate = c1['card-certificate']
        bank_name = c1['bank-name']
        card = c1['card-number']
        bank_certificate = c1['bank-certificate']
        challenge = c1['challenge']
        signature = c1['signature']

        # Use 'with' to be more compact
        myFile = open('card_certificate.txt', 'w')
        myFile.write(card_certificate)
        myFile.close()

        myFile = open('bank_name.txt', 'w')
        myFile.write(bank_name)
        myFile.close()

        myFile = open('card_number.txt', 'w')
        myFile.write(card)
        myFile.close()

        myFile = open('bank_certificate.txt', 'w')
        myFile.write(bank_certificate)
        myFile.close()

        myFile = open('challenge.txt', 'w')
        myFile.write(challenge)
        myFile.close()

        myFile = open('signature.txt', 'w')
        myFile.write(signature)
        myFile.close()

        myFile = open('CA.txt', 'w')
        myFile.write(CA)
        myFile.close()

        #### VERIFICATION ####

        #- The bank's certificate is not signed by the certification authority's certificate
            #openssl verify -trusted CA.txt -untrusted card_certificate.txt bank_certificate.txt
        
        args1 = ['openssl', 'verify', '-trusted', 'CA.txt', 'bank_certificate.txt']
        result1 = subprocess.run(args1, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        test1 = result1.stdout.decode()

        if "OK" in test1:
            print(test1)

        else:
            print(test1)
            print("The bank's certificate is not signed by the certification authority's certificate")
            statuses.append(False)
            continue

        #- The bank's certificate does not have the 'CA bit' that allows it to sign other certificates

        args2 = ['openssl', 'x509', '-text', '-noout', '-in', 'bank_certificate.txt']
        result2 = subprocess.run(args2, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        test2 = result2.stdout.decode()

        #Verifier la presence de "CA:TRUE":
        if "CA:TRUE" in test2:
            print(test2)
        else:
            print(test2)
            print("The bank's certificate does not have the CA bit that allows it to sign other certificates")
            print("OR The bank's certificate is self-signed")
            statuses.append(False)
            continue

        #- The bank's certificate is invalid (the signature does not verify) 
            #openssl verify -CAfile CA.txt bank_certificate.txt 
        
        args3 = ['openssl', 'verify', '-trusted', 'CA.txt', 'bank_certificate.txt']
        result3 = subprocess.run(args3, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        test3 = result3.stdout.decode()

        if "OK" in test3:
            print(test3)
        else:
            print(test3)
            print("The bank's certificate is invalid (the signature does not verify) ")
            statuses.append(False)
            continue

        #- The bank's certificate is self-signed
            #openssl verify -trusted CA.txt bank_certificate.txt
        
        args4 = ['openssl', 'verify', '-trusted', 'CA.txt', '-untrusted', 'card_certificate.txt', 'bank_certificate.txt']
        result4 = subprocess.run(args4, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        test4 = result4.stdout.decode()
        if "OK" in test4:
            print(test4)
        else:
            print(test4)
            print("The bank's certificate is self-signed")
            statuses.append(False)
            continue

        #- The card certificate is not signed by the bank's certificate (but by another one, probably fake)
            #openssl verify -trusted CA.txt -untrusted bank_certificate.txt card_certificate.txt ## OK PRESENT ##
        
        args5 = ['openssl', 'verify', '-trusted', 'CA.txt', '-untrusted', 'bank_certificate.txt', 'card_certificate.txt']    
        result5 = subprocess.run(args5, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        test5 = result5.stdout.decode()
        if "OK" in test5:
            print(test5)
        else:
            print(test5)
            print("The card certificate is not signed by the bank's certificate (but by another one, probably fake)\n")
            print("Or The bank's certificate is self-signed\n")
            print("Or The bank's certificate is invalid (the signature does not verify)\n")
            statuses.append(False)
            continue


        #- The card contains a secret key that does not match its certificate (invalid challenge signature)
        #= The signature of the challenge by the card is invalid.
            #openssl x509 -pubkey -noout -in card_certificate.txt
        
        args6 = ['openssl', 'x509', '-pubkey', '-noout', '-in', 'card_certificate.txt', '-out', 'key']
        subprocess.run(args6)

        # SIGNATURE IN BINARY (for openssl) :
        #args7 = ['base64', '-d', 'signature.txt']
        #out = subprocess.run(args7, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        signa = openssl.decode_base64(str(signature))

        myFile = open('signature.bin', 'wb')
        myFile.write(signa)
        myFile.close()
        
        # Take the output and check the signature of the challenge :
        #openssl dgst -sha256 -verify key -signature signature.bin challenge.txt 

        args8 = ['openssl', 'dgst', '-sha256', '-verify', 'key', '-signature', 'signature.bin', 'challenge.txt']
        result8 = subprocess.run(args8, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        test8 = result8.stdout.decode()

        if "OK" in test8:
            print(test8)
        else:
            print(test8)
            print("The card contains a secret key that does not match its certificate (invalid challenge signature)\n")
            print("Or The signature of the challenge by the card is invalid\n")
            statuses.append(False)
            continue

        #- The bank name and/or card number in the card certificate does not match the transaction data
        #openssl x509 -subject -noout -in card_certificate.txt

        args9 = ['openssl', 'x509', '-subject', '-noout', '-in', 'card_certificate.txt']
        
        result9 = subprocess.run(args9, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        test9 = result9.stdout.decode()

        # Check if we have bank name / num card

        if bank_name and card in test9:
            print(test9)
        else:
            print(test9)
            print("The bank name and/or card number in the card certificate does not match the transaction data")
            statuses.append(False)
            continue

        #- The name in the bank's certificate does not match the bank's name in the transaction
        #openssl x509 -subject -noout -in bank_certificate.txt
        # check if we find bank name in the output 

        args10 = ['openssl', 'x509', '-subject', '-noout', '-in', 'bank_certificate.txt']
        
        result10 = subprocess.run(args10, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        test10 = result10.stdout.decode()

        # Check if we have bank name / num card

        if bank_name in test10:
            print(test10)
        else:
            print(test10)
            print("The name in the bank's certificate does not match the bank's name in the transaction")
            statuses.append(False)
            continue

        statuses.append(True)

    return statuses


"""
RSA with oracle returning MSB
"""
def RSAoracle(n, e, ciphertext, MSB):
        
    C = ciphertext

    stop = 1
    #0
    a = 0
    b = n
    i = 0

    #1
    while True:
        if a == b:
            while True:    
                if pow(a,e,n) == C % n:
                    return a
                else:
                    a+=1
                    
        #2
        if i==1650 or i == 1800 or i==2000 or i==2200 or i==2500:
            RSAoracle(n, e, ciphertext, MSB) # pb 
            
        Cp = (pow(2**i,e)*C) % n
        res = oracle('/most-significant-bit', c=Cp) #oracle to code
        rep = res['MSB']
        print(rep)
        if not rep:
            
            if stop == 1:
                print(stop)
                j = 0
                while(j < 490):
                    b = (a+b) // 2
                    j+=1
                    i+=1
                        
            b = (a+b) // 2
        else:
            a = (a+b) // 2
            
        #3
        i+=1
        print(i)
        stop+=1

    #plain = a
    #tmp = str(plain)
    #print(plain.to_bytes(length=len(tmp), byteorder='big'))



#Pollard RHO 
def F(x, c):
    N = 51236634648386213823681213894281566384322765033598607167583847485631067741913287110004977863159451355891171261857680965509129824255390015985834594319696414376345399368633910213965197696436338621957445623842620520556278798910388353216781247390099032256832702209515002139656635806283376428189606144257666040443670690314197636255023777097429682940477448514907133565165315142412357848658468256707279063420192871172309774863645412939382503025944197866945175604344266854531693196836951934837745241524673099758619783116729051093674460987371493713020950455804949437921990077864852407333443900512613895027456766372388827757231
    return (x*x + c) % N
    
def brent(c=1):
    N = 51236634648386213823681213894281566384322765033598607167583847485631067741913287110004977863159451355891171261857680965509129824255390015985834594319696414376345399368633910213965197696436338621957445623842620520556278798910388353216781247390099032256832702209515002139656635806283376428189606144257666040443670690314197636255023777097429682940477448514907133565165315142412357848658468256707279063420192871172309774863645412939382503025944197866945175604344266854531693196836951934837745241524673099758619783116729051093674460987371493713020950455804949437921990077864852407333443900512613895027456766372388827757231
    
    #init
    trap = 0
    nexttrap = 1
    slow = 2
    fast = slow
    i = 0

    while True:
        
        #+1
        #fast = F(fast, c)
        
        #improved '+1'
        
        slow = F(slow,c) 
        fast = F(F(fast,c),c)

        slow2 = F(slow,c)
        fast2 = F(F(slow,c),c)

        #done?
        g = my_gcd((fast-slow)*(fast2-slow2), N) #use math.gcd if pb
        
        #g = my_gcd(fast-slow, N)
        if g == N:
            c=c+1
            print(c)
            brent(c)
        if g != 1:
            print('yes :')
            print(g)
            return g
        
        #slow MAJ
        if i == nexttrap:
            trap = nexttrap
            slow = fast
            nexttrap = 2*trap
        i+=1 

#Pollard RHO 'p-1'
def RSA_p_1_pollard(B, N, beg):
    #B = int(1e7)
    #N = 117827681420271584017432903522327303325344948050665323956545863

    a = 2
    for i in range(beg,B):    # i == M(B)
        print(i)
        a = pow(a, i, N) 
        p = a-1 % N             # a^M(B) -1 = 0 mod N
        g = my_gcd(p, N) #use math.gcd si pb

        '''
        elif g == N:
            B-=10
            RSA_p_1_pollard(B, N)
        '''
        
        if g == 1 or g == N: #increase i
            continue

        else:
            print('yes')
            print(g)
            return g

#Pollard RHO 'p-1' 2nd phase
def RSA_p_1_phase2(B, C, N):
    #B = int(1e7)
    #C = int(1e9)
    
    # y (Primes avec Pk) 
    L = []
    Pk = 2*3*5*7*11*13*17
    for i in range(1, Pk):
        if my_gcd(i, Pk) == 1:
            L.append(i)

    a = 2
    debut = 1 #int(4e6) #????
    for i in range(debut, C):    # i == M(B)
        print('i =', i)
        time.sleep(1)
        b = pow(a, i, N)
        
        for x in range(B//Pk + 1, C//Pk + 1):
            print('x =', x)
            
            D = pow(b, Pk*x, N)

            X = 1
            #for y in L:
            for z in L:
                print('z =', z)

                #E = pow(b, y, N)
                E = pow(b, z, N)

                #X *= (D*E - 1) % N
                X = ( X*(D - E) ) % N #X *= (D - E) % N

                g = my_gcd(X, N)
                if g !=1:
                    print('yes !!')
                    print(g)
                    return g

#Discret log using Pollard RHO
def Flog(x, q, c=1):
    return (x*x + c) % q

def H(alpha, g, h, q):
    return ( alpha*h*pow(g, Flog(alpha, q), q) ) % q

def discreteLog_rho_pollard(challenge):
    dico = challenge
    p = int(dico['p'])
    g = int(dico['g'])
    h = int(dico['h'])

    i = 1
    x = 0
    alpha = h
    y = Flog(alpha, p)
    beta = H(alpha, g, h, p)

    while alpha != beta:
        print(i)
        x = ( x + Flog(alpha, p) ) % p; alpha = H(alpha, g, h, p)
        y = ( y + Flog(beta, p) ) % p; beta = H(beta, g, h, p)
        #y = ( y + Flog(beta, p) ) % p; beta = H(beta, g, h, p)
        i+=1
    if i < p:
        return ( (x-y)*modinv(i, p) ) %p
    else:
        return None


def discreteLog(challenge):

    dico = challenge
    
    p = int(dico['p'])
    #print(p)
    g = int(dico['g'])
    #print(g)
    h = int(dico['h'])
    
    q = 2**(8*challenge) #ordre < 2**8i et premier
    
    q = ademy.generate_big_prime(2**48, q)

    H = {}
    T = int(math.ceil(math.sqrt(q)) + 1)
    #print(T)

    for i in range(T):
        print(i)
        H[exp_mod(g, i ,p)] = i
    
    gT = exp_mod(g, T, p)
    S = my_inverse(gT, p)

    u = h
    i = 0
    while True:
        try:
            H[u] == H[u]
            x = i*T + H[u]
            print('yes :')
            print(x)
            return x
            
        except KeyError:
            print(i)
            u = (u*S) % p
            i+=1


"""
Exploit DSA signature
using same k twice
"""
def dsa_attack(p, q, g, h, signature1, signature2, message1, message2):
    
    # after parsing signature1 & signature2 : (man asn1parse...)
    r = 95483180317932451814263334791572618975463287486857963143340873671819885984274
    s1 = 59036617777983599705042477901485475224108187316635007267871704808648638709151
    s2 = 33748183497872691790493208335673201846013013435195484997484079916573999053984
    
    sha_msg1 = sha256(message1.encode()).hexdigest()
    sha_msg2 = sha256(message2.encode()).hexdigest()
    sha1 = int(sha_msg1, base=16)
    sha2 = int(sha_msg2, base=16)

    k = ( (sha1 - sha2)*modinv(s1-s2, q) ) % q
    x = ( (k*s1 - sha1)*modinv(r,q) ) % q
    oracle('/sbin/key-hack/dsa/sk', p=p, g=g, q=q, x=x) # produce secret key
    
    
    
    
    ############################################################################
    #                          SOME PROTOCOLS                                  #
    ############################################################################
    
    

def kerberos():
    """
    Kerberos scheme
    """
    dico = #######

    client_TGS_session_key = dico['Client-TGS-session-key']
    TGT = dico['TGT']

    dico = {'username': 'XXXXXXXX', 'timestamp': time.time()}
    jsoned = json.dumps(dico)

    session_key = openssl.decrypt(client_TGS_session_key, 'password')
    authenticator = openssl.encrypt(jsoned, session_key)

    dico = #######

    Client_Server_ticket = dico['Client-Server-ticket']
    Client_Server_session_key = dico['Client-Server-session-key']

    dico = {'username': 'XXXXXXXX', 'timestamp': time.time()}
    jsoned = json.dumps(dico)
    Server_session_key = openssl.decrypt(Client_Server_session_key, session_key)
    authenticator = openssl.encrypt(jsoned, Server_session_key)

    key = Server_session_key

    authenticator = #######
    
    #print(gatewayKerb(dico, key)) # comm with server through gateway (to code)
    
    

from hashlib import sha256

def dh():
    """
    Diffie-Hellman exchange
    """
    
    param = #####
    p = int(param['p'])
    g = int(param['g'])

    CA = #CA certif (sent by a third-party server) 
    myFile = open('CA.txt', 'w')
    myFile.write(CA)
    myFile.close()

    # DH EXCHANGE
    x = rd.randint(p + 1, int(p**2)) % p

    A = pow(g, x, p)
    dico = # server send its B = g**y

    B = dico['B']
    k = dico['k']
    signature = dico['signature']

    myFile = open('signature.txt', 'w')
    myFile.write(signature)
    myFile.close()
    
    ''' Signature into binary '''
    signa = openssl.decode_base64(str(signature))
    myFile = open('signature.bin', 'wb')
    myFile.write(signa)
    myFile.close()

    string = str(A)+','+str(B)+','+str(k)+','+'amycastro'

    myFile = open('string.txt', 'w')
    myFile.write(string)
    myFile.close()

    args = ['openssl', 'x509', '-pubkey', '-noout', '-in', 'CA.txt', '-out', 'key']
    subprocess.run(args)

    #verify signature
    args1 = ['openssl', 'dgst', '-sha256', '-verify', 'key', '-signature', 'signature.bin', 'string.txt']
    result1 = subprocess.run(args1, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    test = result1.stdout.decode()

    if "OK" in test:
        print(test)
    else:
        print(test)
        exit()

    #shared key
    AB = pow(int(B), x, p)
    
    K = sha256(AB.to_bytes(256, byteorder='big')).hexdigest()

    T = str(A)+','+str(B)+','+str(k)+','+'UGLIX'
    confirmation = openssl.sign(T)

    # comm with server through gateway (to code...)
    
    #dico = confirmation 
    #print(gateway(dico, K))

    return K
        
    
    


'''
def obfuscation():
    
    iu = <obfucated_code>


    >>> g = zlib.compress
    >>> h = base64.b85encode
    >>> public_key = h(h(g(g(h(g(h(h(original))))))))  # ENCODAGE (original doit etre Bytes)

    THEN
    
    >>> u = zlib.decompress
    >>> v = base64.b85decode
    >>> exec(v(v(u(v(u(u(v(v(public_key)))))))))


    >>> seed = h(h(g(h(g(h(g(h("exec"))))))))

    >>> v(u(v(u(v(u(v(v(seed))))))))
    "exec"
    >>> eval(v(u(v(u(v(u(v(v(seed)))))))))

    eval(v(u(v(u(v(u(v(v(seed)))))))))(v(v(u(v(u(u(v(v(public_key)))))))))

'''


