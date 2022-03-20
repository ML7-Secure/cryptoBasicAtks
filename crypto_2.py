
import base64
import time


############################################################################################################

def int2str(m):
	size = len(str(m))
	print("".join([chr((m >> j) & 0xff) for j in reversed(range(0, size << 3, 8))]))
	
	
import functools
import operator

def XOR(*seqs):
    return bytearray([functools.reduce(operator.xor, t, 0) for t in zip(*seqs)])


############################################################################################################


from functools import reduce

def chinese_remainder(n, a): # n : list of n_i, a : list of a_i 
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

def find_invpow(x,n):
	"""Finds the integer component of the n'th root of x,
	an integer such that y ** n <= x < (y + 1) ** n.
	"""
	high = 1
	while high ** n <= x:
		high *= 2
	low = high//2
	while low < high:
		mid = (low + high) // 2
		if low < mid and mid**n < x:
			low = mid
		elif high > mid and mid**n > x:
			high = mid
		else:
			return mid
	return mid + 1
	
def inversePow(x) :
    m = find_invpow(x, 3) #exponent e = 3
    size = int(math.log(m, 10) + 1)

    print("".join([chr((m >> j) & 0xff) for j in reversed(range(0, size << 3, 8))])) # int -> str


n1 = 26509383079757729374361821335566463340679677420460038491665579683616976911128418667416887613568919704775424532884175604720478841185657319115097415218930264737142671667012715020750639982192915302212523754250715575158049956313812050644740213053637653414897908716790699633384991768546077429753042059508162119547635927545087499875683564886148889969969214540432692933227448478777068082768081383604266132373744975817568252688015001772473257950818441953331085993886792555910772889184286705784317556474546749076992262139799501013985870270367986856830456522701212348726388043343753526439919423358488017659274368002346575363609
n2 = 24971828782013633869990334770689398604825637354353205744758294421546019612927246003066101492583044750673276034168507657032031629812581982258058507180860205778615864034659554371184515308986659731256567899185519319662707073112114187610268208484273157351575252046784361136735654190453134804266392564720741868081157844532970165859555032055257769817167455902455108762152231267967365577333371828003465133644429499536268097333093364171605297590159220134673107150553296782366104376209453585894385305396516922976972144776285763565046722256633611754857244076000605858464624860318060942796144105640086696239779536663560931058751
n3 = 28704197778369156717701608250108541540731624656902849101916972370340017832959380557120076777240811373834678321790467920505263803618614785844852603569832342496105351172703653664858374803594921115288131310463112487046531928176162619648989774488913510815772024580538104712100693918218511723194255201311606125185106899668696566704145629546016285505116033693030922721606033663310472576216969552243007874199656316379077355287113834598187347321936853260483733812989414306362689440394694525254744849731195657290810261198772136650777669677753742276079289337291970500037836131315371274223109240367950286214700193861680028490149

c1 = int("b70fe1840962335f6b943b033a8247e7b4340495e550f28ff4b93f0655d4294d5bd58777a2dd04a6941c1fadf09a9dedba28b0b9f551a274deff3700b0b733d8741ff14196de87a27b44c89aa14bb0e6f370a9690f060966af35828bd083f257fe4e3b7f57697a484fdfee1670c5d37362ae8d6be993aeebf722c26927662a52df7d19f57d29bbfdff9facb0fc72ba2d0fc9f838900e7bd5b36725bd4a5710b4b92c90b8a555aa1c91e6d6fea81a8e90974775d23d954457ffea7e4e546f1e8108e3bf07b87deaf9d0e68b7f01fb3cece8c7c2683041d0225b7da537c87ec90344194821d9fd7defafd1de83f9e9a903c9495f9482d1d759b6094fe38f251ecc", 16)
c2 = int("44be1f53135747c04dd4fdc946ac57b08f186b5d75e5a2c5aee34290cd1753de218e8555517da831019150299c1839f6915fd287e1c809abb2ec84ffdfd7397c9be2e23974e29b011c6efe39d2f66d669b358ad003a24815c301066558b736af1638f6a3624645db8719e5e981e4caeaba27f3a4c29ad06e18eac859042c2ce319c127fb444d92c41ce95f851b39dc227cc2c01cde0b600f972d1e28bd21269d117c7fbdaafc3a61b91e211798a63074bdd6c6814614cfb1c5b02e63fa60f9feaa3d10db0785d9066aa4c8e327925631128984ee0de140698d52c2e1b7b755a38b955f911d48d43fe8f501aa5a4b79638ce8567df27978fffdc34a80b91aba79", 16)
c3 = int("78688f7588c0c1249b2247bd9a9049ee556a8253e6f85a25e7c2c709c79c3f10177bb9e23138c05cf1d0caddf8b6e67889ecf4b5e4b69583033eecccd6b0e6027d3d7cdfe782f1e57500e5e25506c73866beb000e919fe61d507508cb52b8b2aeff865a34db69a8c48ef0ea4ac961034e7906a73129b040c6dbf7248354641a47e98602532e21271471a64f3eb2d2903b6db5620b40a5479d9838905997f19e04630aea0ad2f1966892dbf8df4741f64501a8630382e9a52fe5fb655f30e5fe57fb0b93a5684d9196d827110ea8bdd161957d1298291a09bd3c49993d2d0e038bad85fdd6cbc46886edbd26ca2189fdf999e9ea5155a0212eb49df9d55e00763", 16)

'''
x = chinese_remainder([n1, n2, n3], [c1, c2, c3])

# Then remove the exponent 'e = 3' 
inversePow(x)
'''

############################################################################################################

from math import gcd
import sys
sys.setrecursionlimit(10000)


def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def modinv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        raise ValueError('Modular inverse does not exist.')
    else:
        return x % m

def attack(c1, c2, e1, e2, N):
    if gcd(e1, e2) != 1:
        raise ValueError("Exponents e1 and e2 must be coprime")
    s1 = modinv(e1,e2)
    s2 = (gcd(e1,e2) - e1 * s1) // e2
    temp = modinv(c2, N)
    m1 = pow(c1,s1,N)
    m2 = pow(temp,-s2,N)
    return (m1 * m2) % N

def RSAsameModulus():
	modulus = 28183968387059668200642458248057716733652295618039278784203296865663977123732851158241552543162802091503168436398822406381888586753727666689169168863858770646662491176891215623973229761305568424652718550468106794209096340712072154620370625046489647128697565767684590029231732968691024504527017933532298688746733494367988869400366356627732145293192972833456049716215060333632670786154088086127675404704950187603964795187133755476501136674627272637654671280359452614480652409319668524012149842793053219760123295630680943514851848255752097045910457800112446086927045427912014279667854724067393312507022174687159370315011
	e1 = 12170559117698521819955250243566569147026579404003724750674488316631871172163721681656153129465829215234649655043612070680218223383110316249292581656542026770927663263905744047892234348851756508444326653142674013285416650449589766213317788238204171579004696716590851385564666309960568315849050359728843495100803148900436473519858032836574890333704455047469701542318792814498275350362420400144493280184124519225556642092278322243437332229018388221763871823027539991934053711208325119643630549591741443177678506164955072013671563462384116764237262087670682551085544813047728152890605294289417654434115406046671167806449
	e2 = 4299862168655641713044275147316638525821909220418066503994418816312720919180329782468341561979479290105636509672401154292306359853983392880872928500114221929950192678246611347056701523839266105514092372199649464533371062336159078601844667559124558128121191723479366463147394398127930150518742797732183112312976873281059501145602543347203341635227531177584241227078451715737311404754202767259520282980446986068094629442422913158738821211778399337774011335977165391600410605330729227702254071934520203165643346794832741338109179740531096746557455604571798752959467303078172371074045847701834909788067501534201155030189 
	ct1 = int("816d5099838a2dc663ed04d4d85af2d7a45ffa679cbbddd65e50e168603cbe2ef5a34a6de2d75d678d5c015f1bb929b97cdb9790a46e835f8c66afebdf8c90524c9c9f992a2fa22d138f1e6dc04fa5b08df97faf536ab1e5d36d638929758dc56e1c52b37be0f18a411aaea549647807503054e950626f32194f9327963be45d08b967de35ad9f21581ced3dd051aa8d364757e273362ea6f530bedaff7895d59f36fb4f384185fac201bb75098e9ff6ae0d16eccafe351211074a2f9d1231eb7b97a9cb5be2f1f2b0e4106c70936703d0ac5ebe8bcdb8b5d077d43cb64ba52c7dd39e413dc3a8773c0781dc0a447e04e3a3e4aac69455f02f02c3f14f8f2122", 16)
	ct2 = int("d8641f6a677c0592de9903e3e49e58bfe7e0e8599c878d8c1149075ffb4f19add068ad5930238c886015a45ec029487fadcb96aacde813630fd385aee9a86a5a0d1aeef2e2b5cf4d2006cd2da70d3443372cb7119873bccff00ad2a17a6cb905623b7c5d6dbd3e4058322defbf8a5a23bf3e98f5acd4e3276daae6ba8e3f2282ba481a48f539769bc26cf53bab9930e7e97af9621ddf4e49483163d9cf489000ccbcf97a563745c8bb1d2bd5594016539557ec0ddea8d37fa92dd50ada33af93defc181c86534bbe9d99e5e62d487bf4f0057b433e1fa88a7d0f9fe694236134eb0dbcfea52595790b72c00ac3ee1853449a93d02cadb90b87b5f4a27a241b13", 16)
	print('Started attack...')
	
	message = attack(ct1, ct2, e1, e2, modulus)
	print('Attack done')
	print('Message:', format(message, 'x')) # in hex, convert to ascii


############################################################################################################


def sdsa_attack():

    p = 32103350193216866699084389784654264256642443788953664365664396449147788929780524052790241213097281642317446848155374333255115863564802320269154328545453275641270417847449673827042781673292749423336356877089738239303460345648240312312424630290356457754840838107685033431207983307460442968141902027222993928891225348588403282519230147996368271939207819439016149159582760855623599288070246080050300317787375511573118906295251011927575480131197195951438278749299335976469535222505388284007147918241899108016479729436941938409338798947559676538970046815141091084305492786275393676369005435381818167682976935328497197922509
    q = 108527423945311937134504112952590721069656646382056042230340968328541269567883
    g = 19070099188261131869578707714015499593262222797949407843771788498975684947408952156208180936247898293240493826997571932855010570708353857810328312294651078393828223890173447678205267479511900167944430829378492812160464495193425371500357569716331270778326755522728036091525618633071751959998812045119924922331361954703991697306768601964993460317274152612947063576717586472978812415887075615099517993601235023187979283182525187295295639144617202798028040811578429110895764248568541863871537701717199102532345065167835100897427438397399967439965133443067101762566626954551279405942564706119203326117127044929245034473387
    #h = 27956889212109587321074342473493143554087986201642221918255234340252001103635453075268303480652137358377283645475322525221271641308593646980775004945992463729954040533533654806467528460158871576286455955600881839203039462065073066684024715715666256521829393643948490942138243200312242064879252231441738658598868250218931198660000848306983080833207023659920197571781024532057578150396409330582709648267312466879021251968640342406408285089568029726211861870186061761663581424453474240885568137242343106223106021522538982066752465316562222158454256392394711353720775302675183983208259651462063157016109030115983305060011 #g^x

    #signature1 = '3045022100BA03418D5CB21043D2FC5E91FA20885591CB3BFA635E0B6D3329778B325EDD0A022058654D237A106D636AD4B3C83028A8DE6296B4EE863E4F394E46E26A014A0CFE'
    #signature2 = '30460221009C1CB5586F42C4E12BDB397EDA2DC066315F8D497BC6B11E06E4CDF028B51A95022100BC09882124EC7E7CA04729C541F2139FEB571C725C2D39E4AC7675E99C889211'
    
    # after parsing signatures : openssl asn1parse -inform DER -in sign.der
    c1 = int("FDD2949CB4C0E6BA28D0C8E9454B9EF5467992481DC942ADD9061CBE862CB668", 16)
    s1 = int("3890C5BA7437ADB19C05942A520F35901F4138256801D20590299B23961599BA", 16)
    
    c2 = int("B4D66776EAEA3C668542D19125AE8DBC90B608A9BDABC7EFB10C7F826A56DF5D", 16)
    s2 = int("B33C07F72676D647C7EBBA84F4E6BD6C5441670291CCD9B56F01A30A6D1D6745", 16)
    
    	
    x = ( (s1-s2) * modinv(c1-c2, q) ) % q # vuln = same random k re-used
    
    print(x) # x = 24590207231390600909863107430259094608295144781632333348330626364935095055144
    

import random as rd
from hashlib import sha256

def signSDSA(x, p, q, g, M, h):
	
	k = rd.randint(1, p)
	r = pow(g, k, p)
	
	toHash = M.encode() + r.to_bytes(2048 // 8, 'big')
	print('toHAsh = ', toHash)
	
	Hash = sha256(toHash).hexdigest()
	print('HAsh = ', Hash)
	
	c = int(str(Hash), 16) #int('BA03418D5CB21043D2FC5E91FA20885591CB3BFA635E0B6D3329778B325EDD0A', 16)
	
	s = ( k + (c*x) ) % q #int('58654D237A106D636AD4B3C83028A8DE6296B4EE863E4F394E46E26A014A0CFE', 16)
	
	print('c = ', hex(c))
	print('s = ', hex(s))
	
	
	print('\n## Verification ##\n')
	
	h_c = pow(h, c, p)
	
	inv_h_c = modinv(h_c, p)
	
	r_verif = ( pow(g,s,p) * inv_h_c ) % p
	#print(r == r_verif)
	
	toHash_verif = M.encode() + r_verif.to_bytes(2048 // 8, 'big')
	Hash_verif = sha256(toHash_verif).hexdigest()
	print( c == int(str(Hash_verif), 16) )
	
	return c, s
'''	
p = 32103350193216866699084389784654264256642443788953664365664396449147788929780524052790241213097281642317446848155374333255115863564802320269154328545453275641270417847449673827042781673292749423336356877089738239303460345648240312312424630290356457754840838107685033431207983307460442968141902027222993928891225348588403282519230147996368271939207819439016149159582760855623599288070246080050300317787375511573118906295251011927575480131197195951438278749299335976469535222505388284007147918241899108016479729436941938409338798947559676538970046815141091084305492786275393676369005435381818167682976935328497197922509
q = 108527423945311937134504112952590721069656646382056042230340968328541269567883
g = 19070099188261131869578707714015499593262222797949407843771788498975684947408952156208180936247898293240493826997571932855010570708353857810328312294651078393828223890173447678205267479511900167944430829378492812160464495193425371500357569716331270778326755522728036091525618633071751959998812045119924922331361954703991697306768601964993460317274152612947063576717586472978812415887075615099517993601235023187979283182525187295295639144617202798028040811578429110895764248568541863871537701717199102532345065167835100897427438397399967439965133443067101762566626954551279405942564706119203326117127044929245034473387

x = 24590207231390600909863107430259094608295144781632333348330626364935095055144

h = 27956889212109587321074342473493143554087986201642221918255234340252001103635453075268303480652137358377283645475322525221271641308593646980775004945992463729954040533533654806467528460158871576286455955600881839203039462065073066684024715715666256521829393643948490942138243200312242064879252231441738658598868250218931198660000848306983080833207023659920197571781024532057578150396409330582709648267312466879021251968640342406408285089568029726211861870186061761663581424453474240885568137242343106223106021522538982066752465316562222158454256392394711353720775302675183983208259651462063157016109030115983305060011

M = str(input("what is the challenge : "))


c, s = signSDSA(x, p, q, g, M, h)

from asn1crypto.core import *

class MySequence(Sequence):
    _fields = [
	('field_one', Integer),
	('field_two', Integer),
    ]

mySeq = MySequence()
mySeq['field_one'] = c
mySeq['field_two'] = s

#print(mySeq['field_one'].dump(force=True))
#print(mySeq['field_two'].dump(force=True))

print(mySeq.dump(force=True).hex()) # signature to der format
'''

############################################################################################################

#from fpylll import *

def franklinReiter(): # Sage

	n = 0x00f02e1f0265e0b18cb449cbe5db30bfb50656239158572bd54922d2a94877b7db6c49ff7969d00f5129e161aea2d26629d2faed05a8e0d7739ee66efd4b989487c5e06fc21ba2c534f6730fc3bf5db5672e5da73e92ee293f1d33358f9d6979259b8353b16e279a88a50c8feed9a028718d9aabba489ff500462fd44d25a94b2023f2b0b59542b7df11e5c32c5a444527641cebb59b0077710a11dbef33f4609e073483808e08c005b3b5bf9e64535c709e5ffde2d3cfcc3cf9264f02ae94898d32e93f4306cc3d829a0af3a602d04628a249305a908ae4d9e4e41e786a73e196888ba76e3292118bbc8b92d621c812b947cff5f9e59583443414d81c9f01de07

    c = 0x12c3156417add22bedf406471fa7b1650679e48d4b8f4cd042399e2ff3cbbe625371f1a9ce8eb5cb9be69ca6065f8c90c1bfdc99628ac332f588b16c8412071237064edeb5c6997b1bbf494534c55fdc7f54ee78d2650459c5063427f187c25cbbbfb919062ff5b3eb75143b3105ae7d7074bb7510256aa6e06a3af3526af342f9315d3b4cf24a972fde1ae8fa5be34dd282f3d56ed1acf968140d9d3015647ff6aacc1d65485c06aac8294bb336ef6b68d47a0e8607b4c8c59f9e3a283634d426c4029c46e6789667c608c610a2a5e234e4413863ffc587d27c3bd37d96f1384e9a5e7fff69c51f8a240d3f01df85c530e9f1b2ea031d45b47dc203f0d5bfb6 

    cp = 0xb523cd55033b7e640dfdec9cd72f070ded607e4e46e173d08f86cc44a4fa417f4dc2b5246dccfe8c625ec54650021346ae9887fa2d267bafdbbd485ac4be1c7bc083fc8f29e2dd2ea38238be05ae3a1a21326d90d205e41a29467155717060682551b2b8933b2e7806e284ca057c0a8f5b14e4b759ce7f57da43c513c5ae5678885dc0d991338f0118f5471488e543a0785ba67e65e425fbd3afa0937392898b1f44dedacb785a8cabff47ae71c7f9967efeb3d1bbad9ab49a69ee3bb37add897d54a4bb81e96addfcc39aad2290733b7f932d8bfd4ab4e36a26d0d17d5fbc5ae231c775e540920b711d7cd4efe65ecfa106d5390dd37952950c0a382dcb3638

    e = 2053

    x = PolynomialRing(ZZ.quo(n*ZZ), 'x').gen()

    p = x^e - c
    q = (x+1)^e - cp

    a = p
    b = q

    msg = 0

    run = True
    while run:
	    r = a % b
	    if r == 0:
		    print('found :', rp) #rp == gcd(p,q)
		    
		    coeffs = rp.coefficients(x)
		    
		    #FinalPoly = A*X + B == X - M => M = B * inv(A) mod N 
		    msg = ( -pow(coeffs[1], -1, n) * coeffs[0] ) % n #<=> tmp = modinv(coeffs[1], n); msg = -tmp * coeffs[0] % n
		    
		    run = False
	    rp = r
	    a,b = b,r

	#int to ascii
	int2str(msg)


############################################################################################################

'''#SVP
from numpy import linalg as LA
import numpy as np

from fpylll import *
'''
def svp(): #fpylll
	a = int("0x9f9c4d766db8cefd", 16)
	p = int("0xba917061a1065abf", 16)

	#hardcode
	l0 = [1, pow(a, 1,p), pow(a, 2,p), pow(a, 3,p), pow(a, 4,p), pow(a, 5,p), pow(a, 6,p)]

	l1 = [0, p, 0, 0, 0, 0, 0]
	l2 = [0, 0, p, 0, 0, 0, 0]
	l3 = [0, 0, 0, p, 0, 0, 0]
	l4 = [0, 0, 0, 0, p, 0, 0]
	l5 = [0, 0, 0, 0, 0, p, 0]
	l6 = [0, 0, 0, 0, 0, 0, p]

	A = IntegerMatrix.from_matrix([l0, l1, l2, l3, l4, l5, l6])

	L = int("2**55.88", 16)
	if (A[0].norm() < L):
		#### fpylll ####
		a = SVP.shortest_vector(A)

		sol = []
		for e in a:
			sol.append(hex(e).lstrip('0x'))

	return sol

#svp()

######################### xA = b ##############################
'''
CVP(A, b) does not work immediatly

CVP(Phi, b||0) avec Phi = A || Id_m renvoie b||x => extraire x
'''

def sis(coeff): #fpylll
	m = 36
	n = 10
	p = 257
	a = 2798646620
	q = 10
	
	B = [141, 67, 161, 111, 175, 177, 128, 34, 38, 204]

	b = []
	for e in B:
		 b.append(e*coeff)
	
	b = tuple(b)
	b = b+(0,0)*18
	
	modulus = 4294967297

	# Initializing A
	A = [ [0 for i in range(m+10)] for j in range(n+36) ]

	# Building A
	for i in range(m):
		for j in range(n):
			A[i][j] = coeff * pow(a, n*i + j, modulus)

		# Concatenate A with identity
		for k in range(10, n+36):
			if (i == k-10):
				A[i][k] = 1	
	
	#'''
	# We work modulus p
	for i in range(36, m+10):
		for j in range(n):
			if (i == j+36):
				A[i][j] = coeff*p
	#'''
	'''	
	print(A)
		
	print(len(A))
	print(len(A[0]))
	'''
	##### LLL #####

	A = IntegerMatrix.from_matrix(A)
	
	#LLL.reduction(A, None, 0.99, 0.51, None, None, 0, 1)
	LLL.reduction(A)
	
	res = CVP.closest_vector(A, b)
	return res

#coeff = 46**5	
#print(sis(coeff))

############################################################################################################


from Crypto.PublicKey import RSA
key_encoded='''-----BEGIN PUBLIC KEY-----
MIIBIDANBgkqhkiG9w0BAQEFAAOCAQ0AMIIBCAKCAQEA7WuPBtWwMNJH6XmU3R0E
8wFpgdveiOW495mjCBfOj3CXxuigul91eeunRud9G9NifQHKt5PRGFcQqcTYgiMQ
Rxu2ASb5W/wNrDjuJsT7Z0LTmsSfDoKOTNr0L0kbEMyO3YR0r4jm7TGl994gt4Qm
Copwr24CBmBkkaMrK6PWzFBlwxg1qQu58zd5YiyEALSkbEokfWrr8CI1hhNTy5zt
lIAncGYzlENiv2+CNzDcfZ9uXCtSqINZ/P0YLb0ddTAJAwpkiKoWDQIccV41Bl1U
LA9vYu9f0f0C1PgOzvavKr9pXL/fri/naQNO8oScmpLQEJcKNznA8377ca/fFqva
TQIBAw==
-----END PUBLIC KEY-----'''


pubkey = RSA.importKey(key_encoded)


import binascii
from pkcs1 import *

def pkcs1Exploit():

	N = pubkey.n
	e = pubkey.e #3
	
	HASH_ID = b'010\r\x06\t`\x86H\x01e\x03\x04\x02\x01\x05\x00\x04 ' #magic string
	"""
	hh = b"3031300d060960864801650304020105000420" #ok
	print(i2osp(int(hh, 16), len(hh) // 2)) #ok
	"""

	#id = "ppti_server_room"
	b0 = '0001'
	b0 = i2osp(int(b0, 16), len(b0) // 2)

	
	nbr = 10
	b1 = nbr*'FF'
	b1 = i2osp(int(b1, 16), len(b1) // 2)
	
	b2 = '00'
	b2 = i2osp(int(b2, 16), len(b2) // 2)
	
	b3 = HASH_ID
	
	M = b"MY WANTED MESSAGE"
	#M = binascii.hexlify(M)
	#M = i2osp(int(M, 16), len(M) // 2)
	h = sha256(M).hexdigest()
	
	b4 = sha256(M).digest() #i2osp(int(h, 16), len(h) // 2)
	
	#### 
	nbr = 256 - len(b0+b1+b2+b3+b4)
	junk = nbr*'FF'#b'junk'
	b5 = i2osp(int(junk, 16), len(junk) // 2)
	####
	
	bloc = b0+b1+b2+b3+b4+b5
	print(len(bloc))
	
	
	x = os2ip(bloc)
	print(x < N)
	
	#Cube root of X
	X = find_invpow(x,e) #power e done in the verification
	
	#Signature
	S = i2osp(X, key_length(N))
	print(S)
	print(S.hex())
	
	print(rsa_pkcs_verify(N, e, M, S))

#pkcs1Exploit()


############################################################################################################

'''
# Montgomery batch inversion (2**24 inversions into 1 inversion and 3 * 2**24 multiplications)

    ## 1.
    
    a = a
    ab = a*b % N
    abc = ab*c % N
    abcd = abc*d % N

    ## 2. X <-- (abcd)^(-1)

    X = (c *) modinv(abcd, N)

    ## 3. All inverses retrieved in 2 multiplications for each :

    inv_d = X*abc % N
    X = X*d % N

    inv_c = X*ab % N
    X = X*c % N

    inv_b = X*a % N
    X = X*b % N

    inv_a = X % N
'''


def multi_inv(values, N):
    partials = [1]
    for i in range(len(values)):
        partials.append( ( (partials[-1] * values[i]) % N ) or 1)
        partials.append( ( (partials[-1] * values[i]) % N ) or 1)
        
    inv = modinv(partials[-1], N)
    outputs = [0] * len(values)
    for i in range(len(values), 0, -1):
        outputs[i-1] = ( (partials[i-1] * inv) % N ) if values[i-1] else 0
        inv = ( (inv * values[i-1]) % N ) or 1
    return outputs
    
    
def rfid():
    e = int(0x10001)
    N =  int('0x1ea982ba8f01d5e03163b0409a554484b8e145af768a8d3e66b84c9723d8604a33bd7c52033def81adfaf49beaa4f0f2b3b92370efb88f07665c5c35afdfd94752eacc4cf24ff3b96954ff391abaf39108df0cf11c26567ac2aa408143038ed11d53172667b95637a7cd3d6bc8972e6a4d7a503730db2af935d3baf8d5a5465d', 16)
    
    size = 2**10 #24
    
    # Construction
    print("[-] Construction...")
    D = {}
    for i in range(1,size):
        x = pow(i,e,N)
        D[x] = i
       
    print("[+] Construction done")
    
    # example
    c = pow(530*999, e, N)
    print('solution to find : ', c)
    
    #strc = '135ef757c5436bd48eb6fb9cba57fcea2a034c27e4f2286f1621f0c1476f04374dd025ac6223ea24ada771dde2b9e082ab6277bea61f09aa96d57c9a4f507fd099f2d587dbb72c23a8cfd34acb3f68a544bc202ec364361ef6549379a2065ff75f52661a4bc1c39ddb400b1d079a1d81cb63bfb3f6cf9753253c0294fc1c933c'

    #c = int(strc,16)
    
    # Test if c can be cut in two int of the same length
    if 2*(c // 2) == c:
        print("[+] Cipher can be cut in two ints of the same length")
        
        
        print("[-] Computing inversions...")
        L = []
        Prod = []
        tmp = 1
        for j in range(1, size): #24
            x = pow(j,e,N)
			
            L.append(x)
            
			## 1
            tmp = (tmp * x) % N
            Prod.append(tmp)
		
		## 2 
        X = ( c * modinv(Prod[-1], N) ) % N  # Prod[-1] == prod of all the x
        
        ## 3
        invs = []
        for j in range(2, size):
            inv = (X*Prod[-j]) % N
            invs.append(inv)
            X = (X * L[-j+1]) % N
        invs.append(X) # inverse of the first number
            
        #invs = multi_inv(L, N)
        print("[+] Inversions done")
        
        # Research
        print("[-] Research...")

        for j in range(size-1):
        
			#y = c * modinv(x, N) % N # invs[j]
            y = invs[j]
            
            #pb with j !!!
            try:
                i = D[y]
                K = pow(i*j, e, N) #K = A*B => RSA(K) = RSA(A) * RSA(B) % N
                print('K = ', i, '*', j)
                
                if (K == c):
                    print('[+] found K :', K)
                    print('K = ', i, '*', j)
                    return i*j

            except KeyError:
                continue
        
    else:
        print("[-] Cipher cannot be cut in two int of the same length")
        return


import hmac

def rfid_ok():
    e = int(0x10001)
    N =  int('0x1ea982ba8f01d5e03163b0409a554484b8e145af768a8d3e66b84c9723d8604a33bd7c52033def81adfaf49beaa4f0f2b3b92370efb88f07665c5c35afdfd94752eacc4cf24ff3b96954ff391abaf39108df0cf11c26567ac2aa408143038ed11d53172667b95637a7cd3d6bc8972e6a4d7a503730db2af935d3baf8d5a5465d', 16)
     	    
    D = {}
    print("[-] Construction...")
    for i in range(1, 2**24):
        if not i % 1000000: # <=> if j % 1000000 == 0:
        	print(i)
        	
        x = pow(i, e, N)
        D[x] = i
        
    print("[+] Construction done")
    
    print("[-] Research...")
    
    while True:
        
        strc = str(input('cipher : '))
        
        c = int(strc,16)
    	
        #c = pow(35001*62980, e, N)
    	
        for j in range(1, 2**24):
            x = pow(j, e, N)
            y = c * modinv(x, N) % N
		    
            if not j % 100000:
            	print(j)
		    
            try:
                i = D[y]
                #tmp = pow(i*j,e,N)
                print('[+] found K :')
                print('[+] K = ', i, '*', j)
                K = i*j

                key = K.to_bytes(6, 'big')
                z = 0x00000000
                z = z.to_bytes(4, 'big')

                a = hmac.new(key, z, 'sha256')
                res = '00000000' + a.hexdigest()

                print(c == pow(K, e, N))

                print(res)
                return
	        
            except KeyError:
                continue

   
#rfid_ok()


############################################################################################################


#implementation of Brent's algorithm for detecting a cycle - credit to Wikipedia

def brent(f, x0):
    power = 1 
    lam = 1
    tort = x0
    hare = f(x0)
    
    while tort != hare:
        if power == lam:
            tort = hare 
            power *= 2
            lam = 0
        hare = f(hare)
        lam += 1

    # Find the position of the first repetition of length λ
    tort = hare = x0
    for i in range(lam):
    # range(lam) produces a list with the values 0, 1, ... , lam-1
        hare = f(hare)
    # The distance between the hare and tortoise is now λ.

    mu = 0
    while tort != hare:
        tort = f(tort)
        hare = f(hare)
        mu += 1

    return #lam, mu #lam : cycle length // mu : cycle start index



#name = 'Louisa'


def F_to_collision(x):
    x = name + str(x)
    #print('key : ', x)

    #print('key_hex : ', x.encode().hex()) #keys in hex at the end

    h = sha256(x.encode()).hexdigest()[:14]
    #print('hash :', h)
    
    return int(h, 16)
    
    
#x0 = 1#rd.randint(1,2**28)

#beg = time.time()

#brent(F_to_collision, x0)

#print('time : ', time.time() - beg)
#print('x0 : ', x0)


#easier and not that long...
def dist_SHA256_Coll(ntries):
    str1 = b'Nametest1 {%d}' 
    str2 = b'Nametest2 {%d}'
    seen = {}
    for n in range(1,ntries):
        h = sha256(str1 % n).hexdigest()[:14]
        seen[h] = n
    for n in range(1,ntries):
        h = sha256(str2 % n).hexdigest()[:14]
        if h in seen:
            print(str1 % seen[h])
            print(str2 % n)

#dist_SHA256_Coll(2**28)

''' These strings have same sha256 56 first bits 
b'Louisatest1 {244615218}'
b'Louisatest2 {166696267}'
'''

############################################################################################################

'''
Hash function
M : bytes encoded string
'''
def SDBM(M):
    h = 0
    i = 0
    
    M = M.hex()
    
    while i < len(M):
        m = int(M[i:i+2], 16)
        h = (65599*h + m) % 2**128
        i+=2
        
    return hex(h)[2:]


#NOTE FOR SDBM :

#For ASCII input
'''
M = """PLAIN RML-Code"""
print(SDBM(M.encode()))
'''

#OR for hex input : 
'''
M = """CIPHER RML-Code"""
M = "".join(M.split('\n'))
M = base64.b16decode(M) #to_bytes
print(SDBM(M))
'''

from aes import *
'''
>>> key = base64.b16decode("00000000000000000000000000000000")
>>> a = AES(key)
>>> plaintext = base64.b16decode("80000000000000000000000000000000")
>>> ciphertext = base64.b16decode("3ad78e726c1ec02b7ebfe92b23d9ec34", casefold=True)
>>> assert a.encrypt(plaintext) == ciphertext
'''

'''
M : Message ASCII
K : Key, hex value (str)
'''
def tarMAC(M, K):
    H = SDBM(M).upper()
    
    key = base64.b16decode(K) #to_bytes
    a = AES(key)
    plaintext = base64.b16decode(H) #to_bytes
    tag = a.encrypt(plaintext)
    return tag.hex()

'''
M = "Les sanglots longs\nDes violons\nDe l'automne\nBlessent mon coeur\nD'une langueur\nMonotone."
K = "00000000000000000000000000000000"
print(tarMAC(M, K))
'''


'''
M : bytes encoded message
K : seed (128 bits)
'''
def TLCG(M, K):
    cipher = ''
    i = 0
    M = M.hex()
    
    x = K
    while i < len(M):
        Y = x // 2**120

        P = int(M[i:i+2], 16)

        C = hex(Y ^ P)[2:].zfill(2)

        cipher += C
        x = (x * 47026247687942121848144207491837523525) % 2**128
        i+=2

    return cipher.upper()

'''
M = "Les sanglots longs\nDes violons\nDe l'automne\nBlessent mon coeur\nD'une langueur\nMonotone."
M = M.encode()
K = 0x71d05909e13748ff733ffccfbfbf40eb
print(TLCG(M, K))
'''



#from fpylll import * 

'''
C : hex string
P : bytes encoded string
'''
def TLCG_attack(P, C):
    P = P.hex()
    C = C.hex()
    
    i = 0
    z = []
    n = len(P) #16 ?
    while i < n:
        p = int(P[i:i+2],16)
        c = int(C[i:i+2],16)
        y = hex(p ^ c)[2:].zfill(2) # y[i] == x[i] // 2**120 ... x == TLCG seed
        #z.append(y)
        z.append( int(y, 16) * 2**120 ) # z = [ y[i] * 2**120 for i in range(n) ]
        i+=2
    
    #build the lattice
    ''' 
    [ 1,      a,   a**2,   a**3, ..., a**(n-1) ]
    [ 0, 2**128,      0,      0, ...,        0 ]
    [ 0,      0, 2**128,      0, ...,        0 ]
    [ 0,      0,      0, 2**128, ...,        0 ]
    [..........................................]
    [ 0,      0,      0,      0, ...,   2**128 ]
    '''
    
    n = 80 # max 85-88
    # Initializing A
    A = [ [0 for i in range(n)] for j in range(n) ]
	
    #first row
    A[0][0] = 1
    a = 47026247687942121848144207491837523525
    for i in range(1, n):
        A[0][i] = a**i

    #diagonal
    for i in range(1, n):
        A[i][i] = 2**128
        

    ##### LLL #####

    A = IntegerMatrix.from_matrix(A)

    LLL.reduction(A)

    res = CVP.closest_vector(A, z)
    return res
        
        
'''
P = ""
P = P.encode()


C = """3d8a065b3ccba48c74c53c4b9d7dbbbcc1b3ba9c8ae689687a31517b3bd79814b133a3b6671124e8
bae01efba766c3ebd9f6908e65000995a99a873cd085bfeada8db8e6565539b1ffb3f703f386b41c
2d37f2bb5b351c""".upper()
C = "".join(C.split('\n'))
C = base64.b16decode(C)

test = TLCG_attack(P, C)
K = test[0]
print(K)
'''

'''
P = P.encode()

C = "".join(C.split('\n'))
C = base64.b16decode(C)

tmp = TLCG_attack(P, C)
K = tmp[0]
print(K)
'''

def TLCG_decypt(C, K):
    plaintext = ''
    i = 0
    C = C.hex()
    
    x = K # seed
    while i < len(C):
        Y = x // 2**120 # MSB

        c = int(C[i:i+2], 16) # cipher

        P = hex(Y ^ c)[2:].zfill(2) # plain

        plaintext += P
        x = (x * 47026247687942121848144207491837523525) % 2**128 # state
        i+=2

    plaintext = int(plaintext,16)
    
    return int2str(plaintext)

'''
C = "".join(C.split('\n'))
C = base64.b16decode(C) # to bytes

K = 36449343409391193238006126174636460803 #obtained with TLCG_attack

print(TLCG_decrypt(C, K))
'''



import string
'''
Search for a suffix to forge a message (RML code) 
which has the same SDBM_hash 'h' than a valid message 
Thus we obtain a correct tarMAC
'''
def SDBM_attack(M_1, M_2, asciNum):

    h = SDBM(M_1.encode())
    print(h)
    h_int = int(h, 16)
    
    n = 21 # n must be >= 16
    
    #building the lattice
    """
     [ 1,  0,  0,  0, ...,  0,  0,  a**(n-1) ]
     [ 0,  1,  0,  0, ...,  0,  0,  a**(n-2) ]
     [ 0,  0,  1,  0, ...,  0,  0,  a**(n-3) ]
     [.......................................]
     [ 0,  0,  0,  0, ...,  1,  0,         a ]
     [ 0,  0,  0,  0, ...,  0,  1,         1 ]
     [ 0,  0,  0,  0, ...,  0,  0,    2**128 ]
    """
    
    # Initializing G
    G = [ [0 for i in range(n+1)] for j in range(n+1) ]

    a = 65599

    #last column
    for i in range(1, n+1):
        G[i-1][n] = a**(n-i)
        
    #diagonal
    for i in range(n):
        G[i][i] = 1
        
    #just before last element
    G[n-1][n] = 1

    #last element
    G[n][n] = 2**128
    
    target = [asciNum for i in range(n+1)] # someties ok with 80 (pay attetion on non-ascii chars in the suffix)...
    
    SDBM_M_2 = int(SDBM(M_2.encode()),16)
    
    verif = ( h_int - SDBM_M_2 * a**n ) % 2**128 # n == len(suffix)
    print('verif :', verif)#hex(verif).lstrip('0x'))

    target[-1] = verif ###### !!!!!!!!!! ########

    target = tuple(target)
    ##### LLL #####
    G = IntegerMatrix.from_matrix(G)

    LLL.reduction(G)

    res = CVP.closest_vector(G, target) # res == (S_0, ..., S_{n-1}, h) 
    print(res)
    
    res = list(res)

    for j in range(128):
        res[-2] = j    
        rsum = 0
        for i in range(n):
            rsum = (rsum + a**(n-1-i) * res[i]) % 2**128
        
        print('sum :', rsum)
        
        if rsum == verif:
            print(j)
            res[-2] = j # found by bruteforce to have 'verif' as the last element 
            break
    

    print(res)
    
    suffix = ''
    for e in res[:-1]:
        suffix +=  chr(e) #avoid int2str function for only one char, only useful for big messages 
    
    print(suffix)
    print('len_suffix :', len(suffix))
    
    # Check if there are unprintable ascii chars
    filtered_suffix = filter(lambda x: x in string.printable, suffix)
    cpt = 0
    for e in filtered_suffix :
        cpt+=1
   
    print('len_filtered_suffix :', cpt)
    
    if cpt != n:
        asciNum += 1
        SDBM_attack(M_1, M_2, asciNum)
        return
    
    else:
        print('Good printable chars, with asciNum =', asciNum)

    
    #SDBM(suffix) =?= (h - SDBM(M_2) * 65599**len(suffix)) % 2**128
    SDBM_suffix = SDBM(suffix.encode())
    print('SDBM_suffix :', SDBM_suffix)
    
    print('verif :', hex(verif).lstrip('0x'))
    
    print('equal ? :', verif == int(SDBM_suffix, 16))
    
    #SDBM((M_2+suffix).encode()) =?= h
    final = M_2+suffix
    
    #print(final)
    SDBM_final = SDBM(final.encode())
    print('SDBM(M2 || suffix) :', SDBM_final)
    print('SDBM(M1) :', h)
    print('equal ? :', h == SDBM_final)

'''
with open('path/to/file.rml', 'r') as f:
        M_1 = f.read()

with open('path/to/file.rml', 'r') as f:
        M_2 = f.read()


M_1 = M_1.rstrip('\n') # '\n' added by f.read()
M_2 = M_2.rstrip('\n') # '\n' added by f.read()
asciNum = 64

SDBM_attack(M_1, M_2, asciNum)
'''

###############################################################################################


"""
AES-128-CTR used as a PRG

K : AES key (bytes)
IV : CTR IV (int)
"""
def aes_128_ctr(Kaes, IV, output_size):
    
    # 1.
    CTR = IV
    
    # 2. (init)
    P = CTR.to_bytes(16, 'big')
        
    a = AES(Kaes)
    #a = AES.new(Kaes, AES.MODE_ECB) #from Crypto.Cipher import AES 
    
    C = a.encrypt(P)
    
    i = 0
    CTR += 1
    
    # 3.
    res = ''
    count = 0
    
    while count < output_size:
        
        if i == 16 :
            P = CTR.to_bytes(16, 'big')
            C = a.encrypt(P)
            i = 0
            CTR += 1
        
        else:
            res += hex(C[i]).lstrip('0x').zfill(2)
            i += 1
            
            count += 1
    
    print('CTR : ', CTR)
    return res
    '''
    pl = "Les sanglots longs\nDes violons\nDe l'automne\nBlessent mon coeur\nD'une langueur\nMonotone."
    res = bytes(XOR(bytes.fromhex(res), pl.encode()))
    return res#.hex()
    '''
    
'''
Kaes = bytes(16)
IV = 0
print(aes_128_ctr(Kaes, IV, 3))
'''


from Crypto.Cipher import AES
from Crypto.Util import Counter
"""
AES-128-CTR encryption
using Cryptodome lib
"""
def aes_ctr_enc(Kaes, IV, data):
    ctr = Counter.new(128, initial_value=IV)
    cipher = AES.new(Kaes, AES.MODE_CTR, counter=ctr)
    ct = cipher.encrypt(data)
    return ct.hex()

'''
data = b"Les sanglots longs\nDes violons\nDe l'automne\nBlessent mon coeur\nD'une langueur\nMonotone."
Kaes = bytes(16)
IV = 0
print(aes_ctr_enc(Kaes, IV, data))
'''
    
"""
AES-128-CTR decryption
using Cryptodome lib
"""
def aes_ctr_dec(Kaes, IV, data):
    ctr = Counter.new(128, initial_value=IV)
    cipher = AES.new(Kaes, AES.MODE_CTR, counter=ctr)
    ct = cipher.decrypt(data)
    return ct#.decode()

'''
data = "2a8c38f49ceb425ce4238e2aea5844403f91f68a9f0d10175f107138ca944f1e66a8b6e901c3d7fd9e46a7b333de9b0b84f0c4df6926364dd79ee69ae1f9cba407777f446e1ff5b447fcd3d9a299e68fa722d6767feb07"
data = bytes.fromhex(data)
Kaes = bytes(16)
IV = 0
print(aes_ctr_dec(Kaes, IV, data)[:-32])
'''

###############################################################################################


########################################### FEISTEL ###########################################


'''
text : bytes (64 bits)
K : bytes (64 bits)
'''
def feistel(text, K, w, debug=0):

    k = []
    for i in range(0,16,4):
       k.append(K[i:i+4])

    L = text[0:16]
    R = text[16:32]
    
    if w == 'Enc':
        
        h0 = sha256(k[0]+R).digest()[:16] # [:16] == 8 bytes == 64 bits...
        L = XOR(L, h0)
        
        h1 = sha256(k[1]+L).digest()[:16]
        R = XOR(R, h1)
        #'''
        h2 = sha256(k[2]+R).digest()[:16]
        L = XOR(L, h2)
        
        h3 = sha256(k[3]+L).digest()[:16]
        R = XOR(R, h3)
        #'''

        #ciphertext = L+R
    
    elif w == 'Dec':
        #'''
        h3 = sha256(k[3]+L).digest()[:16]
        R = XOR(R, h3)

        h2 = sha256(k[2]+R).digest()[:16]
        L = XOR(L, h2)
        #'''
        h1 = sha256(k[1]+L).digest()[:16]
        R = XOR(R, h1)
        
        h0 = sha256(k[0]+R).digest()[:16] # [:16] == 8 bytes
        L = XOR(L, h0)
         
        #plaintext = L+R
        
    else:
        exit("'Enc' to encrypt, 'Dec' to decrypt")
    
    return L,R


''' # Feistel Network ok (Encryption and Decryption)
key = b'0001020304051617'

XL = b'0000000000000000'
YL = b'1111111111111111'
R = b'08090a0b0c0d0e0f'

X = XL+R
Y = YL+R 

X4L,X4R = feistel(X, key, 'Enc')
print(X4L+X4R)
t = bytearray(b'\x98\xadHBFE[(\xed\xe2(?<\xdeo\x7f\x1d\xd8uqV\x97\x9e\xe63w\xc8\xfb\x08F\xad\xbc')
X4L,X4R = feistel(t, key, 'Dec')
print(X4L+X4R)
'''

'''
Last feistel round
'''
def lastFeistelRound(text, keyTest, w):

    L = text[0:16]
    R = text[16:32]
    
    if w == 'Enc': # ???
        h3 = sha256(keyTest+L).digest()[:16]
        R = XOR(R, h3)
        
    elif w == 'Dec': # ok
        h3 = sha256(keyTest+L).digest()[:16]
        R = XOR(R, h3)
    

    return L,R


########################################### Attack ###########################################

"""
1) Request the encryption of 2 messages X = (XL, R) and Y = (YL, R) 
to obtain X4 = (X4L, X4R) and Y4 = (Y4L, Y4R). 

2) For each possibility of key k4 (1 to 2**16)
Reverse a feistel round on X4 and Y4 using the possibility of k4 to obtain X3 and Y3

3) Form the new cipher of the 3 turn distinguisher as Z3 = (Y3L, Y3R xor XL xor YL)

4) Do a new feistel round using the k4 possibility on Z3 to get Z4

5) Pass Z4 to the deciphering oracle to decipher Z4 into Z = (ZL, ZR)

6) If good value of k4 to undo and redo the last round 
   Then property (ZR == Y3L xor X3L xor R) satisfied & X3 and Y3 correct
    
   Else bad value of k4 --back to--> 2)
"""

def feistelAttack():

    key = b'0001020304051617'
    
    #1)
    XL = b'0000000000000000'
    YL = b'1111111111111111'
    R = b'08090a0b0c0d0e0f'

    X = XL+R
    Y = YL+R 

    X4L,X4R = feistel(X, key, 'Enc') # oracle
    X4 = X4L+X4R
    #print('X4 :', X4)

    Y4L,Y4R = feistel(Y, key, 'Enc') # oracle
    Y4 = Y4L+Y4R
    #print('Y4 :', Y4)

    #2)
    '''for trial in range(1, 2**16):'''
    trial = '1617'                          # The correct 'trial' for the 4th key
    '''trial = hex(trial)[2:].zfill(4)'''
    
    keyTest = trial.encode()

    # Reverse Feistel last round
    X3L,X3R = lastFeistelRound(X4, keyTest, 'Dec') # not oracle
    Y3L,Y3R = lastFeistelRound(Y4, keyTest, 'Dec') # not oracle
    '''
    X3 = X3L+X3R 
    print('X3 :', X3)
    Y3 = Y3L+Y3R 
    print('Y3 :', Y3)
    '''

    #3)
    #The new cipher for the 3 round distinguisher
    Z3 = Y3L + XOR(Y3R, XL, YL)
    #print('Z3 : ', Z3)

    #4)
    # re-make the last round
    Z4L,Z4R = lastFeistelRound(Z3, keyTest, 'Enc') # not oracle
    Z4 = Z4L+Z4R

    #5)
    ZL,ZR = feistel(Z4, key, 'Dec') # oracle

    #6) Verification
    '''
    if ZR == XOR(Y3L, X3L, R):
        print(ZR)
        print(XOR(Y3L, X3L, R))
        break
    '''
    print(ZR)
    print(XOR(Y3L, X3L, R))


feistelAttack()



"""
i

tmp = f_i(R_i, k_i)
L_i+1 = R_i
R_i+1 = xor(tmp, L_i)


*i = 0*

h0 = f_0(R_0, k_0)
L_1 = R_0
R_1 = xor(tmp, L_0)

*i = 1*

h1 = f_1(R_1, k_1)
L_2 = R_1
R_2 = xor(tmp, L_1)

*i = 2*

h2 = f_2(R_2, k_2)
L_3 = R_2
R_3 = xor(tmp, L_2)

*i = 3*

h3 = f_3(R_3, k_3)
L_4 = R_3
R_4 = xor(tmp, L_3)

out_L = R_4
out_R = L_4
"""




'''
Test of the Distinguisher for a 2-round Feistel network
'''
def test2roundFeistelDistinguisher():

    key = b'0001020304050607'
    
    XL = b'0000000000000000'
    YL = b'1111111111111111'
    R = b'08090a0b0c0d0e0f'

    X = XL+R
    Y = YL+R 

    X2L,X2R = feistel(X, key, 'Enc')

    Y2L,Y2R = feistel(Y, key, 'Enc')

    # Should be equal for a 2-round Feistel network => ok
    print(XOR(X2L,Y2L)) 
    print(XOR(XL,YL))


#test2roundFeistelDistinguisher() # (have to change the feistel network above)
    
